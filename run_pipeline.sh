#!/usr/bin/env bash
# Exact command sequence for the TREM2 fold-space network figures.
# Companion to METHODS_ARCHIVE.md, which explains every parameter choice.
# Run from a directory containing the five third-party inputs:
#   designs.fasta  design_summary.csv  TREM2.zip  muni_files.zip
#   proteinbase_collection_muni-proteina-complex-auto-research.csv
# plus, for the three-collection extension (section 9):
#   proteinbase_collection_adaptyv-x-muni-hackathon-ai-agents-vs-humans.csv
#   hackathon_boltz2.zip
# Set UP=<dir> if they are not in the current directory.
set -euo pipefail

UP="${UP:-.}"
BIN="${BIN:-$HOME/bin}"
export PATH="$PATH:$BIN/foldseek/bin:$BIN/mmseqs/bin"

# ---------------------------------------------------------------- 0. tooling
# foldseek 9-427df8a  |  mmseqs2 15-6f452  |  py2Dmol @70e9b96
curl -sL https://github.com/steineggerlab/foldseek/releases/download/9-427df8a/foldseek-linux-avx2.tar.gz | tar xz -C "$BIN"
curl -sL https://github.com/soedinglab/MMseqs2/releases/download/15-6f452/mmseqs-linux-avx2.tar.gz     | tar xz -C "$BIN"
pip install biopython gemmi networkx pyarrow statsmodels --break-system-packages -q
git clone --depth 1 https://github.com/sokrypton/py2Dmol.git && pip install -e py2Dmol --break-system-packages -q

# ---------------------------------------------------------------- 1. inputs
# Muni assets (10 Boltz2 complexes + ESMFold + PAE + BLI curves = 60 files).
# The bundled zip already contains them; download_muni.sh re-fetches from ProteinBase
# if you would rather pull them fresh (needs network).
unzip -q -o "$UP/muni_files.zip" -d .          # -> ./muni_files/
# ./download_muni.sh "$UP/proteinbase_collection_muni-proteina-complex-auto-research.csv"

# Anthropic Boltz2 co-folds (binder = chain B) — NOT designed.cif (binder = chain A)
unzip -q -o "$UP/TREM2.zip" '*/insilico/predicted_boltz2_1to1.cif' -d cifs_b2

# ---------------------------------------------------------------- 2. parse
python3 parse_anth.py     # 90 TREM2 records from designs.fasta, stratified by workflow
python3 parse_muni.py     # 10 Muni records; NOTE encoding='utf-8-sig' (file has a BOM)
python3 build_dataset.py --stage seqs
# -> all.fasta (100) and all_binder_pdb/*.pdb (100)
# chains matched to designs BY SEQUENCE, never by chain ID

# ---------------------------------------------------------------- 3. structure
# --alignment-type 1  : TM-align mode, returns a real TM-score
# --exhaustive-search : all-vs-all; without it the prefilter sparsifies the matrix
# -e 10000            : disable E-value filtering, same reason
foldseek easy-search all_binder_pdb all_binder_pdb fs_all/aln.tsv fs_all/tmp \
  --alignment-type 1 \
  --exhaustive-search 1 \
  -e 10000 \
  --max-seqs 2000 \
  --format-output "query,target,alntmscore,lddt,fident"
# expect 10000 rows (100x100)

# residue correspondence for superposition + rendering (CIGAR backtrace)
foldseek easy-search all_binder_pdb all_binder_pdb fs_sup/sup.tsv fs_sup/tmp \
  --alignment-type 1 \
  --exhaustive-search 1 \
  -e 10000 \
  --max-seqs 2000 \
  --format-output "query,target,alntmscore,qstart,qend,tstart,tend,cigar"

# Ca coordinates + alignment table consumed by render_nodes.py
python3 superpose.py                      # -> ca.json, aln_sup.json

# ---------------------------------------------------------------- 4. sequence
# Biopython global BLOSUM62, gap -11/-1, free end gaps, identity / shorter sequence
python3 pairwise.py                       # -> gid.npy, names.json

# assemble every shared artifact the figure scripts consume
python3 build_dataset.py --stage matrices
# -> order_all.json, tm_all.npy (max-symmetrised), sid_all.npy,
#    anth_full.json, archetypes.csv

# independent cross-check (found zero cross-collection / cross-workflow clusters)
mmseqs easy-search  all.fasta all.fasta mm/res.m8 mm/tmp -s 7.5 --max-seqs 1000 -e 1000 \
  --format-output "query,target,fident,alnlen,evalue,bits,qcov,tcov"
for ID in 0.3 0.5 0.7; do
  mmseqs easy-cluster all.fasta "mm/clu$ID" "mm/tmp$ID" --min-seq-id "$ID" -c 0.5 --cov-mode 1
done

# ---------------------------------------------------------------- 5. figures
# argv = TM threshold, % sequence-identity threshold
python3 fig_network.py 0.80 40
#   -> fig_network_binders_tm80_structure.png     (81 nodes, 150 edges)
#   -> fig_network_binders_tm80_sequence.png      (81 nodes,  78 edges, 38 shared)
#   -> fig_network_nonbinders_tm80_structure.png  (19 nodes,  10 edges)
#   -> fig_network_nonbinders_tm80_sequence.png   (19 nodes,   5 edges,  3 shared)

# alternative threshold rendered for comparison
python3 fig_network.py 0.75 40

# ------------------------------------------- 6. three-collection extension (METHODS 9)
# extra input: the Adaptyv x MUNI hackathon collection (ODC-BY) + its 100 Boltz2 CIFs
unzip -q -o "$UP/hackathon_boltz2.zip" -d hack_cifs        # flat, <id>__boltz2_*.cif
python3 build_dataset_hack.py --stage seqs
#   -> all_hack.fasta (200), all_binder_pdb_hack/ (200), epitopes_hack.json
#      epitope residues remapped to 156-numbering; His-tag contacts dropped

for O in "aln.tsv|query,target,alntmscore,lddt,fident|fs_hack"          "sup.tsv|query,target,alntmscore,qstart,qend,tstart,tend,cigar|fs_hack_sup"; do
  IFS='|' read -r OUT FMT DIR <<< "$O"
  foldseek easy-search all_binder_pdb_hack all_binder_pdb_hack "$DIR/$OUT" "$DIR/tmp" \
    --alignment-type 1 --exhaustive-search 1 -e 10000 --max-seqs 4000 \
    --format-output "$FMT"
done                                    # expect 40000 rows each (200x200)

python3 pairwise_hack.py                # -> sid_hack.npy  (19,900 pairs, 8 processes)
python3 build_dataset_hack.py --stage matrices   # -> order_hack.json, tm_hack.npy

# novelty gate emulation: MMseqs2 vs the hackathon corpus, per Anthropic's prompt
mmseqs easy-search q_anth.fasta db_hack.fasta mm_nov/anth.m8 mm_nov/tmp_anth \
  -s 7.5 --max-seqs 500 -e 10000 \
  --format-output "query,target,fident,alnlen,qcov,tcov,evalue,bits"
mmseqs easy-search q_muni.fasta db_hack.fasta mm_nov/muni.m8 mm_nov/tmp_muni \
  -s 7.5 --max-seqs 500 -e 10000 \
  --format-output "query,target,fident,alnlen,qcov,tcov,evalue,bits"

python3 fig_network_hack.py 0.80 40     # 200-node structure + sequence networks
python3 fig_hack_analyses.py            # novelty audit + three-collection epitope map
python3 fig_heatmap_hack.py             # TM / identity heatmaps, two orderings
