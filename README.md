# Comparison of Muni and Anthropic agentic TREM2 binder design results
This repository compares AI agents' designs of TREM2 protein binders that were recently shared by [Muni](https://muni.bio/research/closing-the-autoresearch-loop) and [Anthropic](https://www.anthropic.com/research/Claude-accelerates-protein-design).
All artifacts except those in `/data` are from Claude Opus 5.
The code may not run locally, and external software dependencies are not provided.
[`METHODS_ARCHIVE.md`](METHODS_ARCHIVE.md) provides an overview of the analyses, limitations, and software Claude used.

## Workflow overview

```
STAGE 1 -- two collections (100 designs)         artifacts unsuffixed, figures unchanged
=======================================================================================
designs.fasta --filter target=TREM2--> 90 seqs, 3 workflows -+
design_summary.csv --binder_final----> 72 bind / 18 not -----+
                                                             +--> 100 designs
Muni CSV (utf-8-sig) --binding_strength--> 9 Strong / 1 null -+   (null = expression
muni_files.zip --10 Boltz2 complexes-------------------------+    failure, a measured
                                                                  negative)
TREM2.zip --predicted_boltz2_1to1.cif--> chain matched BY SEQUENCE --> 100 binder PDBs
                                          (chain B here, chain A in designed.cif)
        |
        +--> Foldseek TM-align, exhaustive --> 100x100 TM --max-symmetrized--> tm_all.npy
        |                                  \-> CIGAR backtrace --> superposition
        |
        \--> Biopython BLOSUM62 semi-global --> 100x100 % identity --> sid_all.npy
                                             (cross-checked vs MMseqs2)

tm_all --> Louvain modules --> community-seeded force layout --> node positions (FROZEN)
        \> edges at TM >= 0.80         -+
sid_all \> edges at identity >= 40%    -+--> same positions, two edge layers

structures --> superpose onto global medoid --> py2Dmol best_view (ONE shared frame)
            --> tube render, hue = workflow --> node bitmaps
                                                        |
                              binders x2 panels <-------+-------> non-binders x2 panels
                                          \---> combined x2 panels, hue = binder status


STAGE 2 -- add the hackathon (200 designs)                 artifacts suffixed  _hack
=======================================================================================
hackathon CSV (utf-8-sig) --binding_strength in {Strong,Medium,Weak}--> 37 bind / 63 not
        |                   (agrees with the `binding` boolean 100/100)
        |                   classes: 79 miniprotein . 15 other . 4 nanobody
        |                            . 1 scFv . 1 peptide
        |
        +--> hackathon_boltz2.zip --100 CIFs, binder = chain B--+
        |                                                       +-> + STAGE 1's 100 PDBs
        |                                                       |
        |                                                       \-> 200 binder PDBs
        |
        \--> interface_residues --> TWO TREM2 constructs found --> NUMBERING FIX
                                   92 x 126-mer = TG + TREM2(1-115) + GTKHHHHHH
                                    8 x 156-mer = TREM2, identical to Muni's target
                                          |
                                          +- shift 126-numbering by -2  (verified by
                                          |  residue identity: 0 mismatches)
                                          +- DROP contacts at >=118 (His-tag, 9 contacts)
                                          \-> epitopes_hack.json  -- without this the
                                              hackathon looks like a DIFFERENT epitope
                                              (56/53/54/58 vs 51/53/54/56/57/71)

200 PDBs -> Foldseek exhaustive -> 40,000 rows -> 200x200 TM -> tm_hack.npy
         |                      \-> CIGAR --> superposition (fs_hack_sup)
         \-> Biopython, 19,900 pairs, 8 procs -> 200x200 % id -> sid_hack.npy
                    |
                    \- CROSS-CHECK vs stage 1: TM exact, identity 4,949/4,950
                       (1 pair has 2 co-optimal alignments; below threshold, no effect)

tm_hack --> same layout engine --> node positions (FROZEN across both panels)
         \> edges TM >= 0.80  -> 608         hue = CAMPAIGN, not workflow
sid_hack \> edges id >= 40%   -> 628         "not recorded" class now 432/608 -> recessive
                    |
                    +---> fig_network_hack_tm80_{structure,sequence}
                    +---> fig_heatmap_hack_{bycollection,clustered}
                    |       \- 5 antibody scaffolds = isolated clique
                    |          (10/10 pairs TM >= 0.80; max 0.578 to the other 195)
                    \---> fig_epitope_hack  -- hackathon & Anthropic share top-8 residues

Anthropic prompt ref #11 = the hackathon collection, staged as known-binder corpus
        |
        \--> MMseqs2 emulation of the novelty gate (>60% id over >50% cov)
                    |
                    \---> fig_novelty_audit_hack   Anthropic 51/90 at TM >= 0.80
                         -- gate rejects 0 --      Muni 4/10 . max id 58.3% / 41.7%
                         (gate's TM criterion applied only to target/control chains,
                          never to the binder corpus -> structure was unscreened)
```

## Third-party resources
The `/data` directory contains inputs from third parties.

### muni / Adaptyv
- `proteinbase_collection_muni-proteina-complex-auto-research.csv` from [Proteinbase](https://proteinbase.com/collections/muni-proteina-complex-auto-research)
- `muni_files.zip` files downloaded from the URLs in the above csv
- `proteinbase_collection_adaptyv-x-muni-hackathon-ai-agents-vs-humans.csv` from [Proteinbase](https://proteinbase.com/collections/adaptyv-x-muni-hackathon-ai-agents-vs-humans)
- `hackathon_boltz2.zip` files downloaded from the URLs in the above csv, only the Boltz2 structure .cif files.

[License](https://proteinbase.com/download):
> This work used Proteinbase by Adaptyv Bio under ODC-BY license

### Anthropic
- `design_summary.csv` from [Hugging Face](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/blob/main/data/tables/design_summary.csv)
- `designs.fasta` from [Hugging Face](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/blob/main/data/designs/designs.fasta)
- `TREM2.zip` zipped from [Hugging Face](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/tree/main/data/designs/TREM2)

Additional documentation files were provided to Claude as context.
This uses data release v1.0.

[License](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/blob/main/data/LICENSE.md):
> Data and documentation: CC BY 4.0. Scripts included in the archives: MIT. Third-party material (structure-prediction outputs, reference sequences, vendor report images, reagent names) keeps its own terms; see LICENSE.md in each archive. Please cite as given in [CITATION.cff](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/blob/main/data/CITATION.cff).
