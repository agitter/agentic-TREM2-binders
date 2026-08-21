# Comparison of Muni and Anthropic agentic TREM2 binder design results
This repository compares AI agents' designs of TREM2 protein binders that were recently shared by [Muni](https://muni.bio/research/closing-the-autoresearch-loop) and [Anthropic](https://www.anthropic.com/research/Claude-accelerates-protein-design).
All artifacts except those in `/data` are from Claude Opus 5.
The code may not run locally, and external software dependencies are not provided.
[`METHODS_ARCHIVE.md`](METHODS_ARCHIVE.md) provides an overview of the analyses, limitations, and software Claude used.

## Workflow overview
```
designs.fasta ──filter target=TREM2──► 90 seqs, 3 workflows ─┐
design_summary.csv ──binder_final────► 72 bind / 18 not ─────┤
                                                             ├─► 100 designs
Muni CSV (utf-8-sig) ──binding_strength──► 9 Strong / 1 null ─┤
muni_files.zip ──10 Boltz2 complexes────────────────────────┘

TREM2.zip ──predicted_boltz2_1to1.cif──► chain matched BY SEQUENCE ──► 100 binder PDBs
                                          (chain B here, chain A in designed.cif)
        │
        ├─► Foldseek TM-align, exhaustive ──► 100×100 TM ──max-symmetrized──► tm_all.npy
        │                                  └─► CIGAR backtrace ──► superposition
        │
        └─► Biopython BLOSUM62 semi-global ──► 100×100 % identity ──► sid_all.npy
                                             (cross-checked vs MMseqs2)

tm_all ──► Louvain modules ──► community-seeded force layout ──► node positions (FROZEN)
        └► edges at TM ≥ 0.80          ─┐
sid_all └► edges at identity ≥ 40%     ─┴─► same positions, two edge layers

structures ──► superpose onto global medoid ──► py2Dmol best_view (ONE shared frame)
            ──► tube render, hue = workflow ──► node bitmaps
                                                        │
                              binders ×2 panels ◄───────┴───────► non-binders ×2 panels
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
