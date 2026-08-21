# TREM2 binder comparison — archival methods record

Everything needed to reproduce the figures in this repository.
Every free parameter is listed with the value used, why, and what else was tried.

The record covers two analyses that share a pipeline:

* **Two-collection** (100 designs: Anthropic + Muni) — sections 0–8.
  Artifacts unsuffixed: `tm_all.npy`, `sid_all.npy`, `fig_network_{binders,nonbinders}_tm80_*`.
* **Three-collection** (200 designs: + the Adaptyv × MUNI hackathon) — section 9.
  Artifacts suffixed `_hack`, so nothing in the two-collection analysis is overwritten
  and its figures stay byte-reproducible.

---

## 0. Environment

| Component | Version / build |
|---|---|
| Foldseek | `9-427df8a` (`427df8a6b5d0ef78bee0f98cd3e6faaca18f172d`), linux-avx2 |
| MMseqs2 | `15-6f452` (`6f45232ac8daca14e354ae320a4359056ec524c2`), linux-avx2 |
| py2Dmol | commit `70e9b96` (used for `best_view` only; see §6) |
| Python | 3.12.3 |
| Biopython 1.88 · numpy 2.4.4 · scipy 1.17.1 · pandas 3.0.2 | |
| networkx 3.6.1 · scikit-learn 1.8.0 · matplotlib 3.10.8 · gemmi 0.7.5 | |

**Required inputs (5 files) — sufficient to reproduce every figure:**
`designs.fasta`, `design_summary.csv`, `TREM2.zip` (Anthropic, CC BY 4.0);
`proteinbase_collection_muni-proteina-complex-auto-research.csv`, `muni_files.zip`
(Muni, via ProteinBase — observe the ProteinBase licence).

`muni_files.zip` is the output of `download_muni.sh`; bundling it means the pipeline needs
no network access to ProteinBase. `download_muni.sh` is retained for re-fetching.

**Consulted but not required to run:** `DATA_NOTES.md` and `LOOKUP_TABLES.md` from the
Anthropic release informed the caveats in section 8; no script reads them.

**Verified end to end**: a clean-room run from these five files alone reproduced all four
network figures with matching md5 sums, and `tm_all.npy` / `sid_all.npy` bit-identically.

---

## 1. Sequence extraction

Anthropic: 90 TREM2 records pulled from the 1,440-record FASTA by `target=TREM2`.
Workflow parsed from the header stem `<model>_<campaign>_trem2_rank<NN>`.
Verified: all 90 sequences equal `design_summary.csv.sequence` exactly.

Muni: 10 rows; `evaluations` column is JSON-inside-CSV. **The file carries a UTF-8 BOM**,
so it must be read with `encoding='utf-8-sig'` (pandas strips it silently, `csv` does not).

---

## 2. Structures

```bash
unzip -q TREM2.zip '*/insilico/predicted_boltz2_1to1.cif' -d cifs_b2
./download_muni.sh proteinbase_collection_muni-proteina-complex-auto-research.csv
```

**Choice — which structure file.** `TREM2.zip` ships `designed.cif` (the design model) and
ten `predicted_<pred>_1to1.cif` co-folds. Muni supplies a Boltz2 co-fold. We use
**`predicted_boltz2_1to1.cif`** so both collections come from the same predictor at the same
1:1 stoichiometry. Using `designed.cif` would confound "different design" with
"design model vs prediction". Measured impact: design model vs Boltz2 for the same design is
mean TM 0.964; the all-vs-all matrices correlate at Spearman 0.85, so conclusions are stable
either way, but the cross-collection comparison is only valid with matched predictors.

**Choice — chain selection.** The binder is chain **B** in `predicted_boltz2_1to1.cif` and
chain **A** in `designed.cif`; Muni's binder is chain B (160 aa) against a 156-aa TREM2.
Chains are therefore matched **by exact sequence equality to the design sequence**, never by
chain ID. All 100 resolved uniquely.

---

## 3. Structural similarity — Foldseek

```bash
foldseek easy-search all_binder_pdb all_binder_pdb aln.tsv tmp \
  --alignment-type 1 \
  --exhaustive-search 1 \
  -e 10000 \
  --max-seqs 2000 \
  --format-output "query,target,alntmscore,lddt,fident"
```

| Flag | Value | Why |
|---|---|---|
| `--alignment-type` | `1` | TM-align mode. Gives a true TM-score; the default (2, 3Di+AA) returns a Foldseek score, not TM. |
| `--exhaustive-search` | `1` | Forces all-vs-all. Without it the prefilter drops distant pairs and the matrix is sparse — fatal for a network where absence of an edge is meaningful. |
| `-e` | `10000` | Effectively disables E-value filtering, for the same reason. |
| `--max-seqs` | `2000` | Above the 100 structures, so no per-query truncation. |
| output field | `alntmscore` | The TM-score. `lddt`/`fident` captured but unused in the figures. |

Verified complete: 10,000 rows = 100 × 100, every pair present.

**Choice — symmetrization.** `alntmscore` is direction-dependent (mean |q→t − t→q| = 0.022,
max 0.42) because TM normalises by chain length. We take the **elementwise max**, equivalent
to normalising by the shorter chain, matching the sequence-identity convention in §4.

*Alternatives tested:* mean and min both give 122 binder edges at TM ≥ 0.80 instead of 150.
All conclusions hold under all three (binder–binder > binder–non-binder TM; nearest-neighbour
density p = 0.023–0.024; sequence overlap = 38 in every case). **The edge count is
convention-dependent and should be quoted as such.**

---

## 4. Sequence similarity — Biopython

`Bio.Align.PairwiseAligner`, all 4,950 pairs:

| Parameter | Value | Why |
|---|---|---|
| `mode` | `global` | Full-length comparison; these are complete designed proteins. |
| `substitution_matrix` | BLOSUM62 | Standard for protein identity work. |
| `open_gap_score` | `-11` | EMBOSS/BLASTP default pairing. |
| `extend_gap_score` | `-1` | As above. |
| end gaps | free (score 0) | Semi-global: lengths run 56–160 aa, so terminal-gap penalties would distort short-vs-long pairs. |
| identity denominator | **length of the shorter sequence** | Consistent with the TM max convention. Alignment-length normalisation would inflate identity for pairs with long gapped regions. |

Cross-checked against MMseqs2 (`easy-search -s 7.5 --max-seqs 1000 -e 1000`;
`easy-cluster --min-seq-id {0.3,0.5,0.7} -c 0.5 --cov-mode 1`), which independently found
zero cross-collection and zero cross-workflow clusters.

Six random pairs were re-derived from raw sequences and matched stored values to 0.01%.

---

## 5. Binder / non-binder assignment

**Anthropic: `design_summary.csv → binder_final`, used verbatim. No re-derivation.**
End-to-end check: 354 True release-wide, 262 Mythos Preview / 92 Opus 4.8 — exactly the
README's published figures. TREM2 subset: 72 True, 18 False, 0 null.

**Muni: `binding_strength == "Strong"`** from the collection CSV — Muni's own field
(9 Strong, 1 null). The single null is `brisk-crane-granite`, which the CSV records as
`expressed: False` / `binding: False` in both replicates: a measured expression failure, and
therefore a non-binder on the same basis as the 7 Anthropic designs that failed to express at
Adaptyv and carry `binder_final = False`. Each collection is labelled by its own vendor's
call; no cross-collection rubric was imposed.

---

## 6. Node rendering (`render_nodes.py`)

py2Dmol draws to an HTML5 canvas via JavaScript and has no headless raster path; the Chrome
download CDN was unreachable in this environment. **`py2Dmol.viewer.best_view()` is used
directly** for the canonical orientation; the draw stage is reimplemented on the pipeline in
py2Dmol's `technical_readme.md` (rotate → orthographic project → painter's-algorithm z-sort →
grid shadow map → outline).

| Parameter | Value | Note |
|---|---|---|
| superposition reference | `mythos_preview_single_target_trem2_rank18` | Global medoid (max row-sum of TM). All 100 aligned to it via the Foldseek CIGAR, then Kabsch. |
| shared view frame | one `best_view` from the reference | **Essential**: per-structure `best_view` gives each protein its own axes, so similar folds would look different. |
| spline | Catmull–Rom, 6 subdivisions/residue | Smooth tube. |
| bitmap | 300 × 300 px, tube width 12.0 | |
| shadow grid | 26 × 26, occlusion weight 0.55 | Depth cue. |
| colour | one hue per workflow; N→C light→dark ×0.45, depth shade 0.45–1.0 | Hue carries workflow *and* 3-D read in one channel. |
| draw extent | `xlim/ylim = ±0.53` | Margin so the tube fits inside the node footprint. |

Verified rigid: superposition changes no internal Cα–Cα distance by >3.6e-14 Å;
|det(R) − 1| < 3e-15 (rotation only, no reflection or scaling).

---

## 7. Network construction and layout (`fig_network.py`)

```bash
python3 fig_network.py 0.80 40      # argv: TM threshold, % identity threshold
```

| Parameter | Value | Why / alternatives |
|---|---|---|
| `TM_THR` | **0.80** | TM ≥ 0.5 is the same-fold convention; 0.80 is a strict "near-identical backbone" bar chosen for legibility. Sweep: 0.65→786 edges (Q=0.33), 0.70→481, 0.75→284 (Q=0.43), **0.80→150 (Q=0.62)**, 0.82→119. Both 0.75 and 0.80 were rendered and compared. |
| `ID_THR` | **40 %** | Above the 20–35 % twilight zone, i.e. where a homology search calls two sequences related. Sweep: 25%→1744 edges, 30%→906, 35%→331, **40%→78**, 50%→32. 40 % chosen as the defensible convention, not the most flattering. |
| Louvain `resolution` | 1.0 | Default; no tuning. |
| Louvain / layout `seed` | 7 | Fixed. |
| community graph layout | `spring_layout(k=1.6, iterations=600)` | Gives each module its own territory before node placement. |
| main force loop | 1500 iters, cooling 1.0→0.25 | Edge spring `L = min_d × 1.02`, stiffness ×0.045 scaled by TM; repulsion 0.0045 out to `min_d × 1.7`; collision 0.55; gravity 0.0030. |
| settle pass | 900 iters | Gravity 0.010 pulls detached components into the concavities; edge springs preserve module shape. |
| final de-overlap | ≤800 iters | Hard separation only. |
| `min_d` (node diameter, data units) | 0.108 binders / 0.215 non-binders | Set per panel so 81 and 19 nodes are each legible. |
| image zoom | `2·ring_px/300 · (72/DPI) · 1.16` | `OffsetImage` scales by `zoom × dpi/72`; the `72/DPI` term is required or nodes render 1.74× too large. |
| DPI | 125 | |

**Determinism.** `nx.connected_components()` and `louvain_communities()` return *sets* of
string node names; set iteration order varies with Python's per-process hash randomisation,
which made the layout irreproducible across runs. Fixed by imposing a total order
(`sorted(c)`, then by `(-len, first member)`) before any iteration. Two fresh processes now
produce bit-identical PNGs. No reported statistic was affected — only node placement.

**Layout affects appearance only.** Edges are built from the similarity matrices and are
independent of layout. Verified across all four figures: 0 below-threshold edges drawn,
0 above-threshold edges missing.

---

## 7b. Figures NOT regenerated with the final settings

`fig3_structure.png` was built from `tm.npy` — the **`designed.cif`** matrix (90 Anthropic
designs only), not the Boltz2 matrix used everywhere else. Over the same 90x90 pairs the two
differ by mean 0.034 TM (187 vs 176 edges at TM >= 0.80). It is superseded by panel A of
`fig_master.png` and should not be published alongside the final figures without relabelling.

`fig1_identity_heatmap.png` is sequence-only and unaffected by the structure-source choice,
but predates the final legend wording.

All `fig_network_*` files other than the four `*_tm80_{structure,sequence}.png` are drafts
produced before the determinism fix in section 7 and before the final legend/threshold
choices. They are not reproducible from this repository and should be deleted.

---

## 7c. Backbone non-independence — the analysis treats every design as one node

**The nodes in these networks are not statistically independent, and this analysis does not
correct for it.** Every figure and statistic in this repository counts one node per delivered
*sequence*. Claude's pipeline generates a backbone first and then designs one or more
sequences onto it, so two designs can be different sequences threaded onto the same fold.
Shanehsazzadeh et al. flag this directly ("we count designs per sequence, although sequence
variants of one backbone are not independent"); release-wide they report 809 distinct
backbones behind 1,315 tested designs, and 200/809 (24.7 %) binders when counting one
sequence per backbone.

`provenance_summary.csv → root_backbone_id` records the backbone of each design.

**Scale on TREM2.** The 90 designs come from **61 distinct backbones**. 29 backbones
contributed exactly 2 designs each (never more), so **58 of 90 designs have a sibling**.
Siblings are near-duplicates in both dimensions:

| Within the 29 sibling pairs | median | range |
|---|---|---|
| TM-score | 0.98 | 0.68 – 0.99 (28/29 ≥ 0.80) |
| Sequence identity | 68 % | 49 – 87 % (29/29 ≥ 40 %) |
| Binder length | identical in 29/29 pairs | |

Both members bind in 26 of 29 pairs, against 18.6 expected if the two were independent
draws at the 80 % TREM2 hit rate. Siblings differ mainly in how much in-silico optimization
they received (e.g. `predict-and-redesign x5` vs `x3` on the same backbone).

**Effect on the plotted networks.** Sibling pairs are over-represented in every edge set, and
they dominate the overlap statistic:

| Signal (binder panel) | total | sibling | share |
|---|---|---|---|
| Structural edges (TM ≥ 0.80) | 150 | 25 | 17 % |
| Sequence edges (identity ≥ 40 %) | 78 | 26 | 33 % |
| **Both** (the headline "38") | 38 | **25** | **66 %** |

16 of the 20 highest-TM binder pairs are siblings.

**The headline conclusion survives de-duplication, and strengthens:**

| Fold-sharing pairs that are also sequence-alike | |
|---|---|
| As published (81 nodes) | 38/150 = **25 %** |
| Excluding sibling pairs | 13/125 = **10 %** |
| One design per backbone (55 binder nodes) | 10/54 = **19 %** |

Among genuinely independent backbones only about 1 in 10 fold-sharing pairs shares detectable
sequence similarity, so the structure-without-sequence convergence is not a duplicate artifact
— duplicates were masking its strength. Likewise 125 of the 150 structural edges join designs
from different backbones, so the fold clustering is mostly real.

**Suggested caption wording:** *"The 90 Anthropic designs derive from 61 distinct backbones;
17 % of structural edges join sibling sequences of a shared backbone. Restricting to the 55
independent backbones, 19 % of fold-sharing pairs remain sequence-alike (10 % when sibling
pairs are simply excluded)."*

**Asymmetry.** The Muni collection carries no backbone identifiers, so the same de-duplication
cannot be applied to it. Its 10 designs have their own redundancy — mean pairwise identity
31.7 %, one pair at 95.5 % — arising from related MCTS search configurations rather than a
shared backbone (per the muni autoresearch report). Any de-duplicated comparison is therefore
corrected on the Anthropic side only.

Per-pair data: `backbone_families.csv` (29 pairs with TM, identity, length, outcome).

---

## 8. Known limitations to carry into any caption

0. **Nodes are sequences, not backbones** — see §7c. 58 of the 90 Anthropic designs share a
   backbone with another design, 17 % of structural edges and 66 % of the structure-plus-
   sequence overlap are sibling pairs, and no statistic here is corrected for this.
1. **Cross-group edges are not plotted.** Splitting binders from non-binders drops
   47 of 207 structural edges (23 %) and 27 of 110 sequence edges (25 %). 11 of 19
   non-binders have ≥1 fold link to a binder; with cross edges included, **no non-binder is
   isolated**. The sparse look of the non-binder panel is partly an artifact of the split.
2. **Edge counts are symmetrization-dependent** (§3).
3. **`binder_final` is an adjudicated rubric, not a raw measurement** (`DATA_NOTES` §17).
   Within the 18 Anthropic non-binders: 7 did not express at Adaptyv, 5 were called binders
   by Twist and overruled. This is Anthropic's encoding, adopted unchanged — but it means
   "non-binder" spans designs that were expressed and failed to bind, and designs that
   failed to express. **The two collections are consistent on this point**: Muni's single
   non-binder, `brisk-crane-granite`, records `expressed: False` and `binding: False` in
   both replicates (its null `binding_strength` reflects that there is no affinity to
   report, not that it went untested), so expression failure is a measured negative outcome
   under both vendors' rubrics and is counted the same way on both sides.
4. **Muni is a curated set** (9/10 Strong) versus a complete unfiltered Anthropic campaign.
   Hit rates are not comparable across collections; only sequence/structure comparisons are.
5. Muni's Boltz2 co-folds come from Muni's own autoresearch pipeline (Boltz-2, ranked by
   default ipSAE); Anthropic's are seed-best-of-five by `ipSAE_min`. Same predictor, but not
   an identical seed-selection protocol.
6. **Hackathon affinities** are Adaptyv BLI/SPR, directly comparable to Muni and
   *not* to Anthropic (see item 7).  The hackathon is in-silico pre-filtered
   (top 100 of 141 by Boltz-2 ipSAE), a third selection regime; hit rates are not
   comparable across any of the three collections.
7. **Affinities are not comparable across the Anthropic/Muni collections.** Anthropic's TREM2 plate had
   dissociation too slow to measure, so kinetic and steady-state fits disagreed by >30-fold
   and `kd_nM_final` values are steady-state figures or bounds; six TREM2 binders sit at the
   ~100 pM assay floor. Muni's values are Adaptyv BLI replicate means. Cite the same-plate
   comparison in Shanehsazzadeh et al. instead: the tightest Muni design was re-synthesised
   as a control on Anthropic's TREM2 plate and bound below 100 pM, as did the best design of
   each of the three Claude campaigns.

---

## 9. Three-collection extension — the Adaptyv × MUNI hackathon

### 9.1 Input

`proteinbase_collection_adaptyv-x-muni-hackathon-ai-agents-vs-humans.csv`
(ProteinBase download endpoint; **ODC-BY**, attribution required) plus
`hackathon_boltz2.zip`, the 100 `boltz2_structure_prediction` CIFs fetched with
`download_muni.sh`.  Same six-column schema as the Muni collection and the same UTF-8
BOM, so the same `encoding='utf-8-sig'` applies.

100 designs, 141 submitted of which the top 100 by Boltz-2 ipSAE went to the lab — so
this collection is **in-silico pre-filtered**, a third and distinct selection regime
alongside Muni's curated 10-of-13,916 and Anthropic's unfiltered ranks 1–30.

Verified on receipt: all 100 CIFs parse, exactly 2 chains each, binder is chain **B**
in all 100, and every binder chain equals its CSV sequence exactly.  No design_id
collides with the Muni collection.

### 9.2 Binder definition

`binding_strength ∈ {Strong, Medium, Weak}` = 37 binders, which reproduces the
collection's published headline.  This agrees with the `binding` boolean on **100 of
100** designs, so the two definitions are interchangeable here.  `binding_strength` is
retained as a four-level column for downstream use.  Expression rate 89/100.

Design classes are heterogeneous by intent: Miniprotein 79, Other 15, Nanobody 4,
Peptide 1, scFv 1; lengths 12–245 aa.  All classes are retained.  Note Anthropic's
prompt explicitly places single-domain antibody and VHH formats out of scope, so the
antibody scaffolds are hackathon-only **by construction** and are not evidence about
tool or format choice.

### 9.3 TREM2 construct heterogeneity — a mandatory correction

The hackathon used **two different TREM2 constructs**, both numbered from 1:

| construct | n | composition |
|---|---|---|
| 126-mer | 92 | `TG` + TREM2(1–115) + `GTKHHHHHH` (linker + His6) |
| 156-mer | 8 | TREM2, byte-identical to Muni's target |

Residue *N* in 126-numbering is residue *N−2* in 156-numbering (consistent at 115 of
116 aligned positions; the 9 exceptions are the purification tag).  Epitope residues
are therefore mapped onto 156-numbering with a −2 shift for the 92, and contacts at
126-numbering ≥ 118 are **dropped** as contacts to the tag rather than the antigen
(9 contacts across 4 designs).  Every mapping is verified against the residue identity
in the structure; after the fix, **0 mismatches**.

**This correction is load-bearing.** Pooling raw numbers puts the hackathon consensus
at residues 56/53/54/58, two off from the 51/53/54/56/57/71 hotspots seen in the other
collections, which reads as a *different* epitope.  After correction the top residues
are 54/53/56/51/71/52/26 and hits in the known hotspot set rise from 410 to 628 (+53%).

Across all three collections the TREM2 construct is now four-way heterogeneous
(109 / 110 / 126 / 156).  Structural comparison is unaffected — only the binder chain
is used, matched by sequence.

### 9.4 Predictor provenance — a three-way asymmetry

| collection | co-folded structures |
|---|---|
| Hackathon | Adaptyv ProteinBase Boltz-2 pipeline |
| Muni | Adaptyv ProteinBase Boltz-2 pipeline |
| Anthropic | Anthropic's own pipeline, seed-best-of-five by `ipSAE_min` |

Hackathon↔Muni is therefore a like-for-like comparison; any Anthropic-versus-others
structural difference is **partly confounded with predictor protocol** and must be
reported as such.

### 9.5 What was rerun, and what was not

Rerun on 200 designs: binder-chain extraction, `all_hack.fasta`, Foldseek all-vs-all
(40,000 rows) and the CIGAR run, Biopython pairwise (19,900 pairs, 8-way multiprocessed
in `pairwise_hack.py`), `tm_hack.npy`, `sid_hack.npy`, epitopes, node bitmaps, figures.
Settings are identical to sections 3 and 4 in every respect.

Not rerun: Anthropic provenance and backbone-family analysis, wet-lab analysis, and the
four two-collection network figures.

**Cross-check:** the 100×100 submatrices of the new matrices reproduce the verified
originals — TM exactly (max difference 0.0), sequence identity in 4,949 of 4,950 pairs.
The single exception, `MUNI_deep-gecko-maple` / `MUNI_noble-boar-reed`, has exactly
**two co-optimal alignments** at score 48.0 giving 22.56% or 25.56%; `align()[0]`
returns a different representative.  Both values are far below the 40% edge threshold,
so no figure or statistic changes.  Aligner settings are identical — verified by
re-running the original configuration, which now also returns 25.56%.

### 9.6 Novelty audit (`fig_novelty_audit_hack.png`)

Anthropic's prompt lists the hackathon collection as **reference #11** and stages the
ProteinBase collections as the known-binder corpus ("02 ProteinBase").  Its novelty
gate reads:

> REJECT at >60% identity over >50% coverage to UniRef90 or the binder corpus, OR at
> ≥30% gapped local identity over ≥40 aligned residues OR TM-score ≥0.5 **to any target
> or control chain**

The TM-score criterion is scoped to target and control chains — its stated purpose is
catching target-mimic protomers.  **The known-binder corpus was tested on sequence
only.**  The gate is emulated with MMseqs2 (`-s 7.5`, `fident`/`qcov`/`tcov`).

| | Anthropic (90) | Muni (10) |
|---|---|---|
| TM ≥ 0.50 to some hackathon design | 90 (100%) | 10 (100%) |
| TM ≥ 0.80 | **51 (57%)** | 4 (40%) |
| TM ≥ 0.90 | 18 | 1 |
| max sequence identity to corpus | 58.3% | 41.7% |
| rejected by the >60%/>50% gate | **0** | 0 |

The Muni panel is **descriptive only**: Muni published no equivalent novelty rule, so it
is the same measurement without a stated criterion to audit against.

**Do not read this as imitation of what worked.**  Anthropic's median best TM to a
hackathon *binder* is 0.790 versus 0.797 to a *non-binder* (paired Wilcoxon p = 0.061,
if anything the wrong direction), and the single closest match (TM 0.968) is to an
RFpeptides design with no binding call.  The defensible statement is a shared fold
vocabulary, not selective copying.  Timing is consistent with exposure but does not
establish it — ULID asset timestamps place all hackathon assets before all Muni ones,
and all three campaigns targeted the same epitope on the same protein with overlapping
tools, so convergence and derivation are not separable from this data.

### 9.7 Three-collection figures

`fig_network_hack_tm80_{structure,sequence}.png` — 200 nodes, node hue = campaign,
`min_d` 0.072, `width_in` 20.0; layout and edge semantics otherwise as section 7.
**Edge-colour balance is inverted here:** only Anthropic publishes `root_backbone_id`,
so 432 of 608 structural edges fall in the "backbone not recorded" class.  That class is
therefore drawn recessive (alpha 0.42, 0.72× width, light cyan) so the two informative
classes stay legible — a presentation change only; class membership is unchanged.

`fig_heatmap_hack_{bycollection,clustered}.png` — TM and identity matrices side by side
at matched ordering.  `_bycollection` blocks by campaign with UPGMA within each block;
`_clustered` uses one global UPGMA on TM with campaign in the side bar only.  The 5
antibody-scaffold designs are outlined in pink: they form a perfectly isolated clique
(all 10 pairs TM ≥ 0.80, mean 0.918) with **max TM 0.578 to any of the other 195**.

`fig_epitope_hack.png` — hackathon and Anthropic share an identical top-8 contacted
residue set (26, 51, 52, 53, 54, 56, 57, 71).  Median cross-campaign epitope overlap
(Hackathon–Muni 0.55, Muni–Anthropic 0.50, Hackathon–Anthropic 0.47) meets or exceeds
each collection's internal overlap (Hackathon 0.52, Anthropic 0.46).

---
