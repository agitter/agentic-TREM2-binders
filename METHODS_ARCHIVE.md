# TREM2 binder comparison — archival methods record

Everything needed to reproduce `fig_network_{binders,nonbinders}_tm80_{structure,sequence}.png`.
Every free parameter is listed with the value used, why, and what else was tried.

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
(9 Strong, 1 null). Each collection is labelled by its own vendor's call; no cross-collection
rubric was imposed.

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

## 8. Known limitations to carry into any caption

1. **Cross-group edges are not plotted.** Splitting binders from non-binders drops
   47 of 207 structural edges (23 %) and 27 of 110 sequence edges (25 %). 11 of 19
   non-binders have ≥1 fold link to a binder; with cross edges included, **no non-binder is
   isolated**. The sparse look of the non-binder panel is partly an artifact of the split.
2. **Edge counts are symmetrization-dependent** (§3).
3. **`binder_final` is an adjudicated rubric, not a raw measurement** (`DATA_NOTES` §17).
   Within the 18 Anthropic non-binders: 7 did not express at Adaptyv, 5 were called binders
   by Twist and overruled. This is Anthropic's encoding, adopted unchanged — but it means
   "non-binder" spans true non-binders and expression failures.
4. **Muni is a curated set** (9/10 Strong) versus a complete unfiltered Anthropic campaign.
   Hit rates are not comparable across collections; only sequence/structure comparisons are.
5. Muni's Boltz2 co-folds come from Muni's pipeline at unknown version/settings; Anthropic's
   are seed-best-of-five by `ipsae_min`. Same predictor family, not an identical protocol.
