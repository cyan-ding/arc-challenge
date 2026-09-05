# vcc2026 metrics: what the spec actually says

Extracted from the two official documents shipped in the `cell-eval2` repo under
`docs/vcc2026_metrics/` — `vcc2026-metrics.pdf` (full, with reasoning) and
`vcc2026-metrics-brief.pdf` (definitions and measured values). Both are dated **2026/08/19**.
Measured values in them come from the three official reference bundles (contexts A, B, C) built
with **cell-eval2 0.15.0, competition rule_version 3**. We have **cell-eval2 0.16.0** pinned in
`pyproject.toml`; the CLI surface below was checked against that install.

This file answers the four questions `docs/claude-recipe.md` left open, records the anchor
values, and notes where it corrects our other docs. Where anything here disagrees with
`.cursor/rules/project-rules.mdc`, the rules file wins until we reconcile it — but the two
corrections flagged in §7 are worth folding into the rules file.

## 0. The one thing that resolves most of the confusion: two contracts

The spec is explicit that there are **two different file contracts** and they deliberately differ:

| | upload (`vcc prep`, what we send) | scoring (`cell_eval2`, what the scorer is handed) |
|---|---|---|
| `non-targeting` rows | **forbidden** — rejected | **required**, label must match the reference |
| cells per perturbation | **exactly 400** | unconstrained |
| perturbation column | `target_gene` | `pert_col`, whatever the run declares |

The platform bridges them: it takes the held-out real control cells out of the reference,
concatenates them onto our prediction, and stamps the scorer's perturbation column from
`target_gene`. So every piece of `cell-eval2` documentation that says "your prediction must
include the control label" is describing the *scoring* contract, not what we upload. When we
score locally we reproduce the bridge ourselves: `control_source: real` (the `vcc2026` preset
default) means the real side supplies the controls, and our prediction file carries none.

## 1. Shared pipeline parameters

All from `configs/vcc2026.yaml`, confirmed by loading `EvalConfig.from_preset("vcc2026")` in
0.16.0:

| parameter | value | used by |
|---|---|---|
| `bulk_target_sum` (TS) | 5 × 10⁴ | `pds`, `mse` (pseudobulk) |
| `target_sum` (per cell) | 10⁶ (CPM) | all four DE metrics |
| DE test | Wilcoxon rank-sum, two-sided, normal approximation | DE metrics |
| `filter_gene_min_cpm_cell` | 5 CPM, **read from the reference control cells only** | DE metrics |
| `p_adj_threshold` | 0.05, Benjamini–Hochberg, **per perturbation**, applied after the filter | DE metrics |
| ε (fold change) | 10⁻⁹ | DE metrics |
| `mean_calc` | arithmetic | DE fold change |
| `nan_lfc_policy` | mask | DE |
| `REACH_PURITY_FLOOR` (P₀) | 0.9 | `reach` |
| `min_gate_size` | 10 | `nmae` |
| replicate anchor | 5 disjoint half-splits, base seed 0 | all six |
| `discrimination` | cosine, rank denominator n−1, midrank ties, `exclusion_scope: panel` | `pds` |
| `max_counts_per_cell` | 10⁶ | validation |

Evaluation data per context: 300 target constructs (one per gene), a pooled NTC from 46
constructs, every perturbation **downsampled to exactly 400 cells and to a median depth of
20,000 UMI/cell**, on 18,533 genes.

Fold change on both sides: `lfc = log2((mean_pert + ε) / (mean_ctrl + ε))` where the means are
arithmetic means of the CPM-normalized values over the group's cells. The comparison group for
*our* DE table is the **measured control**, never anything we submit.

The 5 CPM gate reading the control alone is load-bearing: it makes the tested gene universe
identical for every perturbation and for both sides, so the mere presence of a (perturbation,
gene) pair in the table can't leak the sign. An earlier rule that admitted genes above the cutoff
in *either* group leaked direction for ~1 % of the table (100 % of those pairs were positive).
Genes silent in the control and strongly induced by a knockdown are **not scored at all** — a
known cost.

## 2. Per-metric definitions and what they imply for us

### `pds_cosine` — perturbation discrimination

Profile: sum counts over the group's cells, normalize to TS = 5e4, `log1p`. Effect δ = profile −
reference control profile, on both sides. For each perturbation p, cosine distance from our δ̂ₚ to
all 300 measured δ_q. Score = 1 − k/(n−1) where k is the (midrank) position of the true match.
Context value = mean over 300.

**Gene exclusion is the whole panel, not just its own gene.** All 300 target genes are removed
from every vector so all n² distances live in one fixed feature space. The spec explains why: with
own-gene-only exclusion, a submission could spike the other 299 target genes and win every rank
on nothing but the published target list (it scored 0.76–0.80 against a 0.51–0.53 baseline).

Degenerate points: a prediction with one shared profile, or zero effect, scores **exactly 0.5**.
The baseline is exactly 0.500 on all three contexts.

### `expr_mse_unbiased_capped_norm` — expression error

A ratio of two panel-wide sums, so **it has no per-perturbation value**. Numerator per
perturbation: squared distance between our pseudobulk profile and the measured one, minus a
sampling-noise correction; denominator: squared distance from the measured perturbation to the
measured control, minus its own correction. Both exclude the perturbation's own target gene.
`MSE = Σₚ Nₚ / Σₚ Dₚ`.

The correction is a delete-one-cell jackknife on the log-profile, not a closed form. Our credited
correction is capped by (a) the reference's own correction and (b) a factor ρ ≤ 1 that bounds the
total correction by our own across-perturbation spread — so a submission that emits *identical*
cells for every perturbation forfeits the correction entirely and reads above 1. Any submission
with real signal has ρ = 1 and is unaffected.

Fixed points: perfect = 0; emitting the control unchanged ≈ 1 (measured slightly above 1 because
the jackknife runs 0.32 % high). Scaled value clamped to [0, 1]. The four "diagnostics" in the
profile are this metric's numerator, denominator, and two audit columns.

### `de_wilcoxon_direction_fidelity_yield_raw` — direction fidelity ("fid")

**This is the answer to "how does fid scale by yield."** With n_real = genes the reference calls
significant, n_pred = genes *we* call significant that are adjudicable (reference lfc defined and
non-zero), and k = those with matching sign:

```
F_p = k / max(n_pred, n_real)  =  (k / n_pred) · min(1, n_pred / n_real)
```

So it's directional precision times a coverage term capped at 1. Calling *more* genes than the
reference buys nothing (coverage saturates); calling *fewer* costs proportionally. A gene we call
but leave without a direction (null/NaN/zero lfc) stays in n_pred and counts as a miss. Calling
nothing scores 0 whenever the reference found anything. If the reference found nothing, the
perturbation is dropped from the mean (cohort is 295–300, not 300).

**Consequence for us:** fid's baseline is chance, ≈ 0.5, and its span is the narrowest of the six
(0.28–0.33), so raw movement here converts into score fastest. And the spec states plainly that
**a genuine half of the reference scores negative on fid (−0.15 to −0.34)** because a half calls
fewer genes and its raw value (0.41–0.46) falls below the mean-response arm's 0.505. Under-calling
is punished on fid; over-calling is free on fid (but punished on `jac`, see below).

### `de_wilcoxon_direction_reach_raw` — direction reach ("reach")

**This is the answer to "what does reach rank by."** The pool is the reference's significant set
for p (own gene removed). Order that pool by **our** confidence, tie-broken in this exact order:

1. our significant calls first,
2. then ascending predicted **adjusted p-value**,
3. then ascending predicted raw p-value,
4. then descending |predicted lfc|,
5. then gene name.

Walk down counting only adjudicable genes; P(k) = fraction of the first k with correct sign.
`k* = max{k : P(k) ≥ 0.9}` (0 if none), `R_p = k* / n_real`. "Deepest" is literal: purity can dip
below 0.9 and recover, and the deeper crossing counts. A wrong call at position j leaves k* ≥ j−1,
so only a wrong call at position 1 can zero it, and only if the ranking never recovers. The
shallowest depth at which one wrong call is tolerated is 10.

**Consequence for us:** reach is driven by predicted p-values, which come out of the Wilcoxon on
our 400 cells — so it is a direct function of `phi` (our within-perturbation spread) and Δ
magnitude jointly, not of Δ's sign alone. The spec notes the fid/reach split explicitly: right
genes but poor ranking → good fid, bad reach; correct significance ranking but badly calibrated
FDR → bad fid, good reach.

### `de_wilcoxon_sig_jaccard` — significant-set agreement ("jac")

`J_p = |R_p ∩ R̂_p| / |R_p ∪ R̂_p|`, own gene removed from both sets. Empty union is defined as 1
(both agree there was no response), so it returns a value for every perturbation and the cohort is
the full 300. No chance correction; two independent sets of the observed sizes overlap by chance
at ≈ 0.006–0.012.

**Consequence for us:** symmetric — a missed gene and an invented gene cost the same. This is the
metric that punishes over-calling. Its replicate is only 0.38–0.42, meaning most set
disagreement in the real data is the reference's own instability, which the scaling removes.
Largest headroom of the six: exact reproduction of the reference scores **2.48–2.85**.

### `de_wilcoxon_lfc_nmae` — fold-change accuracy ("nmae")

Gate Sₚ = the reference's significant set with finite reference lfc, own gene removed.
`NMAE_p = Σ|lfĉ − lfc| / Σ|lfc|` over the gate. A gate gene we don't carry, or carry as
null/non-finite, is read as **predicted lfc = 0** (filled, not masked — so we can't shrink our
cohort by emitting NaNs). Perturbations with a gate < 10 genes are omitted; that reads the
reference alone, so the cohort (209–261 of 300) is identical for every submission.

Fixed points: perfect = 0; predicting zero fold change everywhere = 1 exactly (term by term).
Because our lfc is computed against the *measured* control, broadcasting the true control lands at
≈ 1.008–1.010, not 1.

**Special replicate estimator.** nmae is the one metric whose anchor is not the plain split-half:
a half calls 21–35 % fewer genes, so its gate would differ. The anchor takes the gate and
denominator from the *full* reference table and only the two fold-change vectors from the halves.
The bundle records which estimator produced each anchor and refuses a mismatch. In the CLI this
is `run --lfc-nmae-ref`, separate from `--anchor`.

Scaled: linear, floored at −6. A raw value of 2.0 already scores −1.58 to −1.75; the floor binds at
raw 4.4–4.8. Effectively bounded above too (raw error can't go below 0): exact reproduction earns
1.58–1.75.

## 3. Reference scaling and the anchor values

`s = (u − b) / (r − b)` per metric per context; context score = equal-weight mean of six;
overall = mean over contexts. Above 1 and below 0 are both reported.

- **Baseline b** is the context's mean perturbation response: an equal-weight average, over the
  constructs passing a **cell-count and knockdown-efficiency filter**, of the per-construct mean
  count vector, assigned identically to every perturbation, then scored as an ordinary submission.
  Note this is computed *from the real perturbed cells* — it is an oracle comparator a submitter
  cannot reproduce from the controls alone, not a floor we could reach by pasting our own
  control mean. (`cell-eval2 baseline --emit dispersed` is the current arm; `tile` is a legacy,
  known-biased arm kept only for reproducing old numbers.)
- **Replicate r**: cells split into two disjoint halves per perturbation *and per control*, one
  half scored as a prediction of the other, five such splits (seeds from base seed 0) averaged.
  **Each half uses its own control cells**, so the two halves are not correlated through a shared
  reference.

Measured on the official A/B/C bundles (cell-eval2 0.15.0, rule_version 3):

| metric | b | r | span \|r−b\| | cohort | exact-reproduction scaled score | clamp |
|---|---|---|---|---|---|---|
| `pds_cosine` | 0.500 | 0.927–0.984 | 0.43–0.48 | 300 | 1.03–1.17 | none |
| `expr_mse_unbiased_capped_norm` | 0.986–0.992 | 0.028–0.045 | 0.95–0.96 | panel ratio | 1.000 | [0, 1] |
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.505–0.522 | 0.795–0.832 | 0.28–0.33 | 295–300 | 1.51–1.72 | none |
| `de_wilcoxon_direction_reach_raw` | 0.047–0.097 | 0.958–0.978 | 0.86–0.93 | 290–300 | 1.02–1.05 | none |
| `de_wilcoxon_sig_jaccard` | 0.021–0.037 | 0.375–0.423 | 0.34–0.39 | 300 | 2.48–2.85 | none |
| `de_wilcoxon_lfc_nmae` | 1.0009–1.0017 | 0.369–0.431 | 0.57–0.63 | 209–261 | 1.58–1.75 | floor −6 |

Floors of the unclamped members (their raw minimum, scaled): `pds` −1.03 to −1.17; `fid` −1.55 to
−1.85 (the deepest); `reach` no worse than −0.12; `jac` −0.05 to −0.11.

**Seed noise.** Standard deviation of r over the five splits, as a fraction of the span the score
divides by: `pds` 0.5–1.3 %, `mse` 0.8–2.2 %, `fid` 1.5–4.8 %, `reach` 0.2–0.9 %, `jac` 0.8–2.2 %,
`nmae` 0.7–2.7 %. The splits are seeded and reproducible, so this is "the scale on which a
difference is too small for the metric to resolve," not a run-to-run error bar. It is a first
estimate of the noise floor for our own CV until we measure it ourselves (action item 12).

## 4. The four open questions, answered

1. **How does `fid` scale by yield?** `F = k / max(n_pred, n_real)`. Precision × coverage capped
   at 1. Over-calling is free here; under-calling is charged linearly. See §2.
2. **What does `reach` rank by?** Our own significance calls first, then ascending predicted
   adjusted p, then raw p, then descending |lfc|, then gene name. So predominantly by our
   p-values, i.e. by `phi` and Δ jointly. See §2.
3. **Published anchors?** Yes — the table in §3, per metric, as ranges over A/B/C. The exact
   per-context values are stamped in the reference bundles (`Anchor set` line in `vcc submit`
   output), not in the PDF.
4. **How are the five half-splits built?** Disjoint halves per perturbation *and per control*,
   each half scored against the other with its own controls (`control_source` forced to the half),
   seeds derived from base seed 0, five splits averaged. nmae uses the special estimator above.

## 5. Tooling that already exists in cell-eval2 0.16.0 (this shrinks action item 10)

Earlier notes assumed we'd have to hand-build the replicate anchor because `compute_ceiling` only
covers `pds_cosine`. That is still true of `--ceiling`, but 0.16.0 ships the competition anchors
directly:

| command | what it produces |
|---|---|
| `cell-eval2 baseline -ar real.h5ad --preset vcc2026 --pert-col target_gene -o baseline/` | the 0 end: mean-response arm built from the real perturbed cells (`--emit dispersed`), scored as a prediction → `baseline_agg.csv` + `baseline_meta.json` |
| `cell-eval2 run -ar real.h5ad --anchor --anchor-splits 5 --anchor-base-seed 0 ...` | the 1 end: 5-split replicate anchor → `anchor_agg.parquet`, `anchor_splits.parquet`, `anchor_meta.json`. One full metric run per split. |
| `cell-eval2 run -ar real.h5ad --lfc-nmae-ref ...` | nmae's special full-gate replicate estimator → `lfc_nmae_ref_agg.csv` |
| `cell-eval2 prep-real-bundle --real real.h5ad --baseline baseline/ -o bundle/ --preset vcc2026 ...` | **both ends plus the manifest binding them**, aggregates only, no matrices. This is the local analogue of Arc's reference bundle. |
| `cell-eval2 run -ap pred.h5ad -ar real.h5ad --preset vcc2026 --pert-col target_gene -o user/` | score our prediction → `agg_results.csv`, `results.csv`, `run_meta.json` |
| `cell-eval2 score --user-agg user/agg_results.csv --real-bundle bundle/ -o scored.csv` | **the competition score**: `from_replicate` takes `avg_score`. A submission whose config/engine/reference don't match the bundle is refused. |

`score --user-agg ... --baseline-agg baseline/baseline_agg.csv` alone gives only the
`from_baseline` column (margin over b, no r). `--anchor` / `--anchor-cache` add `from_replicate`
as a diagnostic that never takes `avg_score`; only `--real-bundle` produces the competition
number.

Practical notes from the install:

- Preset default `pert_col` is `target`, not `target_gene`. Always pass `--pert-col target_gene`
  (or `replace(EvalConfig.from_preset("vcc2026"), pert_col="target_gene")` in Python).
- `--input-type` defaults to `lognorm` on the CLI subcommands; the preset sets `counts`. Passing
  `--preset vcc2026` is what carries the right value — don't hand-assemble flags.
- On a laptop pin `--set de.backend=pdex`. `auto` warns and falls back on a no-CUDA host, but
  **raises** on a host with a visible CUDA device and no `gpudge` — it will not silently switch
  engines.
- Use persistent `--cache-real` / `--cache-pred` directories so the real side (including DE) is
  computed once per reference and reused across every prediction in a sweep.
- Bundles and partials from different `cell-eval2` releases must not be mixed; record the version
  in every run (it is already in `run_meta.json`).

## 6. VCC CLI: what `vcc prep` catches, verified locally (action item 4)

`vcc-cli 0.2.0` installed via `uv tool install vcc-cli` (binary at `~/.local/bin/vcc`; run
`uv tool update-shell` once if it's not on `PATH`). The smoke test below was run before we had the
real controls bundle, against **placeholder** gene and perturbation lists of the right shape
(18,533 symbols, 300 targets) in `/tmp`. Nothing from it is committed. Now that the real
`gene_names.csv` and `pert_counts.csv` are in `data/contracts/`, the same two commands should be
re-run against them.

- `vcc sample -g genes.csv -p perts.csv --h5ad -o sample.h5ad --seed 0` → 360,000 cells ×
  18,533 genes, CSR float32, `obs` = `target_gene`, `context` (A/B/C × 120,000), ~876 MB on
  disk, 17 s.
- `vcc prep sample.h5ad -g genes.csv --perts perts.csv --dry-run` → passes in 1.6 s and reports
  cells, genes, encoding, per-context counts, `targets: verified against the official list`,
  `normalization: counts-preserved`.

Deliberate rejections, with the exact leading line of each message:

| break | message |
|---|---|
| add 0.5 to some entries | `Submission values are FRACTIONAL, but submissions must be RAW INTEGER COUNTS.` |
| relabel 400 cells `non-targeting` | `Submission contains 400 'non-targeting' control cell(s), but controls must NOT be submitted.` |
| drop one cell | `Context A: 1 perturbation(s) have the wrong number of cells.` then `GENE00000 has 399 (expected 400)` |
| move one block from A to B | `Context A: perturbation labels do not match the official list.` — caught **only because** the move made A short one perturbation |

**What `prep` cannot catch:** a clean swap that relabels *all* of context A as B and all of B as
A. Every per-context check still passes and every metric quietly degrades toward the baseline.
This is the check our own `validate.py` (action item 11) must add — e.g. confirm each context's
predicted pseudobulk correlates best with *that* context's own downloaded controls.

Flag gotchas: in `prep`, `-p` is `--pert-col` and the perturbation list is `--perts` (no short
form); in `sample`, `-p` *is* the perturbation list. `--max-nnz` defaults to the scoring
hardware's limit (4.75 B stored entries); dense storage is 1.40× over on its own.

## 7. Corrections to our other docs

1. **`pds_cosine` excludes all 300 panel target genes, not just the perturbation's own.**
   `project-rules.mdc` ("The perturbed gene's own row is excluded from all six metrics") and
   `docs/2026-metric-explanation.md` both state own-gene-only exclusion. The other five metrics
   are own-gene-only; `pds` is panel-wide. The modeling consequence is small (the 300 target
   rows are excluded from the distance either way) but any local re-implementation of `pds` has
   to match it.
2. **The baseline `b` is an oracle, not a reachable floor.** It is the mean of the *real perturbed*
   cells (filtered by cell count and knockdown efficiency), not the mean of the downloaded
   controls. "The context mean scores exactly 0.0" in the rules is correct as a statement about
   the scale's zero point, but a submitter cannot construct that arm from the controls bundle.
   Pasting the control unchanged scores ≈ 0 on `pds` (both read 0.5), ≈ 0 on `nmae`, and
   *slightly below* 0 on `mse` (raw ≈ 1 vs b ≈ 0.99).
3. **The replicate is half-depth data, which is why >1 is common.** Exact reproduction of the
   reference scores 1.03–2.85 depending on the metric. A leaderboard `jac` of 1.5 is not
   suspicious.
4. **`fid`'s baseline is chance (0.5), not "no calls."** A real half-replicate scores *negative*
   on fid. The recipe's assumption that fid rewards sharp confident calls is right, but the
   mechanism is coverage: we must call at least n_real genes per perturbation to stop paying the
   coverage penalty, while `jac` charges every extra call. `phi` sits exactly on that tension.

## 8. Local harness validation on 2025 H1 data (action items 5–7, done 2026-09-04)

### 8.1 The scoring contract, confirmed empirically

Running `cell-eval2 run --preset vcc2026` on a prediction **without** `non-targeting` rows fails
hard: `ValueError: perturbation sets differ: 5 pred vs 6 real; real-only=['non-targeting']`.
So locally we must do what the platform does — concatenate the real control cells onto our
prediction before scoring. `score.py` (action item 9) owns that bridge. The upload file still
carries no controls.

The `vcc2026` run emits ten columns: the six scored metrics plus `expr_distance_unbiased`,
`expr_mse_unbiased`, `expr_mse_unbiased_capped`, `expr_real_mass_ratio`. `results.csv` is long
form (`perturbation, metric, value`); `agg_results.csv` is wide with `count/mean/std/min/max/median`
rows; `expr_mse_unbiased_capped_norm` has `count = NaN` because it is a panel ratio.
`run_meta.json` stamps `cell_eval2_version`, `config_digest`, `resolved_de_backend`, and an
`anchor_semantic_identity` that `score --real-bundle` uses to refuse mismatched runs.

If a target label does not resolve to a gene in `var_names`, the scorer warns and scores that
perturbation **without** target exclusion, which can inflate `pds`. On the full gene space every
2026 panel gene resolves (checked below); `EvalConfig.target_gene_map` exists for the case where a
symbol differs.

### 8.2 Where the 2025 data lives and how it compares to 2026

The Hugging Face repo `arcinstitute/VCC_train` referenced in older notes **no longer exists**.
The processed 2025 files are publicly readable over plain HTTPS from the Arc bucket, no GCP
project or requester-pays needed:

```
https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/
  gene_names.csv                          18,080 genes, headerless
  train/adata_Training.h5ad               15.5 GB, 150 perturbations
  validation/adata_Validation.h5ad         6.9 GB,  50 perturbations, 98,927 cells, 38,176 NTC
  test/adata_Test.h5ad                    12.0 GB, 100 perturbations
  {train,validation,test}/pert_counts_*.csv   with n_cells and median_umi_per_cell columns
```

Local copies: `data/vcc2025/` (ignored by git). `.X` is CSR float32 holding integer counts;
`obs` has `target_gene`, `guide_id`, `batch`; `var_names` equals the 2025 `gene_names.csv` order.
Per-perturbation depth is 161–2,925 cells at a median ~44,000 UMI/cell — about twice the 2026
reference depth of 20,000.

Gene space vs 2026 `data/contracts/gene_names.csv`:

- 18,077 of 18,080 2025 genes are in the 2026 space. The 3 that are not — `HSPA14-1`, `TBCE-1`,
  `TMSB15B-1` — are `var_names_make_unique` suffixes, i.e. 2025 had duplicate symbols that 2026
  collapsed.
- 456 genes are 2026-only, largely `AC…` clone-ID lncRNAs from a newer annotation.
- The order differs, so nothing may be aligned positionally.
- All 300 of the 2026 panel genes are measured in both spaces. **25 of them were also 2025
  targets** (`ACLY`, `ADNP`, `AKT2`, `ANKZF1`, `BRPF1`, `HDAC8`, `HSBP1`, `MED13`, `MED15`,
  `MED25`, `MTA1`, `NFE2L1`, `PLAGL2`, `RNF2`, `SHPRH`, `SIN3B`, `SLIRP`, `SMARCA5`, `STAT6`,
  `TARBP2`, …) — real H1 ground truth for 25 panel genes.

For **scoring against 2025 ground truth no harmonization is needed** — pred and real only have to
share an axis. Harmonization (reindex onto the 2026 order, zero-fill the 456, carry a measured
mask) is only required when 2025 data becomes *training input*. That belongs to the data team.

### 8.3 Landmarks

Local reference: the 2025 validation set downsampled to the challenge shape — 44 perturbations ×
400 cells (six perturbations with < 400 cells dropped: `EIF3H`, `FUBP1`, `MAX`, `SUPV3L1`, `WAC`,
`ZNF598`) plus 8,000 NTC cells, seed 0. Pipeline, all with
`--preset vcc2026 --pert-col target_gene --set de.backend=pdex --cache-real cache/real`:

```
cell-eval2 baseline -ar ref.h5ad ... -o baseline --save-pred pred_baseline.h5ad
cell-eval2 prep-real-bundle --real ref.h5ad --baseline pred_baseline.h5ad -o bundle ... --id <name>
cell-eval2 run -ap pred.h5ad -ar ref.h5ad ... -o user
cell-eval2 score --user-agg user/agg_results.csv --real-bundle bundle -o scored.csv
```

Two flag facts learned the hard way: `prep-real-bundle --baseline` takes the baseline
**prediction h5ad**, not the `baseline/` output directory; and the saved baseline prediction is
fractional (a mean profile) and only passes the counts check because `baseline` sets
`allow_fractional_counts` internally — it cannot be re-run as an ordinary submission, and it
cannot be self-scored against the bundle (`baseline_meta.json` has no anchor identity). The 0 end
is the bundle's zero by construction.

Scaled results (`from_replicate` column — the competition number):

| metric | exact reproduction (real vs real) | spec's exact-reproduction range | control paste |
|---|---|---|---|
| `pds_cosine` | 1.106 | 1.03–1.17 | −0.071 |
| `expr_mse_unbiased_capped_norm` | **1.000000** | 1.000 | 0.000 |
| `de_wilcoxon_direction_fidelity_yield_raw` | 1.163 | 1.51–1.72 | −0.318 |
| `de_wilcoxon_direction_reach_raw` | 1.019 | 1.02–1.05 | −0.104 |
| `de_wilcoxon_sig_jaccard` | 2.281 | 2.48–2.85 | −0.013 |
| `de_wilcoxon_lfc_nmae` | 1.587 | 1.58–1.75 | −0.084 |
| `avg_score` | 1.359 | — | −0.098 |

Both landmarks behave as the spec says: `mse` is exactly 1.0 at exact reproduction and 0
(clamped) for control paste; the five unclamped members all exceed 1 at exact reproduction in the
spec's pattern (`jac` highest, `reach` barely above 1); control paste lands just below 0 on
`nmae` (raw ≈ 1.008, exactly the spec's "broadcasting the true control" value) and near 0 on
`pds`. The harness is wired correctly.

Local anchors vs the published A/B/C ones, for calibration only:

| metric | local b | spec b | local r | spec r |
|---|---|---|---|---|
| `pds_cosine` | 0.515 | 0.500 | 0.954 | 0.927–0.984 |
| `expr_mse_unbiased_capped_norm` | 0.965 | 0.986–0.992 | 0.004 | 0.028–0.045 |
| `fid` | 0.219 | 0.505–0.522 | 0.890 | 0.795–0.832 |
| `reach` | 0.152 | 0.047–0.097 | 0.984 | 0.958–0.978 |
| `jac` | 0.006 | 0.021–0.037 | 0.442 | 0.375–0.423 |
| `nmae` | 0.949 | 1.0009–1.0017 | 0.351 | 0.369–0.431 |

The `r` values sit slightly *above* Arc's on every DE metric, as expected for data at twice the
depth with 8,000 controls. The `b` values differ more — H1's mean response calls fewer genes than
A/B/C's (`fid` b = 0.22 vs 0.51) and actually helps on `nmae` (0.95 vs 1.00). This is exactly the
caveat in action item 14: **local scaled scores rank configurations; they do not predict the
leaderboard.** A panel of 44 also makes `pds` noisier (control paste reads −0.07 instead of the
spec's ≈ 0).

### 8.4 Timing (action item 12, first numbers)

Apple M4 Max, 14 cores, 38 GB, CPU-only `pdex`, 44 perturbations × 400 cells + 8,000 controls,
18,080 genes:

| stage | wall | notes |
|---|---|---|
| `baseline` (build + score the 0 end) | 64 s | 6m50s CPU — `pdex` parallelizes across cores |
| `prep-real-bundle` (5 splits, each a full DE run on halves) | 5 min 8 s | 38 CPU-minutes; one-time per reference |
| `run` on one prediction, real side cached | **38.5 s** | this is the per-config cost in a sweep |

Extrapolating linearly in perturbation count: one 300-perturbation context ≈ 4.4 min per
prediction, so **≈ 13 min for a full three-context prediction on this laptop**, and ≈ 35 min per
context to build a bundle once. A 20-configuration sweep over three contexts is therefore an
afternoon on CPU; a `phi` × `weights` grid in the hundreds wants the GPU box. Real-side cache
for this reference is 56 MB; the bundle is 28 KB (aggregates only, no matrices).

## 9. Still open

- Run the same landmark check on the 150-perturbation training set (closer to the 300 panel) and
  measure the resampling noise floor (action item 12, second number).
- Read the per-context anchor values out of a real `vcc submit` (`Anchor set` stamp) once we have
  a scored submission, and compare with the local bundle to see how far off our local scale sits.
- Fold §8.3's pipeline into `src/vcc/` (action items 8–10): `score()` must perform the control
  bridge of §8.1 and refuse to score without it.
