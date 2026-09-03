# The plan

Two stages, run one after the other:

1. **Predict the effect.** For each of the 300 knocked-down genes, in each of the 3 cell lines, predict how much every one of the 18,533 genes goes up or down. Call that the effect vector, or Δ.
2. **Turn Δ into cells.** The submission has to be 400 individual cells per perturbation, not one average. So take real control cells from that cell line and shift them by Δ.

Only two stages — Part 3 below isn't a third one. It's how we pick the few numbers these two stages leave open: the weights that blend stage 1's competing predictions of Δ, and the spread knob in stage 2. Both get chosen by score, which is why they come last even though they live inside the two stages.

## Part 1: predicting the effect

### What the function looks like

`delta(target_gene, cell_line) -> array of 18,533 numbers`

Internally we'll store it as "difference in average log1p expression, perturbed minus control."

Careful: the scorer doesn't use that space. It uses two others.

- `pds` and `mse` work on pseudobulk normalized to `bulk_target_sum` = 5e4.
- The four DE metrics work on log2 fold changes from counts normalized to 1e6 per cell.

So write an explicit converter into each one, and test both against 2025 data where we already know the right answer. If a conversion is wrong, nothing crashes — the score just quietly drops.

### Step 1: collect what's already known

For every public dataset and cell line, get a log fold change per perturbation plus its standard error (how noisy that number is). Nadig's figshare files already ship log2FC, SE, and p-values for all five Replogle/Nadig datasets, so we don't have to compute pseudobulk ourselves.

Map everything onto our 18,533 genes. Keep a true/false mask of which genes were actually measured in which dataset — "we never measured this gene" and "we measured it and it didn't change" both look like a zero, and they must not get confused.

### Step 2: shrink noisy numbers toward zero

About fifteen lines of code, and the biggest statistical win in the project.

The idea: if an effect was measured in 7 cells it's mostly noise, so pull it most of the way to zero. If it was measured in 500 cells, trust it and leave it alone.

```python
tau2 = np.maximum(0, L.var(axis=0) - (S**2).mean(axis=0))   # per gene
L_shrunk = L * (tau2 / (tau2 + S**2))                        # per (pert, gene)
```

`tau2` is how much a gene genuinely varies across perturbations once measurement noise is subtracted off. The fraction `tau2 / (tau2 + S**2)` is near 1 for a reliable measurement and near 0 for a noisy one.

This is the same logic as `limma`'s moderated statistics, and it's the direct fix for Nadig's finding that reproducibility runs anywhere from 0.16 to 0.90 depending on how many cells were behind the number.

Make the shrinkage strength a knob rather than just using the formula, because the metrics want different things:

- Shrinking harder helps `mse` and `nmae`, which reward getting magnitudes right.
- Shrinking harder hurts `pds`, which ranks our 300 predictions against each other — if everything is pulled toward zero, they all start to look alike.
- And the context mean scores exactly 0.0, so shrinking all the way to the panel average earns nothing at all.

### Step 3: one global number per perturbation

For each target gene, average the shrunk estimates across every cell line where that gene was targeted, weighting by `1/S²` so more precise measurements count more.

If a gene was only ever tested in K562-GenomeWide, this just returns the K562 answer. That means "copy the K562 result over" falls out as the automatic fallback instead of being a special case we bolt on.

### Step 4: find the shared response programs

Stack the global effects into a (perturbations × 18,533) matrix and factor it into 50–100 components. That gives a program dictionary `H` and per-perturbation loadings `W`. Rebuilding a prediction through this basis smooths out noise on its own.

NMF run separately on the up- and down-regulated parts gives more interpretable programs. SVD is faster and fine to start with.

### Step 5: adjust for cell line, in stages

Don't jump to the fancy version. Three tiers, and only move up a tier if leave-one-cell-line-out CV actually improves.

**v1, one scaling number.** `Δ(g,C) = α(g,C) · Δ_global(g)`, where α is a single number per perturbation-and-cell-line pair. It only answers "does this knockdown do more or less in this cell line than average." One parameter is very hard to overfit, and it probably captures most of the real cell-line difference. The most useful feature for predicting α is how strongly the target gene is expressed in this cell line, both in absolute terms and relative to the source lines.

**v2, one scaling number per program.** `Δ(g,C) = Σ_j W_gj · m_j(g,C) · H_j`, so different response programs can scale differently by cell line. Features: how active each program is in C's baseline, how similar C's baseline is to each source line, and target gene expression.

**v3, an unconstrained residual.** Only if there's spare time, which there won't be.

Fit α and m with ridge regression, but shrink toward 1 rather than toward 0 — the default assumption is "same effect everywhere," not "no effect."

### Step 6: validate by holding out whole cell lines

Never hold out random perturbations. That scores well and tells us nothing, because the real test is always a cell line we've never seen. Roughly ten source contexts gives us ten folds. Thin, but honest.

### What to do when data is missing

Three cases, handled differently:

- Measured in several source lines: full weighted and cell-line-corrected estimate.
- Measured only in K562-GenomeWide: copy it over and shrink it, no cell-line correction possible.
- Measured nowhere: borrow from similar genes — ones whose K562 effects correlate with it, or STRING/CORUM partners.

Plus one override that applies to all three: **if the target gene isn't expressed in that cell line, predict roughly no effect.** Cheap, safe, and usually right.

## Part 2: turning Δ into 400 cells

Never build a cell from scratch. We're handed 18,400 real unperturbed cells from each target cell line, and resampling those gives us, free of charge, the right library sizes, the right dropout pattern, the right mix of cell-cycle states, and the right gene-gene correlations. No generative model comes close for the effort.

```python
N_CELLS = 400                                   # flat for every perturbation in 2026
idx     = rng.choice(basal_cells_of_C, size=N_CELLS, replace=True)
X       = basal_log1p[idx]                      # (400, 18533)
eps     = rng.gamma(shape=1/phi**2, scale=phi**2, size=(N_CELLS, 1))  # mean 1, CV = phi
X       = X + eps * delta[None, :]              # per-cell response heterogeneity
X[:, target_gene_idx] = 0.0                     # cosmetic: the target row is dropped from all six metrics
```

`eps` is doing real work. Without it, every cell shifts by exactly the same amount, so our fake perturbed population has exactly the control's spread — biologically wrong and statistically fragile. With it, cells respond by different amounts, which is what actually happens with CRISPRi and is the whole reason Jiang's Mixscale exists. `phi` is the spread knob.

### Don't submit fake control cells

Everyone did this in 2025, PRiMeFlow included (200 of them). The 2026 rules reject any `non-targeting` rows outright. The controls we download are model input only; scoring compares our 400 cells against Arc's own held-out real controls.

That makes `phi` both harder to set and more important. In 2025 we controlled the spread on both sides of the statistical test, so only the ratio mattered. Now one side is real data with real spread, and ours has to match it in absolute terms — otherwise every p-value comes out skewed in the same direction.

### Post-processing

Borrowed from PRiMeFlow, but it has to stop earlier: their version targeted 2025 rules that wanted normalized data, and this one demands raw counts.

1. `expm1` to undo the log.
2. Scale each cell so its total is in line with the controls' ~20,000 median UMI.
3. Round to whole numbers.

Stop there. Don't renormalize after rounding and don't go back to log1p — both reintroduce fractions, and `.X` has to be non-negative, whole, and finite. Rounding is what makes the output look like real count data, and it should be the last thing that touches the matrix.

Then check the limits before packaging: no cell over 1,000,000 total counts, and store CSR sparse, since a dense 360,000 × 18,533 array is on its own 1.4× over the stored-entry cap.

## Part 3: picking the handful of numbers the first two stages left open

Nothing in this part touches the packaged submission. There is exactly one `.h5ad`, built once, from one Δ per (target gene, cell line) pair. We never build several submissions and combine them.

What gets blended is Δ, back in stage 1. Steps 1–5 gave us four different ways to predict the same 18,533-number vector, so we average them into one vector and hand that single vector to the sampler:

```python
candidates = [
    transplant(g, C),        # copy from the most basally-similar source line
    global_avg(g),           # the precision-weighted average from Step 3
    lowrank_gain(g, C),      # low-rank rebuild with the v1 scalar gain
    program_gain(g, C),      # the v2 per-program version, if it survived CV
]
delta = sum(w * c for w, c in zip(weights, candidates))   # still 18,533 numbers
```

That's the whole idea. Four guesses at one quantity, averaged, and `weights` is four numbers that add up to 1. By the time Part 2 runs, the blend has already collapsed into a single Δ and there's nothing left to combine.

Stop at three or four candidates. We only have about ten CV folds, and a blend with more weights than that is fitting the mix to noise.

### Why this is written last if it belongs to stage 1

Because we can't choose `weights` without a score, and getting a score means running all the way through to cells. So the search needs the sampler to already exist. Same for `phi`, the spread knob from Part 2. Those are the only two things left open, and they get chosen together.

### How the search actually runs

For each candidate combination of `weights` and `phi`:

1. Hold out one source cell line — one we have real perturbation data for, so we know the right answer.
2. Run stage 1 with those weights, fit on the remaining lines, to get Δ for the held-out line.
3. Run stage 2 with that `phi` to turn each Δ into 400 cells.
4. Score those cells against the real held-out data with `cell-eval2`.
5. Repeat for every fold and average.

Then keep the combination that scored best. `weights` is three or four numbers and `phi` is one, so this is a coarse grid over three or four dimensions — entirely tractable, and probably worth more than a week of modeling improvements.

Note that this scoring loop runs on our *source* data, not on A/B/C. We have no ground truth for the real validation contexts, so the leaderboard is the only way to score those, and that's capped at two submissions a day.

Score plain equal weights alongside the fitted ones every time, and only take the fitted weights if they win by a clear margin on the held-out folds. With CV this thin, equal weighting wins more often than people expect.

### What the score rewards

The target is the actual 2026 composite: the flat unweighted mean of the scaled scores over 6 metrics × 3 contexts. Every metric is scaled so 0 means "no better than the context mean" and 1 means "as good as repeating the experiment," which means a metric with a low baseline turns raw improvement into points fastest:

- `reach` (baseline ≈ 0.05–0.10) and `jac` (≈ 0.02–0.04) have the most headroom. Focus here.
- `mse` is clamped to `[0, 1]`, so there's a ceiling on what it can ever contribute.
- `nmae` starts at a baseline of ~1.0 with a floor of −6. It's the one metric that can actively drag the total down if our fold changes are badly scaled.

### Why `phi` can't be tuned separately

The metrics want opposite things:

- `fid` and `reach` want sharp, confident, well-ordered DE calls.
- `jac` divides by the union of the two significant-gene sets, so over-calling is punished just as hard as missing genes.
- `nmae` and `mse` want honestly calibrated magnitudes on top of both.

`phi` sits in the middle of all of it, and it has more leverage than in 2025 because four of the six metrics now run through the Wilcoxon test. `phi` sets the within-condition spread, which sets the p-values, which sets how many genes clear BH at α = 0.05. That gene count is simultaneously `jac`'s denominator, `fid`'s yield term, and how deep `reach` gets to search.

- `phi` too small: p-values collapse, we call thousands of genes, directional accuracy drifts toward chance, and the union wrecks `jac`.
- `phi` too large: nothing clears the threshold, so `fid`'s yield and `reach` both go to zero.

Since the best `phi` depends on how big and how noisy Δ is, and the weights determine that, picking one and then the other gives a worse answer than sweeping both at once.

Run the sweep with a perfect (oracle) Δ first, to see what the sampler alone contributes, then with the real Δ.

## Schedule

About nine weeks. The final test set drops October 22 and submissions close November 5. So it's roughly seven weeks building against the A/B/C validation panel, then a two-week final phase on D/E/F.

**September.** Reference table, shrinkage, global effect, the copy-over baseline, and the sampler with `phi` at a rough guess. Get something onto the validation board early — two scored submissions per team per day makes the leaderboard a working signal rather than a scarce resource, so there's no reason to sit on a rough pipeline.

**Early October.** Low-rank basis, scalar gain, leave-one-cell-line-out harness.

**Mid October.** Blend weights and the joint `phi` sweep, finished before the 22nd.

**Last two weeks.** Re-run the fitted pipeline on the new cell lines, plus buffer. No new modeling: the final scores aren't comparable to anything we saw in validation, so tuning then means tuning blind.

## Two questions that are now settled

- Δ has to reach two scoring spaces, not one: pseudobulk at `bulk_target_sum` = 5e4 for `pds` and `mse`, and log2 fold change on CPM for the four DE metrics. Build a converter for each and verify both against 2025 ground truth instead of picking one canonical space.
- The test set gives three cell lines per round, not one, so blend weights can be fit per similarity regime. Just don't read across rounds: validation is A/B/C, the final is D/E/F, the panels differ, and `pds` is a rank within its own panel.

## Two things still to check in the `cell-eval2` source

- Exactly how `fid` scales by yield.
- What `reach` uses to rank our genes: predicted p-value, or predicted effect size.

Both change how `phi` translates into score, so the joint sweep is guesswork until they're pinned down.

---

Want me to start on the reference table and the shrinkage step? That's the piece everything else depends on, and it's self-contained enough to be useful even if the rest of the plan shifts.
