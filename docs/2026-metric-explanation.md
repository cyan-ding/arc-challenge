Good — the rules file has the exact DE parameters, so I can be concrete. Let me build this up from your `.h5ad`.

## Step 1: what you actually hand over

One table, 360,000 rows by 18,533 columns. Each row is one cell, each column one gene, each entry a raw integer count. Every row carries two labels: `target_gene` (which knockdown) and `context` (which cell line).

So really you're handing over **900 little blocks of 400 cells** — 300 knockdowns × 3 contexts. Every metric is computed per block. That's the whole submission.

You do *not* submit control cells. The scorer uses its own held-out real NTC cells. (This is a change from 2025 worth noting, and I'll come back to it.)

## Step 2: the scorer turns each block into exactly two objects

Take one block: your 400 predicted cells for knockdown `ADNP` in context `A`. The scorer computes two things from it, and two matching things from the 400 *real* `ADNP` cells.

**Object 1 — a pseudobulk vector.** Average your 400 cells into a single 18,533-long vector (after normalizing each cell to a total of 5e4). Do the same for the real control cells. Subtract: that's your predicted *effect*, one number per gene. One vector per block.

**Object 2 — a DE table.** Normalize each cell to 1e6 counts (CPM). Then for each gene, run a two-sided Wilcoxon rank-sum test: your 400 perturbed cells against the real held-out control cells. Genes below 5 CPM in the control cells are dropped before testing. Every surviving gene gets a p-value and a log2 fold change. Apply BH correction at α = 0.05 within this block. Now you have, per gene: a log2FC, and a yes/no "significant" flag. The perturbed gene's own row is deleted.

The real data goes through the identical pipeline, producing a reference pseudobulk effect and a reference DE table.

**This is the answer to your Wilcoxon question.** The scorer runs the test on your fake cells exactly as it runs it on real cells. It doesn't ask you for p-values — it derives them from the spread of the 400 cells you invented. So your cell-to-cell variance is a scored quantity, whether you thought about it or not.

## Step 3: the six metrics are six comparisons of those objects

Two use the pseudobulk vectors, four use the DE tables. Let me use a toy example with six genes for the DE ones.

Reference DE says: significant genes are G1, G2, G3, G4, with log2FC of +2.0, +1.0, −1.5, −0.5. G5 and G6 aren't significant, though G5's real log2FC is −0.3.

Your DE says: significant genes are G1, G2, G3, G5, with log2FC of +1.5, +0.5, +1.0, −0.8. Your G4 log2FC is −0.1 but not significant.

**`jac` — do the two significant lists match?** Intersection is {G1, G2, G3} = 3 genes. Union is {G1, G2, G3, G4, G5} = 5. Score 3/5 = 0.6. Missing G4 and falsely calling G5 hurt equally.

**`fid` — of the genes *you* called, do they move the right way?** Your list is {G1, G2, G3, G5}. In the real data G1 is up (you said up ✓), G2 is up ✓, G3 is down but you said up ✗, G5 is down and you said down ✓. Three of four = 0.75. Note G5 counted even though the reference didn't call it significant — only the sign matters. Then it's multiplied by a yield term so that calling one gene right doesn't beat calling forty mostly right.

**`nmae` — are your effect sizes right, on the genes the reference cares about?** Genes are the reference's four. Absolute errors: |2.0−1.5| = 0.5, |1.0−0.5| = 0.5, |−1.5−1.0| = 2.5, |−0.5−(−0.1)| = 0.4. Mean = 0.975. Divide by the mean size of the real effects, (2.0+1.0+1.5+0.5)/4 = 1.25, giving 0.78. Notice that predicting zero fold change everywhere gives exactly 1.0 — which is why the rules say the baseline is ~1.0 and a zero prediction scores 0. Significance doesn't matter here; G4 was included even though you didn't call it.

**`reach` — how deep does your own ranking stay clean?** Take the reference's four significant genes and sort them by *your* confidence. Say your ordering is G1, G3, G2, G4. Walk down: at depth 1 you have G1, direction correct, purity 1.0, above the 0.9 floor ✓. At depth 2 you've added G3, which you got backwards, so purity is 0.5, below the floor ✗. Deepest valid cut is 1 gene out of the reference's 4 → 0.25. One confidently-wrong gene near the top of your ranking wrecks this, even though `fid` gave you 0.75 on nearly the same information.

**`pds` — is your effect vector recognizable as this specific knockdown?** No Wilcoxon. Take your pseudobulk effect for `ADNP`, measure cosine distance to all 300 *real* effect vectors in context A, sort them, and see where real `ADNP` landed. Rank 1 is perfect. A generic prediction lands mid-pack because it's equally close to everything.

**`mse` — is your absolute expression profile calibrated?** Also pseudobulk, no Wilcoxon. Squared difference between your profile and the real one, minus the noise you'd get purely from 400 cells being a finite sample, divided by the size of the real effect. Reported once across the panel, not per knockdown.

## Step 4: every raw number gets rescaled

The raw numbers above aren't your score. Each gets mapped through \(s = (u - b)/(r - b)\), per context, where \(b\) is what a context-mean predictor scores and \(r\) is what a real replicate scores (the real cells split in half and scored against themselves, averaged over five disjoint splits).

So if your raw `jac` is 0.6, the context-mean baseline is 0.03, and a real replicate scores, say, 0.5, you'd get (0.6 − 0.03)/(0.5 − 0.03) ≈ 1.2. **Zero means you did no better than predicting the average perturbation response; one means you did as well as re-running the experiment.** Above 1 is allowed on five of the six.

Your final number is the flat average of 18 values — 6 metrics × 3 contexts.

## What this means mechanically

Four of six metrics come out of that Wilcoxon step, and its only inputs are your 400 cells' means *and* their spread. Make the cells too similar to each other and p-values collapse toward zero, so you call thousands of genes significant, `fid`'s direction rate falls to near-random and `jac` gets crushed by the union. Make them too noisy and nothing clears BH, so `fid`'s yield term and `reach` both go to zero. There's a single variance knob sitting between those two failures, and it's worth more than most modeling improvements.

This exposed two stale assumptions in `cursor-recipe.md`, both since corrected there. It had you emitting synthetic control cells so the grader would compare two clean resamples — but 2026 rejects `non-targeting` rows and scores against held-out real controls, so your 400 cells face real data with real variance and matching it is your problem alone. And it treated DES and AUPRC as the sharpness-rewarding metrics, which no longer exist.

One caveat on the above: the pipeline shape and the DE parameters are solid, but two exact functional forms are not — how `fid`'s yield factor is defined, and what quantity `reach` ranks your genes by (predicted p-value or predicted effect size). Both live in `cell-eval2` and are worth reading directly before tuning against them.