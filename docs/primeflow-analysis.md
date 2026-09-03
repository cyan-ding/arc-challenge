# Important details from paper

*Exported on 8/30/2026 at 22:25:26 PDT from Cursor (3.17.19)*

---

**User**

pls read through [https://arxiv.org/pdf/2604.13986](https://arxiv.org/pdf/2604.13986)  and extract the most important details that I have to know.

---

**Cursor**

I'll read the paper — it looks like it's already been saved locally, so I'll start there.

The paper is long (~107k chars), so I'll read it in sections.

I've read the whole paper. Here's what matters.

## What it is

**PRiMeFlow** (Altos Labs, arXiv 2604.13986v2, May 2026) is a conditional flow-matching model that generates single-cell perturbation responses **directly in log1p gene expression space** — no PCA, no foundation-model embeddings, no VAE latent. The velocity field is a **U-Net** (adapted from OpenAI's guided-diffusion code), conditioned on a compositional one-hot perturbation encoding, a one-hot covariate (cell line), and sinusoidally-encoded time. Source distribution is a plain standard Gaussian with **independent coupling** — deliberately not optimal transport. Trained with the standard CFM loss, 20% condition dropout so classifier-free guidance works at inference.

## The three ablation results that are the actual contribution

These are the parts worth trusting, because they're clean A/B comparisons within one codebase:

1. **Full gene space beats PCA space, decisively.** MMD in gene expression drops from ~1.2 (FM-PCA) to 0.13 (PRiMeFlow) on Srivatsan20, and top-50 DEG recall goes from 0.042 to 0.26. Compressing to PCA is where the perturbation-specific signal dies.
2. **U-Net beats MLP by an enormous margin** — the MLP variant is catastrophically bad (MMD 14 vs 0.13, cosine logFC ≈ 0 vs 0.46). The authors openly admit they don't know why, since gene ordering is arbitrary and a U-Net's spatial inductive bias shouldn't apply. This is the paper's weakest point intellectually but its strongest empirically.
3. **Optimal-transport coupling adds nothing.** FM-PCA-OT (essentially CellFlow) is no better than FM-PCA with a Gaussian source. That's a direct shot at CellFlow's central design choice.

## The metric trade-off you should internalize

PRiMeFlow is the **only** method in their benchmark that achieves non-trivial recall of ground-truth DEGs — CPA, SAMS-VAE, Biolord, and Latent Additive all score between 0 and 0.03, i.e. essentially zero. But on classic pseudobulk metrics (cosine logFC, RMSE of the mean), it **loses** to CPA, SAMS-VAE, and even the trivial Latent Additive baseline. Distribution-fitting and mean-prediction are close to orthogonal skills. Raising CFG weight from 1 to 5 trades distributional fit for DEG recall, so the guidance weight is effectively a dial between the two metric families.

## The VCC section — most relevant to you

They competed on the **2025** H1 hESC challenge (150 train / 50 public test / 100 private test perturbations), not the cross-cell-line 2026 task your `public-data.md` is scoped to. Their private-test numbers:


| model            | DES ↑     | PDS ↑     | MAE ↓     | AUPRC ↑   |
| ---------------- | --------- | --------- | --------- | --------- |
| PRiMeFlow-250pts | 0.228     | 0.747     | **0.086** | 0.266     |
| xBio             | 0.305     | **0.811** | 0.770     | 0.252     |
| Mean Predictors  | 0.305     | 0.741     | 6.723     | **0.582** |
| Shippers         | **0.354** | 0.699     | 0.231     | 0.576     |


Read that honestly: they lose on DES to all three top teams and to a mean predictor, and they lose badly on AUPRC. Their claim to "lowest average rank on the Generalist leaderboard" rests heavily on MAE, where competitors score 0.23 to 6.7 — those teams clearly aren't producing calibrated expression profiles at all. So the headline is less "we won" and more "we're the only entry that's simultaneously decent everywhere and numerically well-scaled."

**Their pretraining atlas is a direct check on your data plan:** 6,975,767 cells, 14 cell lines, 10,944 knockdowns, 89% public. The composition is Replogle 2022, Nadig 2025, Replogle 2020, McFaline-Figueroa 2024, Jiang 2024, the VCC training data, plus proprietary differentiated lines. Two things stand out. First, **they use no Tahoe-100M and no scBaseCount** — an independent vote for your Tier-1/Tier-4 ranking and against chemical-perturbation pretraining. Second, they include a dataset that isn't in your list at all: **Feng et al. 2026, *Cell Genomics*** — a genome-scale CRISPRi map of 7,226 gene knockdowns across **34 iPSC lines from 24 donors**. For a pluripotent-context challenge that's arguably the single most on-target public dataset in existence, and it's a real gap in `public-data.md`.

## Practical recipe details worth stealing

- **Finetuning on test-perturbation identity.** After pretraining, they finetune on a subset containing the 150 H1 training perturbations *plus the 100 private test perturbations as measured in the external datasets* (no H1 ground truth used — legal under covariate transfer). Narrowing that scope helped: the 250-perturbation variant beat the 300-perturbation variant on every metric. Their framing is that pretrained virtual-cell models need context-specific adaptation to unlock predictive power, achievable with zero new ground truth.
- **Asymmetric CFG at inference:** weight **20** for perturbed cells, **5** for synthetic controls. They report this specifically inflates PDS by exaggerating the uniqueness of each perturbation delta. That's metric-gaming, but it's disclosed and it works.
- **Normalization sandwich:** `expm1` → median-normalize to the VCC-recommended `median_umi_per_cell` → **round to integer** → median-normalize again → back to log1p. Then manually zero out the target gene. This is almost certainly where the MAE advantage comes from.
- Dopri5 ODE solver, exact cell counts from `pert_counts_*.csv`, fixed 200 synthetic controls.
- Validation split: 43 of the 150 training perturbations held out (gene list given in full in footnote 6), cut to just 5 during finetuning.

## Caveats

No code release and no acknowledgment section — the hyperparameter tables cover only the three PerturBench datasets, **not the VCC model**, so the headline result isn't reproducible from the paper. Batch sizes of 38–54 with 4-step gradient accumulation on 4 GPUs tell you the U-Net-in-gene-space approach is expensive; hyperparameters for PRiMeFlow itself were hand-picked rather than swept, explicitly due to compute constraints, while the baselines got 60-trial Optuna sweeps. That asymmetry cuts in the baselines' favor, so the ablation conclusions are probably conservative rather than inflated.

Want me to add the Feng 2026 iPSC dataset to `public-data.md` with a proper tier assignment?

---

**User**

can you explain some of the jargon you just used?

---

**Cursor**

Sure. Grouping them by theme, since a lot of the terms only make sense relative to each other.

## The generative modeling stack

**Flow matching** is a way to train a generative model. You want to sample realistic single cells, but you only have examples, not a formula. So you set up a path that gradually morphs a simple distribution (Gaussian noise) into your data distribution over a fake time variable t running from 0 to 1, and you train a network to predict which direction to move at every point along the way. To generate a new cell you start from random noise and follow those directions. It's a cousin of diffusion models — same idea of "learn to walk from noise to data," slightly cleaner math.

The thing the network predicts is a **velocity field**: given a position (a vector of ~18,000 gene expression values) and a time, it returns another vector of the same size saying which way and how fast to move. "Field" just means it's defined at every point in space, like a wind map.

**Conditional flow matching** is the trick that makes this trainable. The mathematically correct target velocity is an intractable integral over the whole dataset. But if you first pick one specific noise sample x_0 and one specific real cell x_1, the velocity along the straight line between them is just x_1 - x_0 — trivial. Regressing against that noisy-but-unbiased target turns out to give the same average gradient as the impossible version. The practical payoff is that training never requires simulating the generation process, only single random snapshots of it.

**Coupling** refers to how you pair noise samples with real cells during training. **Independent coupling** means pairing them at random. **Optimal transport coupling** means solving a matching problem within each minibatch so each noise sample is paired with the real cell it's already closest to, which makes the learned paths straighter and shorter. CellFlow's core design bet is that OT coupling helps; this paper's finding is that it doesn't.

**Classifier-free guidance (CFG)** is borrowed from image generation — it's the same dial as the "guidance scale" slider in Stable Diffusion. During training you randomly hide the conditioning information some fraction of the time (here 20%), so the model learns both "generate a cell given perturbation X" and "generate any cell." At sampling time you extrapolate away from the unconditional prediction: v = v_{\text{uncond}} + w(v_{\text{cond}} - v_{\text{uncond}}). At w=1 you get the plain conditional model. Above that you exaggerate whatever makes perturbation X distinctive — more on-target, but more distorted and less diverse. That's exactly why cranking CFG improved their DEG recall while wrecking their distributional fit, and why they used a huge weight of 20 to inflate the perturbation-discrimination metric.

A **U-Net** is a convolutional architecture from image segmentation: progressively downsample the input, then upsample back, with skip connections between matching levels. It's designed for data where neighbors matter. Applying it to a gene expression vector means treating the ~18,000 genes as a 1D sequence and sliding filters along it — which is odd, because gene order in an AnnData object is arbitrary (usually alphabetical or by genome position). Hence the paper's admission that they don't know why it works.

**log1p space** means the data has been transformed by \log(1+x). Raw counts are wildly skewed (most genes near zero, a few in the thousands), so logging compresses that; the +1 avoids \log 0. `expm1` is the inverse, e^x - 1, which is how they get back to count space for normalization.

**Latent space** or **PCA space** modeling means you first squash cells from ~18,000 genes down to something like 50–100 dimensions, do all the modeling there, then decode back. Much cheaper, but the compression discards information — and the paper's central claim is that the discarded part is exactly the subtle perturbation-specific signal you're being scored on. A **VAE** (variational autoencoder) is one such compression scheme, with a learned encoder and decoder trained to reconstruct the input; CPA and SAMS-VAE are built on this.

**ODE** is an ordinary differential equation — here dx/dt = v(x,t), which you integrate from t=0 to t=1 to generate a cell. You solve it numerically. **Euler** integration takes many small fixed-size steps, the crudest method. **Dopri5** (Dormand–Prince) is an adaptive Runge-Kutta method that adjusts step size on the fly and is much more accurate per unit of compute.

## Evaluation metrics

*The metric names below are the ones PRiMeFlow was scored on, which is the **2025** set. None of DES, MAE, AUPRC, Spearman, or Pearson exists in 2026; see `2026-metric-explanation.md` and `.cursor/rules/project-rules.mdc` for the six that replaced them. PDS survives in modified form as `pds`. The general concepts here — pseudobulk, DEG, MMD, rank metrics, mean predictors — all still apply.*

**Pseudobulk** means averaging expression across all cells in a condition, collapsing thousands of cells into one vector per perturbation. It's what bulk RNA-seq would have given you. It's convenient and it's what most metrics use, but it throws away all cell-to-cell variation — which is the specific thing PRiMeFlow is trying to model. That tension is the paper's main narrative.

**DEG** stands for differentially expressed gene: a gene whose expression changes significantly between perturbed and control cells, usually established with a Wilcoxon rank-sum test plus a multiple-testing correction. **DEG recall** asks: of the true top-50 DEGs, what fraction did the model also identify? Recall only — it doesn't penalize false positives, which is why the DES metric adds a cap.

**MMD** (maximum mean discrepancy) is a distance between two *distributions* given only samples from each. Loosely: run both point clouds through a kernel, compare their means in that transformed space. Zero means statistically indistinguishable. This is the metric that asks "does my generated cloud of 1,000 cells actually look like the real cloud?" rather than "is my average right?"

The **rank versions** of PerturBench metrics are the important subtlety. Instead of asking "how close is my prediction for perturbation A to the truth for A," they ask "among all my predictions, where does my prediction for A rank when compared against the truth for A?" A model that outputs nearly the same profile for every perturbation can score well on the raw metric and terribly on the rank version. Lower is better for rank metrics, since they're a normalized rank position.

For the 2025 challenge's own seven metrics: **DES** (differential expression score) is the overlap between your predicted DEG set and the true one, normalized by the true set's size, with your set truncated so you can't win by declaring every gene differentially expressed. **PDS** (perturbation discrimination score) computes the distance from your prediction to *every* true perturbation's effect and checks where the correct one ranks — 1.0 means your prediction is closest to the right answer, 0.5 is chance. **MAE** is mean absolute error averaged over genes, and it's very sensitive to overall scale, which is why teams submitting z-scored or unnormalized output got MAEs of 0.77 and 6.7. **AUPRC** treats "is this gene really a DEG" as a binary classification problem and measures how well your scores rank true DEGs above the rest across all thresholds; it's the standard choice when positives are rare. **Spearman** correlations operate on ranks (do you get the ordering right), **Pearson** on the actual values (is the linear relationship right), and "delta" or "LFC" means the correlation is computed on the *change* from control rather than on absolute expression.

A **perturbation delta** is pseudobulk(perturbed) minus pseudobulk(control) — it isolates the effect of the perturbation from the baseline state of the cell line.

A **mean predictor** is a baseline that predicts the identical profile for every perturbation, typically the average effect across the training set. It's the key sanity check in this field: if a mean predictor beats your model on a metric, that metric isn't rewarding perturbation specificity. It scored DES 0.305 on the 2025 private test set, higher than PRiMeFlow's 0.228 — which is exactly why 2026 rebuilt its scoring around it. The context mean is now the reference-scaling zero point, so a mean predictor scores precisely 0.0 by construction rather than winning a metric.

## Benchmark task names

**Covariate transfer** is the task shape of both the challenge and most of this paper. A "covariate" is the context variable — cell line, drug treatment, timepoint. You've observed perturbation P in cell lines A and B, and you must predict P in cell line C where it was never measured. Crucially, you *know in advance* which perturbations you'll be asked about, which is what makes the paper's finetuning strategy legal: they trained on those perturbation identities as measured in other cell lines, never on H1 ground truth.

**Combo prediction** is the other task: you've seen gene A knocked down alone and gene B alone, now predict what happens when both are hit simultaneously.

An **ablation study** means removing or swapping one component, retraining, and measuring what you lose. It's how you show a design choice actually earns its keep rather than just coexisting with good results.

## Biology and data

**CRISPRi** is CRISPR interference. A catalytically dead Cas9 (dCas9 — it binds DNA but can't cut) is fused to a repressor domain like KRAB and steered to a gene's promoter by a guide RNA, shutting down transcription. The result is a **knockdown**: expression is reduced, typically partially and reversibly. Contrast with a knockout, where the gene is physically disrupted. This matters because knockdown strength varies from cell to cell and from one dCas9 effector variant to another, which is a real source of label noise.

**Perturb-seq** is a pooled CRISPR screen read out by single-cell RNA sequencing. Each cell receives one guide RNA, and you sequence both the guide (so you know which gene was hit) and the full transcriptome (so you see the consequences). One experiment yields thousands of perturbation-response pairs.

**hESC** is a human embryonic stem cell; H1 is one specific well-characterized line, and it's the cell line of the 2025 challenge. **iPSC** is an induced pluripotent stem cell — an adult cell reprogrammed back to a stem-cell-like state. Both are pluripotent, which is why the Feng 2026 iPSC dataset is unusually on-target for a stem-cell prediction task.

**UMI** stands for unique molecular identifier: a random barcode attached to each mRNA molecule before amplification, so you can count original molecules rather than PCR copies. UMI counts per cell vary several-fold for purely technical reasons (**sequencing depth**), so before comparing cells you rescale each one to a common total — usually the median across the dataset, hence **median normalization** and the `median_umi_per_cell` value the 2025 challenge supplied. In 2026 you infer the equivalent from the provided control cells, whose median is around 20,000 UMI.

A **non-targeting control** (NTC) cell received a guide RNA that matches nothing in the genome. It controls for the stress of the CRISPR machinery itself, so your measured effect is "gene knocked down" rather than "gene knocked down plus lentivirus plus dCas9."

## Engineering odds and ends

**HPO** is hyperparameter optimization, and **Optuna** is a library that automates the search. A "60-trial sweep" means 60 separate training runs with different settings. The relevant point is that the baselines got 60 trials each while PRiMeFlow's settings were chosen by hand, so the comparison is if anything tilted against the paper's own model.

**Gradient accumulation** is a memory workaround. If a batch of 216 cells won't fit on the GPU, you run four batches of 54, sum the gradients, and only then update the weights — mathematically close to one large batch. Seeing batch sizes of 38 to 54 with 4-step accumulation tells you the model is memory-hungry, which is the real cost of working in full gene space.

Happy to go deeper on any of these — the flow matching math and the metric definitions are the two places where a bit more detail would probably change how you'd design against them.