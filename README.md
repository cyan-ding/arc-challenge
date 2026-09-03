# Arc Virtual Cell Challenge 2026

Zero-shot prediction of CRISPRi knockdown responses in unseen cell lines.

Arc ran CRISPRi Perturb-seq on six anonymized cell lines. For three of them per round we get
only non-targeting control profiles plus a list of 300 knockdowns, and must predict
post-perturbation expression with no matched perturbation data from that context.

Submission is one `.h5ad`: 360,000 cells (300 perturbations × 400 cells × 3 contexts),
18,533 genes, raw integer counts, scored by `cell-eval2` on the `vcc2026` preset over six
reference-scaled metrics where 0 = the context mean and 1 = a replicate experiment.

Final test set released Oct 22, 2026 · submissions close Nov 5, 2026.

## Start here

`.cursor/rules/project-rules.mdc` is the authoritative spec — task, submission contract,
metrics, timeline. Read it first; it overrides everything in `docs/`, and it overrides
pretrained knowledge of the 2025 challenge, which differed in nearly every detail.

## Docs

| file | contents |
| --- | --- |
| `docs/claude-recipe.md` | The modeling plan. Two stages: predict a per-(target, context) effect vector Δ over 18,533 genes, then shift real control cells by Δ to produce 400 cells. Covers shrinkage of noisy public log2FC estimates, response-program factorization, context transfer, and which knobs to tune against which metric. |
| `docs/2026-metric-explanation.md` | End-to-end walkthrough of how a submission becomes six scores: the pseudobulk vector and Wilcoxon DE table per 400-cell block, worked toy examples of `jac`, `fid`, `nmae`, `reach`, `pds`, and why predicted cell-to-cell variance is itself scored. |
| `docs/public-data.md` | Training-data survey, tiered by usefulness. Replogle 2022 + Nadig 2025 (~2,000 shared perturbations × 4 cell lines, CRISPRi, plus precomputed log2FC/SE/p-values on figshare), Jiang's crossed six-line grid, Feng 2026 iPSC lines, and the per-perturbation depth caveats that make raw pseudobulk LFCs mostly noise. |
| `docs/arc-atlas-data.md` | Arc Virtual Cell Atlas access notes: the single Requester-Pays GCS bucket that replaced the dead ones, plus VCC 2025 H1 data (exact task match, use as a local validation harness), Tahoe-100M (drugs across 50 lines — the cross-context transfer signal), and scBaseCount. |
| `docs/primeflow-analysis.md` | Analysis of PRiMeFlow (Altos Labs), a conditional flow-matching model that competed in **2025**. Useful ablations: full gene space beats PCA, U-Net beats MLP, optimal-transport coupling adds nothing. Its metric commentary is historical (DES/PDS/MAE/AUPRC) and does not apply to 2026. |
| `docs/primeflow-paper.pdf` | Source PDF for the above. |

## Resources

- [Evaluation](https://virtualcellchallenge.org/evaluation) · [Datasets](https://virtualcellchallenge.org/datasets)
- [CLI guide](https://vcc-cli-wiki.virtualcellchallenge.org/) · [cell-eval2](https://github.com/ArcInstitute/cell-eval2)
