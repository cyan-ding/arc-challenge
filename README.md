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

`.cursor/rules/challenge-rules.mdc` is the authoritative spec — task, submission contract,
metrics, timeline. Read it first; it overrides everything in `docs/`, and it overrides
pretrained knowledge of the 2025 challenge, which differed in nearly every detail.

## Setup

Python 3.12 via `uv` (pinned in `.python-version`; 3.11 also works, 3.13+ does not — see
`pyproject.toml`). `cell-eval2` is pinned exactly because the `vcc2026` profile changed between
releases.

```bash
uv sync --extra dev                          # creates .venv with cell-eval2 0.16.0 + pdex
uv run cell-eval2 --help                     # scorer
uv tool install vcc-cli && vcc --version     # submission CLI (needs `uv tool update-shell` once)
vcc login --token-stdin                      # key from Account -> Credentials on the portal
vcc datasets download controls -d ~/vcc      # A/B/C controls + gene_names.csv + pert_counts.csv
```

After downloading, copy `gene_names.csv` and `pert_counts.csv` into `data/contracts/` and commit
them. They are the gene-space and perturbation contracts every team codes against. Nothing else
under `data/` is tracked.

The 2025 H1 data (our local ground truth) is public over plain HTTPS, no GCP account needed:

```bash
B=https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025
curl -C - -o data/vcc2025/adata_Validation.h5ad $B/validation/adata_Validation.h5ad   # 6.9 GB
curl -C - -o data/vcc2025/adata_Training.h5ad   $B/train/adata_Training.h5ad           # 15.5 GB
```

See `docs/vcc2026-metrics-notes.md` §8 for its gene space vs 2026 and the local scoring bundle.

## Docs

| file | contents |
| --- | --- |
| `docs/lab-notebook.md` | Chronological record of every working session: which action items were closed, what was found, which commit it landed in, and what's still open. Read this to catch up on where the eval stack is. |
| `docs/eval-action-items.md` | The 15-item plan for the evaluation stack, in the order the items should be done. |
| `docs/vcc2026-metrics-notes.md` | What the official metric spec PDFs actually say: exact `fid` and `reach` definitions, the published baseline/replicate anchors per metric, how the five half-splits are built, the upload-vs-scoring contract split that explains the `non-targeting` confusion, which `cell-eval2` 0.16.0 subcommands build the anchors locally, and what `vcc prep` does and doesn't catch. Includes corrections to the other docs. |
| `docs/claude-recipe.md` | The modeling plan. Two stages: predict a per-(target, context) effect vector Δ over 18,533 genes, then shift real control cells by Δ to produce 400 cells. Covers shrinkage of noisy public log2FC estimates, response-program factorization, context transfer, and which knobs to tune against which metric. |
| `docs/2026-metric-explanation.md` | End-to-end walkthrough of how a submission becomes six scores: the pseudobulk vector and Wilcoxon DE table per 400-cell block, worked toy examples of `jac`, `fid`, `nmae`, `reach`, `pds`, and why predicted cell-to-cell variance is itself scored. |
| `docs/public-data.md` | Training-data survey, tiered by usefulness. Replogle 2022 + Nadig 2025 (~2,000 shared perturbations × 4 cell lines, CRISPRi, plus precomputed log2FC/SE/p-values on figshare), Jiang's crossed six-line grid, Feng 2026 iPSC lines, and the per-perturbation depth caveats that make raw pseudobulk LFCs mostly noise. |
| `docs/arc-atlas-data.md` | Arc Virtual Cell Atlas access notes: the single Requester-Pays GCS bucket that replaced the dead ones, plus VCC 2025 H1 data (exact task match, use as a local validation harness), Tahoe-100M (drugs across 50 lines — the cross-context transfer signal), and scBaseCount. |
| `docs/primeflow-analysis.md` | Analysis of PRiMeFlow (Altos Labs), a conditional flow-matching model that competed in **2025**. Useful ablations: full gene space beats PCA, U-Net beats MLP, optimal-transport coupling adds nothing. Its metric commentary is historical (DES/PDS/MAE/AUPRC) and does not apply to 2026. |
| `docs/primeflow-paper.pdf` | Source PDF for the above. |

## Resources

- [Evaluation](https://virtualcellchallenge.org/evaluation) · [Datasets](https://virtualcellchallenge.org/datasets)
- [CLI guide](https://vcc-cli-wiki.virtualcellchallenge.org/) · [cell-eval2](https://github.com/ArcInstitute/cell-eval2)
