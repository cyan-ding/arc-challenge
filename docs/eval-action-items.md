Here's an ordered list, grouped by when to do it. Everything through step 7 is doable this week with nothing but a laptop; it's designed so the very first things you build are the ones the modeling team is blocked on.

## Today: environment and ground truth

1. **Set up the Python environment and pin it.** Python 3.11 or 3.12 — not 3.13, because `cellstream` (the `[scale]` extra you may want later) only ships wheels for 3.11/3.12 and falls back to a Rust build otherwise. Use `uv`, create a `pyproject.toml` at the repo root, and pin `cell-eval2==0.16.0` (current PyPI) plus `pdex` explicitly. `pdex` matters: without it the DE backend falls back to scanpy, which is much slower, and the version pin matters because `vcc2026` profile membership changed between v0.7.x and v0.8.0. Every run should record the version in its metadata.

2. **Install the VCC CLI and pull the controls.** `uv tool install vcc-cli`, `vcc login --token-stdin`, `vcc whoami` (confirm it says ready to submit — registration approval can lag), then `vcc datasets download controls -d ~/vcc`. Unzip `gene_names.csv` and `pert_counts.csv` and **commit those two files to the repo**. They're small and they *are* the gene-space and perturbation contracts the other two teams need to code against.

3. **Read the metrics spec PDFs.** `docs/vcc2026_metrics/vcc2026-metrics.pdf` (full) and `-brief.pdf` (definitions and measured values) in the `cell-eval2` repo. You're reading for four specific answers: how `fid` scales by yield, what quantity `reach` ranks genes by, the published baseline and replicate anchor values for each metric, and exactly how the five disjoint half-splits are constructed. Write those into `docs/vcc2026-metrics-notes.md`. This is the single highest-leverage thing you can do this week — the modeling team's `phi` sweep is guesswork until those are pinned down.

4. **Run `vcc sample` → `vcc prep --dry-run`.** Not to submit, just to see what the validator's output looks like and confirm the toolchain works end to end. `vcc sample -g gene_names.csv -p pert_counts.csv -o sample.vcc`, then `vcc prep` on an h5ad (note `--perts` has no short form in `prep`; `-p` there means `--pert-col`, and the wiki calls out this exact mixup).

## Days 1–2: get the scorer running on real ground truth

5. **Run the `cell-eval2` tutorial end to end.** Clone the repo, use the committed 5 MB H1 subset (`docs/data/H1-VCC-2025-training.h5ad`), and run `cell-eval2 run ... --preset vcc2026 --pert-col target_gene --set de.backend=pdex`. You want to see all ten output columns (six scored, four diagnostics) and understand the shape of `results.csv` (long form, one row per perturbation × metric), `agg_results.csv`, and `run_meta.json`.

6. **Get the full 2025 H1 data locally.** It's on Hugging Face as `arcinstitute/VCC_train` and on the Arc atlas GCS bucket (see `docs/arc-atlas-data.md`; the GCS route is Requester Pays). This is your only local dataset with the exact assay and real ground truth. First thing to check when it loads: whether its gene space is 18,533 and matches `gene_names.csv`. The docs flag that it might not, and if it doesn't, harmonizing it is your job before anyone can score against it.

7. **Validate the two landmarks.** Build two predictions from the H1 data: the context-mean (every perturbation gets the global mean profile) and a split-half replicate (half the real cells scored as a prediction of the other half). Run `cell-eval2 baseline`, `cell-eval2 run`, and `cell-eval2 score --user-agg --baseline-agg`. After scaling, the context-mean should land at ~0 and the replicate at ~1 on every metric. If it doesn't, either your scaling is wrong or your understanding of the anchors is — and finding that out on day two rather than week six is the whole point.

## Week 1: turn it into a package the other teams call

8. **Repo skeleton.** `src/vcc/` package, `pyproject.toml`, `Makefile` with `make score`, `make validate`, `make submit`, and a `tests/` directory. The tutorial's 5 MB h5ad is a fine committed test fixture.

9. **`score.py`.** One function, `score(pred, real, ...) -> DataFrame`, wrapping `compute_metrics` with `EvalConfig.from_preset("vcc2026")`, `pert_col="target_gene"`, `control_source="real"` (the v2 default — matches how the challenge uses held-out controls), `de.backend="pdex"` pinned, and **persistent `cache_real`/`cache_pred` directories** so the real side is computed once and reused across every prediction in a sweep. Return raw per-perturbation results plus per-context aggregates.

10. **`anchors.py` and `scale.py`.** This is the real engineering. The `b` anchor you can get from `cell-eval2 baseline`. The `r` anchor you largely have to build yourself: `compute_ceiling` exists, but it's a Spearman-Brown correction on a hand-maintained list of metrics, and of the six scored in `vcc2026` only `pds_cosine` is on it — the other five return `NaN`. So implement the five-disjoint-half-split replicate directly, then apply `s = (u − b)/(r − b)` per context with the documented clamps (`mse` to `[0, 1]`, `nmae` floored at −6). Validate against step 7.

11. **`validate.py`.** `vcc prep --dry-run` catches format problems, but you want to catch them *before* writing a 3 GB file, and there's one check `prep` can't do: swapped context labels. Add a sanity check that each context's predicted pseudobulk correlates best with its own context's controls. Mirror the rest of the requirements table too — exactly 360,000 cells, 400 per perturbation per context, 18,533 genes in order, whole non-negative finite values, no `non-targeting` rows, ≤1,000,000 counts per cell, CSR sparse, ≤4.75 B stored entries.

12. **Measure two numbers and write them down.** First, wall-clock time to score one 300-perturbation × 400-cell × 18,533-gene block on your CPU, extrapolated to three contexts — that decides whether the joint `weights` × `phi` sweep is feasible on laptops or needs a GPU box (Arc reports ~34 s per 300-perturbation set on an H100). Second, the **noise floor**: take one fixed Δ, resample the 400 cells with five different seeds, score each, and record the standard deviation of all eighteen numbers. That's the smallest improvement anyone on the modeling team is allowed to believe.

## Week 2: make it a shared instrument

13. **Submission log and discipline.** `submissions/log.csv` with entry id, git SHA, config hash, local CV score, and the six per-context scores. A `make submit` target that refuses to run on a dirty tree and puts the SHA in the `-m` message. Agree on one named submitter so two people don't burn both daily slots.

14. **Local-panel helpers and the caveat doc.** A utility to subsample any local reference to 300 perturbations × 400 cells, because `pds` is a within-panel rank and Replogle Essential's 2,393 perturbations is a different problem from the leaderboard's 300. Document plainly that local scaled scores are for *ranking configs*, not predicting the board — the source data is shallower than 400 cells per perturbation, so the local replicate anchor sits lower and inflates scaled scores.

15. **Hand off the contracts.** Give modeling `score()` and the noise floor. Give the data team the exact reference-h5ad format the scorer expects: raw counts, `obs.target_gene`, `obs.context`, `non-targeting` rows present on the *real* side only.

Two things I'd flag rather than assert. I confirmed the two spec PDFs exist and what the README says they contain, but I haven't read them, so the anchor values and the `fid`/`reach` definitions are yours to extract in step 3. And the tutorial and older `cell-eval` material assume predictions contain `non-targeting` cells — that's the 2025 convention. The 2026 CLI wiki is explicit that a submission containing them is rejected, and `control_source="real"` in the scorer is what makes local scoring consistent with that.

I'm in Ask mode so I can't scaffold any of this, but if you switch to Agent mode I can lay down the `pyproject.toml`, the package skeleton, and a first `score.py` against the tutorial fixture so step 5 through 9 have a starting point.