# Lab notebook

Chronological record of what was done, what was found, and what it changed. One entry per
working session, newest at the bottom. Each entry says which action items from
`eval-action-items.md` (or later plans) it closes, the commit(s) it landed in, and anything a
teammate would need to reproduce or build on it. Findings that matter beyond the day get promoted
into the reference docs (`vcc2026-metrics-notes.md`, the rules file); this file is the trail.

Conventions: dates are the session date (US Pacific). "Item N" refers to `eval-action-items.md`.
Anything under `data/` is git-ignored; paths there are on the author's machine unless noted.

---

## 2026-09-03 — Team split, eval plan

**Done.** Team divided into data/engineering (Hemani + 1), modeling (Anisha, Ben + 1), and
evaluation (cyan + 1). Wrote the eval-stack plan that became `eval-action-items.md`: 15 items,
front-loaded so the modeling team's blocker (a local scorer) comes first.

**Decided.** Three interface contracts to freeze before anything else: the gene space
(`gene_names.csv` order), the Δ table (`(target, context) → float32[18533]` plus SE), and the
scorer (`score(pred, real) → 18 raw + 6 scaled`). A throwaway end-to-end submission in week one to
force all three to exist.

**Looked up.** `cell-eval2` README and tutorial. `compute_ceiling` covers only `pds_cosine` among
the six scored metrics — flagged as meaning the replicate anchor might need hand-building
(revised on 09-04, see below).

---

## 2026-09-04 (a) — Items 1–4: environment, CLI, spec, smoke test

**Commits.** `0ac1203`, `f08d3a9`.

**Item 1 — environment.** `pyproject.toml` with `cell-eval2==0.16.0` and `pdex` pinned,
`requires-python >=3.11,<3.13` (`cellstream` wheels stop at 3.12). `uv python pin 3.12`,
`uv sync --extra dev`, `uv.lock` committed. Confirmed `EvalConfig.from_preset("vcc2026")` loads
with `control_source=real`, 5 CPM filter, per-perturbation BH, cosine PDS with
`exclusion_scope=panel`. Preset default `pert_col` is `target`, so `--pert-col target_gene` must
always be passed. `.gitignore` blocks every matrix format and `data/**` except `data/contracts/`.

**Item 2 — VCC CLI.** `uv tool install vcc-cli` → 0.2.0 at `~/.local/bin/vcc`. First attempt hit
a 30 s network timeout on one wheel; `UV_HTTP_TIMEOUT=180` fixed it. Login needs the user's API
key; done by cyan, controls bundle downloaded. `gene_names.csv` (18,533 genes) and
`pert_counts.csv` (300 targets) committed to `data/contracts/`. Both have a header row
(`gene_name`, `target_gene`) even though `vcc sample --help` calls the gene list headerless.

**Item 3 — spec.** Read `vcc2026-metrics.pdf` and `-brief.pdf` (dated 2026-08-19, measured on
cell-eval2 0.15.0, rule_version 3). Wrote `docs/vcc2026-metrics-notes.md`. Answers to the four
open questions from `claude-recipe.md`:

- `fid = k / max(n_pred, n_real)` — precision × coverage capped at 1. Over-calling free,
  under-calling charged.
- `reach` ranks by our significant calls first, then ascending adjusted p, raw p, descending |lfc|,
  gene name. Driven by our p-values, hence by `phi`.
- Published `b`/`r` per metric as ranges over A/B/C (table in the notes §3).
- Five disjoint half-splits per perturbation *and per control*, each half using its own controls,
  base seed 0. `nmae` uses a full-gate estimator.

Two things the spec settles that our docs had wrong: `pds_cosine` excludes **all 300** panel
target genes, not just its own; and the baseline `b` is an oracle built from the *real perturbed*
cells (filtered by knockdown efficiency), not from the controls. Both flagged in the notes §7 for
folding into the rules file. Also: the spec names **two contracts** — upload (`vcc prep`, no
controls, exactly 400 cells) and scoring (`cell_eval2`, controls required) — which is why
`cell-eval2` docs and the CLI wiki appear to contradict each other on `non-targeting` rows.

Revised the 09-03 assumption: 0.16.0 ships `prep-real-bundle`, `run --anchor`,
`run --lfc-nmae-ref`, and `score --real-bundle`, so the replicate anchor does **not** need
hand-building. Item 10 shrinks to "run these and validate."

**Item 4 — smoke test.** With placeholder lists of the right shape (before the real bundle was
downloaded): `vcc sample --h5ad` → 360,000 × 18,533 CSR in 17 s; `vcc prep --dry-run` passes in
1.6 s. Four deliberate breaks each rejected with a clear message (fractional values,
`non-targeting` rows, 399 cells, moved block). A moved block was caught **only** because it left
one context short a perturbation — a clean A↔B relabel of whole contexts passes `prep` silently.
That is the check `validate.py` (item 11) has to add.

**Open.** Re-run the smoke test against the real contract files.

---

## 2026-09-04 (b) — Items 5–7: scorer end to end, 2025 data, landmarks

**Commit.** `4407093`. Full write-up in `vcc2026-metrics-notes.md` §8.

**Item 5 — tutorial fixture.** `cell-eval2 run --preset vcc2026 --pert-col target_gene --set
de.backend=pdex` on the 5 MB `H1-VCC-2025-training.h5ad` (600 cells × 1,000 genes, 5 perts).
All ten columns present. A prediction **without** `non-targeting` rows is refused:
`perturbation sets differ: 5 pred vs 6 real; real-only=['non-targeting']`. So `score()` must
concatenate the real controls onto the prediction before calling the scorer — the same bridge the
platform performs. Only 1 of 5 targets resolved to a gene in the 1,000-gene fixture; the scorer
warns and scores the rest without target exclusion (`EvalConfig.target_gene_map` exists for
symbol mismatches). Not an issue on the full gene space.

**Item 6 — 2025 data.** `arcinstitute/VCC_train` on Hugging Face **no longer exists**. The Arc
bucket serves the processed 2025 files over plain HTTPS, no GCP account:
`https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/`.
Downloaded validation (6.9 GB, 50 perts, 98,927 cells, 38,176 NTC) and training (15.5 GB, 150
perts) to `data/vcc2025/` at ~50 MB/s. A `nohup` background download died when its shell exited;
relaunched through the tool's background mode with `curl -C -`.

Gene space: 18,080 (2025) vs 18,533 (2026); 18,077 shared; the 3 dropped are
`var_names_make_unique` suffixes (`HSPA14-1`, `TBCE-1`, `TMSB15B-1`); 456 are 2026-only, mostly
`AC…` lncRNAs; order differs. All 300 panel genes measured in both. **25 panel genes were 2025
targets** — real H1 ground truth for them. Scoring against 2025 data needs no harmonization;
using it as training input does (data team).

**Item 7 — landmarks.** Built `data/local_ref/vcc2025-val-ref400/ref.h5ad`: validation set at the
challenge shape, 44 perts × 400 cells (6 perts with < 400 cells dropped) + 8,000 NTC, seed 0.
Pipeline: `baseline --save-pred` → `prep-real-bundle --baseline <pred h5ad>` → `run` →
`score --real-bundle`. Two flag lessons: `--baseline` takes the baseline **prediction h5ad**, not
the output dir; and that prediction is fractional, passes only via `allow_fractional_counts`
inside `baseline`, and cannot be self-scored (no anchor identity in `baseline_meta.json`). The 0
end is the bundle's zero by construction.

Results (`from_replicate`): exact reproduction → `mse` **1.000000**, `pds` 1.106, `reach` 1.019,
`fid` 1.163, `nmae` 1.587, `jac` 2.281, avg 1.359 — the spec's pattern. Control paste → `mse`
0.000, `nmae` −0.084 (raw ≈ 1.008 as the spec predicts), `pds` −0.071, others negative, avg
−0.098. **The harness is wired correctly.**

Local anchors vs Arc's: `r` slightly above on every DE metric (deeper data, 8,000 controls); `b`
differs materially (`fid` 0.22 vs 0.51; `nmae` 0.95 vs 1.00). Local scaled scores rank configs;
they don't predict the board (item 14's caveat, now measured).

**Timing (item 12, first half).** M4 Max, 14 cores, CPU `pdex`, 44 perts: `baseline` 64 s;
`prep-real-bundle` 5 min 8 s; **one prediction 38.5 s** with the real side cached. Extrapolates to
≈ 13 min per full 3 × 300 prediction on this laptop, ≈ 35 min per context to build a bundle once.

**Artifacts kept.** `data/local_ref/vcc2025-val-ref400/{ref.h5ad, bundle/, cache/,
scored_exact.csv, scored_ctrlpaste.csv}`.

**Open.** Same landmark check on the 150-pert training set; resampling noise floor (item 12,
second half); fold the pipeline into `src/vcc/` (items 8–10) with the control bridge built into
`score()`.

---

## 2026-09-04 (c) — Housekeeping

Renamed `.cursor/rules/project-rules.mdc` → `challenge-rules.mdc` (content unchanged); updated
the three docs that referenced the old name. `attitude.mdc` now asks that completed action items
be logged here. Created this notebook with the entries above.
