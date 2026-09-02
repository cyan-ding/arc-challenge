Here's the orientation first, because the three links aren't three parallel things.

The [Google Cloud Marketplace listing](https://console.cloud.google.com/marketplace/product/bigquery-public-data/arc-institute?project=gcp-public-data-arc-institute) isn't a dataset — it's the **access layer** for all of them. Tahoe-100M and scBaseCount are the two datasets Arc "bootstrapped" the atlas with, and they're very different in kind. There's also a third folder in that repo that nobody links prominently but which matters most for you: `[virtual-cell-challenge/](https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/virtual-cell-challenge)`, holding last year's data.

Ranked by relevance to the 2026 zero-shot task:


| dataset       | scale                         | perturbation type           | why you'd use it                                |
| ------------- | ----------------------------- | --------------------------- | ----------------------------------------------- |
| VCC 2025 (H1) | ~300k cells, 300 target genes | **CRISPRi knockdown**       | exact task match; your local validation harness |
| Tahoe-100M    | 100.6M cells, 50 cell lines   | small-molecule drugs        | the cross-context transfer signal               |
| scBaseCount   | 502M+ cells, 27 organisms     | mostly none (observational) | pretraining a basal-state representation        |




## The access layer

Everything now lives in one bucket, `gs://arc-institute-virtual-cell-atlas`. Two things to internalize before you write any download code.

**The old buckets are dead.** `gs://arc-ctc-tahoe100/` and `gs://arc-scbasecount` were removed March 31, 2026. Every tutorial, blog post, and Stack Overflow answer written before then points at paths that no longer resolve, including some cells inside Arc's own notebooks. When something 404s, that's usually why.

**It's Requester Pays.** You need a GCP project with billing enabled, and you pass your project ID with every request so egress bills to you. First 2TB/month is free, which is generous but not unlimited when the raw Tahoe matrices are involved. With `gsutil` that's `-u YOUR_PROJECT`; with `gcsfs` it's `GCSFileSystem(project="YOUR_PROJECT", requester_pays=True)`.

The BigQuery side of the listing is genuinely useful: you can SQL the metadata tables to work out which cells you want before downloading a single matrix.

## Tahoe-100M

100,648,790 cells across 1,344 samples: **50 cancer cell lines × ~1,100 small-molecule drugs**, generated on Vevo Therapeutics' Mosaic platform ([paper](https://doi.org/10.1101/2025.02.20.639398), also mirrored on [HuggingFace](https://huggingface.co/datasets/tahoebio/Tahoe-100M)).

The thing to understand is the mismatch and the match. The **perturbation modality is wrong** for your task — these are drugs, not CRISPRi knockdowns, so you can't learn "knock down gene X" directly from it. But it's the best available source of the thing the 2026 challenge actually tests: the same perturbation applied across 50 different cellular contexts. If you want to learn how a response transfers and re-shapes between cell lines, this is the only dataset at that scale with that design.

Structurally it splits into matrices and metadata. The count matrices are per-plate `h5ad` files (`.h5ad.gz` means internal gzip, which `scanpy` handles), and they're large — the tutorial deliberately ships a 2,000-cell sample instead. Note the var dimension is **62,710 genes**, not the challenge's 18,533, so gene-space alignment is real work you'll have to do.

The cell metadata sits separately in `metadata/obs_metadata.parquet`, which is the right entry point: query it first, decide what you need, then pull only those plates. The columns that matter are `drug` and `drugname_drugconc` (the second one carries dose — the same drug at different concentrations is a different perturbation, don't collapse them), `cell_line` as a Cellosaurus ID with `cell_name` as the human-readable version, `pass_filter` for QC (the "Full" setting is stricter on `gene_count` and `tscp_count`), and `pcnt_mito`, `S_score`, `G2M_score`, `phase` as covariates you may want to condition on or regress out.

Two practical notes. The HuggingFace copy stores expression in a **tokenized ragged format** rather than a sparse matrix — `genes` holds integer token IDs and `expressions` holds the aligned counts, and **the first entry of every row is a CLS marker token you must drop**. Map token IDs to symbols via `gene_metadata.parquet`. And HuggingFace also hosts a precomputed `pseudobulk_differential_expression` table, which will save you an enormous amount of compute if pseudobulk log-fold-changes are what you're modeling — which, given how the challenge is scored, they probably are.

## scBaseCount

502M+ cells, 27 organisms, 75 tissues, assembled by AI agents that crawl the Sequence Read Archive, extract metadata, and reprocess everything uniformly through STARsolo ([paper](https://doi.org/10.1101/2025.02.27.640494)).

This one is **mostly observational** — no designed perturbation experiment. Its value for you is pretraining: your only input at test time is the unperturbed profile of a cell line you've never seen, so a model with a strong prior over what basal cell states look like has a real advantage. There is a `perturbation` field in the sample metadata, but it's incidental and heterogeneous, not a curated panel.

Two structural traps here.

**Use the 2026-01-12 release, not the 2025-02-01 one.** The README calls the initial release deprecated; it's a different cell count (230M vs 502M) and different metadata schema.

**Every h5ad carries multiple layers, and picking the wrong one is a silent bug.** `adata.X` is the `Gene`/`Unique` count (exonic reads, unique UMI assignment) — that's the one that corresponds to conventional 10x processing and therefore to the challenge data. But the files also contain `UniqueAndMult-EM` and `UniqueAndMult-Uniform` layers, and the Velocyto feature type gives you `spliced`/`unspliced`/`ambiguous`. The `GeneFull`* variants count intronic reads too and therefore have systematically higher counts. Mixing feature types across your training set shifts the count distribution with nothing throwing an error.

The metadata is LLM-extracted and correspondingly noisy. Arc is upfront about this: disease labels are pulled at **study level** from abstracts, so a "COVID-19" study propagates that label to its healthy controls. They ship a `single_disease_confidence` field (low/medium/high) plus a reasoning string to help you filter. Treat `cell_line`, `tissue`, and `perturbation` as free-text-ish and expect to normalize. One naming quirk: `NRX` accessions are CZI CELLxGENE datasets that were never uploaded to SRA, so they don't follow SRX conventions.

Organized one h5ad per sample (SRX accession), with parquet tables at the sample level. Obviously filter to human first — 27 organisms means most of it isn't relevant.

## Virtual Cell Challenge 2025

~300,000 cells, 300 target genes, **CRISPRi knockdown in H1 human embryonic stem cells**, at `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/`.

Small compared to the other two, and the only one where the assay exactly matches what you're being asked to predict. Its limitation is that it's a **single cell line**, so it teaches you nothing about cross-context transfer — the entire point of this year's challenge. Its enormous advantage is that all three splits are now public *with* ground truth:

```
virtual-cell-challenge/2025/
├── train/
│   ├── adata_Training.h5ad
│   └── pert_counts_Training.csv
├── validation/
└── test/
```

This is what you build your local scoring harness on. You have real perturbed measurements, so you can run `cell-eval2 --preset vcc2026` against them and get honest metric numbers without spending leaderboard submissions. Note the obs schema follows the convention baked into cell-eval's defaults: `target_gene` for the perturbation column, with `non-targeting` marking controls. Worth confirming the gene space matches 2026's 18,533 when you load it — different year, possibly different reference.

## How this composes

The shape of a serious attempt is roughly: pretrain basal-state representation on scBaseCount, learn context-dependent response structure from Tahoe's 50-cell-line × drug grid, learn the CRISPRi-specific response from the 2025 H1 data plus other public Perturb-seq (Replogle is the obvious addition), and validate locally against 2025 ground truth before touching the leaderboard.

The hard unsolved bit, and where the challenge is actually won, is that no single public dataset gives you *genetic* perturbations across *many* cell lines. Tahoe has the contexts but the wrong perturbation type; Perturb-seq datasets have the right type but few contexts. Bridging that gap is the modeling problem.

Want me to write the loaders — GCP auth with requester-pays, a metadata-first query layer over the parquet tables, and a gene-space harmonizer that maps all three onto the challenge's 18,533 symbols? That's the unglamorous part that otherwise eats your first week.