Short answer: yes, decisively. Replogle and Nadig are better than Tahoe or scBaseCount for the core modeling problem, and together with Jiang they close most of the gap I described in your `data.md`. That note says no public dataset gives you genetic perturbations across many cell lines — that's too pessimistic. Replogle and Nadig share a perturbation panel by design, and Jiang is a fully crossed six-line grid. Here's the revised picture.

## Tier 1: Replogle 2022 + Nadig 2025 (treat as one dataset)

This is your primary training data. Five datasets, all **CRISPRi with dCas9-KRAB on 10x Chromium 3′ v3** — the same perturbation modality as the challenge and the same platform family:

| dataset | cell line (tissue) | perturbations | size |
|---|---|---|---|
| K562-GenomeWide | K562 (CML, blood) | 9,866 — all expressed genes | 61.3 GB |
| K562-Essential | K562 | 2,057 essential | 9.9 GB |
| RPE1-Essential | RPE1 (retinal epithelium, **non-cancer**) | 2,393 essential | 8.1 GB |
| Jurkat-Essential | Jurkat (T-cell leukemia) | 2,393 essential | 8.7 GB |
| HepG2-Essential | HepG2 (liver carcinoma) | 2,393 essential | 5.2 GB |

The reason this is the good stuff isn't the scale — it's that **the four Essential datasets target the same panel** (DepMap 20Q1 common essentials), from the same sgRNA library (`dJR092`), processed through the same Cell Ranger pipeline. Nadig ran Jurkat and HepG2 explicitly to extend Replogle's two lines to four for cross-context comparison. That gives you roughly 2,000 perturbations × 4 cell lines from 4 different tissues of origin, fully crossed. That *is* the cross-context genetic perturbation grid, and it's the closest structural analogue to what you're being scored on.

K562-GenomeWide plays a different role: 9,866 targets covering all expressed genes means that for almost any gene in the challenge's panel, you have a measurement of what knocking it down does *somewhere*. Arc's own dataset page flags this ("overlaps with most of the genes used in the Arc VCC datasets"). Think of it as the marginal prior over perturbations, with the four-line grid supplying the context-dependent correction.

**The critical caveat, and it's a big one.** These screens are shallow per perturbation: median 178 cells for K562-GenomeWide, 45 for HepG2. The challenge asks for 400 cells per perturbation. Nadig's central finding is what this does to your labels — the median correlation of raw log2 fold-change estimates between two K562 replicate experiments was **0.16**, but 0.90 once sampling variance was modeled. So if you naively train on per-perturbation pseudobulk LFCs from these files, the large majority of what you're fitting is estimation noise, and you will not find that out from your training loss. Their worst-case example is `EIF4A3`, measured in 7 cells, whose apparent transcriptome-wide effect was almost entirely sampling variance.

Two consequences. Filter or weight by cells-per-perturbation. And strongly consider taking the shortcut: Nadig's [figshare release](https://doi.org/10.6084/m9.figshare.29498366) ships precomputed perturbation × gene matrices of log2FC, **standard errors**, and p-values for all five datasets. Since four of the six challenge metrics are computed from Wilcoxon DE, training against properly-weighted LFCs is arguably more aligned with the objective than training on raw counts — and it lets you skip ~93 GB of downloads while you prototype.

Minor friction worth knowing: the CRISPRi effector differs across experiments (Jurkat used Zim3-dCas9, HepG2 used dCas9-BFP-KRAB), so knockdown strength isn't constant across contexts and is partly confounded with cell line. Reference is 10x GRCh38 2020-A, so gene-space harmonization applies here too. Raw data is on SRA `PRJNA1100571`, processed on GEO `GSE264667`.

## Tier 2: Jiang 2025 — structurally the closest thing to the challenge

[Systematic reconstruction of molecular pathway signatures](https://doi.org/10.1038/s41556-025-01622-z), Satija lab, *Nature Cell Biology*. Note the design sentence: Perturb-seq in "six different cancer cell lines from different tissues of origin" — A549 (lung), MCF7 (breast), HT29 (colon), HAP1, BxPC3 (pancreas), K562. That is nearly verbatim Arc's description of the 2026 challenge design. CRISPRi again (dCas9-KRAB-MeCP2), 2.6M cells, >1,500 perturbations.

The strength is that it's **fully crossed and balanced**: each of five pathway libraries was infected into all six lines, giving 30 samples with identical guides across contexts. For learning or validating "how does the same knockdown differ by context," that's cleaner than anything else available.

Two limits. The panel is narrow — 44 to 61 known regulators per pathway, so on the order of 250 unique target genes, all signaling regulators rather than a broad functional sweep. And every sample is also under cytokine stimulation (IFN-β, IFN-γ, TGF-β, TNF-α, insulin) for 24 hours, which is an extra experimental axis the challenge doesn't have; you'd need to handle it as a covariate rather than pretend it isn't there.

Also: this used **Parse Biosciences Evercode** combinatorial indexing, not 10x droplets. That's a larger platform gap than anything in Tier 1, with different capture efficiency and depth characteristics. Data is Seurat objects on [Zenodo](https://doi.org/10.5281/zenodo.14518762), so there's R-to-AnnData conversion work. Their **Mixscale** framework is worth reading regardless — it addresses per-cell variation in CRISPRi knockdown efficiency, which is a real source of label noise in every dataset on this list.

I'd use Jiang as a held-out validation grid for cross-context transfer rather than as bulk training data.

## Tier 3: Srivatsan 2020 (sci-Plex) — skip

188 compounds across three cancer lines (A549, K562, MCF7), ~650k cells, sci-RNA-seq with nuclear hashing. It's chemical perturbation, which puts it in the same role as Tahoe-100M — and Tahoe strictly dominates it: 100M cells vs 650k, 50 lines vs 3, 1,100 drugs vs 188. The dose-response design is elegant and the paper is a landmark, but the marginal information over Tahoe is close to zero for your purposes.

## Tier 4: McFaline-Figueroa 2024 (sci-Plex-GxE) — read it, don't train on it

Combined genetic *and* chemical screening in glioblastoma lines. Conceptually adjacent in an interesting way — gene × environment interaction is a cousin of the context-dependence you're modeling — and it does contain genetic perturbations. But it's confined to glioblastoma, the design centers on combinations, and disentangling single-perturbation effects from combinatorial ones is extra work for a narrow slice of context diversity. Worth reading for how they model interaction terms.

## Revised ranking

1. **Replogle Essential + Nadig Jurkat/HepG2** — the 4-context shared-panel CRISPRi grid. Your cross-context signal.
2. **Replogle K562-GenomeWide** — per-gene knockdown prior covering nearly the whole target space.
3. **VCC 2025 H1** — local scoring harness with real ground truth and exact metric match. Unchanged.
4. **Jiang 2025** — 6-context balanced validation grid.
5. **Tahoe-100M** — demoted. Still the widest context diversity by far, and useful for pretraining a context encoder, but it's chemical perturbation and it's no longer your only route to context-dependence.
6. **scBaseCount** — basal-state pretraining only.
7. sci-Plex, sci-Plex-GxE — skip for now.

Counting distinct cell lines with CRISPRi Perturb-seq available: K562, RPE1, Jurkat, HepG2, H1, plus Jiang's A549, MCF7, HT29, HAP1, BxPC3. About ten contexts with the right perturbation modality. That's a small N for learning context transfer — it argues for hierarchical or strongly-regularized approaches over anything with a lot of free parameters per context — but it's a real N, not zero.

One thing to check early, because it would reorder this list: whether the 2026 challenge's ~300 target genes are enriched for DepMap common essentials. If they are, the four Essential datasets are directly on-panel and Tier 1 gets even stronger. If the panel is broader, K562-GenomeWide is doing more of the work. You can answer that in about ten lines of code once you've downloaded `pert_counts.csv` and the DepMap common-essential list.

Want me to append this to `data.md` and write that overlap check?