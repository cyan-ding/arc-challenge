Short answer: yes, decisively. Replogle and Nadig are better than Tahoe or scBaseCount for the core modeling problem, and together with Jiang they close most of the gap I described in your `data.md`. That note says no public dataset gives you genetic perturbations across many cell lines — that's too pessimistic. Replogle and Nadig share a perturbation panel by design, Jiang is a fully crossed six-line grid, and Feng 2026 nominally adds 34 more iPSC lines — though its depth per context makes it a different kind of asset than the headline numbers suggest. Here's the revised picture.

## Tier 1: Replogle 2022 + Nadig 2025 (treat as one dataset)

This is your primary training data. Five datasets, all **CRISPRi with dCas9-KRAB on 10x Chromium 3′ v3** — the same perturbation modality as the challenge and the same platform family:


| dataset          | cell line (tissue)                        | perturbations               | size    |
| ---------------- | ----------------------------------------- | --------------------------- | ------- |
| K562-GenomeWide  | K562 (CML, blood)                         | 9,866 — all expressed genes | 61.3 GB |
| K562-Essential   | K562                                      | 2,057 essential             | 9.9 GB  |
| RPE1-Essential   | RPE1 (retinal epithelium, **non-cancer**) | 2,393 essential             | 8.1 GB  |
| Jurkat-Essential | Jurkat (T-cell leukemia)                  | 2,393 essential             | 8.7 GB  |
| HepG2-Essential  | HepG2 (liver carcinoma)                   | 2,393 essential             | 5.2 GB  |


The reason this is the good stuff isn't the scale — it's that **the four Essential datasets target the same panel** (DepMap 20Q1 common essentials), from the same sgRNA library (`dJR092`), processed through the same Cell Ranger pipeline. Nadig ran Jurkat and HepG2 explicitly to extend Replogle's two lines to four for cross-context comparison. That gives you roughly 2,000 perturbations × 4 cell lines from 4 different tissues of origin, fully crossed. That *is* the cross-context genetic perturbation grid, and it's the closest structural analogue to what you're being scored on.

K562-GenomeWide plays a different role: 9,866 targets covering all expressed genes means that for almost any gene in the challenge's panel, you have a measurement of what knocking it down does *somewhere*. Arc's own dataset page flags this ("overlaps with most of the genes used in the Arc VCC datasets"). Think of it as the marginal prior over perturbations, with the four-line grid supplying the context-dependent correction.

**The critical caveat, and it's a big one.** These screens are shallow per perturbation: median 178 cells for K562-GenomeWide, 45 for HepG2. The challenge asks for 400 cells per perturbation. Nadig's central finding is what this does to your labels — the median correlation of raw log2 fold-change estimates between two K562 replicate experiments was **0.16**, but 0.90 once sampling variance was modeled. So if you naively train on per-perturbation pseudobulk LFCs from these files, the large majority of what you're fitting is estimation noise, and you will not find that out from your training loss. Their worst-case example is `EIF4A3`, measured in 7 cells, whose apparent transcriptome-wide effect was almost entirely sampling variance.

Two consequences. Filter or weight by cells-per-perturbation. And strongly consider taking the shortcut: Nadig's [figshare release](https://doi.org/10.6084/m9.figshare.29498366) ships precomputed perturbation × gene matrices of log2FC, **standard errors**, and p-values for all five datasets. Since four of the six challenge metrics are computed from Wilcoxon DE, training against properly-weighted LFCs is arguably more aligned with the objective than training on raw counts — and it lets you skip ~93 GB of downloads while you prototype.

Minor friction worth knowing: the CRISPRi effector differs across experiments (Jurkat used Zim3-dCas9, HepG2 used dCas9-BFP-KRAB), so knockdown strength isn't constant across contexts and is partly confounded with cell line. Reference is 10x GRCh38 2020-A, so gene-space harmonization applies here too. Raw data is on SRA `PRJNA1100571`, processed on GEO `GSE264667`.

## Tier 2: Feng 2026 — genome-scale CRISPRi in 34 iPSC lines, but read the depth numbers first

[A genome-scale single-cell CRISPRi map of trans gene regulation across human pluripotent stem cell lines](https://doi.org/10.1016/j.xgen.2025.101076), Parts / Stegle / Velten labs at the Wellcome Sanger Institute, *Cell Genomics* 6(2):101076. Preprint on [bioRxiv](https://doi.org/10.1101/2024.11.28.625833).

The headline design: **7,226 genes knocked down by CRISPRi across 34 iPSC lines from 26 healthy donors**, using ~20,000 guides (three per gene from the Dolcetto library, plus 40 non-targeting), with ~2M cells sequenced. On paper that quadruples the number of cellular contexts available with the right perturbation modality, and it's why the PRiMeFlow authors put it in their pretraining atlas.

Then read the post-QC numbers, which change the picture completely. After demultiplexing lines by genotype and assigning guides, they retained **219,206 cells carrying both a source line and a targeting guide** — a median of **25 cells per gene** and 8 per guide. Spread over 6,673 usable targets and 34 lines, that's on the order of **one cell per (target, line) pair**. You cannot compute a per-line pseudobulk for a perturbation out of this. The 34-line design is a population-genetics random-effects design, built to ask "how much does perturbation response vary with genetic background," not a covariate-transfer grid you can train on the way you can with the Replogle/Nadig four-line panel.

So its value is two specific things, both real:

**A pluripotent-context perturbation prior.** They publish an aggregate log-fold-change map of **6,673 targets × 6,471 expressed genes**, pooled across all lines. That's the same role K562-GenomeWide plays — broad coverage of what knocking down gene X does — but in a non-cancer, pluripotent context, which is a genuinely complementary axis rather than a redundant one. The target list also overlaps your Tier 1 panel heavily by construction: selection was 2,264 genes with iPSC growth defects plus 4,594 with cancer-line growth defects (1,725 in both), plus 2,093 highly expressed in iPSCs.

**~500,000 control cells spanning 34 lines and 26 donors.** This is the underrated half. The 2026 task hands you the unperturbed profile of a cell line you've never seen and asks for its perturbation response, so your basal-state encoder needs a good prior over what unperturbed profiles look like. Half a million non-targeting and no-guide cells across 26 genetic backgrounds is a better-matched basal panel than anything you'll pull out of scBaseCount, because the assay, the pipeline, and the reference are all internally consistent.

Two further cautions. Their own central QC observation is that **CRISPRi-induced transcriptional change was small relative to global heterogeneity** — technical factors, cell quality, and cell line dominated the variance at 3–6 days post-infection. And at a median of 25 cells per gene, the Nadig sampling-variance warning above applies here harder than to anything else on this list; 25 is shallower than K562-GenomeWide's 178 and even HepG2's 45. Do not train on raw per-perturbation LFCs from this without weighting by standard error.

Data is post-processed count data on [figshare](https://doi.org/10.6084/m9.figshare.27989294): RNA UMI counts, guide UMI counts, and cell metadata, split into three screens (fitness genes, non-fitness genes, targeted). Convenient release — you don't need to reprocess from raw.

## Tier 3: Jiang 2025 — structurally the closest thing to the challenge

[Systematic reconstruction of molecular pathway signatures](https://doi.org/10.1038/s41556-025-01622-z), Satija lab, *Nature Cell Biology*. Note the design sentence: Perturb-seq in "six different cancer cell lines from different tissues of origin" — A549 (lung), MCF7 (breast), HT29 (colon), HAP1, BxPC3 (pancreas), K562. That is nearly verbatim Arc's description of the 2026 challenge design. CRISPRi again (dCas9-KRAB-MeCP2), 2.6M cells, >1,500 perturbations.

The strength is that it's **fully crossed and balanced**: each of five pathway libraries was infected into all six lines, giving 30 samples with identical guides across contexts. For learning or validating "how does the same knockdown differ by context," that's cleaner than anything else available.

Two limits. The panel is narrow — 44 to 61 known regulators per pathway, so on the order of 250 unique target genes, all signaling regulators rather than a broad functional sweep. And every sample is also under cytokine stimulation (IFN-β, IFN-γ, TGF-β, TNF-α, insulin) for 24 hours, which is an extra experimental axis the challenge doesn't have; you'd need to handle it as a covariate rather than pretend it isn't there.

Also: this used **Parse Biosciences Evercode** combinatorial indexing, not 10x droplets. That's a larger platform gap than anything in Tier 1, with different capture efficiency and depth characteristics. Data is Seurat objects on [Zenodo](https://doi.org/10.5281/zenodo.14518762), so there's R-to-AnnData conversion work. Their **Mixscale** framework is worth reading regardless — it addresses per-cell variation in CRISPRi knockdown efficiency, which is a real source of label noise in every dataset on this list.

I'd use Jiang as a held-out validation grid for cross-context transfer rather than as bulk training data.

## Tier 4: Srivatsan 2020 (sci-Plex) — skip

188 compounds across three cancer lines (A549, K562, MCF7), ~650k cells, sci-RNA-seq with nuclear hashing. It's chemical perturbation, which puts it in the same role as Tahoe-100M — and Tahoe strictly dominates it: 100M cells vs 650k, 50 lines vs 3, 1,100 drugs vs 188. The dose-response design is elegant and the paper is a landmark, but the marginal information over Tahoe is close to zero for your purposes.

## Tier 5: McFaline-Figueroa 2024 (sci-Plex-GxE) — read it, don't train on it

Combined genetic *and* chemical screening in glioblastoma lines. Conceptually adjacent in an interesting way — gene × environment interaction is a cousin of the context-dependence you're modeling — and it does contain genetic perturbations. But it's confined to glioblastoma, the design centers on combinations, and disentangling single-perturbation effects from combinatorial ones is extra work for a narrow slice of context diversity. Worth reading for how they model interaction terms.

## Revised ranking

1. **Replogle Essential + Nadig Jurkat/HepG2** — the 4-context shared-panel CRISPRi grid. Your cross-context signal.
2. **Replogle K562-GenomeWide** — per-gene knockdown prior covering nearly the whole target space.
3. **VCC 2025 H1** — local scoring harness with real ground truth and exact metric match. Unchanged.
4. **Feng 2026** — pluripotent-context perturbation prior over 6,673 targets, plus a 500k-cell basal panel across 26 donors. Not a context-transfer grid; don't use it as one.
5. **Jiang 2025** — 6-context balanced validation grid.
6. **Tahoe-100M** — demoted. Still the widest context diversity by far, and useful for pretraining a context encoder, but it's chemical perturbation and it's no longer your only route to context-dependence.
7. **scBaseCount** — basal-state pretraining only, and now partly superseded by Feng's control cells for pluripotent contexts specifically.
8. sci-Plex, sci-Plex-GxE — skip for now.

Counting distinct cell lines with CRISPRi Perturb-seq available: K562, RPE1, Jurkat, HepG2, H1, plus Jiang's A549, MCF7, HT29, HAP1, BxPC3, plus Feng's 34 iPSC lines. That's 44 nominal contexts, and it's tempting to read that as the cross-context problem being solved. It isn't, for two reasons. Only about ten of those are deep enough to estimate a per-context perturbation effect at all. And Feng's 34 are all the same cell *type* — the variation there is donor genetic background within pluripotency, not tissue of origin, which is the axis the challenge actually tests.

So for learning transfer across tissue-distinct contexts, the usable N is still about ten. That still argues for hierarchical or strongly-regularized approaches over anything with many free parameters per context. But Feng's donor axis is a clean, well-powered place to answer a different and useful question: how much response variation should you expect from genetic background alone, holding cell type fixed? That's an empirical floor on the context-transfer problem, and nothing else on this list gives it to you.

One thing to check early, because it would reorder this list: whether the 2026 challenge's ~300 target genes are enriched for DepMap common essentials. If they are, the four Essential datasets are directly on-panel and Tier 1 gets even stronger — and Feng comes along with it, since its target list was selected from cancer-line and iPSC fitness screens and so is enriched the same way. If the panel is broader, K562-GenomeWide and Feng's non-fitness screen are doing more of the work. You can answer that in about ten lines of code once you've downloaded `pert_counts.csv` and the DepMap common-essential list.

Want me to append this to `data.md` and write that overlap check?