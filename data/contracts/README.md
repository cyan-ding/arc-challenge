# data/contracts

The only files under `data/` that are tracked in git. Both come from the VCC controls bundle:

```bash
vcc datasets download controls -d ~/vcc
unzip -o -j ~/vcc/vcc_2026_controls.zip gene_names.csv pert_counts.csv -d data/contracts/
```

- `gene_names.csv` — one column with header `gene_name`, then 18,533 gene symbols (18,534 lines
  total), in the order every matrix in this repo uses. Every loader returns columns in this
  order. Never reorder. Note the `vcc` CLI help calls this file "headerless"; the actual download
  has the header, so skip line 1 when reading it by hand.
- `pert_counts.csv` — one column with header `target_gene`, then 300 target gene symbols (301
  lines total; the 2026 panel). No cell counts: every perturbation is exactly 400 cells.

When the final phase opens (Oct 22) the D/E/F bundle ships its own versions of these files. Keep
them side by side (`gene_names.final.csv`, `pert_counts.final.csv`) rather than overwriting, so
validation-phase code keeps running.
