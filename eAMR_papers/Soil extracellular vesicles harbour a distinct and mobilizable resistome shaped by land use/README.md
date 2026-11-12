# Soil extracellular vesicles harbour a distinct and mobilizable resistome shaped by land use

R code for the statistical analyses and figures of Huang et al., *Soil extracellular vesicles harbour a distinct and
mobilizable resistome shaped by land use*. One script per figure or table. The scripts start from the summary tables of the
bioinformatic workflow (read quality control, removal of viral sequences, Kraken2, ARGs-OAP, mobileOG-db, assembly,
ARG–MGE co-localization and MAGs; see Methods).

## Requirements

- R 4.5.2. Figures were produced on macOS. PDF files are written with the `quartz` device on macOS and with `cairo_pdf`
  elsewhere; the font is Arial.
- CRAN packages (versions used): tidyverse 2.0.0 (ggplot2 4.0.3, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0), patchwork 1.3.2,
  ggh4x 0.3.1, ggtext 0.2.0, ggrepel 0.9.8, cowplot 1.2.0, scales 1.4.0, vegan 2.7-2, permute 0.9-10, ape 5.8-1,
  gggenes 0.6.0. Missing packages are installed by the scripts.

## Input data

The input tables are not part of this repository; raw sequencing data are available as described in the Data availability
statement of the article. `00_prepare_inputs.R` reads

1. the summary tables of the bioinformatic workflow in `data/downstream/`:
   - `sample_rename.tsv`, `summary/00_sample_info.tsv`
   - `summary/01_virus_removal/virus_removal_per_sample.tsv`
   - `summary/02_kraken2_bacteria/bacteria_{phylum,genus,species}_relabund_pct.tsv`
   - `summary/03_argoap/argoap_{type,subtype,gene}_TPM.tsv`, `argoap_gene_TPM_with_risk_rank.tsv`
   - `summary/07_mge_reads/mobileOG_annotation_detected.tsv`, `tpm.id.tsv`
2. the sample information in `data/sample_info/`: `metadata_complete.csv`, `revision/qc_summary_all.tsv`,
   `Sampling_map/sampling_points.csv`,

and writes the tables used by the other scripts to `data/raw/`. In addition, `Fig4.R`, `FigS10_S11.R`, `Fig5.R` and
`FigS9.R` read `data/downstream/summary/04_arg_mge_contigs/` (`gggenes_data.tsv`, `arg_orfs_detail.tsv`),
`summary/05_mag/` (`mag_all_refined.tsv`, `mag_representatives.tsv`, `mag_arg_SARG_annotated.tsv`,
`mag_abundance_relative_abundance.tsv`, `tree/haining_bacteria_MAG.tree`) and `summary/07_mge_reads/tpm.element.tsv`.

## How to run

Open `soil-EV-resistome.Rproj` in RStudio (all paths are relative to the project folder) and run `scripts/run_all.R`, or run
`00_prepare_inputs.R` first and then any figure or table script. Each script writes to `output/<script name>/`: the figure
as PDF, PNG and TIFF and the values shown in it as CSV.

## Scripts

| Script | Content |
|---|---|
| `00_functions.R` | Settings and functions used by several scripts (colours, Soil–EV pairs, paired PERMANOVA, export) |
| `00_prepare_inputs.R` | Input tables in `data/raw/` |
| `Fig1.R` | Fig. 1; Supplementary Data 1–3 |
| `Fig2.R` | Fig. 2; Supplementary Data 4; resistome vs community composition (distance-based redundancy analysis) |
| `Fig3.R` | Fig. 3; Supplementary Data 8 |
| `Fig4.R` | Fig. 4; Supplementary Data 10 |
| `Fig5.R` | Fig. 5; Supplementary Table S4; Supplementary Data 11 and 12 |
| `FigS2.R` … `FigS9.R` | Supplementary Figs. S2–S9 (`FigS9.R` also Supplementary Data 9) |
| `FigS10_S11.R` | Supplementary Figs. S10 (EV) and S11 (bulk soil) |
| `FigS12_S13.R` | Supplementary Figs. S12 (ARG types) and S13 (ARG subtypes); Supplementary Data 5–7 |
| `TableS1.R` … `TableS3.R` | Supplementary Tables S1–S3 |
| `run_all.R` | Runs all scripts in order |

The Supplementary Data workbook was compiled from the CSV files named `Supplementary_Data_*.csv`. Statistical tests are
described at the top of each script and in the Methods of the article.

## Citation

Huang, H., Ghosh, D., Saraiva, J. P., Chatzinotas, A. & Worrich, A. Soil extracellular vesicles harbour a distinct and
mobilizable resistome shaped by land use (manuscript).
