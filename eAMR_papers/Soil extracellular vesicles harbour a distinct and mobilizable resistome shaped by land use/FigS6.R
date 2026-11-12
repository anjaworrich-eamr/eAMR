# =====================================================================
#  FigS6.R - Supplementary Fig. S6: number of detected ARG types (a) and ARGs (b)
# ---------------------------------------------------------------------
#  Detected = TPM > 0; two-sided paired Wilcoxon signed-rank test within each land-use type. Output: output/FigS6/
# =====================================================================
source("scripts/00_functions.R")
use("patchwork")
OUT <- out_folder("FigS6")
SCALE <- 1.8
fig <- tag_panel(richness_panel("Type", OUT), "a") + tag_panel(richness_panel("Gene", OUT), "b")
save_fig(fig, OUT, "FigS6", w = 180 * SCALE / 25.4, h = 64 * SCALE / 25.4)
