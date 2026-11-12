# =====================================================================
#  FigS3.R - Supplementary Fig. S3: Shannon diversity at genus (a) and species (b) level
# ---------------------------------------------------------------------
#  Two-sided paired Wilcoxon signed-rank test, EV vs bulk soil, within each land-use type. Output: output/FigS3/
# =====================================================================
source("scripts/00_functions.R")
use("patchwork")
OUT <- out_folder("FigS3")
SCALE <- 1.8
fig <- tag_panel(shannon_panel("Genus", OUT), "a") + tag_panel(shannon_panel("Species", OUT), "b")
save_fig(fig, OUT, "FigS3", w = 180 * SCALE / 25.4, h = 72 * SCALE / 25.4)
