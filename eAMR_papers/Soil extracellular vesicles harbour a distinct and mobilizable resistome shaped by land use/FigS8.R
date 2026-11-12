# =====================================================================
#  FigS8.R - Supplementary Fig. S8: prevalence and mean abundance of all Rank I ARG subtypes
# ---------------------------------------------------------------------
#  Same layout as Fig. 3b. Output: output/FigS8/
# =====================================================================
source("scripts/00_functions.R")
OUT <- out_folder("FigS8")
SCALE <- 1.7
stat <- rankI_stats()
fig  <- patchwork::wrap_elements(full = rankI_bubble(stat, abund_max = ceiling(max(stat$log_abund))))
save_fig(fig, OUT, "FigS8", w = 125 * SCALE / 25.4, h = 230 * SCALE / 25.4)
