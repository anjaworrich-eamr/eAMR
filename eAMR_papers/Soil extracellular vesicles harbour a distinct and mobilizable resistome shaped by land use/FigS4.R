# =====================================================================
#  FigS4.R - Supplementary Fig. S4: PCoA of bacterial community composition at genus (a) and species (b) level
# ---------------------------------------------------------------------
#  Bray-Curtis dissimilarity; paired PERMANOVA (see 00_functions.R). Output: output/FigS4/
# =====================================================================
source("scripts/00_functions.R")
use("patchwork")
OUT <- out_folder("FigS4")
SCALE <- 1.8
pc <- function(level) pcoa_panel(read.delim(TAXA_FILES[[level]], row.names = 1, check.names = FALSE), OUT, tolower(level))
fig <- tag_panel(pc("Genus"), "a") + tag_panel(pc("Species"), "b")
save_fig(fig, OUT, "FigS4", w = 180 * SCALE / 25.4, h = 105 * SCALE / 25.4)
