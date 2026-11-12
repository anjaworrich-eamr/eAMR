# =====================================================================
#  FigS7.R - Supplementary Fig. S7: number of Rank I ARG subtypes per resistance class
# ---------------------------------------------------------------------
#  Output: output/FigS7/
# =====================================================================
source("scripts/00_functions.R")
OUT <- out_folder("FigS7")
SCALE <- 1.7
theme_this <- theme_bw(base_size = 12, base_family = FONT) +
  theme(plot.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5), strip.text = element_text(face = "bold", size = 12),
        axis.text = element_text(colour = "black", size = 11), axis.text.x = element_text(face = "bold"),
        axis.title = element_text(face = "bold", colour = "black"), legend.title = element_blank(),
        legend.text = element_text(size = 11), legend.background = element_blank(), legend.key = element_blank())
risk <- read.delim(RANK_FILE, check.names = FALSE, stringsAsFactors = FALSE)
risk <- risk[risk$rank == "I", ]
counts <- risk %>% distinct(Subtype, Type) %>% mutate(Class = norm_type(Type)) %>%
  count(Class, name = "n") %>% mutate(Class = fct_reorder(Class, -n))
save_table(counts, OUT, "FigS7_rankI_subtypes_per_class")
p <- ggplot(counts, aes(Class, n)) +
  geom_col(fill = "grey70", colour = "grey30", width = 0.7, linewidth = 0.4) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.8, family = FONT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of Rank I ARG subtypes") +
  theme_this + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
save_fig(patchwork::wrap_elements(full = p), OUT, "FigS7", w = 90 * SCALE / 25.4, h = 66 * SCALE / 25.4)
