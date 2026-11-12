# =====================================================================
#  FigS5.R - Supplementary Fig. S5: genus composition of EV-associated vs bulk-soil DNA in each land-use type
# ---------------------------------------------------------------------
#  Mean relative abundance of each genus (> 0.001% in either fraction) in bulk soil and EVs; Spearman correlation on
#  untransformed means; axes log10 with a pseudocount of 0.001%. Output: output/FigS5/
# =====================================================================
source("scripts/00_functions.R")
use("ggh4x")
OUT <- out_folder("FigS5")
MIN_ABUND <- 0.001; PSEUDO <- 0.001
meta <- read.csv(META_FILE, stringsAsFactors = FALSE) %>% mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)))
gen  <- read.delim(TAXA_FILES[["Genus"]], row.names = 1, check.names = FALSE)[, meta$SampleID]
dat <- map_dfr(names(COL_HABITAT), function(hab) {
  ev_s <- meta$SampleID[meta$Habitat == hab & meta$SampleType == "EV"]
  so_s <- meta$SampleID[meta$Habitat == hab & meta$SampleType == "Soil"]
  tibble(Habitat = hab, Genus = rownames(gen), EV = rowMeans(gen[, ev_s, drop = FALSE]), Soil = rowMeans(gen[, so_s, drop = FALSE])) %>%
    filter(EV > MIN_ABUND | Soil > MIN_ABUND)
}) %>% mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)))
stats <- dat %>% group_by(Habitat) %>%
  summarise(n_genera = n(), rho = unname(cor.test(EV, Soil, method = "spearman", exact = FALSE)$estimate),
            p = cor.test(EV, Soil, method = "spearman", exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(label = sprintf("atop(italic(rho) == '%.2f', %s)", rho, ifelse(p < 0.001, "italic(p) < 0.001", sprintf("italic(p) == '%.3f'", p))))
save_table(select(stats, -label), OUT, "FigS5_correlation")
brks <- c(0.001, 0.01, 0.1, 1, 10, 100); lims <- c(-3.15, 2.2)
p <- ggplot(dat, aes(log10(Soil + PSEUDO), log10(EV + PSEUDO))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_point(aes(colour = Habitat), size = 1.3, alpha = 0.35, stroke = 0) +
  geom_text(data = stats, aes(x = lims[1] + 0.15, y = lims[2] - 0.1, label = label), parse = TRUE, hjust = 0, vjust = 1, size = 3.9, family = FONT) +
  facet_wrap2(~ Habitat, nrow = 1, strip = strip_themed(background_x = elem_list_rect(fill = alpha(COL_HABITAT, 0.35), colour = "black"))) +
  scale_colour_manual(values = COL_HABITAT) +
  scale_x_continuous(breaks = log10(brks), labels = brks) + scale_y_continuous(breaks = log10(brks), labels = brks) +
  coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
  labs(x = "Mean relative abundance in bulk soil (%)", y = "Mean relative abundance in EVs (%)") +
  theme_bw(base_size = 12, base_family = FONT) +
  theme(plot.title = element_blank(), panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        panel.spacing.x = unit(0.9, "lines"), axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black", size = 11),
        axis.text.x = element_text(face = "bold"), axis.title = element_text(face = "bold", colour = "black"),
        strip.text.x = element_text(face = "bold", size = 12), legend.position = "none")
save_fig(p, OUT, "FigS5", w = 180 * 1.8 / 25.4, h = 68 * 1.8 / 25.4)
