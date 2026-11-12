# =====================================================================
#  Fig1.R - Figure 1: bacterial community of bulk soil and EV-associated DNA (phylum level)
# ---------------------------------------------------------------------
#  a  Mean relative abundance of the 10 most abundant phyla per land-use type (dumbbell plot)
#  b  Shannon diversity; two-sided paired Wilcoxon signed-rank test, EV vs bulk soil, within each land-use type
#  c  PCoA (Bray-Curtis) with paired PERMANOVA
#  d  Phyla differing between paired EV and bulk-soil samples (log2 fold-change, bootstrapped 95% CI, heatmap per land use)
#  Also writes Supplementary Data 1-3 (paired EV vs bulk-soil comparison at phylum, genus and species level).
#  Output: output/Fig1/
# =====================================================================
source("scripts/00_functions.R")
use(c("ggh4x", "patchwork"))
OUT <- out_folder("Fig1")
meta <- read.csv(META_FILE, stringsAsFactors = FALSE)

# ---------------- a. Top 10 phyla: group means and individual samples ----------------
lab_taxon <- function(x) {                       # taxon names in bold italic, "sp. xxx" upright, Others bold
  esc <- gsub("'", "\\\\'", x)
  is_sp <- grepl(" sp\\.", esc)
  out <- sprintf("bolditalic('%s')", esc)
  out[is_sp] <- sprintf("bolditalic('%s')~bold('%s')", sub(" sp\\..*$", "", esc[is_sp]),
                        sub("^.*? (sp\\..*)$", "\\1", esc[is_sp], perl = TRUE))
  out[x == "Others"] <- "bold('Others')"
  parse(text = out)
}
p_a <- local({
  NOT_TOP <- "unclassified|unassigned|uncultured|norank|unknown|^Others$"
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey35", linewidth = 0.5),
          axis.text = element_text(colour = "black", size = 11), axis.text.x = element_text(face = "bold"),
          axis.title = element_text(face = "bold", colour = "black"), strip.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold"), legend.background = element_blank(), legend.key = element_blank())
  SAMPLE_ORDER <- c("UR_CHL", "UR_CNH", "UR_CZS", "UR_SJG", "UR_SNB", "UR_SRD", "UR_Park", "UR_Roadside", "UR_Sandy",
    "UR_CHL_EVs", "UR_CNH_EVs", "UR_CZS_EVs", "UR_SJG_EVs", "UR_SNB_EVs", "UR_SRD_EVs", "UR_Park_EVs", "UR_Roadside_EVs", "UR_Sandy_EVs",
    "FO_Forest", "FO_VTE", "FO_VTS", "FO_VTY", "FO_Forest_EVs", "FO_VTE_EVs", "FO_VTS_EVs", "FO_VTY_EVs",
    "FA_soil0_202410", "FA_soil0_202503", "FA_soil20_202410", "FA_soil20_202503", "FA_soil30_202410", "FA_soil30_202503", "FA_Farm",
    "FA_CK1_Soil", "FA_CK2_Soil", "FA_CK3_Soil", "FA_CK4_Soil", "FA_CK5_Soil", "FA_CK6_Soil", "FA_CK7_Soil",
    "FA_EV0_202410", "FA_EV0_202503", "FA_EV20_202410", "FA_EV20_202503", "FA_EV30_202410", "FA_EV30_202503", "FA_Farm_EVs",
    "FA_TCK1_EVs", "FA_TCK2_EVs", "FA_TCK3_EVs", "FA_TCK4_EVs", "FA_TCK5_EVs", "FA_TCK6_EVs", "FA_TCK7_EVs")
  tab <- read.delim(TAXA_FILES[["Phylum"]], row.names = 1, check.names = FALSE)[, intersect(SAMPLE_ORDER, meta$SampleID)]
  cand  <- rownames(tab)[!grepl(NOT_TOP, rownames(tab), ignore.case = TRUE)]
  top10 <- names(sort(rowMeans(tab[cand, ]), decreasing = TRUE))[1:10]
  top   <- rbind(tab[top10, ], Others = colSums(tab[!rownames(tab) %in% top10, , drop = FALSE]))
  y_levels <- rev(c(top10, "Others"))
  df <- top %>% rownames_to_column("Taxon") %>%
    pivot_longer(-Taxon, names_to = "SampleID", values_to = "Abundance") %>%
    left_join(meta, by = "SampleID") %>%
    mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)), SampleType = factor(SampleType, levels = names(SHP_TYPE)),
           Taxon = factor(Taxon, levels = c(top10, "Others")), y = match(as.character(Taxon), y_levels))
  mean_df <- df %>% group_by(Habitat, SampleType, Taxon, y) %>%
    summarise(Mean = mean(Abundance), SD = sd(Abundance), n = n(), .groups = "drop")
  seg_df <- mean_df %>% select(Habitat, SampleType, Taxon, Mean, y) %>%
    pivot_wider(names_from = SampleType, values_from = c(Mean, y))
  save_table(select(mean_df, -y), OUT, "Fig1a_phylum_group_means")

  x_max <- ceiling(max(df$Abundance) / 20) * 20
  xpos  <- function(frac) (frac * sqrt(x_max))^2
  # legend (filled = Soil, open = EV) in a separate row at the bottom of each panel
  leg_pts <- expand_grid(Habitat = factor(names(COL_HABITAT), levels = names(COL_HABITAT)),
                         SampleType = factor(names(SHP_TYPE), levels = names(SHP_TYPE))) %>%
    mutate(x = ifelse(SampleType == "Soil", xpos(0.56), xpos(0.79)), xt = ifelse(SampleType == "Soil", xpos(0.60), xpos(0.83)), y = 0)
  leg_box <- distinct(leg_pts, Habitat) %>% mutate(xmin = xpos(0.50), xmax = xpos(0.99), ymin = -0.38, ymax = 0.38)

  ggplot() +
    geom_rect(data = leg_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey95", colour = NA) +
    geom_point(data = df, aes(Abundance, y, colour = Habitat, shape = SampleType), position = position_identity(),
               fill = "white", alpha = 0.35, size = 1.3, stroke = 0.5) +
    geom_segment(data = seg_df, aes(x = Mean_Soil, xend = Mean_EV, y = y_Soil, yend = y_EV, colour = Habitat), linewidth = 0.7) +
    geom_point(data = filter(mean_df, SampleType == "Soil"), aes(Mean, y, colour = Habitat, shape = SampleType), size = 2.3, stroke = 1.1) +
    geom_point(data = filter(mean_df, SampleType == "EV"), aes(Mean, y, colour = Habitat, shape = SampleType),
               fill = NA, size = 3.4, stroke = 1.1) +
    geom_point(data = leg_pts, aes(x, y, colour = Habitat, shape = SampleType), fill = "white", size = 2.4, stroke = 1) +
    geom_text(data = leg_pts, aes(xt, y, label = SampleType), hjust = 0, size = 3.4, family = FONT) +
    scale_colour_manual(values = COL_HABITAT, guide = "none") +
    scale_shape_manual(values = SHP_TYPE, guide = "none") +
    scale_y_continuous(breaks = seq_along(y_levels), labels = lab_taxon(y_levels), expand = expansion(add = 0.5)) +
    scale_x_sqrt(breaks = c(0, 1, 5, 10, 25, 50, 100), limits = c(0, x_max), expand = expansion(mult = c(0.02, 0.02))) +
    facet_wrap2(~ Habitat, nrow = 1, strip = strip_themed(background_x = elem_list_rect(fill = alpha(COL_HABITAT, 0.35), colour = NA))) +
    labs(x = "Relative abundance (%)", y = NULL) +
    theme_this +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3), panel.spacing.x = unit(0.8, "lines"))
})

# ---------------- b. Shannon diversity ----------------
p_b <- shannon_panel("Phylum", OUT)

# ---------------- c. PCoA ----------------
block_c <- pcoa_panel(read.delim(TAXA_FILES[["Phylum"]], row.names = 1, check.names = FALSE), OUT, "Fig1c_phylum", perm_height = 0.24)

# ---------------- d. Phyla differing between paired EV and bulk-soil samples; Supplementary Data 1-3 ----------------
# Paired tests over all 27 pairs (ev_soil_test in 00_functions.R). If no phylum differs, the 10 most abundant are shown.
st <- ev_soil_setup()
for (k in seq_along(TAXA_FILES)) {
  lv <- ev_soil_test(read.delim(TAXA_FILES[[k]], row.names = 1, check.names = FALSE), st)
  save_table(lv$res %>% arrange(Direction, q), OUT, sprintf("Supplementary_Data_%d_%s", k, tolower(names(TAXA_FILES)[k])))
  if (k == 1) phylum <- lv
}
show <- ev_soil_shown(phylum, st, n_show = 15)
save_table(show, OUT, "Fig1d_phylum_plotted")
lab_italic <- function(x) {
  sp <- grepl(" sp\\.", x)
  ifelse(sp, sprintf("<i>%s</i> %s", sub(" sp\\..*$", "", x), sub("^.*? (sp\\..*)$", "\\1", x, perl = TRUE)), sprintf("<i>%s</i>", x))
}
p_d <- ev_soil_plot(show, lab_italic, heat_width = 1.5, xtitle_size = 10)

# ---------------- Assemble (180 mm wide at print size; canvas x 1.8) ----------------
SCALE <- 1.8; PRINT_W <- 180; ROW_H <- c(top = 60, bottom = 105)
p_a <- p_a + guides(x = guide_axis(check.overlap = TRUE)) + theme(plot.margin = margin(5, 14, 5, 5))
p_b <- p_b + theme(plot.margin = margin(5, 5, 5, 22))
top    <- tag_panel(p_a, "a") + tag_panel(p_b, "b") + plot_layout(widths = c(122, 58))
bottom <- tag_panel(block_c, "c") + tag_panel(p_d, "d") + plot_layout(widths = c(62, 118))
fig    <- wrap_elements(full = top) / wrap_elements(full = bottom) + plot_layout(heights = ROW_H)
save_fig(fig, OUT, "Fig1", w = PRINT_W * SCALE / 25.4, h = sum(ROW_H) * SCALE / 25.4)
message("Figure 1 written to ", OUT)
