# =====================================================================
#  Fig4.R - Figure 4: ARG-MGE co-localization on the assembled contigs
# ---------------------------------------------------------------------
#  a  Element type of the core mobility genes on the co-localized contigs; centre: number of co-localized contigs and their
#     share of all ARG-carrying contigs of the group (definition in 00_functions.R, read_coloc)
#  b  Share of ARG-carrying contigs that are co-localized, per sample; paired Wilcoxon signed-rank test within each land use
#  Also writes the paired test across all 27 pairs (text) and Supplementary Data 10 (ARG-carrying contigs with an MGE gene).
#  Output: output/Fig4/
# =====================================================================
source("scripts/00_functions.R")
use(c("ggh4x", "patchwork"))
OUT <- out_folder("Fig4")
co  <- read_coloc()
arg_ctg <- co$arg_ctg; mge <- co$mge

# ---------------- Supplementary Data 10 ----------------
save_table(arg_ctg %>% filter(C0) %>%
             left_join(co$arg_genes %>% group_by(molecule) %>%
                         summarise(ARGs = paste(unique(sub("^.*?__", "", arg_subtype, perl = TRUE)), collapse = "; "), .groups = "drop"),
                       by = "molecule") %>%
             left_join(mge %>% group_by(molecule) %>%
                         summarise(MGE_genes = paste(unique(gene_label), collapse = "; "), MGE_core = paste(unique(gene_label[core]), collapse = "; "),
                                   MGE_types = paste(unique(Type), collapse = "; "), min_gap_bp = min(gap_to_ARG), .groups = "drop"), by = "molecule") %>%
             arrange(SampleType, Habitat, sample) %>%
             select(sample, Habitat, SampleType, contig_id, length, core_mobility_gene = C1, ARGs, MGE_genes, MGE_core, MGE_types, min_gap_bp),
           OUT, "Supplementary_Data_10_ARG_MGE_contigs")

# ---------------- Counts per land use and fraction ----------------
counts <- arg_ctg %>%
  mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)), SampleType = factor(SampleType, levels = c("Soil", "EV"))) %>%
  group_by(Habitat, SampleType) %>%
  summarise(n_samples = n_distinct(sample), ARG_contigs = n(), coloc = sum(C1), coloc_pct = 100 * mean(C1), .groups = "drop") %>%
  bind_rows(arg_ctg %>% mutate(SampleType = factor(SampleType, levels = c("Soil", "EV"))) %>% group_by(SampleType) %>%
              summarise(Habitat = factor("All"), n_samples = n_distinct(sample), ARG_contigs = n(), coloc = sum(C1),
                        coloc_pct = 100 * mean(C1), .groups = "drop"))
save_table(counts, OUT, "Fig4_colocalized_contigs_counts")

# ---------------- Share of co-localized contigs per sample; paired test across all 27 pairs ----------------
share_full <- co$pairs %>% pivot_longer(c(Soil, EV), names_to = "SampleType", values_to = "sample") %>%
  left_join(arg_ctg %>% group_by(sample) %>% summarise(n_ARG_contigs = n(), share = 100 * mean(C1), .groups = "drop"), by = "sample")
save_table(share_full, OUT, "Fig4b_share_per_sample")   # samples without ARG-carrying contigs have no share and drop out of the test
sh <- share_full %>% select(Pair, Habitat, SampleType, share) %>% pivot_wider(names_from = SampleType, values_from = share) %>%
  filter(!is.na(EV), !is.na(Soil))
save_table(tibble(pairs_tested = nrow(sh), pairs_EV_higher = sum(sh$EV > sh$Soil), pairs_Soil_higher = sum(sh$EV < sh$Soil),
                  pairs_tied = sum(sh$EV == sh$Soil), median_share_EV = median(sh$EV), median_share_Soil = median(sh$Soil),
                  Wilcoxon_p = suppressWarnings(wilcox.test(sh$EV, sh$Soil, paired = TRUE, exact = FALSE))$p.value),
           OUT, "Fig4_paired_test_all_pairs")

# ---------------- a. Element types of the core mobility genes (donut charts) ----------------
p_a <- local({
  MGE_LEVELS <- c("Plasmid", "IS", "IGE", "Others")
  MGE_LABELS <- c(Plasmid = "Plasmid-associated genes", IS = "Insertion sequences (IS)", IGE = "Integrative elements (IGE)", Others = "Others")
  COL_MGE    <- c(Plasmid = "#5B4B8A", IS = "#E0A030", IGE = "#A0522D", Others = "grey78")
  LABEL_MIN_PCT <- 4
  grid6 <- expand_grid(Habitat = factor(names(COL_HABITAT), levels = names(COL_HABITAT)),
                       SampleType = factor(c("Soil", "EV"), levels = c("Soil", "EV")))
  tab_a <- mge %>% filter(molecule %in% arg_ctg$molecule[arg_ctg$C1], core) %>%
    left_join(select(co$meta, sample = SampleID, Habitat), by = "sample") %>%
    mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)), SampleType = factor(SampleType, levels = c("Soil", "EV"))) %>%
    count(Habitat, SampleType, Type, name = "n_genes") %>%
    group_by(Habitat, SampleType) %>% mutate(pct = 100 * n_genes / sum(n_genes)) %>% ungroup()
  save_table(tab_a, OUT, "Fig4a_MGE_types")
  centre <- grid6 %>%
    left_join(filter(counts, Habitat != "All") %>% mutate(Habitat = factor(as.character(Habitat), levels = names(COL_HABITAT))) %>%
                select(Habitat, SampleType, n = coloc, pct = coloc_pct), by = c("Habitat", "SampleType")) %>%
    mutate(n = replace_na(n, 0L), label = ifelse(is.na(pct), "n = 0", sprintf("n = %d\n(%.1f%%)", n, pct)))
  df <- tab_a %>% mutate(Type = factor(Type, levels = MGE_LEVELS)) %>% arrange(Habitat, SampleType, Type) %>%
    group_by(Habitat, SampleType) %>% mutate(ymax = pmin(cumsum(pct), 100), ymin = ymax - pct, ymid = (ymin + ymax) / 2) %>% ungroup()
  empty <- anti_join(grid6, distinct(df, Habitat, SampleType), by = c("Habitat", "SampleType"))
  ggplot() +
    { if (nrow(empty) > 0) geom_rect(data = empty, aes(xmin = 2, xmax = 3.3, ymin = 0, ymax = 100), fill = "grey95", colour = "grey80", linewidth = 0.3) } +
    geom_rect(data = df, aes(xmin = 2, xmax = 3.3, ymin = ymin, ymax = ymax, fill = Type), colour = "white", linewidth = 0.5) +
    geom_text(data = filter(df, pct >= LABEL_MIN_PCT), aes(x = 3.95, y = ymid, label = sprintf("%.0f%%", pct), colour = Type),
              size = 3.3, fontface = "bold", family = FONT, show.legend = FALSE) +
    geom_text(data = centre, aes(x = 0, y = 0, label = label), size = 3.4, family = FONT, lineheight = 0.95) +
    scale_fill_manual(values = COL_MGE, labels = MGE_LABELS, name = "MGE type", drop = FALSE) +
    scale_colour_manual(values = c(COL_MGE[1:3], Others = "grey45"), guide = "none") +
    scale_x_continuous(limits = c(0, 4.3), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    coord_polar(theta = "y") +
    facet_grid2(SampleType ~ Habitat, switch = "y",
                strip = strip_themed(background_x = elem_list_rect(fill = alpha(COL_HABITAT, 0.35), colour = NA),
                                     background_y = elem_list_rect(fill = c("grey88", "white"), colour = c(NA, "grey60")))) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(colour = NA))) +
    theme_void(base_size = 12, base_family = FONT) +
    theme(strip.text.x = element_text(face = "bold", size = 12, margin = margin(4, 0, 4, 0)),
          strip.text.y.left = element_text(face = "bold", size = 11.5, angle = 90, margin = margin(0, 4, 0, 4)),
          panel.spacing = unit(0.3, "lines"), legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10.5), plot.margin = margin(5, 5, 5, 5))
})

# ---------------- b. Share of co-localized contigs per sample ----------------
p_b <- local({
  ALPHA_BOX <- c(Soil = 0.30, EV = 0); POINT_SPREAD <- 0.075; Y_MAX <- 30
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), panel.spacing.x = unit(0.8, "lines"),
          strip.text = element_text(face = "bold", size = 12), axis.text = element_text(colour = "black", size = 11),
          axis.text.x = element_text(face = "bold"), axis.title = element_text(face = "bold", colour = "black"))
  df <- share_full %>% rename(v = share) %>% filter(!is.na(v)) %>%
    mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)), SampleType = factor(SampleType, levels = names(SHP_TYPE)))
  fmt_p <- function(p) ifelse(p < 0.001, "italic(p) < 0.001", sprintf("italic(p) == '%.2f'", p))
  tests <- df %>% select(Habitat, Pair, SampleType, v) %>% pivot_wider(names_from = SampleType, values_from = v) %>%
    filter(!is.na(Soil), !is.na(EV)) %>% group_by(Habitat) %>%
    summarise(n = n(), p = if (n() >= 3) suppressWarnings(wilcox.test(EV, Soil, paired = TRUE, exact = FALSE))$p.value else NA_real_,
              .groups = "drop") %>% filter(!is.na(p)) %>% mutate(label = fmt_p(p))
  save_table(select(tests, -label), OUT, "Fig4b_paired_tests")
  brk <- tests %>% left_join(df %>% group_by(Habitat) %>% summarise(top = max(v), .groups = "drop"), by = "Habitat") %>%
    mutate(y = pmin(top + Y_MAX * 0.08, Y_MAX * 0.9), y_tick = y - Y_MAX * 0.02, y_text = y + Y_MAX * 0.015)
  df <- df %>% group_by(Habitat, SampleType, bin = round(v / (Y_MAX * 0.02))) %>%
    mutate(x = as.numeric(SampleType) + (row_number() - (n() + 1) / 2) * POINT_SPREAD) %>% ungroup()
  ggplot(df, aes(SampleType, v)) +
    geom_boxplot(aes(fill = Habitat, alpha = SampleType, colour = Habitat), width = 0.5, outlier.shape = NA, linewidth = 0.5) +
    geom_point(aes(x = x, colour = Habitat, shape = SampleType), fill = "white", size = 2.3, stroke = 0.9) +
    geom_segment(data = brk, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = brk, aes(x = 1, xend = 1, y = y, yend = y_tick), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = brk, aes(x = 2, xend = 2, y = y, yend = y_tick), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = brk, aes(x = 1.5, y = y_text, label = label), parse = TRUE, vjust = 0, size = 3.6, family = FONT, inherit.aes = FALSE) +
    scale_colour_manual(values = COL_HABITAT, guide = "none") + scale_fill_manual(values = COL_HABITAT, guide = "none") +
    scale_alpha_manual(values = ALPHA_BOX, guide = "none") + scale_shape_manual(values = SHP_TYPE, guide = "none") +
    scale_y_continuous(limits = c(0, Y_MAX), breaks = seq(0, Y_MAX, 10), expand = expansion(mult = c(0.03, 0))) +
    facet_wrap2(~ Habitat, nrow = 1, strip = strip_themed(background_x = elem_list_rect(fill = alpha(COL_HABITAT, 0.35), colour = NA))) +
    labs(x = NULL, y = "ARG\u2013MGE co-localized contigs (%)") + theme_this
})

# ---------------- Assemble (180 x 85 mm at print size; canvas x 1.8) ----------------
SCALE <- 1.8; PRINT_W <- 180; FIG_H <- 85; B_FILL <- 0.82
col_b <- tag_panel(p_b, "b") / plot_spacer() + plot_layout(heights = c(B_FILL, 1 - B_FILL))
fig   <- tag_panel(p_a, "a") + wrap_elements(full = col_b) + plot_layout(widths = c(98, 82))
save_fig(fig, OUT, "Fig4", w = PRINT_W * SCALE / 25.4, h = FIG_H * SCALE / 25.4)
message("Figure 4 written to ", OUT)
