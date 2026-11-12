# =====================================================================
#  Fig3.R - Figure 3: high-risk (Rank I) ARGs
# ---------------------------------------------------------------------
#  a  Number of Rank I ARG subtypes per sample; two-sided paired Wilcoxon signed-rank test within each land-use type
#  b  Prevalence and mean abundance of the Rank I subtypes detected in >= 10% of all samples
#  Rank I = health-risk ranking of Zhang et al. (2021), Nat. Commun. 12, 4765 (annotation in the ARGs-OAP risk table).
#  Also writes Supplementary Data 8 (prevalence and abundance of all Rank I subtypes per land use and fraction).
#  Output: output/Fig3/
# =====================================================================
source("scripts/00_functions.R")
use(c("ggh4x", "patchwork"))
OUT  <- out_folder("Fig3")
meta <- read.csv(META_FILE, stringsAsFactors = FALSE)

# ---------------- a. Number of Rank I ARG subtypes ----------------
p_a <- local({
  ALPHA_BOX <- c(Soil = 0.30, EV = 0); POINT_SPREAD <- 0.075
  ylim <- c(0, 25); y_brk <- seq(0, 25, 5)
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), panel.spacing.x = unit(0.8, "lines"),
          strip.text = element_text(face = "bold", size = 12), axis.text = element_text(colour = "black", size = 11),
          axis.text.x = element_text(face = "bold"), axis.title = element_text(face = "bold", colour = "black"),
          legend.title = element_blank(), legend.text = element_text(size = 11), legend.background = element_blank(),
          legend.key = element_blank(), legend.position = "right", legend.justification = c(0, 1))
  pairs <- tibble(EV = meta$SampleID[meta$SampleType == "EV"]) %>%
    mutate(Soil = vapply(EV, ev_to_soil, character(1)), Pair = paste0("P", row_number()),
           Habitat = meta$Habitat[match(EV, meta$SampleID)])
  risk <- read.delim(RANK_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  risk <- risk[risk$rank == "I", ]
  n <- colSums(rowsum(as.matrix(risk[, meta$SampleID]), risk$Subtype) > 0)   # variants summed per subtype; detected = TPM > 0
  df <- pairs %>% pivot_longer(c(Soil, EV), names_to = "SampleType", values_to = "SampleID") %>%
    left_join(tibble(SampleID = names(n), n_detected = unname(n)), by = "SampleID") %>%
    mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)), SampleType = factor(SampleType, levels = names(SHP_TYPE)))
  save_table(df %>% select(Habitat, Pair, SampleType, SampleID, n_detected) %>% arrange(Habitat, Pair, SampleType),
             OUT, "Fig3a_rankI_subtype_richness_values")
  fmt_p <- function(p) ifelse(is.na(p), "'n.a.'", ifelse(p < 0.001, "italic(p) < 0.001", sprintf("italic(p) == '%.2f'", p)))
  tests <- df %>% select(Habitat, Pair, SampleType, n_detected) %>%
    pivot_wider(names_from = SampleType, values_from = n_detected) %>% group_by(Habitat) %>%
    summarise(n_pairs = n(), median_Soil = median(Soil), median_EV = median(EV), n_EV_higher = sum(EV > Soil),
              V = unname(suppressWarnings(wilcox.test(EV, Soil, paired = TRUE))$statistic),
              p = suppressWarnings(wilcox.test(EV, Soil, paired = TRUE))$p.value, .groups = "drop") %>%
    mutate(label = fmt_p(p))
  save_table(tests, OUT, "Fig3a_rankI_subtype_paired_tests")
  brk <- tests %>% left_join(df %>% group_by(Habitat) %>% summarise(top = max(n_detected), .groups = "drop"), by = "Habitat") %>%
    mutate(y = top + ylim[2] * 0.08)
  df <- df %>% group_by(Habitat, SampleType, n_detected) %>%       # spread points with identical counts horizontally
    mutate(x = as.numeric(SampleType) + (row_number() - (n() + 1) / 2) * POINT_SPREAD) %>% ungroup()
  ggplot(df, aes(SampleType, n_detected)) +
    geom_boxplot(aes(fill = Habitat, alpha = SampleType, colour = Habitat), width = 0.6, outlier.shape = NA, linewidth = 0.5) +
    geom_point(aes(x = x, colour = Habitat, shape = SampleType), fill = "white", size = 2.3, stroke = 0.9) +
    geom_segment(data = brk, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = brk, aes(x = 1, xend = 1, y = y, yend = y - ylim[2] * 0.02), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = brk, aes(x = 2, xend = 2, y = y, yend = y - ylim[2] * 0.02), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = brk, aes(x = 1.5, y = y + ylim[2] * 0.015, label = label), parse = TRUE,
              vjust = 0, size = 3.6, family = FONT, inherit.aes = FALSE) +
    scale_colour_manual(values = COL_HABITAT, guide = "none") +
    scale_fill_manual(values = COL_HABITAT, guide = "none") +
    scale_alpha_manual(values = ALPHA_BOX, guide = "none") +
    scale_shape_manual(values = SHP_TYPE, guide = "none") +
    scale_y_continuous(limits = ylim, breaks = y_brk, expand = expansion(mult = c(0.03, 0))) +
    facet_wrap2(~ Habitat, nrow = 1, strip = strip_themed(background_x = elem_list_rect(fill = alpha(COL_HABITAT, 0.35), colour = NA))) +
    labs(x = NULL, y = "Number of Rank I ARG subtypes") +
    theme_this
})

# ---------------- b. Prevalence and abundance of Rank I subtypes; Supplementary Data 8 ----------------
stat <- rankI_stats()
save_table(stat %>% select(Class, Gene, Subtype, Habitat, SampleType, n, n_detected, prevalence, mean_TPM, mean_TPM_detected, prev_all) %>%
             arrange(Class, desc(prev_all), Habitat, SampleType), OUT, "Supplementary_Data_8_rankI_prevalence_abundance")
p_b <- rankI_bubble(filter(stat, prev_all >= 10), abund_max = ceiling(max(stat$log_abund)))

# ---------------- Assemble (180 x 110 mm at print size; canvas x 1.8) ----------------
SCALE <- 1.8; PRINT_W <- 180; FIG_H <- 110; A_FILL <- 0.72
col_a <- tag_panel(p_a, "a") / plot_spacer() + plot_layout(heights = c(A_FILL, 1 - A_FILL))
fig   <- wrap_elements(full = col_a) + tag_panel(p_b, "b") + plot_layout(widths = c(62, 118))
save_fig(fig, OUT, "Fig3", w = PRINT_W * SCALE / 25.4, h = FIG_H * SCALE / 25.4)
message("Figure 3 written to ", OUT)
