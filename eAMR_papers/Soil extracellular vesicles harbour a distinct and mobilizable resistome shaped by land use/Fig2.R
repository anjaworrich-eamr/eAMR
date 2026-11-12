# =====================================================================
#  Fig2.R - Figure 2: resistome of bulk soil and EV-associated DNA
# ---------------------------------------------------------------------
#  a  Number of detected ARG subtypes (TPM > 0); two-sided paired Wilcoxon signed-rank test within each land-use type
#  b  ARG subtypes detected only in bulk soil, in both fractions or only in EVs (all samples and per land-use type)
#  c  PCoA of ARG subtype composition (Bray-Curtis) with paired PERMANOVA
#  d  Relative abundance of the 12 most abundant ARG types
#  Also writes Supplementary Data 4 (detection, prevalence and abundance of every ARG subtype) and the distance-based
#  redundancy analysis of the resistome conditioned on community composition (cited in the text).
#  Output: output/Fig2/
# =====================================================================
source("scripts/00_functions.R")
use(c("ggh4x", "patchwork", "vegan", "permute"))
OUT  <- out_folder("Fig2")
meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
tpm_sub <- read.delim(ARG_FILES[["Subtype"]], row.names = 1, check.names = FALSE)

# ---------------- a. Number of detected ARG subtypes ----------------
p_a <- richness_panel("Subtype", OUT)

# ---------------- b. Shared and fraction-specific ARG subtypes ----------------
p_b <- local({
  HABITATS <- c("UR", "FO", "FA")
  CATS     <- c("Soil-only", "Shared", "EV-only")
  COL_CAT  <- c("Soil-only" = "#6FA3DA", "Shared" = "#D9D9D9", "EV-only" = "#F4A3A3")
  theme_bars <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = 0.5), axis.ticks.y = element_blank(),
          axis.text = element_text(colour = "black", size = 11), axis.text.y = element_text(face = "bold"),
          axis.title = element_text(face = "bold", colour = "black"), legend.title = element_blank(),
          legend.text = element_text(size = 11), legend.position = "bottom", legend.key.size = unit(0.9, "lines"))
  # detected = TPM > 0 in at least one sample of the fraction
  detected_in <- function(ids) rownames(tpm_sub)[rowSums(tpm_sub[, ids, drop = FALSE] > 0) > 0]
  res <- lapply(setNames(c("All", HABITATS), c("All", HABITATS)), function(sc) {
    keep <- if (sc == "All") rep(TRUE, nrow(meta)) else meta$Habitat == sc
    s <- detected_in(meta$SampleID[keep & meta$SampleType == "Soil"]); e <- detected_in(meta$SampleID[keep & meta$SampleType == "EV"])
    list(`Soil-only` = setdiff(s, e), Shared = intersect(s, e), `EV-only` = setdiff(e, s))
  })
  cnt <- imap_dfr(res, function(o, sc) tibble(Scope = sc, Category = CATS, n = lengths(o)[CATS])) %>%
    group_by(Scope) %>% mutate(total = sum(n), pct = 100 * n / total) %>% ungroup()
  save_table(cnt, OUT, "Fig2b_overlap_counts")

  circle <- function(x0, r = 1, n = 400) { t <- seq(0, 2 * pi, length.out = n); tibble(x = x0 + r * cos(t), y = r * sin(t)) }
  a <- 0.55; r <- 1
  cS <- circle(-a); cE <- circle(a)
  lens <- bind_rows(filter(cS, (x - a)^2 + y^2 <= r^2), filter(cE, (x + a)^2 + y^2 <= r^2)) %>% mutate(ang = atan2(y, x)) %>% arrange(ang)
  lab <- cnt %>% filter(Scope == "All") %>% mutate(x = c(-a - 0.5, 0, a + 0.5))
  p_venn <- ggplot() +
    geom_polygon(data = cE, aes(x, y), fill = COL_CAT[["EV-only"]], colour = NA) +
    geom_polygon(data = cS, aes(x, y), fill = COL_CAT[["Soil-only"]], colour = NA) +
    geom_polygon(data = lens, aes(x, y), fill = COL_CAT[["Shared"]], colour = NA) +
    geom_text(data = lab, aes(x, 0.28, label = Category), fontface = "bold", size = 4, family = FONT) +
    geom_text(data = lab, aes(x, 0.02, label = n), fontface = "bold", size = 5.2, family = FONT) +
    geom_text(data = lab, aes(x, -0.24, label = sprintf("(%.1f%%)", pct)), size = 3.9, family = FONT) +
    coord_fixed(clip = "off") + theme_void(base_family = FONT)
  bars <- cnt %>% filter(Scope %in% HABITATS) %>%
    mutate(Scope = factor(Scope, levels = rev(HABITATS)), Category = factor(Category, levels = rev(CATS)))
  p_bars <- ggplot(bars, aes(pct, Scope, fill = Category)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.1f", pct)), position = position_stack(vjust = 0.5), size = 3.5, family = FONT) +
    scale_fill_manual(values = COL_CAT, breaks = CATS) +
    scale_x_continuous(limits = c(0, 100.01), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(x = "Proportion of ARG subtypes (%)", y = NULL) + theme_bars
  p_venn / p_bars + plot_layout(heights = c(1.6, 1))
})

# Supplementary Data 4: every ARG subtype detected in bulk soil, in EVs or in both, with prevalence (share of the 27 samples
# of a fraction with TPM > 0), mean TPM over all 27 samples and abundance percentile among the subtypes detected in that fraction
local({
  ids_ev <- meta$SampleID[meta$SampleType == "EV"]; ids_so <- meta$SampleID[meta$SampleType == "Soil"]
  tab <- tpm_sub[rowSums(tpm_sub[, c(ids_ev, ids_so)] > 0) > 0, ]
  pct <- function(x, det) { out <- rep(NA_real_, length(x)); out[det] <- 100 * rank(x[det]) / sum(det); out }
  sub_tab <- tibble(Subtype = rownames(tab), Prev_Soil = rowMeans(tab[, ids_so] > 0), Prev_EV = rowMeans(tab[, ids_ev] > 0),
                    MeanTPM_Soil = rowMeans(tab[, ids_so]), MeanTPM_EV = rowMeans(tab[, ids_ev])) %>%
    mutate(Category = case_when(Prev_Soil > 0 & Prev_EV > 0 ~ "Shared", Prev_EV > 0 ~ "EV_only", TRUE ~ "Soil_only"),
           Pctl_Soil = pct(MeanTPM_Soil, Prev_Soil > 0), Pctl_EV = pct(MeanTPM_EV, Prev_EV > 0))
  save_table(sub_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 4))) %>% arrange(Category, desc(MeanTPM_Soil)),
             OUT, "Supplementary_Data_4_ARG_subtype_detection")
})

# ---------------- c. PCoA of ARG subtype composition ----------------
block_c <- pcoa_panel(tpm_sub, OUT, "Fig2c_ARG_subtype", perm_height = 0.36)

# ---------------- d. ARG types: mean relative abundance per land use and fraction ----------------
p_d <- local({
  TOP_N <- 12
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), panel.spacing.x = unit(0.9, "lines"),
          strip.text = element_text(face = "bold", size = 12), axis.text = element_text(colour = "black", size = 11),
          axis.text.x = element_text(face = "bold"), axis.title = element_text(face = "bold", colour = "black"),
          legend.title = element_blank(), legend.text = element_text(size = 11), legend.background = element_blank(),
          legend.key = element_blank())
  tpm <- read.delim(ARG_FILES[["Type"]], row.names = 1, check.names = FALSE)
  tpm <- as.matrix(tpm[, meta$SampleID])
  top <- names(sort(rowSums(tpm), decreasing = TRUE))[seq_len(TOP_N)]
  grp <- rbind(tpm[top, , drop = FALSE], Others = colSums(tpm[!rownames(tpm) %in% top, , drop = FALSE]))
  rel <- sweep(grp, 2, colSums(tpm), "/") * 100
  rownames(rel) <- norm_type(rownames(rel))
  class_levels <- c(norm_type(top), "Others")
  long <- as.data.frame(rel) %>% rownames_to_column("Class") %>%
    pivot_longer(-Class, names_to = "SampleID", values_to = "RelAbund") %>%
    left_join(meta, by = "SampleID") %>%
    mutate(Class = factor(Class, levels = class_levels), Habitat = factor(Habitat, levels = names(COL_HABITAT)),
           SampleType = factor(SampleType, levels = c("Soil", "EV")), Group = factor(Group, levels = GROUP_ORDER))
  means <- long %>% group_by(Habitat, SampleType, Class) %>%
    summarise(n = n(), Mean = mean(RelAbund), SD = sd(RelAbund), SE = SD / sqrt(n), .groups = "drop")
  save_table(means, OUT, "Fig2d_ARG_type_group_means")

  # paired Wilcoxon test per ARG type within each land use (cited in the text, not shown in the figure)
  ev_ids <- meta$SampleID[meta$SampleType == "EV"]; so_ids <- vapply(ev_ids, ev_to_soil, character(1))
  hab <- meta$Habitat[match(ev_ids, meta$SampleID)]
  tests <- expand_grid(Habitat = names(COL_HABITAT), Class = class_levels) %>%
    mutate(res = map2(Habitat, Class, function(h, cl) {
      e <- rel[cl, ev_ids[hab == h]]; s <- rel[cl, so_ids[hab == h]]
      tibble(n_pairs = length(e), mean_Soil = mean(s), mean_EV = mean(e), p = suppressWarnings(wilcox.test(e, s, paired = TRUE))$p.value)
    })) %>% unnest(res) %>% group_by(Habitat) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup()
  save_table(tests, OUT, "Fig2d_ARG_type_paired_tests")

  x_max <- ceiling(max(means$Mean) * 1.25 / 10) * 10
  plot_df <- means %>%
    mutate(Class = factor(Class, levels = rev(class_levels)),
           fill  = ifelse(SampleType == "Soil", COL_HABITAT[as.character(Habitat)], "white"),
           edge  = COL_HABITAT[as.character(Habitat)], x_lab = Mean + x_max * 0.025)
  n_h <- length(COL_HABITAT)
  strip_fill <- c(alpha(COL_HABITAT, 0.35), rep("white", 2 * n_h))
  strip_txt  <- c(rep(list(element_text(face = "bold", size = 12)), n_h), rep(list(element_text(size = 10.5)), 2 * n_h))
  ggplot(plot_df, aes(Mean, Class)) +
    geom_col(aes(fill = I(fill), colour = I(edge)), width = 0.68, linewidth = 0.5) +
    geom_text(aes(x = x_lab, label = sprintf("%.1f", Mean)), hjust = 0, size = 2.9, family = FONT) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max - 10, 10), expand = c(0, 0)) +
    facet_nested(~ Habitat + SampleType, nest_line = FALSE,
                 strip = strip_nested(background_x = elem_list_rect(fill = strip_fill, colour = NA), text_x = strip_txt)) +
    labs(x = "Relative abundance (%)", y = NULL) +
    theme_this + theme(axis.text.x = element_text(size = 9.5, face = "bold"))
})

# ---------------- Resistome vs community composition (text) ----------------
# Mantel test (Spearman) between ARG and genus Bray-Curtis distances, whole Soil-EV pairs permuted; distance-based redundancy
# analysis of ARG composition ~ fraction, unconstrained and conditioned on the genus PCoA axes that explain >= 80% of the
# genus variation, permutations within pairs.
local({
  set.seed(42)
  genus <- read.delim(TAXA_FILES[["Genus"]], row.names = 1, check.names = FALSE)
  ids2  <- intersect(meta$SampleID, intersect(colnames(tpm_sub), colnames(genus)))
  d_arg <- vegdist(t(as.matrix(tpm_sub[, ids2])), method = "bray")
  d_gen <- vegdist(t(as.matrix(genus[, ids2])), method = "bray")
  pairs <- tibble(EV = meta$SampleID[meta$SampleType == "EV"]) %>% mutate(Soil = vapply(EV, ev_to_soil, character(1)))
  m2 <- meta[match(ids2, meta$SampleID), ] %>%
    mutate(Pair = ifelse(SampleType == "EV", SampleID, pairs$EV[match(SampleID, pairs$Soil)]),
           Fraction = factor(ifelse(SampleType == "Soil", "Soil", "EV"), levels = c("Soil", "EV")))
  perm_units <- how(plots = Plots(strata = factor(m2$Pair), type = "free"), within = Within(type = "none"), nperm = 999)
  mt <- mantel(d_arg, d_gen, method = "spearman", permutations = perm_units)
  pco <- cmdscale(d_gen, k = 20, eig = TRUE)
  ev_pos <- pco$eig[pco$eig > 0]; cum <- cumsum(ev_pos) / sum(ev_pos)
  K <- which(cum >= 0.8)[1]
  G <- as.data.frame(pco$points[, 1:K]); colnames(G) <- paste0("G", 1:K)
  m2 <- bind_cols(m2, G)
  f_cond <- as.formula(paste("d_arg ~ Fraction + Condition(", paste(colnames(G), collapse = " + "), ")"))
  db_c <- capscale(f_cond, data = m2)
  db_u <- capscale(d_arg ~ Fraction, data = m2)
  perm_pair <- how(blocks = factor(m2$Pair), nperm = 999)
  a_c <- anova(db_c, permutations = perm_pair); a_u <- anova(db_u, permutations = perm_pair)
  tot <- db_u$tot.chi
  save_table(tibble(
    Test = c("Mantel (Spearman): genus vs ARG Bray-Curtis, all 54 samples (whole pairs permuted)",
             "dbRDA: ARG ~ fraction (pairs as blocks)",
             sprintf("dbRDA: ARG ~ fraction | genus composition (%d PCoA axes, %.0f%% of genus variation)", K, 100 * cum[K]),
             "dbRDA: share of ARG variation explained by genus composition (conditioning axes)"),
    Statistic = c(sprintf("r = %.2f", mt$statistic),
                  sprintf("F = %.2f; %.1f%% of ARG variation", a_u$F[1], 100 * db_u$CCA$tot.chi / tot),
                  sprintf("F = %.2f; %.1f%% of ARG variation", a_c$F[1], 100 * db_c$CCA$tot.chi / tot),
                  sprintf("%.1f%% of ARG variation", 100 * db_c$pCCA$tot.chi / tot)),
    p = c(mt$signif, a_u$`Pr(>F)`[1], a_c$`Pr(>F)`[1], NA)), OUT, "Resistome_vs_community")
})

# ---------------- Assemble (180 mm wide at print size; canvas x 1.8) ----------------
SCALE <- 1.8; PRINT_W <- 180; ROW_H <- c(top = 84, bottom = 72)
top <- tag_panel(p_a, "a") + tag_panel(p_b, "b") + tag_panel(block_c, "c") + plot_layout(widths = c(66, 54, 60))
fig <- wrap_elements(full = top) / tag_panel(p_d, "d") + plot_layout(heights = ROW_H)
save_fig(fig, OUT, "Fig2", w = PRINT_W * SCALE / 25.4, h = sum(ROW_H) * SCALE / 25.4)
message("Figure 2 written to ", OUT)
