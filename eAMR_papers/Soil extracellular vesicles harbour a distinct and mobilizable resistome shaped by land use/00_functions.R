# =====================================================================
#  00_functions.R - settings and helper functions shared by all figure and table scripts
# ---------------------------------------------------------------------
#  Sourced at the top of every script. Run the scripts from the project folder (open soil-EV-resistome.Rproj).
# =====================================================================

if (!file.exists("soil-EV-resistome.Rproj"))
  stop("Open soil-EV-resistome.Rproj (or set the working directory to the project folder) before running this script")

# Install missing CRAN packages and attach them
use <- function(pkgs) {
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(for (p in pkgs) library(p, character.only = TRUE))
}
use("tidyverse")

# ---------------- Input files (written by 00_prepare_inputs.R) ----------------
META_FILE  <- file.path("data", "raw", "metadata_complete.csv")
TAXA_FILES <- c(Phylum  = file.path("data", "raw", "Microbial_community", "bacteria_phylum_percent.tsv"),
                Genus   = file.path("data", "raw", "Microbial_community", "bacteria_genus_percent.tsv"),
                Species = file.path("data", "raw", "Microbial_community", "bacteria_species_percent.tsv"))
ARG_FILES  <- c(Type    = file.path("data", "raw", "ARG_OAP", "merged_tpm.type.txt"),
                Subtype = file.path("data", "raw", "ARG_OAP", "merged_tpm.subtype.txt"),
                Gene    = file.path("data", "raw", "ARG_OAP", "merged_tpm.gene.txt"))
RANK_FILE  <- file.path("data", "raw", "ARG_OAP", "merged_tpm_gene_high_risk_rankI_II.txt")
MGE_TYPE_FILE <- file.path("data", "raw", "MGE", "mge_per_sample_type_TPM.tsv")
MGE_GENE_FILE <- file.path("data", "raw", "MGE", "mge_no_phage_TPM_matrix.tsv")
# summary tables of the bioinformatic workflow (contigs and MAGs), read directly by Fig4, Fig5, FigS10_S11
DS <- if (dir.exists(file.path("data", "downstream"))) file.path("data", "downstream") else file.path("..", "downstream_results")

# ---------------- Colours and symbols used in all figures ----------------
FONT        <- "Arial"
COL_HABITAT <- c(UR = "#C8553D", FO = "#3E8E5E", FA = "#3A6EA5")   # land use: urban, forest, farmland
SHP_TYPE    <- c(Soil = 16, EV = 21)                                 # bulk soil filled, EV open
GROUP_ORDER <- c("UR_Soil", "UR_EVs", "FO_Soil", "FO_EVs", "FA_Soil", "FA_EVs")

# ---------------- Soil-EV pairs ----------------
# Each EV sample was isolated from one bulk-soil sample (27 pairs: UR 9, FO 4, FA 14). EV sample name -> soil sample name:
ev_to_soil <- function(ev) {
  if (grepl("^FA_TCK", ev)) return(paste0("FA_CK", sub("_EVs$", "", sub("^FA_TCK", "", ev)), "_Soil"))
  if (grepl("^FA_EV[0-9].*_20", ev)) return(sub("^FA_EV", "FA_soil", ev))
  sub("_EVs$", "", ev)
}

# ---------------- Output ----------------
out_folder <- function(id) {
  d <- file.path("output", id); dir.create(d, recursive = TRUE, showWarnings = FALSE); d
}
save_table <- function(df, dir, name) write.csv(df, file.path(dir, paste0(name, ".csv")), row.names = FALSE)

# PDF (quartz on macOS, cairo elsewhere), PNG at 300 dpi and TIFF at 600 dpi
save_fig <- function(p, dir, name, w, h) {
  pdf_dev <- if (capabilities("aqua")) {
    function(filename, width, height, ...) grDevices::quartz(file = filename, width = width, height = height, type = "pdf")
  } else {
    function(filename, width, height, ...) grDevices::cairo_pdf(filename, width = width, height = height)
  }
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = pdf_dev)
  ggsave(file.path(dir, paste0(name, ".png")),  p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(dir, paste0(name, ".tiff")), p, width = w, height = h, dpi = 600, compression = "lzw")
}

# Panel with a bold tag (a, b, ...) for the multi-panel figures
tag_panel <- function(p, tag, size = 18) {
  patchwork::wrap_elements(full = p) + labs(tag = tag) +
    theme(plot.tag = element_text(face = "bold", size = size, family = FONT), plot.tag.position = c(0.01, 0.99))
}

# Paired PERMANOVA (Bray-Curtis, 999 permutations). Land use, constant within a pair, is tested by permuting whole pairs;
# fraction and the land use x fraction interaction with pair as a term in the model and permutations within pairs.
# R2 values are those of ~ Habitat * SampleType. m needs the columns Habitat, SampleType and Pair.
paired_permanova <- function(d, m) {
  units  <- permute::how(plots = permute::Plots(strata = m$Pair, type = "free"), within = permute::Within(type = "none"), nperm = 999)
  within <- permute::how(blocks = m$Pair, nperm = 999)
  set.seed(42); a <- as.data.frame(vegan::adonis2(d ~ Habitat * SampleType, data = m, permutations = units, by = "terms"))
  set.seed(42); b <- as.data.frame(vegan::adonis2(d ~ Pair + SampleType + Habitat:SampleType, data = m, permutations = within, by = "terms"))
  ib <- grep(":", rownames(b), value = TRUE)
  stopifnot(isTRUE(all.equal(a[c("SampleType", "Habitat:SampleType"), "R2"], b[c("SampleType", ib), "R2"])))
  a[c("SampleType", "Habitat:SampleType"), c("F", "Pr(>F)")] <- b[c("SampleType", ib), c("F", "Pr(>F)")]
  a
}

# Metadata with ordered factors and the Soil-EV pair of each sample
read_meta <- function() {
  meta <- read.csv(META_FILE, stringsAsFactors = FALSE) %>%
    mutate(Group      = factor(Group, levels = GROUP_ORDER),
           Habitat    = factor(Habitat, levels = names(COL_HABITAT)),
           SampleType = factor(SampleType, levels = names(SHP_TYPE)))
  evs <- meta$SampleID[meta$SampleType == "EV"]
  pair_of <- setNames(rep(evs, 2), c(evs, vapply(evs, ev_to_soil, character(1))))
  meta$Pair <- factor(unname(pair_of[meta$SampleID]))
  stopifnot(!any(is.na(meta$Pair)), all(table(meta$Pair) == 2))
  meta
}

# =====================================================================
#  Analyses shown in more than one figure
# =====================================================================

# ---------------- Shannon diversity (Fig. 1b: phylum; Supplementary Fig. S3: genus, species) ----------------
# Box plots of bulk soil and EVs per land-use type; brackets: two-sided paired Wilcoxon signed-rank test, EV vs bulk soil
shannon_panel <- function(level, dir) {
  use("vegan")
  ALPHA_TYPE <- c(Soil = 0.35, EV = 0)
  BOX_W <- 0.55; DODGE <- 0.66
  Y_LIM <- list(Phylum  = list(lim = c(0, 2.2), brk = seq(0, 2, 0.5)),
                Genus   = list(lim = c(0, 7),   brk = 0:7),
                Species = list(lim = c(0, 10),  brk = seq(0, 9, 3)))
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey35", linewidth = 0.5),
          axis.text = element_text(colour = "black", size = 11), axis.text.x = element_text(face = "bold"),
          axis.title = element_text(face = "bold", colour = "black"), strip.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold"), legend.background = element_blank(), legend.key = element_blank())

  meta <- read_meta()
  tab  <- read.delim(TAXA_FILES[[level]], row.names = 1, check.names = FALSE)
  stopifnot(setequal(meta$SampleID, colnames(tab)))
  prop <- sweep(tab, 2, colSums(tab), "/")
  sh   <- vegan::diversity(t(prop), index = "shannon")
  df   <- meta %>% mutate(Level = level, Shannon = unname(sh[SampleID]))
  save_table(df %>% select(SampleID, Pair, Group, Habitat, SampleType, Shannon), dir, sprintf("Shannon_%s_values", level))

  # paired Wilcoxon test (EV vs bulk soil) within each land-use type
  fmt_p <- function(p) ifelse(p >= 0.05, "ns", ifelse(p < 0.001, "italic(p) < 0.001", sprintf("italic(p) == %.3f", p)))
  w <- df %>% select(Pair, Habitat, SampleType, Shannon) %>% pivot_wider(names_from = SampleType, values_from = Shannon)
  tests <- w %>% group_by(Habitat) %>%
    summarise(n_pairs = n(), median_Soil = median(Soil), median_EV = median(EV), n_EV_lower = sum(EV < Soil),
              V = unname(suppressWarnings(wilcox.test(EV, Soil, paired = TRUE))$statistic),
              p = suppressWarnings(wilcox.test(EV, Soil, paired = TRUE))$p.value, .groups = "drop") %>%
    mutate(label = fmt_p(p))
  save_table(tests %>% select(-label), dir, sprintf("Shannon_%s_paired_tests", level))
  ymax <- Y_LIM[[level]]$lim[2]
  lab_df <- df %>% group_by(Habitat) %>% summarise(top = max(Shannon), .groups = "drop") %>%
    left_join(tests, by = "Habitat") %>%
    mutate(x = as.numeric(Habitat), y = top + 0.06 * ymax, ypos = y + 0.015 * ymax)

  off <- DODGE / 4; tick <- 0.02 * ymax; dodge <- position_dodge(width = DODGE)
  ggplot(df, aes(Habitat, Shannon)) +
    geom_boxplot(aes(colour = Habitat, fill = Habitat, alpha = SampleType, group = interaction(Habitat, SampleType)),
                 position = dodge, width = BOX_W, outlier.shape = NA, linewidth = 0.6) +
    geom_point(aes(colour = Habitat, shape = SampleType, group = SampleType),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = DODGE, seed = 1),
               fill = "white", size = 1.8, stroke = 0.8, show.legend = FALSE) +
    geom_segment(data = lab_df, aes(x = x - off, xend = x + off, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = lab_df, aes(x = x - off, xend = x - off, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = lab_df, aes(x = x + off, xend = x + off, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = lab_df, aes(x = x, y = ypos, label = label), parse = TRUE, vjust = 0, size = 3.6,
              family = FONT, inherit.aes = FALSE) +
    scale_colour_manual(values = COL_HABITAT, guide = "none") +
    scale_fill_manual(values = COL_HABITAT, guide = "none") +
    scale_alpha_manual(values = ALPHA_TYPE, name = NULL,
                       guide = guide_legend(override.aes = list(fill = "grey45", colour = "grey20"))) +
    scale_shape_manual(values = SHP_TYPE, guide = "none") +
    scale_y_continuous(limits = Y_LIM[[level]]$lim, breaks = Y_LIM[[level]]$brk, expand = c(0, 0)) +
    labs(x = "Land-use type", y = "Shannon diversity") +
    theme_this +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          legend.position = "right", legend.justification = c(0, 1), legend.direction = "vertical",
          legend.key.size = unit(1.1, "lines"), legend.text = element_text(size = 11))
}

# ---------------- PCoA with paired PERMANOVA (Figs. 1c, 2c; Supplementary Figs. S4, S9b) ----------------
# tab: features x samples. Returns the plot with the PERMANOVA results written below it.
pcoa_panel <- function(tab, dir, name, perm_height = 0.3) {
  use(c("vegan", "ape", "permute", "patchwork"))
  ELLIPSE_ALPHA <- 0.12
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black", size = 11),
          axis.text.x = element_text(face = "bold"), axis.title = element_text(face = "bold", colour = "black"),
          legend.title = element_text(face = "bold"), legend.text = element_text(size = 11),
          legend.background = element_blank(), legend.key = element_blank())

  meta <- read_meta()
  ids  <- intersect(colnames(tab), meta$SampleID)
  m    <- meta[match(ids, meta$SampleID), ]
  d    <- vegdist(t(tab[, ids]), method = "bray")
  pc   <- pcoa(d)
  expl <- pc$values$Relative_eig[1:2] * 100
  coords <- tibble(SampleID = ids, PC1 = pc$vectors[, 1], PC2 = pc$vectors[, 2]) %>% left_join(m, by = "SampleID")
  ad2  <- paired_permanova(d, m)
  save_table(coords %>% select(SampleID, Group, Habitat, SampleType, PC1, PC2), dir, paste0(name, "_PCoA_coordinates"))
  save_table(as.data.frame(ad2) %>% rownames_to_column("Term"), dir, paste0(name, "_PERMANOVA"))

  p <- ggplot(coords, aes(PC1, PC2, colour = Habitat)) +
    stat_ellipse(aes(fill = Habitat), geom = "polygon", type = "t", level = 0.95,
                 alpha = ELLIPSE_ALPHA, colour = NA, show.legend = FALSE) +
    geom_point(aes(shape = SampleType), fill = "white", size = 2.8, stroke = 1) +
    coord_cartesian(clip = "off") +
    scale_colour_manual(values = COL_HABITAT, name = "Land-use type", guide = guide_legend(order = 1)) +
    scale_fill_manual(values = COL_HABITAT, guide = "none") +
    scale_shape_manual(values = SHP_TYPE, name = "Fraction",
                       guide = guide_legend(order = 2, override.aes = list(colour = "grey20"))) +
    labs(x = sprintf("PC1 (%.1f%%)", expl[1]), y = sprintf("PC2 (%.1f%%)", expl[2])) +
    theme_this +
    guides(colour = guide_legend(order = 1, ncol = 1, title.position = "top"),
           shape  = guide_legend(order = 2, ncol = 1, title.position = "top", override.aes = list(colour = "grey20"))) +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.justification = "left",
          legend.box.just = "top", plot.margin = margin(5, 8, 0, 5), aspect.ratio = 1)

  fp  <- function(p) if (p < 0.001) "italic(p) < 0.001" else sprintf("italic(p) == '%.3f'", p)
  lab <- c("'PERMANOVA'",
           sprintf("'Land use:'~italic(R)^2 == '%.3f'*','~%s",    ad2["Habitat", "R2"],            fp(ad2["Habitat", "Pr(>F)"])),
           sprintf("'Sample type:'~italic(R)^2 == '%.3f'*','~%s", ad2["SampleType", "R2"],         fp(ad2["SampleType", "Pr(>F)"])),
           sprintf("'Interaction:'~italic(R)^2 == '%.3f'*','~%s", ad2["Habitat:SampleType", "R2"], fp(ad2["Habitat:SampleType", "Pr(>F)"])))
  p_perm <- ggplot(tibble(y = 4:1, lab = lab), aes(0, y, label = lab)) +
    geom_text(parse = TRUE, hjust = 0, size = 3.6, family = FONT) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, 4.6), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(0, 5, 5, 12))
  p / p_perm + plot_layout(heights = c(1, perm_height))
}

# ---------------- Number of detected ARGs per sample (Fig. 2a: subtypes; Supplementary Fig. S6: types, genes) ----------------
# Detected = TPM > 0. Brackets: two-sided paired Wilcoxon signed-rank test, EV vs bulk soil, within each land-use type
richness_panel <- function(level, dir) {
  use("ggh4x")
  ALPHA_BOX <- c(Soil = 0.30, EV = 0)
  Y_LAB <- c(Subtype = "Number of detected ARG subtypes", Type = "Number of detected ARG types", Gene = "Number of detected ARGs")
  Y_LIM <- list(Subtype = list(lim = c(0, 500),  brk = seq(0, 500, 100)),
                Type    = list(lim = c(0, 30),   brk = seq(0, 30, 10)),
                Gene    = list(lim = c(0, 2500), brk = seq(0, 2500, 500)))
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), panel.spacing.x = unit(0.8, "lines"),
          strip.text = element_text(face = "bold", size = 12), axis.text = element_text(colour = "black", size = 11),
          axis.text.x = element_text(face = "bold"), axis.title = element_text(face = "bold", colour = "black"),
          legend.title = element_blank(), legend.text = element_text(size = 11), legend.background = element_blank(),
          legend.key = element_blank(), legend.position = "right", legend.justification = c(0, 1))

  meta  <- read.csv(META_FILE, stringsAsFactors = FALSE)
  pairs <- tibble(EV = meta$SampleID[meta$SampleType == "EV"]) %>%
    mutate(Soil = vapply(EV, ev_to_soil, character(1)), Pair = paste0("P", row_number()),
           Habitat = meta$Habitat[match(EV, meta$SampleID)])
  tab <- read.delim(ARG_FILES[[level]], row.names = 1, check.names = FALSE)
  stopifnot(all(meta$SampleID %in% colnames(tab)))
  n <- colSums(tab[, meta$SampleID] > 0)
  df <- pairs %>% pivot_longer(c(Soil, EV), names_to = "SampleType", values_to = "SampleID") %>%
    left_join(tibble(SampleID = names(n), n_detected = unname(n)), by = "SampleID") %>%
    mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)), SampleType = factor(SampleType, levels = names(SHP_TYPE)))
  save_table(df %>% select(Habitat, Pair, SampleType, SampleID, n_detected) %>% arrange(Habitat, Pair, SampleType),
             dir, sprintf("ARG_%s_richness_values", level))

  fmt_p <- function(p) ifelse(p >= 0.05, "ns", ifelse(p < 0.001, "italic(p) < 0.001", sprintf("italic(p) == %.3f", p)))
  tests <- df %>% select(Habitat, Pair, SampleType, n_detected) %>%
    pivot_wider(names_from = SampleType, values_from = n_detected) %>% group_by(Habitat) %>%
    summarise(n_pairs = n(), median_Soil = median(Soil), median_EV = median(EV), n_EV_higher = sum(EV > Soil),
              V = unname(suppressWarnings(wilcox.test(EV, Soil, paired = TRUE))$statistic),
              p = suppressWarnings(wilcox.test(EV, Soil, paired = TRUE))$p.value, .groups = "drop") %>%
    mutate(label = fmt_p(p))
  save_table(tests, dir, sprintf("ARG_%s_richness_paired_tests", level))

  ylim <- Y_LIM[[level]]$lim
  brk  <- tests %>% mutate(y = ylim[2] * 0.90)
  ggplot(df, aes(SampleType, n_detected)) +
    geom_boxplot(aes(fill = Habitat, alpha = SampleType, colour = Habitat), width = 0.45, outlier.shape = NA, linewidth = 0.5) +
    geom_point(aes(colour = Habitat, shape = SampleType), fill = "white", size = 2.3, stroke = 0.9) +
    geom_segment(data = brk, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = brk, aes(x = 1, xend = 1, y = y, yend = y - ylim[2] * 0.02), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = brk, aes(x = 2, xend = 2, y = y, yend = y - ylim[2] * 0.02), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = brk, aes(x = 1.5, y = y + ylim[2] * 0.015, label = label), parse = TRUE,
              vjust = 0, size = 3.6, family = FONT, inherit.aes = FALSE) +
    scale_colour_manual(values = COL_HABITAT, guide = "none") +
    scale_fill_manual(values = COL_HABITAT, guide = "none") +
    scale_alpha_manual(values = ALPHA_BOX, guide = "none") +
    scale_shape_manual(values = SHP_TYPE, guide = "none") +
    scale_y_continuous(limits = ylim, breaks = Y_LIM[[level]]$brk, expand = c(0, 0)) +
    facet_wrap2(~ Habitat, nrow = 1,
                strip = strip_themed(background_x = elem_list_rect(fill = alpha(COL_HABITAT, 0.35), colour = NA))) +
    labs(x = NULL, y = Y_LAB[[level]]) +
    theme_this
}

# ---------------- ARG class names ----------------
norm_type <- function(x) {
  x <- gsub("macrolide.*lincosamide.*streptogramin", "MLS", x, ignore.case = TRUE, perl = TRUE)
  x <- gsub("beta[ _-]?lactam", "Beta-lactam", x, ignore.case = TRUE, perl = TRUE)
  x <- gsub("_", " ", x)
  special <- x %in% c("MLS", "Others", "Beta-lactam")
  x[!special] <- paste0(toupper(substr(x[!special], 1, 1)), tolower(substr(x[!special], 2, nchar(x[!special]))))
  x
}

# ---------------- High-risk (Rank I) ARG subtypes: prevalence and abundance per group ----------------
# Fig. 3b (subtypes detected in >= 10% of all samples) and Supplementary Fig. S8 (all Rank I subtypes)
rankI_stats <- function() {
  meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
  risk <- read.delim(RANK_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(all(meta$SampleID %in% colnames(risk)))
  risk <- risk[risk$rank == "I", ]
  tpm <- rowsum(as.matrix(risk[, meta$SampleID]), risk$Subtype)       # subtype x sample
  sub_info <- tibble(Subtype = rownames(tpm)) %>%
    mutate(Class = norm_type(sub("__.*", "", Subtype)), Gene = sub(".*__", "", Subtype), prev_all = rowMeans(tpm > 0) * 100)
  sub_info$mean_all <- rowMeans(tpm)[sub_info$Subtype]
  as.data.frame(tpm) %>% rownames_to_column("Subtype") %>%
    pivot_longer(-Subtype, names_to = "SampleID", values_to = "TPM") %>%
    left_join(meta %>% select(SampleID, Habitat, SampleType), by = "SampleID") %>%
    group_by(Subtype, Habitat, SampleType) %>%
    summarise(n = n(), n_detected = sum(TPM > 0), prevalence = 100 * n_detected / n, mean_TPM = mean(TPM),
              mean_TPM_detected = ifelse(n_detected > 0, sum(TPM) / n_detected, NA_real_), .groups = "drop") %>%
    mutate(log_abund = log10(mean_TPM + 1)) %>%
    left_join(sub_info, by = "Subtype") %>%
    mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)), SampleType = factor(SampleType, levels = c("Soil", "EV")))
}

rankI_bubble <- function(df, abund_max) {
  use(c("ggh4x", "patchwork"))
  COL_ABUND   <- c("#FDE7E4", "#F4A6A0", "#D7453F", "#8B0A1A")
  COL_ZERO    <- "grey80"
  PREV_LABELS <- c("0", ">0\u201320", "20\u201340", "40\u201360", "60\u201380", ">80")
  PREV_SIZES  <- setNames(c(0.6, 2.2, 3.4, 4.6, 5.8, 7), PREV_LABELS)
  CLASS_GAP   <- 0.7
  COL_CLASS <- c(Multidrug = "#3F5E7D", Tetracycline = "#5B8FA8", Polymyxin = "#A0522D", Quinolone = "#6F9A8D",
                 Sulfonamide = "#9A8F2A", "Beta-lactam" = "#A05B6F", Aminoglycoside = "#C77F3A", MLS = "#6E5A8C",
                 Chloramphenicol = "#7A8590", Florfenicol = "#8E4A38")
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.4), panel.spacing.x = unit(0.15, "lines"),
          strip.text.x = element_text(face = "bold", size = 11.5), axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), legend.title = element_text(face = "bold", size = 10.5),
          legend.text = element_text(size = 10), legend.background = element_blank(), legend.key = element_blank(),
          legend.position = "bottom", legend.box = "horizontal")

  # rows: classes by their highest prevalence, subtypes within a class by prevalence
  rows <- df %>% distinct(Class, Gene, prev_all, mean_all) %>%
    group_by(Class) %>% mutate(cls_top = max(prev_all), cls_n = n()) %>% ungroup() %>%
    arrange(desc(cls_top), desc(cls_n), Class, desc(prev_all), desc(mean_all)) %>%
    mutate(cls_id = cumsum(Class != lag(Class, default = first(Class))), y = -(row_number() + cls_id * CLASS_GAP))
  df <- df %>% left_join(select(rows, Gene, y), by = "Gene") %>%
    mutate(prev_bin = cut(prevalence, breaks = c(-Inf, 0, 20, 40, 60, 80, 100), labels = PREV_LABELS))
  cls <- rows %>% group_by(Class) %>%
    summarise(ymin = min(y) - 0.38, ymax = max(y) + 0.38, ymid = (min(y) + max(y)) / 2, .groups = "drop") %>%
    mutate(col = ifelse(Class %in% names(COL_CLASS), COL_CLASS[as.character(Class)], "grey35"))
  y_lim <- c(min(rows$y) - 0.7, max(rows$y) + 0.7)

  # left: ARG class (coloured bar and name) and subtype
  p_left <- ggplot() +
    geom_segment(data = cls, aes(x = 0, xend = 0, y = ymin, yend = ymax, colour = I(col)), linewidth = 1.6) +
    geom_text(data = cls, aes(x = 0.035, y = ymid, label = Class, colour = I(col)),
              hjust = 0, fontface = "bold", size = 3.5, family = FONT) +
    geom_text(data = rows, aes(x = 1, y = y, label = Gene), hjust = 1, fontface = "italic", size = 3.5, family = FONT) +
    annotate("text", x = c(0, 1), y = y_lim[2] + 0.35, label = c("ARG class", "ARG subtype"),
             hjust = c(0, 1), vjust = 0, fontface = "bold", size = 3.7, family = FONT) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(-0.01, 1), ylim = y_lim, clip = "off") +
    theme_void(base_size = 12, base_family = FONT) + theme(plot.margin = margin(5, 6, 5, 5))

  # right: bubbles per land use and fraction
  n_h <- length(COL_HABITAT)
  strip_fill <- c(alpha(COL_HABITAT, 0.40), as.vector(rbind(alpha(COL_HABITAT, 0.22), alpha(COL_HABITAT, 0.10))))
  p_main <- ggplot(df, aes(1, y)) +
    geom_point(aes(size = prev_bin, colour = ifelse(n_detected > 0, log_abund, NA_real_)), show.legend = TRUE) +
    scale_size_manual(values = PREV_SIZES, drop = FALSE, name = "Prevalence (% samples)") +
    scale_colour_gradientn(colours = COL_ABUND, limits = c(0, abund_max), na.value = COL_ZERO,
                           name = expression(bold("Abundance (" * log[10] * "(TPM + 1))"))) +
    scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    facet_nested(~ Habitat + SampleType, nest_line = FALSE,
                 strip = strip_nested(background_x = elem_list_rect(fill = strip_fill, colour = NA),
                                      text_x = c(rep(list(element_text(face = "bold", size = 11.5)), n_h),
                                                 rep(list(element_text(face = "plain", size = 10.5)), 2 * n_h)))) +
    guides(size   = guide_legend(order = 1, nrow = 1, title.position = "top", label.position = "bottom",
                                 override.aes = list(colour = c(COL_ZERO, rep("grey35", 5)))),
           colour = guide_colourbar(order = 2, title.position = "top", barwidth = unit(9, "lines"),
                                    barheight = unit(0.7, "lines"))) +
    theme_this
  p_left + p_main + plot_layout(widths = c(1.15, 2))
}

# ---------------- ARG-MGE co-localization on the assembled contigs (Fig. 4; Supplementary Figs. S10, S11) ----------------
# An ARG-carrying contig is co-localized when another ORF on it is a core mobility gene: a named transposition, integration,
# recombination, plasmid replication or copy-control, mobilization or conjugation gene or ICE regulator (CORE_RE), or a
# mobileOG-db hit without a gene name in the integration/excision or transfer category. Hits of the phage category and ORFs
# that are also ARGs are not counted as MGE genes. C0 = any MGE gene (listed in Supplementary Data 10); C1 = core mobility gene.
CORE_RE     <- "^(tnp|tni|tns|ins[A-Z]|int$|intI|gin$|bin[0-9]*$|tIS|ltrA|rep[A-Z]?$|mob[A-Z]?$|pre$|tra[A-Z]|trb[A-Z]|vir[BD]|copR|immR)"
UNNAMED_RE  <- "^$|^NA:|^[A-Z]{1,5}_[0-9]+;?$"           # no gene name: empty, "NA:Keyword" or a locus tag only
UNNAMED_CAT <- c("integration/excision", "transfer")
GG_FILE     <- file.path(DS, "summary", "04_arg_mge_contigs", "gggenes_data.tsv")    # all genes on contigs with an ARG and an MGE hit
ARG_ORF_FILE <- file.path(DS, "summary", "04_arg_mge_contigs", "arg_orfs_detail.tsv") # all ARG ORFs on all contigs
el_type <- function(e) case_when(grepl("IS", e) ~ "IS", grepl("ICE|IME|CIME", e) ~ "IGE", grepl("plasmid", e) ~ "Plasmid", TRUE ~ "Others")

read_coloc <- function() {
  meta  <- read.csv(META_FILE, stringsAsFactors = FALSE)
  pairs <- tibble(EV = meta$SampleID[meta$SampleType == "EV"]) %>%
    mutate(Soil = vapply(EV, ev_to_soil, character(1)), Pair = paste0("P", row_number()), Habitat = meta$Habitat[match(EV, meta$SampleID)])
  pair_of <- setNames(c(pairs$Pair, pairs$Pair), c(pairs$EV, pairs$Soil))
  gg <- read.delim(GG_FILE, stringsAsFactors = FALSE, na.strings = c("", "NA")) %>%
    mutate(is_ARG = as.logical(is_ARG), is_MGE = as.logical(is_MGE),
           mge_name = replace_na(mge_name, ""), element_type = replace_na(element_type, ""))
  mge <- gg %>% filter(is_MGE, !is_ARG, mge_major_category != "phage") %>%
    mutate(unnamed = grepl(UNNAMED_RE, mge_name),
           core = grepl(CORE_RE, mge_name) | (unnamed & mge_major_category %in% UNNAMED_CAT),
           Type = ifelse(unnamed, "Others", el_type(element_type)),
           gene_label = ifelse(unnamed, paste0("unspecified ", mge_major_category, " (", mge_mobileOG_id, ")"), mge_name))
  arg_genes <- gg %>% filter(is_ARG)
  mge <- mge %>% rowwise() %>%                          # distance (bp) between each MGE gene and the nearest ARG
    mutate(gap_to_ARG = { a <- arg_genes[arg_genes$molecule == molecule, ]; min(pmax(0, pmax(a$start, start) - pmin(a$end, end))) }) %>%
    ungroup()
  arg_ctg <- read.delim(ARG_ORF_FILE, stringsAsFactors = FALSE) %>%
    distinct(sample, contig_id) %>%
    mutate(molecule = paste0(sample, "__", contig_id), length = as.numeric(sub(".*_length_([0-9]+)_.*", "\\1", contig_id))) %>%
    left_join(select(meta, sample = SampleID, Habitat, SampleType), by = "sample") %>%
    mutate(Pair = pair_of[sample],
           C0 = molecule %in% mge$molecule[mge$molecule %in% arg_genes$molecule],
           C1 = molecule %in% mge$molecule[mge$core & mge$molecule %in% arg_genes$molecule])
  stopifnot(!anyNA(arg_ctg$Pair), !anyNA(arg_ctg$length))
  list(meta = meta, pairs = pairs, gg = gg, mge = mge, arg_genes = arg_genes, arg_ctg = arg_ctg)
}

# ---------------- EV vs bulk soil over all 27 pairs (Fig. 1d, Supplementary Data 1-3; Supplementary Figs. S12, S13, Data 5, 6) ----
# One-sided paired Wilcoxon signed-rank tests in each direction, Benjamini-Hochberg per direction (q < 0.05); features with a
# mean relative abundance > 0.01% in either fraction are tested. log2FC = log2((mean EV + 1e-4) / (mean soil + 1e-4)), in %.
# 95% CI of log2FC: 2,000 bootstrap resamples of the pairs within each land-use type. Values per land use are descriptive.
ev_soil_setup <- function(n_boot = 2000) {
  meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
  ev_ids <- meta$SampleID[meta$SampleType == "EV"]
  hab <- factor(meta$Habitat[match(ev_ids, meta$SampleID)], levels = c("UR", "FO", "FA"))
  set.seed(42)
  boot <- replicate(n_boot, unlist(lapply(split(seq_along(hab), hab), function(ix) ix[sample.int(length(ix), replace = TRUE)])),
                    simplify = FALSE)
  list(ev_ids = ev_ids, so_ids = vapply(ev_ids, ev_to_soil, character(1)), hab = hab, boot = boot)
}
lfc <- function(e, s) log2((mean(e) + 1e-4) / (mean(s) + 1e-4))

ev_soil_test <- function(tab, st) {             # tab: features x samples, relative abundance (%)
  hab <- st$hab
  E <- as.matrix(tab[, st$ev_ids]); S <- as.matrix(tab[, st$so_ids])
  keep <- rowMeans(E) > 0.01 | rowMeans(S) > 0.01
  E <- E[keep, , drop = FALSE]; S <- S[keep, , drop = FALSE]
  p_ev   <- sapply(seq_len(nrow(E)), function(i) suppressWarnings(wilcox.test(E[i, ], S[i, ], paired = TRUE, alternative = "greater")$p.value))
  p_soil <- sapply(seq_len(nrow(E)), function(i) suppressWarnings(wilcox.test(E[i, ], S[i, ], paired = TRUE, alternative = "less")$p.value))
  res <- tibble(Taxon = rownames(E), mean_EV = rowMeans(E), mean_Soil = rowMeans(S),
                p_EV = p_ev, q_EV = p.adjust(p_ev, "BH"), p_Soil = p_soil, q_Soil = p.adjust(p_soil, "BH")) %>%
    mutate(Direction = case_when(q_EV < 0.05 ~ "EV", q_Soil < 0.05 ~ "Soil", TRUE ~ "ns"),
           q = pmin(q_EV, q_Soil),
           log2FC = map2_dbl(Taxon, Taxon, ~ lfc(E[.x, ], S[.x, ])),
           UR = map_dbl(Taxon, ~ lfc(E[.x, hab == "UR"], S[.x, hab == "UR"])),
           FO = map_dbl(Taxon, ~ lfc(E[.x, hab == "FO"], S[.x, hab == "FO"])),
           FA = map_dbl(Taxon, ~ lfc(E[.x, hab == "FA"], S[.x, hab == "FA"])))
  list(res = res, E = E, S = S)
}

# Features shown: all significant ones, at most n_show per direction (lowest q); the 10 most abundant if none is significant
ev_soil_shown <- function(lv, st, n_show) {
  res <- lv$res
  n_ev <- sum(res$Direction == "EV"); n_so <- sum(res$Direction == "Soil")
  show <- if (n_ev + n_so == 0) {
    res %>% slice_max(mean_EV + mean_Soil, n = 10)
  } else if (n_ev > n_show || n_so > n_show) {
    bind_rows(res %>% filter(Direction == "EV")   %>% arrange(q, desc(abs(log2FC))) %>% head(n_show),
              res %>% filter(Direction == "Soil") %>% arrange(q, desc(abs(log2FC))) %>% head(n_show))
  } else res %>% filter(Direction != "ns")
  ci <- map_dfr(show$Taxon, function(g) {
    bt <- vapply(st$boot, function(ix) lfc(lv$E[g, ix], lv$S[g, ix]), numeric(1))
    tibble(lo = unname(quantile(bt, 0.025)), hi = unname(quantile(bt, 0.975)))
  })
  bind_cols(show, ci) %>%
    mutate(n_same_direction = (sign(UR) == sign(log2FC)) + (sign(FO) == sign(log2FC)) + (sign(FA) == sign(log2FC)))
}

# Forest plot (log2FC with CI and q) and heatmap of log2FC per land use. lab_fun(names) returns markdown labels.
ev_soil_plot <- function(df, lab_fun, heat_width = 1.1, xtitle_size = NULL) {
  use(c("ggtext", "patchwork"))
  COL_DIR <- c(EV = "#542788", Soil = "#B35806", ns = "grey45")   # purple: enriched in EVs; orange-brown: enriched in bulk soil
  HABITATS <- c("UR", "FO", "FA")
  theme_this <- theme_bw(base_size = 12, base_family = FONT) +
    theme(plot.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          axis.text = element_text(colour = "black", size = 11), axis.text.x = element_text(face = "bold"),
          axis.title = element_text(face = "bold", colour = "black"), legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 11), legend.background = element_blank(), legend.key = element_blank())
  fmt_q <- function(q) {
    out <- character(length(q))
    for (i in seq_along(q)) {
      if (q[i] >= 0.05) { out[i] <- "'ns'"; next }
      e <- floor(log10(q[i])); m <- round(q[i] / 10^e, 1)
      if (m >= 10) { m <- 1; e <- e + 1 }
      out[i] <- sprintf("'%.1f' %%*%% 10^%d", m, e)
    }
    out
  }
  df <- df %>% arrange(factor(Direction, levels = c("EV", "ns", "Soil")), log2FC) %>%
    mutate(Taxon = factor(Taxon, levels = rev(Taxon)), q_lab = fmt_q(q))
  n_g   <- nrow(df)
  y_lab <- setNames(sprintf("<span style='color:%s'>%s</span>", COL_DIR[df$Direction], lab_fun(as.character(df$Taxon))),
                    as.character(df$Taxon))
  rng   <- range(c(df$lo, df$hi, 0))
  x_lim <- c(floor(rng[1]), ceiling(rng[2]))
  xq    <- x_lim[2] + diff(x_lim) * 0.10
  x_lab <- paste0("log<sub>2</sub>(EV / Soil)<br>",
                  "<span style='color:", COL_DIR["Soil"], "'>&larr; Enriched in bulk soil</span>",
                  "<span style='color:white'>MMMM</span>",
                  "<span style='color:", COL_DIR["EV"], "'>Enriched in EVs &rarr;</span>")
  p_forest <- ggplot(df, aes(log2FC, Taxon, colour = Direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.3, linewidth = 0.55) +
    geom_point(size = 2.4) +
    geom_text(aes(x = xq, label = q_lab), parse = TRUE, hjust = 0, size = 3.1, family = FONT, colour = "black") +
    annotate("text", x = xq, y = n_g + 0.9, label = "italic(q)~value", parse = TRUE, hjust = 0, size = 3.4, family = FONT) +
    scale_colour_manual(values = COL_DIR, guide = "none") +
    scale_y_discrete(labels = y_lab, expand = expansion(add = 0.5)) +
    coord_cartesian(xlim = x_lim, clip = "off") +
    labs(x = x_lab, y = NULL) +
    theme_this +
    theme(axis.text.y = element_markdown(size = 10.5), axis.title.x = element_markdown(size = xtitle_size),
          plot.margin = margin(15, 75, 5, 5))
  hab_long <- df %>% select(Taxon, all_of(HABITATS)) %>%
    pivot_longer(-Taxon, names_to = "Habitat", values_to = "lfc") %>% mutate(Habitat = factor(Habitat, levels = HABITATS))
  lim_h <- max(abs(hab_long$lfc))
  p_heat <- ggplot(hab_long, aes(Habitat, Taxon, fill = lfc)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    scale_fill_gradient2(low = COL_DIR[["Soil"]], mid = "white", high = COL_DIR[["EV"]], limits = c(-lim_h, lim_h),
                         name = "log<sub>2</sub>(EV/Soil)") +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(x = NULL, y = NULL) +
    theme_this +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title = element_markdown(),
          legend.position = "right", legend.justification = c(0, 1))
  p_forest + p_heat + plot_layout(widths = c(3, heat_width))
}
