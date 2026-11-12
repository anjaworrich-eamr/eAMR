# =====================================================================
#  FigS12_S13.R - Supplementary Figs. S12 (ARG types) and S13 (ARG subtypes) enriched in EV-associated or bulk-soil DNA
# ---------------------------------------------------------------------
#  Relative abundance = TPM as a percentage of the total ARG TPM of the sample. Paired tests over all 27 pairs as for Fig. 1d
#  (ev_soil_test in 00_functions.R); all EV-enriched entries and as many bulk-soil-enriched entries (lowest q) are shown,
#  at most 21 per direction. Also writes Supplementary Data 5 and 6 (all tested types and subtypes) and Supplementary Data 7
#  (ARG subtypes differing between EV and bulk soil within each land-use type).
#  Output: output/FigS12_S13/
# =====================================================================
source("scripts/00_functions.R")
OUT <- out_folder("FigS12_S13")
st  <- ev_soil_setup()

type_name <- function(x) {
  x <- dplyr::recode(x, `macrolide-lincosamide-streptogramin` = "MLS", beta_lactam = "beta-lactam",
                     other_peptide_antibiotics = "other peptide antibiotics", pleuromutilin_tiamulin = "pleuromutilin")
  x <- gsub("_", " ", x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
lab_arg <- function(x) ifelse(grepl("__", x),                # type: "Rifamycin"; subtype: "<i>rpoB2</i> (rifamycin)"
                              sprintf("<i>%s</i> (%s)", sub(".*__", "", x), { t <- type_name(sub("__.*", "", x)); ifelse(t == "MLS", t, tolower(t)) }),
                              type_name(x))

for (level in c("Type", "Subtype")) {
  tab <- read.delim(ARG_FILES[[level]], row.names = 1, check.names = FALSE)
  tab <- sweep(as.matrix(tab), 2, colSums(tab), "/") * 100      # % of the total ARG TPM of each sample
  lv  <- ev_soil_test(tab, st)
  save_table(lv$res %>% arrange(Direction, q), OUT,
             if (level == "Type") "Supplementary_Data_5_ARG_types" else "Supplementary_Data_6_ARG_subtypes")
  show <- ev_soil_shown(lv, st, n_show = 21)
  save_table(show, OUT, sprintf("ARG_%s_plotted", level))
  p <- ev_soil_plot(show, lab_arg)
  save_fig(p, OUT, if (level == "Type") "FigS12_ARG_types" else "FigS13_ARG_subtypes",
           w = 7.5 + 0.06 * max(nchar(show$Taxon)), h = 0.27 * nrow(show) + 2)
}

# ---------------- Supplementary Data 7: within each land-use type ----------------
# One-sided paired Wilcoxon signed-rank tests on TPM in each direction; subtypes with identical values in all pairs are not
# tested; Benjamini-Hochberg across the subtypes of one land use and direction (FDR); log2FC = log2((mean EV + 1)/(mean soil + 1)).
local({
  meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
  tpm  <- as.matrix(read.delim(ARG_FILES[["Subtype"]], row.names = 1, check.names = FALSE))
  pairs <- tibble(EV = meta$SampleID[meta$SampleType == "EV"]) %>%
    mutate(Soil = vapply(EV, ev_to_soil, character(1)), Habitat = meta$Habitat[match(EV, meta$SampleID)])
  res <- list()
  for (h in c("UR", "FO", "FA")) {
    d <- filter(pairs, Habitat == h)
    E <- tpm[, d$EV, drop = FALSE]; S <- tpm[, d$Soil, drop = FALSE]
    keep <- rownames(tpm)[rowSums(E != S) > 0]
    l2fc <- log2((rowMeans(E[keep, , drop = FALSE]) + 1) / (rowMeans(S[keep, , drop = FALSE]) + 1))
    for (alt in c("greater", "less")) {
      pv <- vapply(keep, function(g) suppressWarnings(wilcox.test(E[g, ], S[g, ], paired = TRUE, alternative = alt)$p.value), numeric(1))
      res[[paste(h, alt)]] <- tibble(Land_use = h, Direction = ifelse(alt == "greater", "EV > soil", "soil > EV"),
        n_pairs = nrow(d), Subtype = keep, Type = sub("__.*", "", keep), Gene = sub(".*__", "", keep),
        Mean_EV = rowMeans(E[keep, , drop = FALSE]), Mean_soil = rowMeans(S[keep, , drop = FALSE]),
        n_EV_higher = rowSums(E[keep, , drop = FALSE] > S[keep, , drop = FALSE]), log2FC = l2fc, p = pv, FDR = p.adjust(pv, "BH")) %>%
        mutate(Significant = FDR < 0.05 & if (alt == "greater") log2FC > 0 else log2FC < 0)
    }
  }
  all_tests <- bind_rows(res) %>% mutate(across(c(Mean_EV, Mean_soil), ~ round(.x, 2)), log2FC = round(log2FC, 3),
                                         p = signif(p, 3), FDR = signif(FDR, 3))
  save_table(all_tests %>% filter(Significant) %>% arrange(Land_use, Direction, FDR), OUT, "Supplementary_Data_7_ARG_subtypes_within_land_use")
  save_table(all_tests %>% group_by(Land_use, Direction, n_pairs) %>% summarise(n_tested = n(), n_FDR = sum(Significant), .groups = "drop"),
             OUT, "Supplementary_Data_7_summary")
})
