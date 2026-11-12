# =====================================================================
#  TableS2.R - Supplementary Table S2: viral sequences in bulk-soil and EV metagenomes
# ---------------------------------------------------------------------
#  Per sample: reads classified as viral by Kraken2 (% of read pairs), viral and provirus-containing contigs identified by
#  geNomad (% of contigs >= 500 bp) and read pairs removed before the downstream analyses (%); range, mean and median per
#  fraction; paired comparison of EV and bulk soil. Output: output/TableS2/
# =====================================================================
source("scripts/00_functions.R")
OUT  <- out_folder("TableS2")
meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
vir  <- read.delim(file.path("data", "raw", "Sequencing", "virus_removal_per_sample.tsv"), stringsAsFactors = FALSE)
pairs <- tibble(EV = meta$SampleID[meta$SampleType == "EV"]) %>% mutate(Soil = vapply(EV, ev_to_soil, character(1)))
val <- function(s, col) vir[[col]][match(s, vir$sample)]
per <- bind_rows(tibble(Fraction = "Soil", Sample = pairs$Soil), tibble(Fraction = "EV", Sample = pairs$EV)) %>%
  mutate(Viral_reads_pct = val(Sample, "kraken_viral_pct"), Viral_contigs_pct = val(Sample, "viral_contig_pct"),
         Reads_removed_pct = val(Sample, "removed_pct"))
stopifnot(!anyNA(per))
summ <- per %>% group_by(Fraction) %>%
  summarise(across(c(Viral_reads_pct, Viral_contigs_pct, Reads_removed_pct),
                   list(min = min, max = max, mean = mean, median = median)), .groups = "drop")
save_table(per, OUT, "Table_S2_viral_per_sample")
save_table(summ, OUT, "Table_S2_summary")
# EV vs paired bulk soil
pv <- pairs %>% mutate(across(c(EV, Soil), list(reads = ~ val(.x, "kraken_viral_pct"), contigs = ~ val(.x, "viral_contig_pct"))))
save_table(tibble(Measure = c("Viral reads", "Viral contigs"),
                  Pairs_EV_higher = c(sum(pv$EV_reads > pv$Soil_reads), sum(pv$EV_contigs > pv$Soil_contigs)),
                  n_pairs = nrow(pv),
                  p_paired_Wilcoxon = c(wilcox.test(pv$EV_reads, pv$Soil_reads, paired = TRUE)$p.value,
                                        wilcox.test(pv$EV_contigs, pv$Soil_contigs, paired = TRUE)$p.value)),
           OUT, "Table_S2_paired_test")
