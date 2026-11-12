# =====================================================================
#  TableS3.R - Supplementary Table S3: land use, fraction and sampling region (paired PERMANOVA)
# ---------------------------------------------------------------------
#  Bacterial genera, ARG subtypes and MGE genes; Bray-Curtis, 999 permutations, seed 42.
#  a  All 27 pairs: ~ land use x fraction (M0) and ~ country + land use x fraction (M1, country first). Country and land use
#     (constant within a pair) are tested by permuting whole pairs (land use after country: within country); fraction and the
#     interaction with pair as a term in the model and permutations within pairs.
#  b  Within Germany, China and Xiamen alone: ~ land use + fraction (marginal), whole pairs permuted for land use; land use
#     within the EV and the bulk-soil fraction separately.
#  Output: output/TableS3/
# =====================================================================
source("scripts/00_functions.R")
use(c("vegan", "permute"))
OUT  <- out_folder("TableS3")
meta <- read_meta() %>%
  mutate(Region  = factor(case_when(Origin == "Germany" ~ "Leipzig-Halle", Habitat == "FA" ~ "Ningbo", TRUE ~ "Xiamen"),
                          levels = c("Leipzig-Halle", "Xiamen", "Ningbo")),
         Country = factor(ifelse(Origin == "Germany", "Germany", "China"), levels = c("Germany", "China")))
g   <- read.delim(MGE_GENE_FILE, check.names = FALSE, stringsAsFactors = FALSE)
mge <- rowsum(as.matrix(g[, intersect(meta$SampleID, colnames(g))]), group = g$sseqid_gene_name); rm(g)
mge <- as.data.frame(sweep(mge, 2, colSums(mge), "/") * 1e6)
TABS <- list(Genus = read.delim(TAXA_FILES[["Genus"]], row.names = 1, check.names = FALSE),
             ARG_subtype = read.delim(ARG_FILES[["Subtype"]], row.names = 1, check.names = FALSE),
             MGE_gene = mge)

pair_units <- function(m, blocks = NULL) {
  if (is.null(blocks)) how(plots = Plots(strata = m$Pair, type = "free"), within = Within(type = "none"), nperm = 999)
  else how(blocks = blocks, plots = Plots(strata = m$Pair, type = "free"), within = Within(type = "none"), nperm = 999)
}
within <- function(m) how(blocks = m$Pair, nperm = 999)
ad <- function(f, m, perm, by = "terms") { set.seed(42); as.data.frame(adonis2(f, data = m, permutations = perm, by = by)) }
row_of <- function(a, term) tibble(Df = a[term, "Df"], R2 = round(a[term, "R2"], 3), F = round(a[term, "F"], 2), p = a[term, "Pr(>F)"])
subsets <- list("Germany (Leipzig-Halle)" = quote(Country == "Germany"), "China (Xiamen + Ningbo)" = quote(Country == "China"),
                "China, Xiamen only (UR vs FO)" = quote(Region == "Xiamen"))
models <- list(); wc <- list()
for (ds in names(TABS)) {
  tab <- TABS[[ds]]; ids <- intersect(meta$SampleID, colnames(tab)); m <- meta[match(ids, meta$SampleID), ]
  dm <- as.matrix(vegdist(t(as.matrix(tab[, ids])), method = "bray")); d <- as.dist(dm)
  free <- ad(d ~ Habitat * SampleType, m, 999)                                  # R2 of the interaction
  a0 <- ad(d ~ Habitat * SampleType, m, pair_units(m))
  b  <- ad(d ~ Pair + SampleType + Habitat:SampleType, m, within(m)); ib <- grep(":", rownames(b), value = TRUE)
  m0 <- bind_rows(row_of(a0, "Habitat") %>% mutate(Term = "Land use", Permutation = "whole pairs"),
                  row_of(b, "SampleType") %>% mutate(Term = "Fraction", Permutation = "within pairs (pair in model)"),
                  row_of(b, ib) %>% mutate(Term = "Land use x fraction", Permutation = "within pairs (pair in model)",
                                           R2 = round(free["Habitat:SampleType", "R2"], 3)))
  a1c <- ad(d ~ Country + Habitat * SampleType, m, pair_units(m))
  a1h <- ad(d ~ Country + Habitat * SampleType, m, pair_units(m, blocks = m$Country))
  models[[ds]] <- bind_rows(
    m0 %>% mutate(Model = "M0 ~ land use x fraction"),
    bind_rows(row_of(a1c, "Country") %>% mutate(Term = "Country", Permutation = "whole pairs"),
              row_of(a1h, "Habitat") %>% mutate(Term = "Land use", Permutation = "whole pairs, within country"),
              m0 %>% filter(Term != "Land use")) %>% mutate(Model = "M1 ~ country + land use x fraction")) %>%
    mutate(Dataset = ds, .before = 1)
  for (s in names(subsets)) {
    k <- eval(subsets[[s]], m); mm <- droplevels(m[k, ]); dd <- as.dist(dm[k, k])
    lu <- ad(dd ~ Habitat + SampleType, mm, pair_units(mm), by = "margin")
    fr <- ad(dd ~ Pair + SampleType, mm, within(mm))
    one <- function(frac) { k2 <- mm$SampleType == frac; m3 <- droplevels(mm[k2, ]); d3 <- as.dist(as.matrix(dd)[k2, k2])
      set.seed(42); as.data.frame(adonis2(d3 ~ Habitat, data = m3, permutations = 999)) }   # one sample per pair
    ev <- one("EV"); so <- one("Soil")
    wc[[paste(ds, s)]] <- tibble(Dataset = ds, Subset = s, n_pairs = nlevels(mm$Pair),
      Land_use_groups = paste(names(table(mm$Habitat)), table(mm$Habitat) / 2, collapse = ", "),
      LandUse_R2 = round(lu["Habitat", "R2"], 3), LandUse_F = round(lu["Habitat", "F"], 2), LandUse_p = lu["Habitat", "Pr(>F)"],
      Fraction_R2 = round(fr["SampleType", "R2"], 3), Fraction_p = fr["SampleType", "Pr(>F)"],
      LandUse_R2_EV_only = round(ev[1, "R2"], 3), LandUse_p_EV_only = ev[1, "Pr(>F)"],
      LandUse_R2_soil_only = round(so[1, "R2"], 3), LandUse_p_soil_only = so[1, "Pr(>F)"])
  }
}
save_table(bind_rows(models), OUT, "Table_S3a_all_pairs")
save_table(bind_rows(wc), OUT, "Table_S3b_within_country")
