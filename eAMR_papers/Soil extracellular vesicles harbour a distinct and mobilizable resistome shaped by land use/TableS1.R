# =====================================================================
#  TableS1.R - Supplementary Table S1: Soil-EV sample pairs
# ---------------------------------------------------------------------
#  One row per pair: land use, region, data source, mean coordinates of the pooled sampling points, and clean and
#  virus-free read numbers. Sampling dates were added from the field records. Output: output/TableS1/
# =====================================================================
source("scripts/00_functions.R")
OUT  <- out_folder("TableS1")
meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
pts  <- readr::read_csv(file.path("data", "raw", "Sampling_map", "sampling_points.csv"), show_col_types = FALSE,
                        locale = readr::locale(encoding = "UTF-8")) %>% as.data.frame()
qc   <- read.delim(file.path("data", "raw", "Sequencing", "qc_summary_all.tsv"), check.names = FALSE)
reads <- setNames(qc$Clean_reads, qc$Sample); final <- setNames(qc$Final_reads, qc$Sample)

# sampling points of each soil sample (the Leipzig-Halle samples are composites of three locations)
pt_key <- c(UR_Park = "^Park_soil_", UR_Roadside = "^Roadside_soil_", UR_Sandy = "^Sandy_soil_",
            FO_Forest = "^Forest_soil_", FA_Farm = "^Farm_soil_",
            # Xiamen sites: SampleID in sampling_points.csv is the Chinese site name, written here as Unicode escapes
            UR_CHL = "\u6e56\u91cc\u516c\u56ed",   # Huli Park
            UR_CNH = "\u5357\u6e56\u516c\u56ed",   # Nanhu Park
            UR_CZS = "\u4e2d\u5c71\u516c\u56ed",   # Zhongshan Park
            UR_SJG = "\u5609\u5e9a\u516c\u56ed",   # Jiageng Park
            UR_SNB = "\u5b81\u5b9d\u516c\u56ed",   # Ningbao Park
            UR_SRD = "\u65e5\u4e1c\u516c\u56ed",   # Ridong Park
            FO_VTY = "\u5929\u7afa\u5c71\u4e00\u53f7\u70b9",   # Tianzhushan site 1
            FO_VTE = "\u5929\u7afa\u5c71\u4e8c\u53f7\u70b9",   # Tianzhushan site 2
            FO_VTS = "\u5929\u7afa\u5c71\u4e09\u53f7\u70b9")   # Tianzhushan site 3
bl_pattern <- "^FA_soil(0|20|30)_(202410|202503)$"           # Bad Lauchstaedt field experiment (manure 0, 20, 30)
points_for <- function(s) {
  if (s %in% names(pt_key)) {
    k <- pt_key[[s]]
    return(pts[if (startsWith(k, "^")) grepl(k, pts$SampleID) else pts$SampleID == k, ])
  }
  if (grepl(bl_pattern, s)) return(pts[pts$SampleID == paste0("manure", sub(bl_pattern, "\\1", s), "_", substr(sub(bl_pattern, "\\2", s), 3, 6)), ])
  if (grepl("^FA_CK", s)) return(pts[pts$SampleID == sub("^FA_(CK[0-9]+)_Soil$", "\\1", s), ])
  pts[0, ]
}

tab <- bind_rows(lapply(meta$SampleID[meta$SampleType == "EV"], function(ev) {
  so <- ev_to_soil(ev); m <- meta[meta$SampleID == so, ]; p <- points_for(so)
  tibble(Soil_sample = so, EV_sample = ev, Land_use = m$Habitat,
         Region = if (m$Origin == "Germany") "Leipzig-Halle, Germany" else if (m$Habitat == "FA") "Ningbo, China" else "Xiamen, China",
         Source = if (m$Origin == "Germany") "This study" else if (m$Habitat == "FA") "Qin et al. (2026)" else "Zhu et al. (2022)",
         n_locations = nrow(p), Latitude = round(mean(p$Latitude), 4), Longitude = round(mean(p$Longitude), 4),
         Clean_reads_soil_M = round(reads[[so]] / 1e6, 1), Clean_reads_EV_M = round(reads[[ev]] / 1e6, 1),
         Virus_free_reads_soil_M = round(final[[so]] / 1e6, 1), Virus_free_reads_EV_M = round(final[[ev]] / 1e6, 1))
})) %>% arrange(factor(Land_use, levels = c("UR", "FO", "FA")), Region, Soil_sample)
save_table(tab, OUT, "Table_S1_samples")
