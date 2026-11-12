# =====================================================================
#  00_prepare_inputs.R - input tables for all figure and table scripts
# ---------------------------------------------------------------------
#  Reads the summary tables of the bioinformatic workflow (virus-free reads and contigs) and the sample information,
#  and writes the tables used by the other scripts to data/raw/.
#  Input : data/downstream/ (or ../downstream_results/)    summary tables of the workflow, see README
#          data/sample_info/ (or ../R_code_for_upload/data/raw/)   metadata_complete.csv, revision/qc_summary_all.tsv,
#                                                            Sampling_map/sampling_points.csv
#  Notes : Kraken2: percentages of bacterial reads classified at a rank are re-scaled to 100% per sample; taxa with the
#          same name (different taxid) are summed.
#          MGE: read-based mobileOG-db TPM; entries from phage databases are removed and each sample is re-scaled to
#          1e6 TPM. Entries without a gene name are kept. MGE types from the source databases: IS if any source is IS,
#          integrative element if ICE / AICE / IME / CIME, plasmid if plasmid, otherwise (ACLAME only) other.
#          Read numbers of the quality-control table are renamed with sample_rename.tsv.
# =====================================================================
source("scripts/00_functions.R")
INFO <- if (dir.exists(file.path("data", "sample_info"))) file.path("data", "sample_info") else file.path("..", "R_code_for_upload", "data", "raw")
stopifnot(dir.exists(DS), dir.exists(INFO))
RAW <- file.path("data", "raw")
for (d in c("Microbial_community", "ARG_OAP", "MGE", "Sequencing", "Sampling_map")) dir.create(file.path(RAW, d), recursive = TRUE, showWarnings = FALSE)
rd <- function(...) read.delim(file.path(...), check.names = FALSE, stringsAsFactors = FALSE)
wr <- function(df, ...) write.table(df, file.path(...), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------- Metadata: names and sample types must agree with the workflow ----------------
meta <- read.csv(file.path(INFO, "metadata_complete.csv"), stringsAsFactors = FALSE)
info <- rd(DS, "summary", "00_sample_info.tsv")
chk  <- inner_join(meta, info, by = c("SampleID" = "sample"), suffix = c("", ".server"))
stopifnot(nrow(chk) == 54, nrow(meta) == 54, all(chk$SampleType == chk$SampleType.server),
          all(chk$Habitat == chk$Habitat.server), all(chk$Origin == chk$Origin.server))
write.csv(meta, file.path(RAW, "metadata_complete.csv"), row.names = FALSE)
SAMPLES <- meta$SampleID

# ---------------- Bacterial composition (Kraken2) ----------------
taxa_table <- function(level) {
  x <- rd(DS, "summary", "02_kraken2_bacteria", paste0("bacteria_", level, "_relabund_pct.tsv"))
  names(x)[1] <- "Taxon"
  x$Taxon <- trimws(sub("\\s*\\(taxid [0-9]+\\)$", "", x$Taxon))
  m <- rowsum(as.matrix(x[, SAMPLES]), x$Taxon)
  m <- sweep(m, 2, colSums(m), "/") * 100
  m <- m[rowSums(m) > 0, , drop = FALSE]
  stopifnot(all(abs(colSums(m) - 100) < 1e-6))
  tibble(Taxon = rownames(m)) %>% bind_cols(as_tibble(m))
}
for (lv in c("phylum", "genus", "species")) wr(taxa_table(lv), RAW, "Microbial_community", paste0("bacteria_", lv, "_percent.tsv"))

# ---------------- ARGs (ARGs-OAP, TPM; Rank I / II risk table) ----------------
for (lv in c("type", "subtype", "gene")) {
  x <- rd(DS, "summary", "03_argoap", paste0("argoap_", lv, "_TPM.tsv"))
  names(x)[1] <- lv
  wr(x[, c(lv, SAMPLES)], RAW, "ARG_OAP", paste0("merged_tpm.", lv, ".txt"))
}
risk <- rd(DS, "summary", "03_argoap", "argoap_gene_TPM_with_risk_rank.tsv") %>%
  filter(Rank %in% c("I", "II")) %>%
  transmute(variant = gene, rank = Rank, Type = SARG_type, Subtype = SARG_subtype, Mechanism.group = Mechanism_group, across(all_of(SAMPLES)))
wr(risk, RAW, "ARG_OAP", "merged_tpm_gene_high_risk_rankI_II.txt")

# ---------------- Read numbers and removal of viral sequences ----------------
ren <- rd(DS, "sample_rename.tsv")
qc  <- rd(INFO, "revision", "qc_summary_all.tsv") %>%
  left_join(ren %>% select(raw_sample, corrected_sample), by = c("Sample" = "raw_sample")) %>%
  mutate(Sample = corrected_sample) %>% select(-corrected_sample)
vir <- rd(DS, "summary", "01_virus_removal", "virus_removal_per_sample.tsv")
qc  <- qc %>% left_join(vir %>% transmute(Sample = sample, Final_reads = 2 * final_pairs), by = "Sample")
stopifnot(!anyNA(qc$Sample), !anyNA(qc$Final_reads), setequal(qc$Sample, SAMPLES))
wr(qc, RAW, "Sequencing", "qc_summary_all.tsv")
wr(vir, RAW, "Sequencing", "virus_removal_per_sample.tsv")

# ---------------- MGEs (mobileOG-db, read-based) ----------------
MGE_DIR <- file.path(DS, "summary", "07_mge_reads")
mge_ann <- rd(MGE_DIR, "mobileOG_annotation_detected.tsv")[, c("id", "element", "category", "name")]
mge_tpm <- rd(MGE_DIR, "tpm.id.tsv")
names(mge_tpm)[1] <- "id"
stopifnot(setequal(SAMPLES, names(mge_tpm)[-1]), all(mge_tpm$id %in% mge_ann$id))
mge_ann <- mge_ann[match(mge_tpm$id, mge_ann$id), ]
el_parts <- strsplit(mge_ann$element, "+", fixed = TRUE)
is_phage <- vapply(el_parts, function(x) "phage" %in% x, logical(1))
mge_class <- vapply(el_parts, function(x) {
  if ("IS" %in% x) "Insertion Sequence (IS)"
  else if (any(x %in% c("ICE", "AICE", "IME", "CIME"))) "Integrative Element (IGE)"
  else if ("plasmid" %in% x) "Plasmid" else "Other"
}, character(1))
m <- as.matrix(mge_tpm[!is_phage, SAMPLES])
m <- sweep(m, 2, colSums(m), "/") * 1e6
keep <- rowSums(m) > 0
mge_np <- data.frame(id = mge_tpm$id[!is_phage][keep], sseqid_gene_name = mge_ann$name[!is_phage][keep],
                     `Major mobileOG Category` = mge_ann$category[!is_phage][keep], MGE_Class = mge_class[!is_phage][keep],
                     m[keep, ], check.names = FALSE)
stopifnot(all(abs(colSums(mge_np[, SAMPLES]) - 1e6) < 1e-3))
names(mge_np)[1] <- ""
wr(mge_np, RAW, "MGE", "mge_no_phage_TPM_matrix.tsv")
mge_type <- rowsum(as.matrix(mge_np[, SAMPLES]), mge_np$MGE_Class)[c("Insertion Sequence (IS)", "Integrative Element (IGE)", "Other", "Plasmid"), ]
wr(data.frame(sample = SAMPLES, t(mge_type)[SAMPLES, ], check.names = FALSE), RAW, "MGE", "mge_per_sample_type_TPM.tsv")

# ---------------- Sampling points ----------------
invisible(file.copy(file.path(INFO, "Sampling_map", "sampling_points.csv"), file.path(RAW, "Sampling_map"), overwrite = TRUE))
message("Inputs written to data/raw/")
