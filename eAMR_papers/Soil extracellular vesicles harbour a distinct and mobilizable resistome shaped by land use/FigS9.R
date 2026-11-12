# =====================================================================
#  FigS9.R - Supplementary Fig. S9: overall MGE pool (read-based, mobileOG-db, phage genes removed)
# ---------------------------------------------------------------------
#  a  Mean share of MGE types per land use and fraction (types from the mobileOG-db source databases)
#  b  PCoA of MGE gene composition (Bray-Curtis; TPM summed per gene name, entries without a name removed) with paired PERMANOVA
#  Also writes Supplementary Data 9 (MGE type shares per sample and the share of phage genes removed).
#  Output: output/FigS9/
# =====================================================================
source("scripts/00_functions.R")
use(c("ggh4x", "patchwork"))
OUT  <- out_folder("FigS9")
meta <- read_meta()
MGE_LEVELS <- c("Plasmid", "IS", "IGE", "Others")
MGE_LABELS <- c(Plasmid = "Plasmid-associated genes", IS = "Insertion sequences (IS)", IGE = "Integrative elements (IGE)", Others = "Others")
COL_MGE    <- c(Plasmid = "#5B4B8A", IS = "#E0A030", IGE = "#A0522D", Others = "grey78")
TYPE_MAP   <- c("Plasmid" = "Plasmid", "Insertion Sequence (IS)" = "IS", "Integrative Element (IGE)" = "IGE", "Other" = "Others")

# ---------------- a. MGE types ----------------
share <- read.delim(MGE_TYPE_FILE, check.names = FALSE) %>% pivot_longer(-sample, names_to = "Class", values_to = "TPM") %>%
  mutate(Type = factor(TYPE_MAP[Class], levels = MGE_LEVELS)) %>%
  group_by(sample) %>% mutate(pct = 100 * TPM / sum(TPM)) %>% ungroup() %>%
  left_join(select(meta, sample = SampleID, Habitat, SampleType, Pair), by = "sample")
stopifnot(!anyNA(share$Type), !anyNA(share$Habitat))

# Supplementary Data 9: share of each MGE type per sample, and the share of TPM from phage databases that was removed
el <- read.delim(file.path(DS, "summary", "07_mge_reads", "tpm.element.tsv"), row.names = 1, check.names = FALSE)
phage <- colSums(el[vapply(strsplit(rownames(el), "+", fixed = TRUE), function(x) "phage" %in% x, logical(1)), , drop = FALSE]) / colSums(el) * 100
save_table(share %>% select(sample, Habitat, SampleType, Type, pct) %>% pivot_wider(names_from = Type, values_from = pct) %>%
             mutate(phage_genes_removed_pct = unname(phage[sample])) %>%
             select(sample, Habitat, SampleType, phage_genes_removed_pct, all_of(MGE_LEVELS)) %>% arrange(SampleType, Habitat, sample),
           OUT, "Supplementary_Data_9_MGE_pool")

tab_a <- share %>% group_by(Habitat, SampleType, Type) %>%
  summarise(n_samples = n(), median_pct = median(pct), pct = mean(pct), .groups = "drop")
save_table(tab_a, OUT, "FigS9a_MGE_type_share_by_group")
centre <- tab_a %>% distinct(Habitat, SampleType, n_samples) %>% mutate(label = sprintf("n = %d", n_samples))
p_a <- tab_a %>% arrange(Habitat, SampleType, Type) %>% group_by(Habitat, SampleType) %>%
  mutate(ymax = pmin(cumsum(pct), 100), ymin = ymax - pct, ymid = (ymin + ymax) / 2) %>% ungroup() %>%
  ggplot() +
  geom_rect(aes(xmin = 2, xmax = 3.3, ymin = ymin, ymax = ymax, fill = Type), colour = "white", linewidth = 0.5) +
  geom_text(data = function(d) filter(d, pct >= 4), aes(x = 3.95, y = ymid, label = sprintf("%.0f%%", pct), colour = Type),
            size = 3.3, fontface = "bold", family = FONT, show.legend = FALSE) +
  geom_text(data = centre, aes(x = 0, y = 0, label = label), size = 3.4, family = FONT) +
  scale_fill_manual(values = COL_MGE, labels = MGE_LABELS, name = "MGE type", drop = FALSE) +
  scale_colour_manual(values = c(COL_MGE[1:3], Others = "grey45"), guide = "none") +
  scale_x_continuous(limits = c(0, 4.3), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
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

# ---------------- b. PCoA of MGE gene composition ----------------
gene_raw <- read.delim(MGE_GENE_FILE, check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- intersect(colnames(gene_raw), meta$SampleID)
gene_raw <- gene_raw[gene_raw$sseqid_gene_name != "NA:Keyword", ]
gene_tab <- rowsum(as.matrix(gene_raw[, gene_ids]), group = gene_raw$sseqid_gene_name)
gene_tab <- as.data.frame(sweep(gene_tab, 2, colSums(gene_tab), "/") * 1e6)
block_b  <- pcoa_panel(gene_tab, OUT, "FigS9b_MGE_gene")

SCALE <- 1.8
fig <- tag_panel(p_a, "a") + tag_panel(block_b, "b") + plot_layout(widths = c(0.56, 0.44))
save_fig(fig, OUT, "FigS9", w = 180 * SCALE / 25.4, h = 92 * SCALE / 25.4)
