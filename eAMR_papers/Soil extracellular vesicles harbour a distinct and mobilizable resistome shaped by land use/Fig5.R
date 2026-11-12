# =====================================================================
#  Fig5.R - Figure 5: ARG-carrying metagenome-assembled genomes (MAGs)
# ---------------------------------------------------------------------
#  MAGs: MetaWRAP binning and refinement (completeness >= 50%, contamination <= 10%) of the virus-free contigs, dereplicated
#  at 95% ANI (dRep, 162 representatives), GTDB-Tk taxonomy (R226), relative abundance with CoverM, ARGs with SARG.
#  a  Phylogeny (GTDB-Tk de novo tree, pruned to the ARG-carrying MAGs), phylum, number of ARG subtypes per class and total
#  b  Mean relative abundance (% of reads) of these MAGs per land use and fraction (zeros included)
#  Also writes Supplementary Table S4 (ARG-carrying MAGs), Supplementary Data 11 (all MAGs) and 12 (ARGs in the
#  representative MAGs) and the numbers cited in the text.
#  Output: output/Fig5/
# =====================================================================
source("scripts/00_functions.R")
use(c("ggh4x", "patchwork", "ape"))
OUT <- out_folder("Fig5")
MAG_DIR <- file.path(DS, "summary", "05_mag")
COL_PHYLUM <- c(Pseudomonadota = "#7A6098", Bacteroidota = "#E0A030", Actinomycetota = "#A0522D",
                Patescibacteriota = "#2F4F4F", Acidobacteriota = "#C9A0DC", Thermoproteota = "#6B8E23",
                Methylomirabilota = "#8C8C8C", Planctomycetota = "#D2B48C", Chloroflexota = "#B0C4DE", Others = "grey80")
COL_ARG <- c(Multidrug = "#3F5E7D", Bacitracin = "#D49560", `Beta-lactam` = "#A05B6F", Polymyxin = "#5E8B7E",
             MLS = "#6E5A8C", Quinolone = "#A8C4B8", Tetracycline = "#8AAEC0", Sulfonamide = "#B8AE3F",
             Aminoglycoside = "#D4A8A8", Vancomycin = "#A88FBE", Chloramphenicol = "#7A8590", Fosfomycin = "#B85C8E",
             Rifamycin = "#A8916A", Novobiocin = "#2E6B73", Others = "grey70")
TYPE_MAP <- c(multidrug = "Multidrug", beta_lactam = "Beta-lactam", `macrolide-lincosamide-streptogramin` = "MLS",
              polymyxin = "Polymyxin", quinolone = "Quinolone", bacitracin = "Bacitracin", tetracycline = "Tetracycline",
              sulfonamide = "Sulfonamide", aminoglycoside = "Aminoglycoside", chloramphenicol = "Chloramphenicol",
              fosfomycin = "Fosfomycin", rifamycin = "Rifamycin", vancomycin = "Vancomycin", novobiocin = "Novobiocin")
AB_BREAKS <- c(-Inf, 0, 0.01, 0.1, 0.5, 1, Inf)
AB_LABELS <- c("0", "< 0.01", "0.01\u20130.1", "0.1\u20130.5", "0.5\u20131", "> 1")
AB_SIZES  <- setNames(c(0.6, 1.6, 2.8, 4.2, 5.6, 7.2), AB_LABELS)
CELL_ALPHA <- 0.55
theme_base <- theme_bw(base_size = 12, base_family = FONT) +
  theme(plot.title = element_blank(), panel.grid = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), legend.title = element_text(face = "bold", size = 10.5),
        legend.text = element_text(size = 10), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(4.5, "mm"), legend.position = "bottom", legend.box = "horizontal",
        strip.background = element_rect(colour = NA), plot.margin = margin(4, 2, 4, 2))
STRIP_TXT <- element_text(face = "bold", size = 8.5, colour = "black")
no_prefix <- function(x) sub("^haining_", "", x)        # project prefix of the MAG names on the server

# ---------------- Data ----------------
meta <- read.csv(META_FILE, stringsAsFactors = FALSE) %>%
  mutate(Habitat = factor(Habitat, levels = names(COL_HABITAT)), SampleType = factor(SampleType, levels = c("Soil", "EV")))
reps <- read.delim(file.path(MAG_DIR, "mag_representatives.tsv"), stringsAsFactors = FALSE, na.strings = c("", "NA")) %>%
  mutate(n_ARG_ORF_SARG = replace_na(n_ARG_ORF_SARG, 0), n_unique_ARG_SARG = replace_na(n_unique_ARG_SARG, 0),
         label = case_when(!is.na(species) ~ species, !is.na(genus) ~ paste(genus, "sp."), !is.na(family) ~ paste(family, "sp."),
                           !is.na(order) ~ paste(order, "sp."), TRUE ~ paste(phylum, "MAG"))) %>%
  group_by(label) %>% mutate(label = if (n() > 1) paste0(label, " (", source_sample, " ", bin, ")") else label) %>% ungroup()
argg  <- read.delim(file.path(MAG_DIR, "mag_arg_SARG_annotated.tsv"), stringsAsFactors = FALSE) %>% filter(MAG %in% reps$MAG)
relab <- read.delim(file.path(MAG_DIR, "mag_abundance_relative_abundance.tsv"), stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(MAG %in% reps$MAG)
ab_long <- relab %>% pivot_longer(-MAG, names_to = "SampleID", values_to = "relab") %>%
  inner_join(select(meta, SampleID, Habitat, SampleType), by = "SampleID")
ab_grp <- ab_long %>% group_by(MAG, Habitat, SampleType) %>%
  summarise(mean_relab = mean(relab), n_detected = sum(relab > 0), n_samples = n(), .groups = "drop") %>%
  mutate(bin = cut(mean_relab, breaks = AB_BREAKS, labels = AB_LABELS),
         col = ifelse(mean_relab == 0, "grey70", COL_HABITAT[as.character(Habitat)]),
         fill = ifelse(mean_relab == 0, "grey70", ifelse(SampleType == "Soil", COL_HABITAT[as.character(Habitat)], "white")))

# ---------------- Supplementary Data 11 and 12 ----------------
rep_cols <- reps %>% left_join(ab_long %>% group_by(MAG) %>%
                                 summarise(n_Soil_samples_detected = sum(relab > 0 & SampleType == "Soil"),
                                           n_EV_samples_detected = sum(relab > 0 & SampleType == "EV"), .groups = "drop"), by = "MAG") %>%
  transmute(representative = MAG, name = label, phylum, class, order, family, genus, species, found_in,
            n_Soil_samples_detected, n_EV_samples_detected, n_ARG_ORFs_SARG = n_ARG_ORF_SARG,
            n_ARG_subtypes_SARG = n_unique_ARG_SARG, ARG_genes_SARG, n_ARG_ORFs_CARD = n_ARG_ORF_CARD)
all_mags <- read.delim(file.path(MAG_DIR, "mag_all_refined.tsv"), stringsAsFactors = FALSE) %>%
  left_join(rep_cols, by = "representative") %>%
  transmute(MAG = no_prefix(MAG), source_sample = sample, land_use = Habitat, fraction = SampleType, completeness, contamination, GC,
            size_Mb = round(size_bp / 1e6, 2), N50, quality, checkm_lineage, cluster_representative = ifelse(is_representative, "yes", ""),
            name, phylum, class, order, family, genus, species, found_in, n_Soil_samples_detected, n_EV_samples_detected,
            n_ARG_ORFs_SARG, n_ARG_subtypes_SARG, ARG_genes_SARG, n_ARG_ORFs_CARD, rep = representative) %>%
  arrange(phylum, class, order, family, genus, rep, desc(cluster_representative)) %>% select(-rep)
save_table(all_mags, OUT, "Supplementary_Data_11_MAGs")
save_table(argg %>% left_join(select(reps, MAG, MAG_name = label, phylum, source_sample, Habitat, SampleType), by = "MAG") %>%
             transmute(MAG = no_prefix(MAG), MAG_name, phylum, source_sample, land_use = Habitat, source_fraction = SampleType,
                       ORF_ID = no_prefix(ORF_ID), ARG_type = Type, ARG_subtype = sub("^.*?__", "", Subtype, perl = TRUE), Mechanism_group,
                       identity_pct = pident, query_coverage_pct = qcovhsp, evalue) %>% arrange(MAG_name, ARG_type, ARG_subtype),
           OUT, "Supplementary_Data_12_MAG_ARGs")

# ---------------- a. ARG-carrying MAGs ----------------
arg_cnt <- argg %>% mutate(Class = unname(coalesce(TYPE_MAP[tolower(Type)], "Others"))) %>%
  group_by(MAG, Class) %>% summarise(n = n_distinct(Subtype), .groups = "drop")
mags <- reps %>% filter(n_ARG_ORF_SARG > 0) %>%
  left_join(arg_cnt %>% group_by(MAG) %>% summarise(total = sum(n), .groups = "drop"), by = "MAG") %>%
  arrange(desc(total), desc(n_ARG_ORF_SARG), label)
lab_vec <- setNames(mags$label, mags$MAG)

# rows follow the GTDB-Tk de novo tree (120 bacterial marker genes, FastTree) pruned to these MAGs
tree_a <- read.tree(file.path(MAG_DIR, "tree", "haining_bacteria_MAG.tree"))
stopifnot(all(mags$MAG %in% tree_a$tip.label))
tree_a <- ladderize(keep.tip(tree_a, mags$MAG), right = FALSE)
nt <- length(tree_a$tip.label)
tseq <- tree_a$edge[tree_a$edge[, 2] <= nt, 2]
nx <- node.depth.edgelength(tree_a); ny <- numeric(nt + tree_a$Nnode); ny[tseq] <- seq_len(nt)
po_a <- reorder(tree_a, "postorder")$edge
for (nd in unique(po_a[, 1])) ny[nd] <- mean(range(ny[po_a[po_a[, 1] == nd, 2]]))
MAG_LEVELS <- tree_a$tip.label[tseq]                       # bottom to top
tr_h <- tibble(x = nx[tree_a$edge[, 1]], xend = nx[tree_a$edge[, 2]], y = ny[tree_a$edge[, 2]])
tr_v <- tibble(parent = tree_a$edge[, 1], y = ny[tree_a$edge[, 2]]) %>% group_by(parent) %>%
  summarise(ymin = min(y), ymax = max(y), .groups = "drop") %>% mutate(x = nx[parent])
tr_tip <- tibble(x = nx[seq_len(nt)], y = ny[seq_len(nt)], x_tip = max(nx))
SCALE_LEN <- 0.2                                           # scale bar: substitutions per site
p_tree <- ggplot() +
  geom_segment(data = tr_tip, aes(x = x, xend = x_tip, y = y, yend = y), colour = "grey75", linewidth = 0.25, linetype = "dotted") +
  geom_segment(data = tr_h, aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.35) +
  geom_segment(data = tr_v, aes(x = x, xend = x, y = ymin, yend = ymax), linewidth = 0.35) +
  annotate("segment", x = 0, xend = SCALE_LEN, y = 0.72, yend = 0.72, linewidth = 0.5) +
  annotate("text", x = SCALE_LEN * 1.2, y = 0.72, label = SCALE_LEN, family = FONT, size = 2.7, hjust = 0, vjust = 0.5) +
  facet_wrap(~ "Phylogeny") + scale_y_continuous(limits = c(0.5, nt + 0.5), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.02))) + theme_base +
  theme(axis.text = element_blank(), panel.border = element_blank(), strip.background = element_blank(), strip.text.x = STRIP_TXT, strip.clip = "off")
mags$MAG_f <- factor(mags$MAG, levels = MAG_LEVELS)
class_order <- arg_cnt %>% group_by(Class) %>% summarise(n = sum(n)) %>% arrange(Class == "Others", desc(n)) %>% pull(Class)
grid_df <- expand_grid(MAG = mags$MAG, Class = class_order) %>% left_join(arg_cnt, by = c("MAG", "Class")) %>%
  mutate(MAG = factor(MAG, levels = MAG_LEVELS), Class = factor(Class, levels = class_order))
phy_in <- unique(mags$phylum); phy_cols <- COL_PHYLUM[intersect(names(COL_PHYLUM), phy_in)]
if (length(setdiff(phy_in, names(COL_PHYLUM)))) phy_cols[setdiff(phy_in, names(COL_PHYLUM))] <- "grey80"
p_phy <- ggplot(mutate(mags, hd = "Phylum"), aes(1, MAG_f, fill = phylum)) + geom_tile(width = 1, height = 0.9) +
  facet_wrap(~ hd) + scale_fill_manual(values = phy_cols, name = "Phylum", guide = guide_legend(order = 1, ncol = 1, title.position = "top")) +
  scale_y_discrete(labels = lab_vec, expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + theme_base +
  theme(axis.text.y = element_text(face = "italic", colour = "black", size = 10.5), axis.text.x = element_blank(),
        strip.text.x = STRIP_TXT, strip.background = element_blank(), strip.clip = "off", panel.border = element_blank())
p_heat <- ggplot(grid_df, aes(1, MAG)) +
  geom_tile(fill = "grey96", colour = "white", linewidth = 0.8) +
  geom_tile(data = filter(grid_df, !is.na(n)), aes(fill = Class), alpha = CELL_ALPHA, colour = "white", linewidth = 0.8) +
  geom_text(data = filter(grid_df, !is.na(n)), aes(label = n), family = FONT, fontface = "bold", size = 3.6) +
  facet_wrap2(~ Class, nrow = 1, strip = strip_themed(background_x = elem_list_rect(fill = alpha(COL_ARG[class_order], CELL_ALPHA), colour = NA))) +
  scale_fill_manual(values = COL_ARG[class_order], name = "ARG class", drop = FALSE,
                    guide = guide_legend(order = 2, nrow = 2, byrow = TRUE, title.position = "top", override.aes = list(alpha = CELL_ALPHA))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(subtitle = "ARG classes (number of ARG subtypes)") + theme_base +
  theme(axis.text = element_blank(), plot.subtitle = element_text(face = "bold", size = 10.5, hjust = 0.5, margin = margin(b = 3)),
        strip.text.x = STRIP_TXT, panel.spacing.x = unit(0, "lines"), panel.border = element_blank())
p_bar <- ggplot(mutate(mags, hd = "Total ARG count"), aes(total, MAG_f)) + geom_col(width = 0.62, fill = "grey45") +
  geom_text(aes(label = total), hjust = -0.35, family = FONT, fontface = "bold", size = 3.6) + facet_wrap(~ hd) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45))) + scale_y_discrete(expand = c(0, 0)) + theme_base +
  theme(axis.text = element_blank(), strip.text.x = STRIP_TXT, strip.background = element_rect(fill = "grey90", colour = NA),
        panel.border = element_blank())
# Supplementary Table S4: ARG-carrying MAGs and their ARG subtypes per class
save_table(grid_df %>% filter(!is.na(n)) %>% mutate(MAG_name = lab_vec[as.character(MAG)]) %>%
             pivot_wider(id_cols = c(MAG, MAG_name), names_from = Class, values_from = n, values_fill = 0) %>%
             left_join(select(mags, MAG, phylum, genus, species, source_sample, Habitat, SampleType, Total = total,
                              n_ARG_ORFs = n_ARG_ORF_SARG), by = "MAG") %>% arrange(desc(Total)) %>% mutate(MAG = no_prefix(MAG)),
           OUT, "Table_S4_ARG_carrying_MAGs")

# ---------------- b. Mean relative abundance of the ARG-carrying MAGs ----------------
ab <- ab_grp %>% filter(MAG %in% mags$MAG) %>% mutate(MAG = factor(MAG, levels = MAG_LEVELS))
save_table(ab %>% mutate(MAG_name = lab_vec[as.character(MAG)], MAG = no_prefix(as.character(MAG))) %>%
             select(MAG, MAG_name, Habitat, SampleType, mean_relab, n_detected, n_samples), OUT, "Fig5b_MAG_group_abundance")
strip_fill <- c(alpha(COL_HABITAT, 0.35), rep("grey95", 6))
p_bub <- ggplot(ab, aes(1, MAG)) +
  geom_point(aes(size = bin, colour = col, fill = fill), shape = 21, stroke = 0.8) +
  facet_nested(~ Habitat + SampleType, nest_line = FALSE, strip = strip_nested(background_x = elem_list_rect(fill = strip_fill, colour = NA))) +
  scale_colour_identity() + scale_fill_identity() +
  scale_size_manual(values = AB_SIZES, drop = FALSE, name = "Mean relative abundance (% of reads)",
                    guide = guide_legend(nrow = 1, label.position = "bottom", title.position = "top", order = 3,
                                         override.aes = list(colour = c("grey70", rep("grey25", 5)), fill = c("grey70", rep("grey25", 5))))) +
  scale_y_discrete(limits = MAG_LEVELS, labels = NULL, expand = c(0, 0.5)) + scale_x_continuous(expand = expansion(add = 0.5)) +
  theme_base +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.4), panel.spacing.x = unit(0.15, "lines"),
        strip.text.x = element_text(face = "bold", size = 10.5))

# ---------------- Numbers cited in the text ----------------
sum_car <- ab_long %>% filter(MAG %in% mags$MAG) %>% group_by(SampleID, Habitat, SampleType) %>% summarise(relab = sum(relab), .groups = "drop")
fa <- sum_car %>% filter(Habitat == "FA")
writeLines(c(
  sprintf("Representative MAGs: %d (high quality %d, medium %d); from EV assemblies %d, from bulk-soil assemblies %d",
          nrow(reps), sum(grepl("^high", reps$quality)), sum(grepl("^medium", reps$quality)), sum(reps$SampleType == "EV"), sum(reps$SampleType == "Soil")),
  sprintf("ARG-carrying representatives: %d with %d ARG ORFs; phyla: %s", nrow(mags), sum(mags$n_ARG_ORF_SARG),
          paste(names(table(mags$phylum)), table(mags$phylum), collapse = ", ")),
  sprintf("Pseudomonadota share of ARG subtype occurrences: %.0f%%", 100 * sum(mags$total[mags$phylum == "Pseudomonadota"]) / sum(mags$total)),
  sprintf("Summed relative abundance of ARG-carrying MAGs in farmland: median bulk soil %.3f%%, EV %.3f%% (14 pairs)",
          median(fa$relab[fa$SampleType == "Soil"]), median(fa$relab[fa$SampleType == "EV"]))),
  file.path(OUT, "Fig5_numbers_in_text.txt"))

# ---------------- Assemble (180 mm wide at print size; canvas x 2.0) ----------------
SCALE <- 2.0; PRINT_W <- 180; FIG_H <- 45 + 3.6 * nrow(mags)
p_bub <- p_bub + labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 20, family = FONT), plot.tag.position = c(-0.06, 1.0))
inner <- p_tree + p_phy + p_heat + p_bar + plot_spacer() + p_bub +
  plot_layout(widths = c(2.4, 0.9, length(class_order) * 1.6, 2.6, 0.5, 6.5), guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "left")
fig <- wrap_elements(full = inner) + labs(tag = "a") +
  theme(plot.tag = element_text(face = "bold", size = 20, family = FONT), plot.tag.position = c(0.01, 0.985))
save_fig(fig, OUT, "Fig5", w = PRINT_W * SCALE / 25.4, h = FIG_H * SCALE / 25.4)
message("Figure 5 written to ", OUT)
