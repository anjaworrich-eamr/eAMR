# =====================================================================
#  FigS10_S11.R - Supplementary Figs. S10 (EV) and S11 (bulk soil): gene maps of the ARG-MGE co-localized contigs
# ---------------------------------------------------------------------
#  Co-localized contigs as in Fig. 4 (read_coloc in 00_functions.R). Long contigs are shown from 2 kb before the first to
#  2 kb after the last ARG or core mobility gene. Figures with more than 20 contigs are drawn in three columns.
#  Output: output/FigS10_S11/
# =====================================================================
source("scripts/00_functions.R")
use(c("gggenes", "ggrepel", "cowplot", "scales"))
OUT <- out_folder("FigS10_S11")
co  <- read_coloc()
COL_MGE <- c(Plasmid = "#5B4B8A", IGE = "#A0522D", IS = "#E0A030", Others = "#4D4D4D")
COL_ARG <- c(Multidrug = "#3F5E7D", Tetracycline = "#8AAEC0", Polymyxin = "#5E8B7E", Quinolone = "#A8C4B8",
             Bacitracin = "#D49560", Novobiocin = "#2E6B73", Sulfonamide = "#B8AE3F", Rifamycin = "#A8916A",
             `Beta-lactam` = "#A05B6F", Aminoglycoside = "#D4A8A8", MLS = "#6E5A8C", Vancomycin = "#A88FBE",
             Chloramphenicol = "#7A8590", Trimethoprim = "#2E7D32", Fosfomycin = "#B85C8E", Others = "grey70")
CLASS_MAP <- c(multidrug = "Multidrug", tetracycline = "Tetracycline", polymyxin = "Polymyxin", quinolone = "Quinolone",
               bacitracin = "Bacitracin", novobiocin = "Novobiocin", sulfonamide = "Sulfonamide", rifamycin = "Rifamycin",
               beta_lactam = "Beta-lactam", aminoglycoside = "Aminoglycoside", `macrolide-lincosamide-streptogramin` = "MLS",
               vancomycin = "Vancomycin", chloramphenicol = "Chloramphenicol", trimethoprim = "Trimethoprim", fosfomycin = "Fosfomycin")
FLANK <- 2000; N_COL <- 2; MAX_PER_PAGE <- 20; FIG_WIDTH <- 12.5; ROW_HEIGHT <- 1.55

ctg <- co$arg_ctg %>% filter(C1)
gg  <- co$gg %>% filter(molecule %in% ctg$molecule) %>%
  mutate(type = case_when(is_ARG ~ "ARG", is_MGE & mge_major_category != "phage" ~ "MGE", TRUE ~ "Other"),
         unnamed = type == "MGE" & grepl(UNNAMED_RE, mge_name),
         core = type == "MGE" & (grepl(CORE_RE, mge_name) | (unnamed & mge_major_category %in% UNNAMED_CAT)),
         fill_group = case_when(type == "ARG" ~ unname(coalesce(CLASS_MAP[arg_type], "Others")),
                                unnamed ~ "Others", type == "MGE" ~ el_type(element_type), TRUE ~ "Other ORF"),
         label = case_when(type == "ARG" ~ sub("^.*?__", "", arg_subtype, perl = TRUE), unnamed ~ "unspecified",
                           type == "MGE" ~ sub(";$", "", mge_name), TRUE ~ ""),
         forward = strand == "forward") %>%
  left_join(select(co$meta, sample = SampleID, Habitat), by = "sample")
win <- gg %>% filter(type == "ARG" | core) %>% group_by(molecule) %>%
  summarise(lo = max(0, min(start) - FLANK), hi = min(contig_length[1], max(end) + FLANK), .groups = "drop")
gg <- gg %>% left_join(win, by = "molecule") %>% filter(end >= lo, start <= hi) %>% mutate(start = pmax(start, lo), end = pmin(end, hi))

make_one <- function(mol, d_all, cols) {
  d  <- filter(d_all, molecule == mol)
  dp <- filter(d, type != "Other")
  lab <- dp %>% mutate(x_mid = (start + end) / 2) %>% arrange(x_mid)
  tiers <- if (nrow(lab) >= 8) c(0.05, 0.11, 0.17) else if (nrow(lab) >= 4) c(0.05, 0.12) else 0.05
  lab$nudge <- tiers[(seq_len(nrow(lab)) - 1) %% length(tiers) + 1]
  ttl <- sprintf("%s  |  %s (%s bp%s)", d$sample[1], sub("_length_.*", "", d$contig_id[1]), scales::comma(d$contig_length[1]),
                 ifelse(d$lo[1] > 0 | d$hi[1] < d$contig_length[1], sprintf("; shown %s-%s", scales::comma(d$lo[1]), scales::comma(d$hi[1])), ""))
  ggplot() +
    annotate("segment", x = d$lo[1], xend = d$hi[1], y = 1, yend = 1, colour = "black", linewidth = 0.7) +
    geom_gene_arrow(data = filter(d, type == "Other"), aes(xmin = start, xmax = end, y = 1, forward = forward),
                    fill = "white", colour = "grey60", linewidth = 0.3, arrowhead_height = unit(5, "mm"),
                    arrowhead_width = unit(2.5, "mm"), arrow_body_height = unit(4, "mm")) +
    geom_gene_arrow(data = dp, aes(xmin = start, xmax = end, y = 1, fill = fill_group, forward = forward), colour = NA,
                    arrowhead_height = unit(7, "mm"), arrowhead_width = unit(3, "mm"), arrow_body_height = unit(6, "mm")) +
    geom_text_repel(data = lab, aes(x = x_mid, y = 1, label = label), nudge_y = lab$nudge, size = 3.3, fontface = "italic",
                    family = FONT, colour = "grey10", segment.size = 0.25, segment.colour = "grey60", min.segment.length = 0,
                    box.padding = 0.35, force = 3, max.overlaps = Inf, seed = 1) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_x_continuous(labels = scales::comma, limits = c(d$lo[1], d$hi[1]), expand = expansion(mult = c(0.01, 0.03))) +
    scale_y_continuous(limits = c(0.96, max(lab$nudge) + 1.06), expand = c(0, 0)) +
    labs(title = ttl, x = NULL, y = NULL) +
    theme_classic(base_size = 10, base_family = FONT) +
    theme(plot.title = element_text(size = 10.5, face = "bold", colour = COL_HABITAT[[d$Habitat[1]]], margin = margin(b = 4)),
          aspect.ratio = 0.17, axis.text.x = element_text(size = 9.5, colour = "black"), axis.line.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",
          plot.margin = margin(t = 10, r = 4, b = 2, l = 4))
}
legend_grob <- function(cats, cols, ttl, ncol) {
  d <- data.frame(x = seq_along(cats), grp = factor(cats, levels = cats))
  p <- ggplot(d, aes(x, 1, fill = grp)) + geom_tile() +
    scale_fill_manual(values = cols[cats], name = ttl) +
    guides(fill = guide_legend(title.position = "top", ncol = ncol, byrow = TRUE)) +
    theme_void(base_family = FONT) +
    theme(legend.position = "left", legend.title = element_text(size = 11.5, face = "bold"),
          legend.text = element_text(size = 10), legend.key.size = unit(4.5, "mm"))
  g <- ggplotGrob(p); g$grobs[[which(g$layout$name %in% c("guide-box-left", "guide-box"))[1]]]
}

for (st in c("EV", "Soil")) {
  d <- filter(gg, SampleType == st)
  ord <- ctg %>% filter(SampleType == st) %>% mutate(Habitat = factor(Habitat, levels = c("UR", "FO", "FA"))) %>%
    arrange(Habitat, sample, desc(length)) %>% pull(molecule)
  arg_cls <- sort(unique(d$fill_group[d$type == "ARG"])); mge_cls <- intersect(names(COL_MGE), unique(d$fill_group[d$type == "MGE"]))
  cols  <- c(COL_MGE[mge_cls], COL_ARG[arg_cls])
  plots <- lapply(ord, make_one, d_all = d, cols = cols)
  leg   <- plot_grid(legend_grob(mge_cls, cols, "MGE type", 1), legend_grob(arg_cls, cols, "ARG class", 2), nrow = 1, rel_widths = c(1, 2))
  n_col <- if (length(plots) > MAX_PER_PAGE) 3 else N_COL
  items <- c(plots[1], list(leg), plots[-1])
  fig <- plot_grid(plotlist = items, ncol = n_col) + theme(plot.background = element_rect(fill = "white", colour = NA))
  name <- if (st == "EV") "FigS10_EV" else "FigS11_soil"
  save_fig(fig, OUT, name, w = FIG_WIDTH * n_col / N_COL, h = ceiling(length(items) / n_col) * ROW_HEIGHT)
  message(sprintf("%s: %d contigs", st, length(ord)))
}
