# =====================================================================
#  FigS2.R - Supplementary Fig. S2: DNase treatment of vesicles, genomic DNA and plasmid DNA
# ---------------------------------------------------------------------
#  dsDNA concentration (Qubit dsDNA HS) before and after Turbo DNase treatment; genomic and plasmid DNA of E. coli pKJK5
#  in triplicate, vesicles (~10^8 particles) measured once. Data from Tskhay et al. (2026), entered below.
#  Output: output/FigS2/
# =====================================================================
source("scripts/00_functions.R")
use("ggtext")
OUT <- out_folder("FigS2")
COL_COND  <- c("-DNase" = "grey75", "+DNase" = "white")
LAB_GROUP <- c(Vesicles = "Vesicles (~10<sup>8</sup> particles)", Genomic = "Genomic DNA", Plasmid = "Plasmid DNA")
raw <- tribble(
  ~Group,     ~Condition, ~Value,
  "Vesicles", "-DNase",   1.26,
  "Vesicles", "+DNase",   0.26,
  "Genomic",  "-DNase",   4.50,  "Genomic", "-DNase", 4.24,  "Genomic", "-DNase", 4.90,
  "Genomic",  "+DNase",   0.252, "Genomic", "+DNase", 0.286, "Genomic", "+DNase", 0.256,
  "Plasmid",  "-DNase",   0.700, "Plasmid", "-DNase", 0.748, "Plasmid", "-DNase", 0.630,
  "Plasmid",  "+DNase",   0.170, "Plasmid", "+DNase", 0.138, "Plasmid", "+DNase", 0.134) %>%
  mutate(Group = factor(Group, levels = names(LAB_GROUP)), Condition = factor(Condition, levels = names(COL_COND)))
summ <- raw %>% group_by(Group, Condition) %>%
  summarise(Mean = mean(Value), SD = if (n() > 1) sd(Value) else NA_real_, n = n(), .groups = "drop")
save_table(summ, OUT, "FigS2_DNase_summary")
p <- ggplot(summ, aes(Condition, Mean)) +
  geom_col(aes(fill = Condition), colour = "black", width = 0.65, linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(Mean - SD, 0), ymax = Mean + SD), width = 0.25, linewidth = 0.5, na.rm = TRUE) +
  geom_point(data = raw, aes(Condition, Value), size = 1.6, colour = "black", position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  facet_wrap(~ Group, nrow = 1, labeller = as_labeller(LAB_GROUP)) +
  scale_fill_manual(values = COL_COND) +
  scale_y_continuous(limits = c(0, 5.5), breaks = seq(0, 5, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "DNA concentration (ng/\u00b5L)") +
  theme_bw(base_size = 12, base_family = FONT) +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.6),
        strip.text = element_markdown(face = "bold", size = 10.5), axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 11), axis.text.x = element_text(face = "bold"),
        axis.title = element_text(face = "bold", colour = "black"), legend.position = "none", plot.margin = margin(5, 5, 5, 5))
save_fig(patchwork::wrap_elements(full = p), OUT, "FigS2", w = 110 * 1.8 / 25.4, h = 60 * 1.8 / 25.4)
