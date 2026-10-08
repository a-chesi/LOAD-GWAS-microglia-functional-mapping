# Figure 2C luciferase analysis and plot
# Input: Figure2C_luciferase_input_data.csv
# Columns: Condition, FoldChange

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "Figure2C_luciferase_input_data.csv"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

out_dir <- dirname(normalizePath(input_file))

required_packages <- c("ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Please install required R packages: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)
library(dplyr)

dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("Condition", "FoldChange")
if (!all(required_cols %in% names(dat))) {
  stop("Input must contain columns: ", paste(required_cols, collapse = ", "))
}

condition_order <- c("Empty", "PromOnly", "Risk", "Protective")
if (!all(condition_order %in% unique(dat$Condition))) {
  stop("Expected conditions: ", paste(condition_order, collapse = ", "))
}

dat <- dat %>%
  filter(Condition %in% condition_order) %>%
  mutate(Condition = factor(Condition, levels = condition_order))

# One-way ANOVA followed by Tukey HSD, matching the manuscript analysis.
fit <- aov(FoldChange ~ Condition, data = dat)
anova_tab <- summary(fit)[[1]]
anova_out <- data.frame(
  Test = "One-way ANOVA",
  Df_between = anova_tab["Condition", "Df"],
  Df_within = anova_tab["Residuals", "Df"],
  F = anova_tab["Condition", "F value"],
  P_value = anova_tab["Condition", "Pr(>F)"]
)

tukey <- TukeyHSD(fit, "Condition")$Condition
tukey_out <- data.frame(
  Comparison = rownames(tukey),
  Difference = tukey[, "diff"],
  CI_low = tukey[, "lwr"],
  CI_high = tukey[, "upr"],
  P_adjusted = tukey[, "p adj"],
  row.names = NULL
)

summary_out <- dat %>%
  group_by(Condition) %>%
  summarise(
    N = n(),
    Mean = mean(FoldChange),
    SD = sd(FoldChange),
    SEM = SD / sqrt(N),
    .groups = "drop"
  )

write.csv(summary_out, file.path(out_dir, "Figure_2C_luciferase_summary.csv"), row.names = FALSE)
write.csv(anova_out, file.path(out_dir, "Figure_2C_luciferase_ANOVA.csv"), row.names = FALSE)
write.csv(tukey_out, file.path(out_dir, "Figure_2C_luciferase_TukeyHSD.csv"), row.names = FALSE)

# Manuscript-style bar plot: mean +/- SD, with the three comparisons
# displayed in Figure 2C annotated.
plot_dat <- summary_out

# Extract adjusted P values for the three comparisons displayed in the manuscript.
get_p <- function(g1, g2) {
  candidates <- c(paste0(g2, "-", g1), paste0(g1, "-", g2))
  hit <- which(tukey_out$Comparison %in% candidates)
  if (length(hit) != 1) return(NA_real_)
  tukey_out$P_adjusted[hit]
}

p_prom_risk <- get_p("PromOnly", "Risk")
p_prom_prot <- get_p("PromOnly", "Protective")
p_risk_prot <- get_p("Risk", "Protective")

p_to_label <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.05) return("*")
  "ns"
}

p <- ggplot(plot_dat, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_col(width = 0.72, color = NA) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15, linewidth = 0.6) +
  scale_fill_manual(values = c(Empty = "grey80", PromOnly = "#F9DF62", Risk = "#62BCE0", Protective = "#F47A5D")) +
  scale_x_discrete(labels = c(
    Empty = "Empty",
    PromOnly = "PromOnly",
    Risk = "Risk",
    Protective = "Protective"
  )) +
  labs(x = "Condition", y = "Fold Change to Promoter Only") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# Brackets are drawn manually to avoid an additional plotting dependency.
y1 <- max(plot_dat$Mean + plot_dat$SD) * 1.08
y2 <- max(plot_dat$Mean + plot_dat$SD) * 1.24
y3 <- max(plot_dat$Mean + plot_dat$SD) * 1.40
bracket <- function(x1, x2, y, label) {
  list(
    annotate("segment", x = x1, xend = x2, y = y, yend = y),
    annotate("segment", x = x1, xend = x1, y = y, yend = y - 0.05),
    annotate("segment", x = x2, xend = x2, y = y, yend = y - 0.05),
    annotate("text", x = (x1 + x2)/2, y = y + 0.04, label = label, size = 3.5)
  )
}

p <- p + bracket(2, 3, y1, p_to_label(p_prom_risk)) +
         bracket(2, 4, y2, p_to_label(p_prom_prot)) +
         bracket(3, 4, y3, p_to_label(p_risk_prot)) +
         coord_cartesian(ylim = c(0, y3 + 0.18), clip = "off")

print(p)

ggsave(file.path(out_dir, "Figure_2C_luciferase.pdf"), p, width = 4.0, height = 3.1)
ggsave(file.path(out_dir, "Figure_2C_luciferase.png"), p, width = 4.0, height = 3.1, dpi = 300)

cat("\nFigure 2C analysis complete.\n")
cat("ANOVA P =", format(anova_out$P_value, digits = 5), "\n")
cat("PromOnly vs Risk Tukey-adjusted P =", format(p_prom_risk, digits = 5), "\n")
cat("PromOnly vs Protective Tukey-adjusted P =", format(p_prom_prot, digits = 5), "\n")
cat("Risk vs Protective Tukey-adjusted P =", format(p_risk_prot, digits = 5), "\n")
cat("Outputs saved to:", out_dir, "\n")
