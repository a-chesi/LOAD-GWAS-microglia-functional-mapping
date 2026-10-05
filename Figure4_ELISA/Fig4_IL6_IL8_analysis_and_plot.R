# Figure 4: IL-6 and IL-8 secretion after enhancer deletion
# ----------------------------------------------------------
# This script reproduces the statistical analyses underlying Figure 4.
#
# Input workbook:
#   Figure4_ELISA_input_data.xlsx
#
# Required sheet:
#   ELISA_data
#
# Each row is one biological-replicate ELISA concentration after averaging
# the technical triplicate measurements. Biological replicates are averaged
# within each clone, and clone means are used as the statistical analysis unit
# (3 Mock clones and 3 enhancer-KO clones per cytokine).
#
# Usage:
#   Rscript Fig4_IL6_IL8_analysis_and_plot.R
# or
#   Rscript Fig4_IL6_IL8_analysis_and_plot.R path/to/Figure4_ELISA_input_data.xlsx

required_pkgs <- c("readxl", "dplyr", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(readxl)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Input
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "Figure4_ELISA_input_data.xlsx"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

input_file <- normalizePath(input_file, mustWork = TRUE)
output_dir <- dirname(input_file)

cat("\nInput file: ", input_file, "\n", sep = "")
cat("Output directory: ", output_dir, "\n", sep = "")

dat <- read_excel(input_file, sheet = "ELISA_data")
required_cols <- c(
  "Cytokine", "Genotype", "Clone",
  "Biological_replicate", "Concentration_pg_ml"
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Input file is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

dat <- dat %>%
  select(all_of(required_cols)) %>%
  filter(!is.na(Cytokine), !is.na(Genotype), !is.na(Clone), !is.na(Concentration_pg_ml))

# ------------------------------------------------------------
# 2. Average biological replicates within each clone
# ------------------------------------------------------------

clone_means <- dat %>%
  group_by(Cytokine, Genotype, Clone) %>%
  summarise(
    n_biological_replicates = n(),
    Concentration_pg_ml = mean(Concentration_pg_ml, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Genotype = factor(Genotype, levels = c("Mock", "Enhancer KO")))

if (any(clone_means$n_biological_replicates != 3)) {
  warning("At least one clone does not contain exactly 3 biological replicates.")
}

# ------------------------------------------------------------
# 3. Two-tailed Student's t-tests on clone means
# ------------------------------------------------------------

run_test <- function(cytokine_name) {
  x <- clone_means %>% filter(Cytokine == cytokine_name)
  mock <- x %>% filter(Genotype == "Mock") %>% pull(Concentration_pg_ml)
  ko <- x %>% filter(Genotype == "Enhancer KO") %>% pull(Concentration_pg_ml)

  if (length(mock) != 3 || length(ko) != 3) {
    stop(cytokine_name, ": expected 3 clone means per genotype.")
  }

  tt <- t.test(
    mock, ko,
    alternative = "two.sided",
    var.equal = TRUE,
    conf.level = 0.95
  )

  data.frame(
    Cytokine = cytokine_name,
    n_clones_Mock = length(mock),
    n_clones_Enhancer_KO = length(ko),
    mean_Mock = mean(mock),
    mean_Enhancer_KO = mean(ko),
    mean_difference_KO_minus_Mock = mean(ko) - mean(mock),
    t = unname(tt$statistic),
    df = unname(tt$parameter),
    p_value = tt$p.value
  )
}

stats_out <- bind_rows(run_test("IL-6"), run_test("IL-8"))
cat("\n--- FIGURE 4 PRIMARY STATISTICS ---\n")
print(stats_out)

stats_file <- file.path(output_dir, "Figure_4_IL6_IL8_statistics.csv")
clone_means_file <- file.path(output_dir, "Figure_4_clone_means.csv")

write.csv(stats_out, stats_file, row.names = FALSE)
write.csv(clone_means, clone_means_file, row.names = FALSE)

cat("\nSaved statistics: ", stats_file, "\n", sep = "")
cat("Saved clone means: ", clone_means_file, "\n", sep = "")

# ------------------------------------------------------------
# 4. Generate Figure 4 panels from clone-level values
# ------------------------------------------------------------

make_panel <- function(cytokine_name, p_value, fill_color) {
  plot_dat <- clone_means %>% filter(Cytokine == cytokine_name)
  significance_label <- if (p_value < 0.05) "*" else "n.s."
  y_max <- max(plot_dat$Concentration_pg_ml, na.rm = TRUE)

  ggplot(plot_dat, aes(x = Genotype, y = Concentration_pg_ml, fill = Genotype)) +
    geom_boxplot(
      width = 0.55,
      linewidth = 0.7,
      alpha = 0.75
    ) +
    scale_fill_manual(
      values = c(
        "Mock" = fill_color,
        "Enhancer KO" = fill_color
      )
    ) +
    scale_x_discrete(labels = c("Mock" = "Mock", "Enhancer KO" = "KO")) +
    annotate(
      "text",
      x = 1.5,
      y = y_max * 1.10,
      label = significance_label,
      size = 4
    ) +
    labs(
      title = cytokine_name,
      x = "Genotype",
      y = paste0(cytokine_name, " concentration (pg/ml)")
    ) +
    coord_cartesian(ylim = c(0, y_max * 1.20)) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
}

# Build plots
p_IL6 <- make_panel(
  "IL-6",
  stats_out$p_value[stats_out$Cytokine == "IL-6"],
  "#FC4E07"
)

p_IL8 <- make_panel(
  "IL-8",
  stats_out$p_value[stats_out$Cytokine == "IL-8"],
  "cyan"
)

# Display plots in RStudio / an interactive R session
print(p_IL6)
print(p_IL8)

# Save plots
# (The filenames are explicit so no separate output_stub object is needed.)
ggsave(
  file.path(output_dir, "Figure_4A_IL6.pdf"),
  plot = p_IL6,
  width = 4,
  height = 4,
  units = "in"
)
ggsave(
  file.path(output_dir, "Figure_4A_IL6.png"),
  plot = p_IL6,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(output_dir, "Figure_4B_IL8.pdf"),
  plot = p_IL8,
  width = 4,
  height = 4,
  units = "in"
)
ggsave(
  file.path(output_dir, "Figure_4B_IL8.png"),
  plot = p_IL8,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)

cat("\nExpected manuscript P values:\n")
cat("  IL-6: ", format.pval(stats_out$p_value[stats_out$Cytokine == "IL-6"], digits = 2), "\n", sep = "")
cat("  IL-8: ", format.pval(stats_out$p_value[stats_out$Cytokine == "IL-8"], digits = 2), "\n", sep = "")

cat("\nAll output files were written to: ", output_dir, "\n", sep = "")
