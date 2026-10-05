# Figure 2D: RTFDC1 genotype-expression analysis
# ------------------------------------------------
# This script reproduces the genotype-expression analysis and Figure 2D.
#
# Expected input:
#   An Excel file with at least two columns:
#     Genotype   (values used here: "GG" and "GA")
#     TPM        (RTFDC1 expression in TPM)
#
# The individual-level FreshMicro input data are not redistributed with this
# repository. Authorized users should obtain the source data through the
# original data-access mechanism described in the manuscript.
#
# Usage:
#   Rscript Fig2D_genotype_expression_analysis_and_plot.R
#
# Optional:
#   Rscript Fig2D_genotype_expression_analysis_and_plot.R path/to/input.xlsx
#
# If no path is supplied, the script looks for:
#   RTFDC1_genotype_expression_input.xlsx

required_pkgs <- c("readxl", "dplyr", "ggplot2", "effectsize", "pwr")
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
library(effectsize)
library(pwr)

# ------------------------------------------------------------
# 1. Input
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) {
  args[1]
} else {
  "RTFDC1_genotype_expression_input.xlsx"
}

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nProvide the path as the first command-line argument or place ",
    "'RTFDC1_genotype_expression_input.xlsx' in the working directory."
  )
}

input_file <- normalizePath(input_file, mustWork = TRUE)
output_dir <- dirname(input_file)

cat("\nInput file: ", input_file, "\n", sep = "")
cat("Output directory: ", output_dir, "\n", sep = "")

dat <- read_excel(input_file, sheet = 1)

required_cols <- c("Genotype", "TPM")
missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop(
    "Input file is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

dat <- dat %>%
  select(Genotype, TPM) %>%
  filter(!is.na(Genotype), !is.na(TPM))

# ------------------------------------------------------------
# 2. Restrict to GG and GA
# ------------------------------------------------------------

GG <- dat %>%
  filter(Genotype == "GG") %>%
  pull(TPM)

GA <- dat %>%
  filter(Genotype == "GA") %>%
  pull(TPM)

if (length(GG) < 2 || length(GA) < 2) {
  stop("At least two GG and two GA observations are required.")
}

# ------------------------------------------------------------
# 3. Single-pass +/-3 SD filtering within genotype
# ------------------------------------------------------------

GG_lower <- mean(GG) - 3 * sd(GG)
GG_upper <- mean(GG) + 3 * sd(GG)

GA_lower <- mean(GA) - 3 * sd(GA)
GA_upper <- mean(GA) + 3 * sd(GA)

GG_used <- GG[GG >= GG_lower & GG <= GG_upper]
GA_used <- GA[GA >= GA_lower & GA <= GA_upper]

cat("\n--- SAMPLE COUNTS ---\n")
cat("GG before filtering:", length(GG), "\n")
cat("GA before filtering:", length(GA), "\n")
cat("GG after filtering: ", length(GG_used), "\n")
cat("GA after filtering: ", length(GA_used), "\n")

cat("\n--- REMOVED VALUES ---\n")
cat("GG removed:\n")
print(GG[!(GG %in% GG_used)])
cat("GA removed:\n")
print(GA[!(GA %in% GA_used)])

# ------------------------------------------------------------
# 4. Descriptive statistics
# ------------------------------------------------------------

mean_GG <- mean(GG_used)
mean_GA <- mean(GA_used)
sd_GG <- sd(GG_used)
sd_GA <- sd(GA_used)
mean_diff <- mean_GG - mean_GA

cat("\n--- DESCRIPTIVE STATISTICS ---\n")
cat("GG n:", length(GG_used), "\n")
cat("GG mean:", mean_GG, "\n")
cat("GG SD:", sd_GG, "\n")
cat("GA n:", length(GA_used), "\n")
cat("GA mean:", mean_GA, "\n")
cat("GA SD:", sd_GA, "\n")
cat("Mean difference (GG - GA):", mean_diff, "\n")

# ------------------------------------------------------------
# 5. Two-tailed Student's t-test (equal variances)
# ------------------------------------------------------------

tt <- t.test(
  GG_used,
  GA_used,
  alternative = "two.sided",
  var.equal = TRUE,
  conf.level = 0.95
)

cat("\n--- TWO-TAILED STUDENT'S T-TEST ---\n")
print(tt)

# ------------------------------------------------------------
# 6. Cohen's d with 95% CI
# ------------------------------------------------------------

d_result <- effectsize::cohens_d(
  GG_used,
  GA_used,
  pooled_sd = TRUE,
  ci = 0.95
)

cat("\n--- COHEN'S d ---\n")
print(d_result)

# ------------------------------------------------------------
# 7. Post-hoc power based on observed effect size
# ------------------------------------------------------------

pooled_sd <- sqrt(
  ((length(GG_used) - 1) * var(GG_used) +
     (length(GA_used) - 1) * var(GA_used)) /
    (length(GG_used) + length(GA_used) - 2)
)

d_observed <- abs((mean_GG - mean_GA) / pooled_sd)

power_result <- pwr::pwr.t2n.test(
  n1 = length(GG_used),
  n2 = length(GA_used),
  d = d_observed,
  sig.level = 0.05,
  alternative = "two.sided"
)

cat("\n--- POST-HOC POWER ---\n")
print(power_result)

# ------------------------------------------------------------
# 8. Save concise statistics table
# ------------------------------------------------------------

d_value <- unname(d_result$Cohens_d[1])
d_ci_low <- unname(d_result$CI_low[1])
d_ci_high <- unname(d_result$CI_high[1])

stats_out <- data.frame(
  n_GG = length(GG_used),
  n_GA = length(GA_used),
  mean_GG = mean_GG,
  mean_GA = mean_GA,
  mean_difference_GG_minus_GA = mean_diff,
  t_test_p = tt$p.value,
  mean_difference_CI_low = unname(tt$conf.int[1]),
  mean_difference_CI_high = unname(tt$conf.int[2]),
  cohens_d = d_value,
  cohens_d_CI_low = d_ci_low,
  cohens_d_CI_high = d_ci_high,
  posthoc_power = power_result$power
)

stats_file <- file.path(output_dir, "Figure_2D_genotype_expression_statistics.csv")
write.csv(
  stats_out,
  stats_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Figure 2D
# ------------------------------------------------------------

plot_df <- data.frame(
  Genotype = factor(
    c(rep("GG", length(GG_used)), rep("GA", length(GA_used))),
    levels = c("GG", "GA")
  ),
  RTFDC1 = c(GG_used, GA_used)
)

p_label <- paste0(
  "P = ",
  format.pval(tt$p.value, digits = 2, eps = 0.001)
)

p <- ggplot(plot_df, aes(x = Genotype, y = RTFDC1, fill = Genotype)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.7,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.7,
    alpha = 0.7,
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4,
    color = "#C0392B"
  ) +
  scale_fill_manual(
    values = c(
      "GG" = "#BFD7EA",
      "GA" = "#F2D6A2"
    )
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(plot_df$RTFDC1) + 10,
    label = p_label,
    size = 4
  ) +
  labs(
    x = "Genotype",
    y = expression(italic("RTFDC1") ~ "(TPM)")
  ) +
  coord_cartesian(
    ylim = c(0, max(plot_df$RTFDC1) + 18)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7)
  )

# Display plot in RStudio / an interactive R session
print(p)

pdf_file <- file.path(output_dir, "Figure_2D_genotype_expression.pdf")
png_file <- file.path(output_dir, "Figure_2D_genotype_expression.png")

ggsave(
  pdf_file,
  plot = p,
  width = 4.5,
  height = 4.5,
  units = "in"
)

ggsave(
  png_file,
  plot = p,
  width = 4.5,
  height = 4.5,
  units = "in",
  dpi = 600
)

cat("\nSaved:\n")
cat("  ", stats_file, "\n", sep = "")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")
