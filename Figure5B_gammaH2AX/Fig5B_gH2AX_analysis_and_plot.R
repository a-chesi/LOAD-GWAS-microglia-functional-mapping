# Figure 5B: gamma-H2AX recovery analysis and plotting
#
# Input:
#   well_averages.xlsx (Sheet1)
#
# Usage:
#   Rscript Fig5B_gH2AX_analysis_and_plot.R
# or:
#   Rscript Fig5B_gH2AX_analysis_and_plot.R path/to/well_averages.xlsx
#
# This script reproduces the statistical analyses and quantitative plot
# for Figure 5B from well-level mean CTCFI values.
#
# Statistical approach:
#   1) Enhancer KO vs Mock:
#      - average technical replicate wells within clone x experimental day x timepoint
#      - log2-transform CTCFI
#      - linear mixed-effects model:
#          log2(CTCFI) ~ Group * Timepoint + experimental day + (1 | clone)
#   2) RTFDC1 overexpression vs LentiMock:
#      - same approach as above
#   3) siRTFDC1 vs NT:
#      - average replicate wells within independent experiment x condition x timepoint
#      - calculate log2 recovery from 0 h to 24 h
#      - paired t-test across independent experiments
#
# The plotted values are raw CTCFI; P values are calculated from the
# log2-transformed analyses described above.

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "lme4",
  "lmerTest",
  "emmeans"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)

# -----------------------------------------------------------------------------
# 1. Read and validate input data
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "well_averages.xlsx"

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nProvide the path as the first command-line argument or place ",
    "'well_averages.xlsx' in the working directory."
  )
}

input_file <- normalizePath(input_file, mustWork = TRUE)
output_dir <- dirname(input_file)

cat("\nInput file: ", input_file, "\n", sep = "")
cat("Output directory: ", output_dir, "\n", sep = "")

dat <- read_excel(input_file, sheet = "Sheet1")

required_cols <- c(
  "Timepoint", "Group", "CondClone", "Clone", "Rep",
  "Well", "WellID", "mean_CTCFI", "n_nuclei"
)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

expected_groups <- c("Mock", "KO", "LentiMock", "OE", "NT", "siRTFDC1")
missing_groups <- setdiff(expected_groups, unique(dat$Group))
if (length(missing_groups) > 0) {
  stop("Missing expected groups: ", paste(missing_groups, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 2. Enhancer KO vs Mock: primary mixed-effects analysis
# -----------------------------------------------------------------------------

ko_dat <- dat %>%
  filter(
    Group %in% c("Mock", "KO"),
    Timepoint %in% c("Hour0", "Hour24")
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Mock", "KO")),
    Timepoint = factor(Timepoint, levels = c("Hour0", "Hour24")),
    Rep = factor(Rep),
    CondClone = factor(CondClone)
  ) %>%
  group_by(Group, CondClone, Rep, Timepoint) %>%
  summarise(
    CTCFI = mean(mean_CTCFI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(log2_CTCFI = log2(CTCFI))

ko_model <- lmer(
  log2_CTCFI ~ Group * Timepoint + Rep + (1 | CondClone),
  data = ko_dat
)

ko_anova <- anova(ko_model)
ko_interaction_p <- ko_anova["Group:Timepoint", "Pr(>F)"]
ko_interaction_F <- ko_anova["Group:Timepoint", "F value"]
ko_interaction_df <- ko_anova["Group:Timepoint", "DenDF"]

# Planned group contrasts at 0 h and 24 h
ko_emm <- emmeans(ko_model, ~ Group | Timepoint)
ko_contrasts <- as.data.frame(pairs(ko_emm, adjust = "none"))

# Raw-scale recovery percentages reported for interpretation
ko_raw_means <- ko_dat %>%
  group_by(Group, Timepoint) %>%
  summarise(mean_CTCFI = mean(CTCFI), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = mean_CTCFI) %>%
  mutate(recovery_percent = 100 * (Hour0 - Hour24) / Hour0)

# -----------------------------------------------------------------------------
# 3. RTFDC1 overexpression vs LentiMock: primary mixed-effects analysis
# -----------------------------------------------------------------------------

oe_dat <- dat %>%
  filter(
    Group %in% c("LentiMock", "OE"),
    Timepoint %in% c("Hour0", "Hour24")
  ) %>%
  mutate(
    Group = factor(Group, levels = c("LentiMock", "OE")),
    Timepoint = factor(Timepoint, levels = c("Hour0", "Hour24")),
    Rep = factor(Rep),
    CondClone = factor(CondClone)
  ) %>%
  group_by(Group, CondClone, Rep, Timepoint) %>%
  summarise(
    CTCFI = mean(mean_CTCFI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(log2_CTCFI = log2(CTCFI))

oe_model <- lmer(
  log2_CTCFI ~ Group * Timepoint + Rep + (1 | CondClone),
  data = oe_dat
)

oe_anova <- anova(oe_model)
oe_interaction_p <- oe_anova["Group:Timepoint", "Pr(>F)"]
oe_interaction_F <- oe_anova["Group:Timepoint", "F value"]
oe_interaction_df <- oe_anova["Group:Timepoint", "DenDF"]

# Planned group contrasts at 0 h and 24 h
oe_emm <- emmeans(oe_model, ~ Group | Timepoint)
oe_contrasts <- as.data.frame(pairs(oe_emm, adjust = "none"))

# Raw-scale recovery percentages reported for interpretation
oe_raw_means <- oe_dat %>%
  group_by(Group, Timepoint) %>%
  summarise(mean_CTCFI = mean(CTCFI), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = mean_CTCFI) %>%
  mutate(recovery_percent = 100 * (Hour0 - Hour24) / Hour0)

# -----------------------------------------------------------------------------
# 4. siRTFDC1 vs NT: paired experiment-level recovery analysis
# -----------------------------------------------------------------------------

sirna_dat <- dat %>%
  filter(
    Group %in% c("NT", "siRTFDC1"),
    Timepoint %in% c("Hour0", "Hour24")
  ) %>%
  mutate(
    Group = factor(Group, levels = c("NT", "siRTFDC1")),
    Timepoint = factor(Timepoint, levels = c("Hour0", "Hour24")),
    Rep = factor(Rep)
  ) %>%
  group_by(Group, Rep, Timepoint) %>%
  summarise(
    CTCFI = mean(mean_CTCFI, na.rm = TRUE),
    .groups = "drop"
  )

sirna_recovery <- sirna_dat %>%
  pivot_wider(names_from = Timepoint, values_from = CTCFI) %>%
  mutate(
    recovery_percent = 100 * (Hour0 - Hour24) / Hour0,
    log2_recovery = log2(Hour0) - log2(Hour24)
  )

sirna_paired <- sirna_recovery %>%
  select(Rep, Group, recovery_percent, log2_recovery) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(recovery_percent, log2_recovery)
  )

sirna_test <- t.test(
  sirna_paired$log2_recovery_siRTFDC1,
  sirna_paired$log2_recovery_NT,
  paired = TRUE
)

sirna_p <- sirna_test$p.value

# -----------------------------------------------------------------------------
# 5. Print statistical results
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FIGURE 5B STATISTICAL RESULTS\n")
cat("============================================================\n\n")

cat("Enhancer KO vs Mock\n")
cat("Group x time interaction: F =", round(ko_interaction_F, 2),
    ", denominator df =", round(ko_interaction_df, 2),
    ", P =", signif(ko_interaction_p, 4), "\n")
print(ko_contrasts)
cat("Raw recovery percentages:\n")
print(ko_raw_means)

cat("\nRTFDC1 overexpression vs LentiMock\n")
cat("Group x time interaction: F =", round(oe_interaction_F, 2),
    ", denominator df =", round(oe_interaction_df, 2),
    ", P =", signif(oe_interaction_p, 4), "\n")
print(oe_contrasts)
cat("Raw recovery percentages:\n")
print(oe_raw_means)

cat("\nsiRTFDC1 vs NT\n")
cat("Paired t-test on log2 recovery: P =", signif(sirna_p, 4), "\n")
cat("Raw recovery percentages by independent experiment:\n")
print(sirna_recovery %>% select(Group, Rep, recovery_percent))

# Save a compact machine-readable summary of the three primary tests
stats_summary <- data.frame(
  comparison = c(
    "Enhancer KO vs Mock",
    "RTFDC1 OE vs LentiMock",
    "siRTFDC1 vs NT"
  ),
  test = c(
    "Group x time interaction from linear mixed-effects model",
    "Group x time interaction from linear mixed-effects model",
    "Paired t-test on log2 recovery"
  ),
  p_value = c(
    ko_interaction_p,
    oe_interaction_p,
    sirna_p
  )
)

stats_file <- file.path(output_dir, "Fig5B_primary_statistics.csv")
write.csv(
  stats_summary,
  stats_file,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Build the analysis-level data plotted in Figure 5B
# -----------------------------------------------------------------------------

# For enhancer KO and overexpression, each plotted observation is a
# clone x experimental-day mean after averaging technical replicate wells.
# For siRNA, each plotted observation is an independent-experiment mean after
# averaging replicate wells.

plot_ko <- dat %>%
  filter(Group %in% c("Mock", "KO")) %>%
  group_by(Group, CondClone, Rep, Timepoint) %>%
  summarise(GH2AX = mean(mean_CTCFI, na.rm = TRUE), .groups = "drop") %>%
  transmute(
    Experiment = "RTFDC1 enhancer KO",
    PlotGroup = if_else(Group == "Mock", "Control", "RTFDC1 perturbation"),
    Timepoint,
    GH2AX
  )

plot_oe <- dat %>%
  filter(Group %in% c("LentiMock", "OE")) %>%
  group_by(Group, CondClone, Rep, Timepoint) %>%
  summarise(GH2AX = mean(mean_CTCFI, na.rm = TRUE), .groups = "drop") %>%
  transmute(
    Experiment = "RTFDC1 overexpression",
    PlotGroup = if_else(Group == "LentiMock", "Control", "RTFDC1 perturbation"),
    Timepoint,
    GH2AX
  )

plot_si <- dat %>%
  filter(Group %in% c("NT", "siRTFDC1")) %>%
  group_by(Group, Rep, Timepoint) %>%
  summarise(GH2AX = mean(mean_CTCFI, na.rm = TRUE), .groups = "drop") %>%
  transmute(
    Experiment = "RTFDC1 knockdown",
    PlotGroup = if_else(Group == "NT", "Control", "RTFDC1 perturbation"),
    Timepoint,
    GH2AX
  )

fig_dat <- bind_rows(plot_ko, plot_oe, plot_si) %>%
  mutate(
    Timepoint = recode(
      Timepoint,
      "NoTreatment" = "Untreated",
      "Hour0" = "0 h",
      "Hour24" = "24 h"
    ),
    Timepoint = factor(Timepoint, levels = c("Untreated", "0 h", "24 h")),
    Experiment = factor(
      Experiment,
      levels = c(
        "RTFDC1 enhancer KO",
        "RTFDC1 overexpression",
        "RTFDC1 knockdown"
      )
    ),
    PlotGroup = factor(
      PlotGroup,
      levels = c("Control", "RTFDC1 perturbation")
    )
  )

# -----------------------------------------------------------------------------
# 7. Figure 5B
# -----------------------------------------------------------------------------

format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "f", digits = 4)
  } else {
    formatC(p, format = "f", digits = 3)
  }
}

stats_plot <- data.frame(
  Experiment = factor(
    c(
      "RTFDC1 enhancer KO",
      "RTFDC1 overexpression",
      "RTFDC1 knockdown"
    ),
    levels = levels(fig_dat$Experiment)
  ),
  x1 = c(2, 2, 2),
  x2 = c(3, 3, 3),
  label = c(
    paste0("P = ", format_p(ko_interaction_p)),
    paste0("P = ", format_p(oe_interaction_p)),
    paste0("P = ", format_p(sirna_p))
  )
)

panel_max <- fig_dat %>%
  group_by(Experiment) %>%
  summarise(ymax = max(GH2AX, na.rm = TRUE), .groups = "drop")

stats_plot <- stats_plot %>%
  left_join(panel_max, by = "Experiment") %>%
  mutate(
    y_bracket = ymax + 900,
    y_text = ymax + 1800
  )

y_upper <- max(stats_plot$y_text, na.rm = TRUE) + 1000

p_fig5b <- ggplot(
  fig_dat,
  aes(x = Timepoint, y = GH2AX, fill = PlotGroup)
) +
  geom_boxplot(
    aes(group = interaction(Timepoint, PlotGroup)),
    position = position_dodge(width = 0.68),
    width = 0.54,
    linewidth = 0.55,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_point(
    aes(group = PlotGroup),
    position = position_jitterdodge(
      jitter.width = 0.045,
      jitter.height = 0,
      dodge.width = 0.68,
      seed = 123
    ),
    shape = 21,
    size = 1.65,
    stroke = 0.25,
    color = "black",
    alpha = 0.80
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "grey80",
      "RTFDC1 perturbation" = "#4C78A8"
    )
  ) +
  geom_segment(
    data = stats_plot,
    aes(x = x1, xend = x2, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  geom_segment(
    data = stats_plot,
    aes(x = x1, xend = x1, y = y_bracket, yend = y_bracket - 350),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  geom_segment(
    data = stats_plot,
    aes(x = x2, xend = x2, y = y_bracket, yend = y_bracket - 350),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  geom_text(
    data = stats_plot,
    aes(x = (x1 + x2) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 2.6
  ) +
  facet_wrap(~ Experiment, nrow = 1) +
  labs(
    x = NULL,
    y = expression(gamma * "H2AX mean CTCFI"),
    fill = NULL
  ) +
  coord_cartesian(ylim = c(0, y_upper), clip = "on") +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    panel.spacing.x = grid::unit(0.55, "lines")
  )

print(p_fig5b)

pdf_file <- file.path(output_dir, "Fig5B_gH2AX.pdf")
png_file <- file.path(output_dir, "Fig5B_gH2AX.png")

ggsave(
  pdf_file,
  plot = p_fig5b,
  width = 7.0,
  height = 3.1,
  units = "in"
)

ggsave(
  png_file,
  plot = p_fig5b,
  width = 7.0,
  height = 3.1,
  units = "in",
  dpi = 600
)

cat("\nSaved:\n")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")
cat("  ", stats_file, "\n", sep = "")
cat("\nR session information:\n")
print(sessionInfo())

