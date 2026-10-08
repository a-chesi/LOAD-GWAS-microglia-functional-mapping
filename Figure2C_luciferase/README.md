# Figure 2C — Luciferase reporter analysis

This folder contains the input data and R code used to reproduce the statistical analysis underlying Figure 2C.

## Contents

- `Figure2C_luciferase_input_data.csv` — fold-change values for the luciferase reporter experiments.
- `Fig2C_luciferase_analysis_and_plot.R` — performs the one-way ANOVA, Tukey HSD post-hoc comparisons, and generates a manuscript-style Figure 2C plot.

## Input data

The CSV file contains two columns:

| Column | Description |
|---|---|
| `Condition` | Reporter condition: `Empty`, `PromOnly`, `Risk`, or `Protective` |
| `FoldChange` | Luciferase activity normalized to the RTFDC1 promoter-only condition |

There are seven independent assays per condition (N = 7).

`Risk` corresponds to the rs6024870 major risk allele (G), and `Protective` corresponds to the minor protective allele (A).

## Statistical analysis

The script reproduces the analysis described in the manuscript:

1. Performs a one-way ANOVA across the four reporter conditions.
2. Performs Tukey's HSD post-hoc pairwise comparisons.
3. Reports group means, SD, and SEM.
4. Generates the Figure 2C bar plot using mean ± SD, matching the manuscript figure.

For the data used in the manuscript, the rs6024870 risk-allele enhancer construct shows approximately 2.3-fold luciferase activity relative to the promoter-only control. The risk-versus-promoter comparison is significant by Tukey HSD (`P < 0.001`), the protective-versus-promoter comparison is also significant, and the risk-versus-protective allele comparison is not significant.

## Running the analysis

Place `Figure2C_luciferase_input_data.csv` and `Fig2C_luciferase_analysis_and_plot.R` in the same directory and run:

```bash
Rscript Fig2C_luciferase_analysis_and_plot.R
```

Alternatively, provide the input CSV path as the first command-line argument.

## Output

Running the script creates:

- `Figure_2C_luciferase_summary.csv`
- `Figure_2C_luciferase_ANOVA.csv`
- `Figure_2C_luciferase_TukeyHSD.csv`
- `Figure_2C_luciferase.pdf`
- `Figure_2C_luciferase.png`

When run interactively (for example, by sourcing the script in RStudio), the plot is also displayed in the Plots pane.

## R packages

The script requires:

- `dplyr`
- `ggplot2`
