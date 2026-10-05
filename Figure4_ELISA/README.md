# Figure 4 — IL-6 and IL-8 ELISA analysis

This folder contains the input data and R code used to reproduce the statistical analyses underlying Figure 4.

## Contents

- `Figure4_ELISA_input_data.xlsx` — cleaned ELISA input data used for Figure 4 only.
- `Fig4_IL6_IL8_analysis_and_plot.R` — calculates clone-level means, performs the statistical tests, and generates the Figure 4 panels.

## Input data

The workbook contains one sheet, `ELISA_data`, with the following columns:

| Column | Description |
|---|---|
| `Cytokine` | IL-6 or IL-8 |
| `Genotype` | Mock or Enhancer KO |
| `Clone` | Individual HMC3 clone |
| `Biological_replicate` | Biological replicate number within clone |
| `Concentration_pg_ml` | ELISA concentration in pg/ml |

Each concentration is the value obtained after averaging the technical triplicate ELISA measurements for that biological replicate.

The dataset contains three biological replicates for each of three independent clones per genotype for both cytokines.

## Statistical analysis

For each cytokine:

1. The three biological-replicate values are averaged within each clone.
2. The resulting clone means are used as the statistical analysis units (N = 3 clones per genotype).
3. Mock and enhancer-KO clone means are compared using a two-tailed Student's two-sample t-test assuming equal variances.

For the data used in the manuscript, the analysis gives:

- IL-6: `P = 0.13`
- IL-8: `P = 0.03`

The plotted panels use the clone-level means and match the manuscript convention of `n.s.` for IL-6 and `*` for IL-8; exact P values are retained in the CSV statistical output.

## Output

Running the script creates the following files in the same directory as the input workbook:

- `Figure_4_IL6_IL8_statistics.csv`
- `Figure_4_clone_means.csv`
- `Figure_4A_IL6.pdf` and `.png`
- `Figure_4B_IL8.pdf` and `.png`

When run interactively (for example, by sourcing the script in RStudio), both panels are also displayed in the Plots pane.

## R packages

The script requires:

- `readxl`
- `dplyr`
- `ggplot2`

## Reproducibility

Place `Figure4_ELISA_input_data.xlsx` and `Fig4_IL6_IL8_analysis_and_plot.R` in the same directory and run:

```bash
Rscript Fig4_IL6_IL8_analysis_and_plot.R
```

The script will calculate the clone-level analysis values directly from the biological-replicate input data before performing the statistical tests.
