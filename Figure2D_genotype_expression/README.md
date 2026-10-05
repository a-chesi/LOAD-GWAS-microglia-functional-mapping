# Figure 2D — RTFDC1 genotype-expression analysis

This folder contains the R code used to reproduce the genotype-expression analysis and plot shown in Figure 2D.

## Contents

- `Fig2D_genotype_expression_analysis_and_plot.R` — performs the statistical analysis and generates Figure 2D.
- The individual-level input dataset is **not included** in this repository because it derives from the FreshMicro study and is subject to the source dataset's data-use requirements.

## Input data

The script expects an Excel file containing at least these two columns:

| Column | Description |
|---|---|
| `Genotype` | rs6024870 genotype; the analysis uses `GG` and `GA` individuals |
| `TPM` | RTFDC1 expression in TPM |

By default, the script looks for:

`RTFDC1_genotype_expression_input.xlsx`

Alternatively, provide the input file path on the command line:

```bash
Rscript Fig2D_genotype_expression_analysis_and_plot.R path/to/input.xlsx
```

Authorized users should obtain the underlying FreshMicro individual-level genotype and expression data through the original data-access mechanism described in the manuscript / AD Knowledge Portal rather than from this repository.

## Analysis

The script reproduces the analysis described in the manuscript:

1. Restricts the analysis to `GG` and `GA` individuals.
2. Removes TPM values more than 3 SD from the mean **within each genotype group**, in a single filtering pass.
3. Compares RTFDC1 expression between genotype groups using a two-tailed Student's t-test assuming equal variances.
4. Calculates Cohen's d with a 95% confidence interval.
5. Calculates post-hoc power from the observed effect size and sample sizes at alpha = 0.05.
6. Generates the Figure 2D boxplot.

For the dataset used in the manuscript, the filtering leaves 68 GG and 21 GA individuals.

## Output

Running the script creates the following files in the same directory as the input workbook:

- `Figure_2D_genotype_expression_statistics.csv`
- `Figure_2D_genotype_expression.pdf`
- `Figure_2D_genotype_expression.png`

When run interactively (for example, by sourcing the script in RStudio), the plot is also displayed in the Plots pane.

## R packages

The script requires:

- `readxl`
- `dplyr`
- `ggplot2`
- `effectsize`
- `pwr`

## Data availability

The individual-level FreshMicro input data are not redistributed here. This repository provides the analysis code needed to reproduce Figure 2D after obtaining authorized access to the underlying source data.
