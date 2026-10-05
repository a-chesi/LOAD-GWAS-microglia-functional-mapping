# LOAD GWAS microglia functional mapping

This repository contains R code and, where redistribution is permitted, input data supporting manuscript-specific statistical analyses and quantitative figures for:

**Burton et al. — “Functional mapping of LOAD GWAS loci in microglia implicates inflammation and DNA damage response at the CASS4 locus.”**

## Repository contents

- `Figure2D_genotype_expression/` — rs6024870 genotype–RTFDC1 expression analysis and Figure 2D plotting code.
- `Figure4_ELISA/` — IL-6 and IL-8 ELISA analysis for mock and enhancer-KO HMC3 clones.
- `Figure5B_gammaH2AX/` — γH2AX recovery analyses for enhancer deletion, RTFDC1 overexpression, and RTFDC1 knockdown.

Each folder contains a README describing the input data, statistical analysis, required R packages, and expected output files.

## Running the analyses

The analyses are organized as self-contained figure folders. From the relevant folder, run the corresponding R script with `Rscript`, or source the script in RStudio. Plots are displayed in interactive R sessions and saved as PDF and PNG files. Tabular statistical outputs are saved as CSV files.

For Figure 4 and Figure 5B, the input Excel workbooks are included. For Figure 2D, the individual-level FreshMicro input data are not redistributed because they are subject to the source dataset's data-access requirements; the analysis script can be run after obtaining the underlying data through the original data-access mechanism described in the manuscript.

## Software

Analyses were performed in R. Required packages are listed in the README within each figure folder.

## Data availability

Input data generated in this study are included here where appropriate for the analyses represented in this repository. External individual-level data subject to access restrictions are not redistributed.
