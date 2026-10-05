# Figure 5B — γH2AX recovery analysis

This folder contains the R code and input data used to reproduce the statistical analyses and plot shown in Figure 5B.

## Contents

- `Fig5B_gH2AX_analysis_and_plot.R` — performs the statistical analyses and generates Figure 5B.
- `well_averages.xlsx` — well-level γH2AX corrected total cell fluorescence intensity (CTCFI) values used as input.

## Input data

The script reads `well_averages.xlsx`, which contains well-level mean γH2AX CTCFI values for the following experimental groups:

- Mock
- Enhancer KO
- LentiMock
- RTFDC1 overexpression (OE)
- Non-targeting siRNA control (NT)
- RTFDC1 siRNA (siRTFDC1)

Measurements are available for untreated cells and at 0 h and 24 h following etoposide treatment.

## Analysis

The script reproduces the analyses described in the manuscript.

### Enhancer KO versus Mock

Technical replicate wells are first averaged within each clone, experimental day, and timepoint. CTCFI values are log2-transformed and analyzed using a linear mixed-effects model with:

- condition
- time
- condition × time interaction
- experimental day as a fixed blocking factor
- clone as a random effect

The condition × time interaction is the primary test of differential γH2AX recovery.

### RTFDC1 overexpression versus LentiMock

The same analysis is performed for RTFDC1-overexpression and LentiMock control cells: technical replicate wells are averaged within clone, experimental day, and timepoint, followed by a linear mixed-effects model on log2-transformed CTCFI values.

### RTFDC1 knockdown versus non-targeting control

For the siRNA experiment, replicate wells are averaged within each independent experiment and timepoint. Recovery is calculated as the change in log2-transformed CTCFI from 0 h to 24 h and compared between NT and siRTFDC1 conditions using a paired t-test.

## Figure

Figure 5B displays raw CTCFI values for the analysis-level experimental units. Statistical P values shown in the figure are derived from the log2-transformed analyses described above.

For the dataset used in the manuscript, the primary statistical results are:

- Enhancer KO versus Mock, condition × time: `P = 0.0054`
- RTFDC1 overexpression versus LentiMock, condition × time: `P = 0.0010`
- siRTFDC1 versus NT recovery: `P = 0.078`

## Running the analysis

Place `well_averages.xlsx` and `Fig5B_gH2AX_analysis_and_plot.R` in the same directory and run:

```bash
Rscript Fig5B_gH2AX_analysis_and_plot.R
```

Alternatively, provide the input workbook path as the first command-line argument.

## Output

Running the script creates the following files in the same directory as the input workbook:

- `Fig5B_gH2AX.pdf`
- `Fig5B_gH2AX.png`
- `Fig5B_primary_statistics.csv`

When run interactively (for example, by sourcing the script in RStudio), the plot is also displayed in the Plots pane.

## R packages

The script requires:

- `readxl`
- `dplyr`
- `tidyr`
- `ggplot2`
- `lme4`
- `lmerTest`
- `emmeans`

## Reproducibility

The provided input spreadsheet contains the well-level values used in the analysis. Running the supplied R script from the same directory reproduces the statistical analyses and Figure 5B.
