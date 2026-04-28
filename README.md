# Geostatistical Analysis & Spatial Statistics in R

An end-to-end spatial statistics project covering kriging-based 
spatial prediction, point process analysis, spatio-temporal 
simulation, and monitoring network optimisation — implemented 
entirely in R.

## Project Overview

This project demonstrates a complete geostatistical workflow 
applied to two real datasets and one simulated process, covering 
the full pipeline from exploratory spatial analysis through to 
model validation and practical network design recommendations.

## Analyses Covered

### 1. Kriging & Variogram Modelling
- Empirical variogram computation from the `head` dataset (n=29)
- Gaussian and Exponential model fitting via OLS — model 
  selected by AIC
- Directional variogram analysis to assess spatial anisotropy
- Ordinary kriging for point and grid prediction with 
  uncertainty quantification

### 2. Model Validation
- Leave-One-Out Cross-Validation (LOOCV) with residual mapping
- 5-Fold Cross-Validation — RMSE: 33.87, MAE and bias reported
- Spatial visualisation of prediction errors

### 3. Monitoring Network Optimisation
- Informational value assessment of each observation point
- Identification of high-value points (retain) vs low-value 
  points (candidates for removal)
- Practical recommendations for network design under resource 
  constraints

### 4. Spatial Point Process Analysis
- `hyytiala` forest dataset — four tree species (aspen, birch, 
  pine, rowan)
- Kernel density estimation and intensity surface mapping
- Quadrat test for Complete Spatial Randomness (CSR)
- Species-level intensity comparison
- Cross-type pair correlation function (birch vs rowan)
- F and G summary statistics for clustering assessment

### 5. Spatio-Temporal Poisson Process Simulation
- Non-homogeneous spatio-temporal Poisson process with intensity 
  λ(x,y,t) = 100(x+y)²(5−t)²
- Thinning algorithm implementation
- Empirical vs theoretical intensity comparison
- Subwindow simulations for local process behaviour

## Tech Stack

| Tool | Purpose |
|---|---|
| R | Primary analysis language |
| geoR | Variogram modelling and kriging |
| gstat | Cross-validation and spatial prediction |
| spatstat | Point process analysis |
| sp | Spatial data structures |
| RMarkdown | Reproducible reporting |

## Files

| File | Description |
|---|---|
| `geostatistical_analysis.Rmd` | Full analysis — run in RStudio to reproduce all results and plots |

## How to Run

1. Clone this repository
2. Open `geostatistical_analysis.Rmd` in RStudio
3. Install required packages if needed:
```r
install.packages(c("geoR", "gstat", "spatstat", "sp", "lattice"))
```
4. Click **Knit** to render the full HTML report with all plots

## Key Findings

- Gaussian model outperformed Exponential (AIC: 35M vs 205M) 
  for the `head` dataset
- 5-fold CV RMSE of 33.87 against data SD of 82.78 — reasonable 
  given small sample size (n=29)
- Birch and rowan show inter-species inhibition at short distances 
  (g(r) < 1), suggesting local competition
- Aspen exhibited strongest spatial clustering (4× intensity 
  variation across the plot)
- Spatio-temporal simulation successfully replicated theoretical 
  intensity decay pattern

## Context

Completed as part of the **Master of Data Science (AI 
Specialization)** at La Trobe University, Melbourne, Australia 
(2025).

---

*Author: Nafisa Ahmad | 
[LinkedIn](https://www.linkedin.com/in/nafisa-ahmad-957474199) | 
[GitHub](https://github.com/namoumita)*
