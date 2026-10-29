---
title: "README"
author: "Ligia"
date: "2026-10-28"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Maternal mortality risk among Black women in Brazil: a nationwide spatial consensus approach

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Reproducibility](https://img.shields.io/badge/Reproducibility-Verified-success.svg)](#)

This repository contains the full analytical codebase and computational pipeline for the epidemiological study: **"Maternal mortality risk among Black women in Brazil: a nationwide spatial consensus approach"**. 

The project integrates spatial autocorrelation, flexible spatial scan statistics (`rflexscan`), and Bayesian spatial smoothing via a Besag-York-Mollié (`CARBayes::S.CARbym`) model across all 5,570 Brazilian municipalities. The final output generates an integrated multi-criteria typology to prioritize municipal surveillance and targeted maternal health interventions.

---

## Repository Structure

```text
.
├── scripts/
│   ├── 01_fetch_sim_domat.R       # Automated retrieval of raw public datasets
│   ├── 02_fetch_sinasc.R          # Data harmonization and standardized variable coding
│   ├── 03_data_preprocessing.R     # Case filtering, cause classification, denominator setup
│   ├── 04_mmr_trend.R    # Standardized Mortality Ratio (SMR) & Global/Local Moran's I
│   ├── 05_comparison.R  # Chi-square test, z-scores, and Poisson rate ratios
│   ├── 06_flexible_spatial_scan.R  # Irregular cluster detection using rflexscan
│   └── 07_bayesian_spatial_model.R # CARBYM Bayesian model, MCMC diagnostics, integrated map
├── data_processed/                 # Cleaned Parquet tables and spatial shapefiles (git-ignored)
│   ├── domat_clean_2019_2024.parquet
│   ├── sinasc_clean_2019_2024.parquet
│   ├── tab_cluster.parquet
│   └── risk_black_maternal_mortality.shp
├── figures/                        
│   └── supp_figure_mcmc_diagnostics.png # figure mcmc diagnostic
│   ├── README.md                        # Documentation & execution instructions
├── .gitignore                      # Git exclusion rules for large datasets
└── maternal_mortality_spatial.Rproj# RStudio project configuration file
```
Prerequisites & Installation
To run this pipeline, ensure you have R (version ≥ 4.1.0) and RStudio installed.

Required R Packages
The pipeline uses pacman for automated package management. Run the following code in your R console to set up all dependencies:

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse,
  arrow,
  sf,
  geobr,
  SpatialEpi,
  spdep,
  rflexscan,
  CARBayes,
  coda,
  lattice,
  ggplot2,
  patchwork
)

Execution Pipeline
Execute the R scripts sequentially in the scripts/ directory:

Order	Script	Main Objective	Output Deliverable
01	01_data_acquisition.R	Download raw SIM (mortality) & SINASC (births) records.	data_raw/
02	02_data_cleaning.R	Standardize variables and filter maternal codes (ICD-10 O00–O99).	domat_clean_2019_2024.parquet
03	03_data_preprocessing.R	Map racial self-identification (Black/Parda vs. White).	sinasc_clean_2019_2024.parquet
04	04_exploratory_spatial.R	Compute crude SMR and spatial contiguity matrices.	Moran's I Statistics
05	05_cause_disparity_tests.R	Chi-Square, adjusted residuals (z-scores), and Poisson Rate Ratios.	Statistical summary tables
06	06_flexible_spatial_scan.R	Detect non-circular spatial clusters via FlexScan algorithm.	tab_cluster.parquet
07	07_bayesian_spatial_model.R	Fit Poisson CARBYM model, run MCMC convergence, export final map.	risk_black_maternal_mortality.shp, Figures

Key Methodological Highlights
1. Bayesian Spatial Smoothing (CARBYM)
To resolve the small-area instability inherent to maternal death rates in small municipalities, we fitted a Poisson Besag-York-Mollié (BYM) model using intrinsic conditional autoregressive (iCAR) priors:

$$\log(\mu_i) = \log(E_i) + \beta_0 + \phi_i + \theta_i$$

Where:$E_i$: Expected cases based on Black live birth population.$\phi_i$: Spatial structured random effect (Queen contiguity spatial dependence).$\theta_i$: Unstructured heterogeneity random effect.

2. Contiguity Matrix & Island Fix
Spatial contiguity across all 5,570 municipalities was enforced by dynamically linking isolated coastal/island municipalities (e.g., Ilha de Santa Catarina, Fernando de Noronha) to their nearest geometric neighbor, maintaining a perfectly symmetric W adjacency matrix.

3. MCMC Convergence Diagnostics
Convergence of posterior distributions was visually verified via trace/density plots and formally evaluated using Geweke's Z-score diagnostic (
Z=−1.49,p=0.137), indicating that post-burn-in iterations (20,000 burn-in; 100,000 retained) reached stationary equilibrium.

Integrated Risk Typology
The final integrated surveillance map synthesizes posterior probabilities of elevated risk (P($\RR_i$ >1∣Data)≥0.95) and FlexScan spatial scan results into five operational categories:

Confirmed Critical Risk: High Bayesian posterior probability (≥0.95) AND confirmed FlexScan cluster (Y>0).

Persistent Regional Risk: High Bayesian posterior probability (≥0.95) without FlexScan clustering (Y>0).

Local Excess Alert: Significant FlexScan cluster without regional Bayesian elevation (Y>0).

Epidemiological Silence (Investigate): High modeled spatial risk despite 0 recorded deaths (Y=0), signaling potential underreporting.

Expected Risk Level: Baseline risk within expected statistical boundaries.

Data Transparency & Ethics
Data sources: Brazilian Mortality Information System (SIM) and Live Births Information System (SINASC), provided by the Ministry of Health (SVSA/CGIAE).

Data access: Publicly available, anonymized administrative data downloaded via DATASUS. No institutional IRB approval was required under Brazilian National Health Council Resolution No. 510/2016.

Citation
If you use this codebase or methodology in your research, please cite:

Snippet de código
@article{maternal_mortality_black_brazil_2026,
  title = {Maternal mortality risk among Black women in Brazil: a nationwide spatial consensus approach},
  author = {Author Names},
  journal = {Journal Name},
  year = {2026},
  volume = {XX},
  number = {X},
  pages = {XXX--XXX},
  doi = {10.XXXX/XXXXX}
}
License
This repository is licensed under the MIT License — see the LICENSE file for details.