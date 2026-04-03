# Replication files for *Debiased Machine Learning for Conformal Prediction of Counterfactual Outcomes Under Runtime Confounding*

## Requirements

The following R packages are required:

- `tidyverse`, `SuperLearner`, `ranger`, `grf`, `mvtnorm`, `randomForest`
- `foreach`, `doParallel`
- `ggridges`, `ggtext` (for plotting)

```r
install.packages(c("tidyverse", "SuperLearner", "ranger", "grf", "mvtnorm",
                    "randomForest", "foreach", "doParallel", "ggridges", "ggtext"))
```

## Main code files

- `R/estimation.R`: Core estimation functions — DML conformal procedure (`est_r_dml`) and weighted conformal prediction (`wcp_rc`) with runtime confounding adjustment weights
- `R/conformal-runtime-sims.R`: Fully synthetic simulation study
- `R/clean-acic.R`: Semi-synthetic data application with the 2018 ACIC data
- `R/conformal-plots.R`: Plots output from `conformal-runtime-sims.R`
- `R/plotting-acic.R`: Plots output from `clean-acic.R`

## Reproducing results

All scripts should be run from within the `R/` directory. Run in this order:

1. **Simulations:** `Rscript conformal-runtime-sims.R` — outputs results CSV in `../output/`
2. **ACIC application:** `Rscript clean-acic.R` — outputs results CSV in `../output/`
3. **Plots:** `Rscript conformal-plots.R` and `Rscript plotting-acic.R` — creates PDFs in `../figures/`

Note: The simulation and ACIC scripts use `doParallel` with 40–50 cores by default. Adjust `registerDoParallel(cores = ...)` for your machine.

## Other files

- `R/archive-funcs.R`, `R/demo-plots.R`, `R/test-surr-eff.R`, `R/conformal-runtime-test.R`: Exploratory/archived scripts not used in the final paper
- `output/`: Simulation results (CSV files with timestamps)
- `figures/`: Generated figures (PDF)
- `data/`: Data for the semi-synthetic application (`synthetic_data.csv`)
