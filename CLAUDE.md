# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an epidemiological modeling project in R that assesses the impact of an influenza pandemic on agricultural workers and food production in the United States. It combines a household-structured compartmental disease model with spatial and demographic data to simulate disease spread and agricultural impact.

## Running the Analysis

**Prerequisites:**
- R 4.5.0 with renv for dependency management
- Census API key set in `~/.Renviron` as `CENSUS_API_KEY`
- NAWS data file at `data/naws_all.sas7bdat` (restricted-access; not distributed with this repository)

**First-time setup:**
```r
renv::restore()                  # Install dependencies
source('code/setup_check.R')     # Validate environment
check_setup()                    # Will report any missing components
```

**Run full analysis from R (working directory must be project root):**
```r
source('code/run_analysis.R')
```

This orchestrates, in order:
1. `setup.R` - load dependencies and configuration
2. `summarystats.R` - descriptive statistics
3. `calibrate_model.R` - calibrate beta for each R0 / SAR / fold / gamma combination
4. `parameters.R` - build the sensitivity parameter sets (`pars_list`)
5. `run_regional_batch.R` - regional simulations for every parameter set (27 runs)
6. `sensitivity_analysis.R` - summaries and figures across sensitivity dimensions
7. `decomposition_analysis.R` - decompose the A-vs-C disparity into contributing factors
8. `simulate.R` (×3) - county-level simulations under each `adjust_hhvars` method
9. `crop_calendars.R` - crop production loss analysis (California commodities)
10. `generate_si_tables.R` and `compute_reff.R` - SI tables S3/S4/S6 and R_eff under vaccination
11. `plot_main_figures.R` - main publication figures

**Run individual components:**
```r
source('code/setup.R')                # Load all dependencies and config (do this first)
source('code/calibrate_model.R')      # National-level calibration only
source('code/parameters.R')           # Build pars_list (requires calibration)
source('code/run_regional_batch.R')   # Regional simulations for all of pars_list
source('code/simulate.R')             # County-level simulation (baseline params)
source('code/sensitivity_analysis.R') # Compare results across sensitivity dimensions
source('code/generate_si_tables.R')   # Rebuild SI tables from output/ CSVs
```

## Architecture

### Configuration (`code/config.R`)
Central configuration file containing all paths, parameters, and settings:
- `paths` - File paths for input data and output directories
- `data_settings` - ACS/Census year settings, NAWS date range
- `default_pars` - Baseline epidemiological and vaccination parameters
- `comorbidity_pars` - Obesity effect on the symptomatic fraction (`p_symp`) for each subpopulation
- `sim_settings` - Simulation options (time steps, parallelization, establishment threshold)
- `sensitivity_values` - One-at-a-time sensitivity ranges for every dimension (single source of truth; `parameters.R` reads from here)
- `region_map` - NAWS region name/abbreviation mapping

### Setup Files
- **`setup.R`** - Consolidated dependency loader; source this once at session start
- **`setup_check.R`** - Environment validation; run `check_setup()` to verify prerequisites
- **`utils.R`** - Helper functions with roxygen documentation (household state tables, crowding adjustment, `calculate_tau()`, `compute_p_symp()`, `compute_reff()`, etc.)

### Disease Model (`code/epimodels.R`)
Defines `household_model_twopop_crowding` with the `odin` package: a household-structured SIR model (House & Keeling 2009) with two populations (general community C and agricultural workers A) and household crowding. States are indexed by (susceptible, infected, recovered, household size, crowding status).

Vaccination reduces susceptibility: the group-specific multiplier `vax_mult = 1 - vax_eff * vax_cov` scales **both** the within-household and between-household force of infection on susceptibles.

### Data Pipeline
- **`import_acs.R`** - Downloads county-level ACS data via tidycensus API (household size, crowding)
- **`import_naws.R`** - Processes National Agricultural Workers Survey (SAS file at `data/naws_all.sas7bdat`)
- **`calibrate_model.R`** - Calibrates beta directly for each target R0 (and each SAR / fold / gamma value that affects calibration); outputs `calibrated_betas_df` and `output/calibrated_betas.csv`
- **`run_regional_batch.R`** - Runs regional simulations for all parameter sets in parallel
- **`simulate.R`** - Runs county-level simulations
- **`export_app_data.R`** - One-time export of small CSVs to `app/data/` for the Shiny app (`app/`)

### Parameter Sets (`code/parameters.R`)
`create_parset()` is the single factory for all parameter sets. The one-at-a-time sensitivity analysis varies one dimension at a time from baseline.

**Transmission dimensions (each produces regional simulation runs; baseline in bold):**
1. **R0:** 1.5, **2.0**, 3.0
2. **Mixing (eps):** 1/4, **1/3**, 1/2, 2/3, 3/4, 1 (eps = 0 is complete assortativity, eps = 1 is proportional mixing; assortativity η = 1 − eps, so baseline η = 2/3)
3. **SAR in crowded households:** 20%, 30%, **40%**, 50%, 60%
4. **Crowding fold difference:** 1, **2**, 3
5. **Gamma (recovery rate):** 1/3 (~3-day infectious period), **1/5** (~5-day), 1/10 (~10-day)
6. **Seed target:** C only, **both**, A only (which subpopulations receive initial infections)
7. **Ag worker vaccine coverage (`vax_cov_A`):** 0.2, **0.4**, 0.6, 0.8
8. **Community vaccine coverage (`vax_cov_C`):** 0.3, 0.4, **0.5**, 0.6
9. **Vaccine efficacy (`vax_eff`):** 0.2, 0.4, **0.6**, 0.8

**Comorbidity dimensions (applied post hoc to `p_symp`; reuse the baseline `r0_2` simulation):**
10. **Obesity OR for symptomatic disease:** 1, **1.5**, 3
11. **Ag worker obesity prevalence:** 0.40, 0.50, **0.55**, 0.60, 0.70

**Key Parameters:**

*Primary SAR parameters (defined in config.R):*
- `sar_uncrowded` - Baseline secondary attack rate for uncrowded households (20%)
- `sar_crowded` - Secondary attack rate for crowded households (40%)

*Derived transmission parameters (computed in parameters.R using helpers from utils.R):*
- `tau` - Within-household transmission rate, derived from `sar_uncrowded` using: `tau = SAR * gamma / (1 - SAR)`
- `tau_boost` - Additional transmission rate for crowded households, computed as the difference between tau_crowded and tau
- `beta` - Between-household transmission rate, calibrated by `calibrate_model.R`
- `vax_mult_C`, `vax_mult_A` - Vaccination susceptibility multipliers passed to the odin model

*Other epidemiological parameters:*
- `gamma` - Recovery rate (1/5 = 5-day infectious period)
- `eps` - Mixing between populations
- `crowding_fold_diff` - How much more likely large households are to be crowded

**Parameter Set Naming Convention** (`parset_name = {sens_type}_{round(sens_value, 4)}`):
- `r0_1.5`, `r0_2`, `r0_3` - R0 sensitivity runs (`r0_2` is the baseline run)
- `eps_0.25`, `eps_0.5`, etc. - Mixing sensitivity runs
- `sar_0.2`, `sar_0.3`, etc. - Crowded household SAR sensitivity runs
- `fold_1`, `fold_3` - Crowding fold difference sensitivity runs
- `gamma_0.3333`, `gamma_0.1` - Infectious period sensitivity runs
- `seed_1`, `seed_3` - Seed target sensitivity runs (1=C only, 3=A only; baseline=2=both)
- `vax_A_0.2`, `vax_C_0.3`, `vax_eff_0.8`, etc. - Vaccination sensitivity runs

### Sensitivity Analysis (`code/sensitivity_analysis.R`)
Loads all regional simulation outputs and generates comparative summaries:
- Calculates summary statistics (peak prevalence, attack rates, timing)
- Computes differential metrics (agricultural workers vs community)
- Creates comparison figures across sensitivity dimensions

### Output Files
- `output/regional_sim_{parset_name}.csv` - Regional sensitivity runs (e.g., `regional_sim_r0_2.csv`)
- `output/county_sim_{none,multiplicative,additive}.csv` - County-level daily trajectories (large; git-ignored)
- `output/calibrated_betas.csv` - Calibrated beta for each calibration combination
- `output/sensitivity_summary.csv` - Summary statistics for all parameter sets
- `output/sensitivity_differential.csv` - Differential metrics (A vs C) for all parameter sets
- `output/decomposition_results_{infections,cases}.csv` - Disparity decomposition
- `output/crop_impact_*.csv` - Crop production loss results
- `output/reff_table.csv` - R_eff under baseline vaccination for each R0 scenario
- `output/table_s3.html`, `table_s4.html`, `table_s6.html` - Formatted SI tables
- `figures/` - All figures (PDF and PNG)

## Key Conventions

- All scripts assume working directory is project root
- Source `setup.R` at session start to load all dependencies and configuration
- All configurable values (paths, parameters, sensitivity ranges, years) are in `config.R`
- Simulations use `future.apply::future_lapply()` for parallelization
- Regional mapping defined in `config.R` (`region_map` tibble)
- Maximum household size capped at 7; baseline crowding fold-difference = 2
- Output files include parameter metadata (`parset`, `parset_name`, `sens_type`, `sens_value`)
