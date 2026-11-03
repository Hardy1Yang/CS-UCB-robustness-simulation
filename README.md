# Code for the robustness study

This repository contains the R and C++ code that produced every simulation result in the manuscript
*Robustness in 2PL Adaptive Testing: Confidence-Sequence Upper Confidence Bound Item Selection Under Operational Constraints*.
It generates the item banks, runs the 2 × 2 × 2 factorial simulation (11,904,000 simulated adaptive tests), summarizes the results, and rebuilds the simulation tables and figures of the manuscript.
A last step checks the reproduced values against the archived results and against the tables as printed.

## Quick start

Open a terminal in this folder and run:

```
Rscript run_all.R --quick
```

This takes about 2 minutes on a recent laptop. It ends with `RESULT: all ... checks passed.` when the reproduction matches the archived results.

| Command | What it does | Time |
|---|---|---|
| `Rscript run_all.R --use-archived` | Rebuilds Tables 3–5 and B.1 and Figures 1 and B.1–B.3 from the summaries in `archived/analysis/`, then compares the tables with the manuscript | seconds |
| `Rscript run_all.R --quick` | Generates all 60 item banks and runs two short simulations of 24,800 tests each (see [Reproducibility](#reproducibility)), compares them with the archived results, then rebuilds the tables and figures as above | about 2 minutes |
| `Rscript run_all.R --cores=9` | Full reproduction: banks, the whole simulation, analysis, tables, figures, comparison | about 11 hours with 9 workers |

In RStudio, set the working directory to this folder and run `source("run_all.R")`. That runs the `--quick` version.

## Requirements

- R ≥ 4.3. The reported results were produced with R 4.5.1; `environment.txt` lists the environment of the verification run.
- A C++ compiler for Rcpp: Xcode Command Line Tools on macOS, Rtools on Windows, or `g++` on Linux.
- R packages:
  ```r
  install.packages(c("Rcpp", "doRNG", "doParallel", "foreach", "catR", "tidyverse", "scales"))
  ```
- Disk: a full run writes about 13 GB of intermediate files to `output/robustness_cpp/intermediate/`.
- Memory: the original run used a 10-core Apple M4 with 24 GB. It started with 9 workers and was later reduced to 5 (`parallel_settings$n_cores`) after a worker ran out of memory on the K = 900 banks. Use `--cores` to set the number of workers.

**macOS note.** With some recent macOS SDKs, `Rcpp::sourceCpp()` fails with an "unknown architecture" error from `tapi`.
Pointing `SDKROOT` to an older SDK of the Command Line Tools fixes it, for example:
`SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk Rscript run_all.R --quick`.
The scripts stop with an error if any C++ file fails to compile. They never fall back to slower R code under the same method name.

## What reproduces what

| Manuscript | Produced by | Output file(s) |
|---|---|---|
| Table 2: operational settings | `config/simulation_settings.R` (all design constants) | — |
| Table 3: methods under perfect item parameters | `scripts/04_tables.R` | `output/tables/robustness_method_perfect.tex` |
| Table 4: methods under item-parameter error | `scripts/04_tables.R` | `output/tables/robustness_method_error.tex` |
| Table 5: pooled factorial main effects | `scripts/04_tables.R` | `output/tables/robustness_factorial_effects.tex` |
| Table B.1: performance by calibration-error level | `scripts/04_tables.R` | `output/tables/error_level_summary.tex` |
| Figure 1: method-level RMSE | `scripts/05_figures.R` | `output/figures/robustness_method_rmse_comparison.{png,pdf}` |
| Figure B.1: pooled main effects | `scripts/05_figures.R` | `output/figures/robustness_main_effects_rmse.{png,pdf}` |
| Figure B.2: RMSE by bank sparsity | `scripts/05_figures.R` | `output/figures/robustness_sparsity_interaction.{png,pdf}` |
| Figure B.3: coverage by calibration-error level | `scripts/05_figures.R` | `output/figures/robustness_coverage_by_error.{png,pdf}` |
| Numbers quoted in the Abstract and Sections 6–8 | `scripts/03_analyze.R` | `output/robustness/analysis/*.csv`, `robustness_analysis.rds` |

Table 1 (the CS-UCB variants) contains no simulation results.
`expected/` holds the four simulation tables exactly as printed in the manuscript. `scripts/06_verify.R` compares them with `output/tables/`.

## Pipeline

| Step | Script | Output |
|---|---|---|
| 0 | `scripts/00_setup.R`: loads packages, settings and functions, compiles `src/*.cpp`; every script sources it | — |
| 1 | `scripts/01_generate_item_banks.R`: 15 banks (K = 300, 600, 900 × sparsity 3%, 5%, 10%, 30%, 50%) and 45 calibration-error banks (× σ_err = 0.05, 0.10, 0.20) | `output/item_banks/` |
| 2 | `scripts/02_run_simulation.R [--cores=N] [--quick]`: the factorial simulation under perfect and under noisy item parameters | `output/robustness/` |
| 3 | `scripts/03_analyze.R`: summaries by method, factor, sparsity and error level; ANOVA; exploratory plots | `output/robustness/analysis/`, `output/visualization/` |
| 4 | `scripts/04_tables.R [--use-archived]`: LaTeX tables | `output/tables/` |
| 5 | `scripts/05_figures.R [--use-archived]`: figures | `output/figures/` |
| 6 | `scripts/06_verify.R`: comparison with `archived/` and `expected/` | `output/verification_report.txt` |

Each script runs on its own from this folder, for example `Rscript scripts/02_run_simulation.R --cores=4`.
Step 2 must start in a fresh R process (see [Reproducibility](#reproducibility)); `run_all.R` runs every step in its own process.

### Simulation design

All constants are in `config/simulation_settings.R`.

- Item banks: 2PL with log a ~ N(0, 0.3²), clipped to [0.5, 2.5]. Exactly round(Kρ) items are extreme: b is drawn from N(0, 1.5²) until |b| > 2. The other items have b ~ N(0, 0.8²), clipped to [−2, 2].
- Calibration error: the algorithm sees a + ε_a clipped to [0.1, 3] and b + ε_b clipped to [−4, 4], with ε_a ~ N(0, σ_err²) and ε_b ~ N(0, (2σ_err)²). The noise is redrawn for every replication. Responses always come from the true parameters.
- Examinees: 31 ability levels from −3 to 3 in steps of 0.2, 100 replications each, T = 50 items per test.
- Factorial design (2 × 2 × 2): initialization (5-item staircase or censored MLE), BAS stratification (yes or no), Sympson–Hetter exposure control with maximum rate 0.25 (yes or no).
- Eight selection methods:

| Code name | Manuscript |
|---|---|
| `infomax` | InfoMax |
| `csucb_theory` | CS-UCB Theory |
| `csucb_obsinfo` | CS-UCB ObsInfo |
| `csucb_bernstein` | CS-UCB Bernstein |
| `csucb_mixture` | CS-UCB Mixture |
| `mpwi` | MPWI |
| `kl_information` | KL Information |
| `progressive` | Progressive |

This gives 15 banks × 8 cells × 8 methods × 3,100 tests = 2,976,000 tests under perfect parameters, and 45 × 8 × 8 × 3,100 = 8,928,000 tests under calibration error.

### Code layout

| Folder | Contents |
|---|---|
| `config/` | `packages.R` (packages, global seed), `simulation_settings.R` (design constants) |
| `utils/` | R functions. `robustness_functions.R` runs one adaptive test; `exposure_control.R` holds Sympson–Hetter selection and calibration; `ci_functions.R` holds the confidence sequences; `traditional_wrappers.R` holds MPWI, KL and Progressive |
| `src/` | C++ (Rcpp): MLE, item information, CS-UCB indices, MPWI/KL/Progressive selection, Sympson–Hetter calibration |
| `archived/` | Results of the original run, used by `--use-archived`, `--quick` and the comparison (see below) |
| `expected/` | The four simulation tables as printed in the manuscript |

## Output format

`output/robustness/factorial_perfect_params.rds` has one row per simulated test:

- `K`, `sparsity`, `initialization`, `stratified`, `exposure_control`, `factorial_label`, `selection_method`, `theta_true`, `rep`: the condition and replication.
- `theta_final`, `bias`, `abs_error`: the censored MLE after the last item and its error.
- `ci_coverage`, `ci_width`: whether the terminal Wald interval covers the true ability, and its width.
- `n_sparse_used`, `sparse_ratio`: how many sparse items the test administered.
- `T_actual`: the number of items administered (50).
- `cum_regret`: cumulative regret.
- `mean_selection_time`: mean time per item selection in seconds. It depends on the computer and is excluded from every comparison.

`output/robustness/factorial_param_error.rds` has the same layout plus `sigma_err`. Ability estimates, errors, coverage and widths come in two versions: `_noisy` uses the noisy parameters the algorithm saw, and `_oracle` rescores the same responses with the true parameters.

`robustness_results_full.rds` combines both tables with the run configuration. The per-test results of the original run are not included here (about 440 MB, plus 580 MB of item-level traces). The full run regenerates them.

## Reproducibility

Random numbers:

- Each bank has its own seed: `42 + K + 100ρ` for the item parameters, and `99 + K + 100ρ + 1000σ_err` for the noisy parameters of the error banks.
- `config/packages.R` calls `set.seed(20260113)`. Each simulated test then draws from its own L'Ecuyer-CMRG stream, created by `doRNG` from the session's RNG state. The results therefore do not depend on the number of workers.
- Under calibration error, each test also resets the seed to `20260113 + K + 1000ρ + 10000σ_err + 100f + 1009r`, where f is the factorial cell and r the replication. This seed fixes the noisy parameters and the responses. It is shared by all ability levels and all methods with the same f and r.
- The Sympson–Hetter calibration of each bank (step 2) uses the session's RNG state. That state depends on every earlier step of the same R process.

Because of the last two points, a run reproduces the archived results only if it starts in a fresh R process and runs uninterrupted. If a run is restarted, banks already saved in `output/robustness_cpp/intermediate/` are reused. The banks computed after the restart then get different random streams. Do not set `parallel_settings$use_parallel` to `FALSE`: the sequential branch seeds differently.

What has been checked (macOS arm64, R 4.5.1; see `environment.txt`):

- All 60 item banks generated by step 1 are identical to `archived/robustness_banks.rds` and `archived/error_banks.rds`.
- `--quick` runs:
  - bank K300_rho3 with factorial cell 1 (staircase, BAS, exposure control) under perfect parameters;
  - bank K300_rho3_err5 with cell 5 (staircase, BAS, no exposure control) under calibration error.

  Each is 8 methods × 3,100 tests. Both reproduce the original results in `archived/quick_reference.rds` exactly, row by row.
- Step 3 applied to the original per-test results reproduces every summary in `archived/analysis/` exactly.
- Steps 4 and 5 applied to those summaries reproduce the four tables in `expected/` character for character, and the four figure PNGs of the manuscript byte for byte.

What a full rerun is expected to reproduce (not checked, because it takes about 11 hours):

- The original run was interrupted while computing the 40th of the 45 calibration-error banks, and once more later. Each restart reused the banks already saved.
- A fresh, uninterrupted rerun should reproduce exactly:
  - all 15 banks under perfect parameters;
  - all error banks without exposure control;
  - the exposure-control cells of the first 39 error banks.
- The exposure-control cells of the last six error banks (K900 with ρ = 30% or 50%) used Sympson–Hetter calibrations computed after the restarts, and those calibrations were not saved. A rerun calibrates these banks again, so these 595,200 of the 11,904,000 tests will differ at the level of Monte Carlo error.
- `archived/calibration_params.rds` holds the calibrations of the 15 perfect-parameter banks, and step 6 checks them.

## Provenance

The code was extracted from the authors' dissertation project, in the state that produced the June 2026 run reported in the manuscript.
The selection, estimation, interval and exposure-control code is unchanged. Only the following was changed:

- File names and output folders (`processed_data/` and `visualization/` became `output/`).
- The package list: unused packages (mstR, mirtCAT, mirt, ltm and others) were removed, and a missing package now stops the run.
- The precision-study banks, functions and settings that the robustness code does not use were removed. Each bank has its own seed, so the removal does not change the robustness banks.
- The options `--cores` and `--quick` were added to step 2.
- Step 4 contains the robustness part of the original table script, with the captions printed in the manuscript.
- Step 5 contains the robustness part of the original figure script, without in-image titles. Figure B.2 there is drawn from the same run as the tables.
- Steps 6 and `run_all.R` are new.

## License

The code is released under the MIT License (see `LICENSE`).

## Citation

If you use this code, please cite the article. The citation will be added here after publication.
