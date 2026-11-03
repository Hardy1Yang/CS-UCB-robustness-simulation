# =============================================================================
# 04_tables.R - LaTeX Tables of the Manuscript
# =============================================================================
# Description: Writes the four simulation tables of the manuscript from the
#              summaries of scripts/03_analyze.R:
#                robustness_method_perfect.tex     Table 3
#                robustness_method_error.tex       Table 4
#                robustness_factorial_effects.tex  Table 5
#                error_level_summary.tex           Table B.1
#              The table code is the robustness block of the original
#              10_dissertation_tables.R; only the captions were changed to the
#              ones printed in the manuscript. Table 2 (operational settings)
#              is descriptive and is not generated; its values are the
#              settings in config/simulation_settings.R.
# Usage: Rscript scripts/04_tables.R [--use-archived]
#   --use-archived  read the summaries in archived/analysis/ instead of
#                   output/robustness/analysis/
# Output: output/tables/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if ("--use-archived" %in% args) "archived/analysis" else "output/robustness/analysis"
out_dir <- "output/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- file.path(in_dir, c("method_summary.csv", "error_summary.csv",
                                      "robustness_analysis.rds"))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required inputs for the tables (run scripts/03_analyze.R, or use --use-archived):\n",
    paste(" -", missing_files, collapse = "\n")
  )
}
cat("Reading summaries from:", in_dir, "\n")

method_label <- function(x) {
  recode(
    x,
    infomax = "InfoMax",
    csucb_bernstein = "CS-UCB Bernstein",
    csucb_theory = "CS-UCB Theory",
    csucb_obsinfo = "CS-UCB ObsInfo",
    csucb_mixture = "CS-UCB Mixture",
    progressive = "Progressive",
    kl_information = "KL Information",
    mpwi = "MPWI"
  )
}

fmt_num <- function(x, digits = 3) {
  sprintf(paste0("%.", digits, "f"), x)
}

write_latex_table <- function(filename, caption, label, column_spec, headers,
                              rows, size = "small", extra_before_tabular = NULL) {
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\", size),
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}")
  )

  if (!is.null(extra_before_tabular)) {
    lines <- c(lines, extra_before_tabular)
  }

  lines <- c(
    lines,
    paste0("\\begin{tabular}{", column_spec, "}"),
    "\\toprule",
    paste0(paste(headers, collapse = " & "), " \\\\"),
    "\\midrule"
  )

  if (length(rows) > 0) {
    lines <- c(lines, paste0(rows, " \\\\"))
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  out_path <- file.path(out_dir, filename)
  writeLines(lines, out_path)
  cat("Wrote:", out_path, "\n")
}

rows_from_df <- function(df) {
  df[] <- lapply(df, as.character)
  apply(df, 1, function(x) paste(x, collapse = " & "))
}

robustness_rds <- readRDS(file.path(in_dir, "robustness_analysis.rds"))
method_summary <- read_csv(
  file.path(in_dir, "method_summary.csv"),
  show_col_types = FALSE
)
error_summary <- read_csv(
  file.path(in_dir, "error_summary.csv"),
  show_col_types = FALSE
)

# ---------------------------------------------------------------------------
# Robustness tables
# ---------------------------------------------------------------------------
method_summary_rows <- method_summary %>%
  arrange(RMSE) %>%
  transmute(
    Method = method_label(selection_method),
    RMSE = fmt_num(RMSE, 3),
    `Coverage (\\%)` = fmt_num(coverage, 2),
    `Mean width` = fmt_num(mean_width, 3)
  ) %>%
  rows_from_df()

write_latex_table(
  filename = "robustness_method_perfect.tex",
  caption = "Method Comparison Under Perfect Item-Parameter Knowledge",
  label = "tab:robustness-method-perfect",
  column_spec = "lccc",
  headers = c("Method", "RMSE", "Coverage (\\%)", "Mean width"),
  rows = method_summary_rows
)

method_error_rows <- robustness_rds$method_error %>%
  arrange(RMSE_noisy) %>%
  transmute(
    Method = method_label(selection_method),
    `RMSE (noisy)` = fmt_num(RMSE_noisy, 3),
    `RMSE (oracle)` = fmt_num(RMSE_oracle, 3),
    `Coverage (noisy, \\%)` = fmt_num(coverage_noisy, 2)
  ) %>%
  rows_from_df()

write_latex_table(
  filename = "robustness_method_error.tex",
  caption = "Method Comparison Under Item-Parameter Error",
  label = "tab:robustness-method-error",
  column_spec = "lccc",
  headers = c("Method", "RMSE (noisy)", "RMSE (oracle)", "Coverage (noisy, \\%)"),
  rows = method_error_rows
)

factorial_effects <- bind_rows(
  robustness_rds$init_effect %>%
    transmute(
      factor = "Initialization",
      level = recode(initialization, censored_mle = "Censored MLE", staircase = "Staircase"),
      RMSE,
      coverage
    ),
  robustness_rds$strat_effect %>%
    transmute(
      factor = "Stratification",
      level = if_else(stratified, "BAS", "No BAS"),
      RMSE,
      coverage
    ),
  robustness_rds$exp_effect %>%
    transmute(
      factor = "Exposure control",
      level = if_else(exposure_control, "Yes", "No"),
      RMSE,
      coverage
    )
)

factorial_rows <- factorial_effects %>%
  transmute(
    Factor = factor,
    Level = level,
    RMSE = fmt_num(RMSE, 3),
    `Coverage (\\%)` = fmt_num(coverage, 2)
  ) %>%
  rows_from_df()

write_latex_table(
  filename = "robustness_factorial_effects.tex",
  caption = "Pooled Factorial Main Effects",
  label = "tab:robustness-factorial-effects",
  column_spec = "llcc",
  headers = c("Factor", "Level", "RMSE", "Coverage (\\%)"),
  rows = factorial_rows
)

error_summary_rows <- error_summary %>%
  arrange(sigma_err) %>%
  transmute(
    sigma = fmt_num(sigma_err, 2),
    `RMSE (noisy)` = fmt_num(RMSE_noisy, 3),
    `RMSE (oracle)` = fmt_num(RMSE_oracle, 3),
    `Coverage (noisy, \\%)` = fmt_num(coverage_noisy, 2),
    `Coverage (oracle, \\%)` = fmt_num(coverage_oracle, 2)
  ) %>%
  rows_from_df()

write_latex_table(
  filename = "error_level_summary.tex",
  caption = "Aggregate Robustness Performance by Calibration-Error Level",
  label = "tab:error-level-summary",
  column_spec = "lcccc",
  headers = c("$\\sigma_{\\mathrm{err}}$", "RMSE (noisy)", "RMSE (oracle)",
              "Coverage (noisy, \\%)", "Coverage (oracle, \\%)"),
  rows = error_summary_rows
)
