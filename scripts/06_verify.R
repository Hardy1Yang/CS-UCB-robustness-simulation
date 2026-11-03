# =============================================================================
# 06_verify.R - Comparison With the Archived Results and the Manuscript
# =============================================================================
# Description: Runs every check whose inputs exist in output/:
#   1. item banks      output/item_banks/            vs archived/*_banks.rds
#   2. quick check     output/quick/                 vs archived/quick_reference.rds
#   3. full run        output/robustness/            vs archived/calibration_params.rds
#                      output/robustness/analysis/   vs archived/analysis/
#   4. tables          output/tables/                vs expected/ (as printed in the manuscript)
#   Results 1-3 are compared exactly (identical()). Two things are excluded:
#   selection times, which measure the computer, not the method, and the
#   output folder names stored in the run configuration.
# Usage: Rscript scripts/06_verify.R
# Output: output/verification_report.txt
# =============================================================================

report <- character(0)
n_fail <- 0L
n_checks <- 0L

note <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  report <<- c(report, line)
}

check <- function(label, ok, detail = "") {
  n_checks <<- n_checks + 1L
  if (!ok) n_fail <<- n_fail + 1L
  note(sprintf("[%s] %s%s", if (ok) "SAME" else "DIFFERS", label,
               if (nzchar(detail)) paste0("  (", detail, ")") else ""))
}

drop_machine_specific <- function(d) {
  if (is.data.frame(d)) {
    d$mean_selection_time <- NULL
    rownames(d) <- NULL
  } else if (is.list(d)) {
    d$output_paths <- NULL
  }
  d
}

same_data <- function(x, y) {
  identical(drop_machine_specific(x), drop_machine_specific(y))
}

describe_difference <- function(x, y) {
  out <- all.equal(x, y)
  if (isTRUE(out)) "equal up to floating-point tolerance" else paste(head(out, 3), collapse = "; ")
}

note("Verification report, ", format(Sys.Date()))
note("")

# --- 1. Item banks ----------------------------------------------------------
for (f in c("robustness_banks.rds", "error_banks.rds")) {
  out <- file.path("output/item_banks", f)
  if (file.exists(out)) {
    check(paste("item banks:", f), identical(readRDS(out), readRDS(file.path("archived", f))))
  }
}

# --- 2. Quick check ---------------------------------------------------------
quick_files <- c(perfect = "output/quick/factorial_perfect_params.rds",
                 error = "output/quick/factorial_param_error.rds")
if (all(file.exists(quick_files))) {
  ref <- readRDS("archived/quick_reference.rds")
  for (cond in names(quick_files)) {
    new <- readRDS(quick_files[[cond]])
    ok <- same_data(new, ref[[cond]])
    check(sprintf("quick check, %s parameters: %d tests",
                  if (cond == "perfect") "perfect" else "noisy", nrow(new)),
          ok, if (ok) "" else describe_difference(new, ref[[cond]]))
  }
}

# --- 3. Full run ------------------------------------------------------------
if (file.exists("output/robustness/calibration_params.rds")) {
  new <- readRDS("output/robustness/calibration_params.rds")
  ref <- readRDS("archived/calibration_params.rds")
  same_bank <- vapply(names(ref), function(b) identical(new[[b]], ref[[b]]), logical(1))
  check(sprintf("full run: Sympson-Hetter calibration of the %d perfect-parameter banks", length(ref)),
        all(same_bank),
        if (all(same_bank)) "" else paste("differs for", paste(names(ref)[!same_bank], collapse = ", ")))
}
if (file.exists("output/robustness/analysis/robustness_analysis.rds")) {
  new <- readRDS("output/robustness/analysis/robustness_analysis.rds")
  ref <- readRDS("archived/analysis/robustness_analysis.rds")
  for (component in names(ref)) {
    new_c <- drop_machine_specific(new[[component]])
    ref_c <- drop_machine_specific(ref[[component]])
    ok <- identical(new_c, ref_c)
    check(paste("full run: summary", component), ok,
          if (ok) "" else describe_difference(new_c, ref_c))
  }
  for (f in c("method_summary.csv", "error_summary.csv", "factorial_summary.csv", "sparsity_summary.csv")) {
    new_csv <- read.csv(file.path("output/robustness/analysis", f))
    ref_csv <- read.csv(file.path("archived/analysis", f))
    check(paste("full run:", f), same_data(new_csv, ref_csv))
  }
}

# --- 4. Tables --------------------------------------------------------------
expected <- list.files("expected", pattern = "\\.tex$")
for (f in expected) {
  out <- file.path("output/tables", f)
  if (!file.exists(out)) {
    check(paste("table:", f), FALSE, "not produced; run scripts/04_tables.R")
    next
  }
  new_lines <- readLines(out)
  ref_lines <- readLines(file.path("expected", f))
  ok <- identical(new_lines, ref_lines)
  check(paste("table:", f, "matches the manuscript"), ok,
        if (ok) "" else paste("first differing line:",
                              which(new_lines[seq_len(min(length(new_lines), length(ref_lines)))] !=
                                      ref_lines[seq_len(min(length(new_lines), length(ref_lines)))])[1]))
}

note("")
if (n_checks == 0L) {
  note("RESULT: nothing to check; run the pipeline first (see README.md).")
} else if (n_fail == 0L) {
  note(sprintf("RESULT: all %d checks passed.", n_checks))
} else {
  note(sprintf("RESULT: %d of %d checks differ.", n_fail, n_checks))
}

dir.create("output", showWarnings = FALSE)
writeLines(report, "output/verification_report.txt")
cat("Report written to output/verification_report.txt\n")
quit(status = if (n_checks > 0L && n_fail == 0L) 0L else 1L)
