# =============================================================================
# run_all.R - Reproduce the Robustness Study
# =============================================================================
# Usage (from this folder):
#   Rscript run_all.R [--cores=N]   full reproduction: simulation, analysis,
#                                   tables, figures, verification (about 8 hours
#                                   with 9 workers)
#   Rscript run_all.R --quick       item banks, two short simulation checks,
#                                   tables and figures from archived/ (minutes)
#   Rscript run_all.R --use-archived  tables and figures from archived/ only
#                                   (seconds)
# In an interactive session, source("run_all.R") runs the --quick version.
#
# Each step runs in its own R process. The simulation must start in a fresh
# process, because its random streams derive from the session's RNG state.
# =============================================================================

args <- if (interactive()) "--quick" else commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
use_archived <- "--use-archived" %in% args
cores <- grep("^--cores=", args, value = TRUE)
unknown <- setdiff(args, c("--quick", "--use-archived", cores))
if (length(unknown) > 0) stop("unknown option(s): ", paste(unknown, collapse = " "))
if (quick && use_archived) stop("use either --quick or --use-archived, not both")

if (!file.exists("scripts/00_setup.R")) {
  stop("run this script from the folder that contains run_all.R")
}

rscript <- file.path(R.home("bin"), "Rscript")
run_step <- function(script, ...) {
  step_args <- c(...)
  cat("\n>>> Rscript", script, step_args, "\n")
  started <- Sys.time()
  status <- system2(rscript, c(script, step_args))
  cat(sprintf(">>> %s finished in %.1f minutes (exit status %d)\n",
              script, as.numeric(difftime(Sys.time(), started, units = "mins")), status))
  status
}
stop_on_error <- function(status, script) {
  if (status != 0) stop(script, " failed; see the messages above")
}

if (use_archived) {
  stop_on_error(run_step("scripts/04_tables.R", "--use-archived"), "04_tables.R")
  stop_on_error(run_step("scripts/05_figures.R", "--use-archived"), "05_figures.R")
} else if (quick) {
  stop_on_error(run_step("scripts/01_generate_item_banks.R"), "01_generate_item_banks.R")
  stop_on_error(run_step("scripts/02_run_simulation.R", "--quick", cores), "02_run_simulation.R")
  stop_on_error(run_step("scripts/04_tables.R", "--use-archived"), "04_tables.R")
  stop_on_error(run_step("scripts/05_figures.R", "--use-archived"), "05_figures.R")
} else {
  stop_on_error(run_step("scripts/01_generate_item_banks.R"), "01_generate_item_banks.R")
  stop_on_error(run_step("scripts/02_run_simulation.R", cores), "02_run_simulation.R")
  stop_on_error(run_step("scripts/03_analyze.R"), "03_analyze.R")
  stop_on_error(run_step("scripts/04_tables.R"), "04_tables.R")
  stop_on_error(run_step("scripts/05_figures.R"), "05_figures.R")
}

verify_status <- run_step("scripts/06_verify.R")
if (!interactive()) quit(status = verify_status)
