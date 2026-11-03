# =============================================================================
# Helper Functions for Simulation Scripts
# =============================================================================
# Description: Utility functions for printing, saving, and loading results
# Source: Extracted from script/00_setup.R
# =============================================================================

#' Print Section Header
#'
#' @param title Section title string
print_section <- function(title) {
  width <- 60
  line <- paste(rep("=", width), collapse = "")
  cat("\n", line, "\n", sep = "")
  cat(" ", title, "\n", sep = "")
  cat(line, "\n\n", sep = "")
}

#' Print Progress Bar
#'
#' @param current Current iteration
#' @param total Total iterations
#' @param prefix Prefix label (default "Progress")
print_progress <- function(current, total, prefix = "Progress") {
  pct <- round(100 * current / total)
  bar_width <- 40
  filled <- round(bar_width * current / total)
  bar <- paste0(
    "[",
    paste(rep("=", filled), collapse = ""),
    paste(rep(" ", bar_width - filled), collapse = ""),
    "]"
  )
  cat("\r", prefix, ": ", bar, " ", pct, "% (", current, "/", total, ")", sep = "")
  if (current == total) cat("\n")
}

#' Save Results with Timestamp
#'
#' @param obj R object to save
#' @param name Base name for file
#' @param dir Directory to save in (default "processed_data")
#' @return File path of saved file
save_with_timestamp <- function(obj, name, dir = "processed_data") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0(name, "_", timestamp, ".rds")
  filepath <- file.path(dir, filename)
  saveRDS(obj, filepath)
  cat("Saved:", filepath, "\n")
  return(filepath)
}

#' Load Most Recent Results
#'
#' @param pattern File name pattern to match
#' @param dir Directory to search in (default "processed_data")
#' @return Loaded R object
load_most_recent <- function(pattern, dir = "processed_data") {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No files matching pattern: ", pattern)
  }
  # Sort by modification time (most recent first)
  files <- files[order(file.mtime(files), decreasing = TRUE)]
  cat("Loading:", files[1], "\n")
  readRDS(files[1])
}

message("=== Helper functions loaded ===")
