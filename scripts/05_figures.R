# =============================================================================
# 05_figures.R - Figures of the Manuscript
# =============================================================================
# Description: Draws Figure 1 and Figures B.1-B.3 of the manuscript from the
#              summaries of scripts/03_analyze.R:
#                robustness_method_rmse_comparison  Figure 1
#                robustness_main_effects_rmse       Figure B.1
#                robustness_sparsity_interaction    Figure B.2
#                robustness_coverage_by_error       Figure B.3
#              The plotting code is the robustness block of the original
#              09_dissertation_figures.R, without the in-image titles, and with
#              Figure B.2 drawn from the same simulation run as the tables
#              (robustness_analysis.rds$sparsity_method).
# Usage: Rscript scripts/05_figures.R [--use-archived]
#   --use-archived  read the summaries in archived/analysis/ instead of
#                   output/robustness/analysis/
# Output: output/figures/ (.png at 600 dpi and .pdf)
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(forcats)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
sim <- if ("--use-archived" %in% args) "archived/analysis" else "output/robustness/analysis"
out_dir <- "output/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- file.path(sim, c("method_summary.csv", "robustness_analysis.rds", "error_summary.csv"))
missing <- inputs[!file.exists(inputs)]
if (length(missing)) {
  stop("Missing required inputs for the figures (run scripts/03_analyze.R, or use --use-archived): ",
       paste(missing, collapse = ", "))
}
cat("Reading summaries from:", sim, "\n")

family_palette <- c(
  "InfoMax" = "#16324F",
  "CS-UCB" = "#2A7F62",
  "Classical baseline" = "#8C5E3C"
)

theme_dissertation <- function() {
  theme_minimal(base_size = 12, base_family = "serif") +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "#555555", size = 10.5),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#DDDDDD", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "#EAEAEA", linewidth = 0.3),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}

save_plot <- function(plot, filename, width, height) {
  ggsave(file.path(out_dir, filename), plot = plot, width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(out_dir, sub("\\.png$", ".pdf", filename)), plot = plot, width = width, height = height, bg = "white")
  cat("wrote", file.path(out_dir, filename), "(+ .pdf)\n")
}

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

method_family <- function(label) {
  case_when(
    label == "InfoMax" ~ "InfoMax",
    grepl("^CS-UCB", label) ~ "CS-UCB",
    TRUE ~ "Classical baseline"
  )
}

method_summary <- read_csv(file.path(sim, "method_summary.csv"), show_col_types = FALSE)
robustness_rds <- readRDS(file.path(sim, "robustness_analysis.rds"))
method_error <- robustness_rds$method_error
error_summary <- read_csv(file.path(sim, "error_summary.csv"), show_col_types = FALSE)

# --- Figure 1: method-level RMSE, perfect calibration vs parameter error ---
robust_methods <- method_summary %>%
  transmute(method = method_label(selection_method), condition = "Perfect calibration", rmse = RMSE, coverage = coverage, width = mean_width) %>%
  bind_rows(
    method_error %>%
      transmute(method = method_label(selection_method), condition = "Parameter error", rmse = RMSE_noisy, coverage = coverage_noisy, width = NA_real_)
  )

method_order <- robust_methods %>%
  group_by(method) %>%
  summarise(reference = min(rmse), .groups = "drop") %>%
  arrange(reference) %>%
  pull(method)

robust_wide <- robust_methods %>%
  select(method, condition, rmse) %>%
  pivot_wider(names_from = condition, values_from = rmse) %>%
  mutate(method = factor(method, levels = rev(method_order)))

robust_methods_labeled <- robust_methods %>%
  mutate(method = factor(method, levels = rev(method_order)), family = method_family(as.character(method)))
robust_wide_labeled <- robust_wide %>%
  mutate(family = method_family(as.character(method)))

robustness_plot <- ggplot() +
  geom_segment(
    data = robust_wide_labeled,
    aes(y = method, yend = method, x = `Perfect calibration`, xend = `Parameter error`, color = family),
    linewidth = 1.0,
    alpha = 0.3
  ) +
  geom_point(data = robust_methods_labeled, aes(x = rmse, y = method, shape = condition, color = family), size = 3.2) +
  scale_color_manual(values = family_palette, guide = guide_legend(order = 1)) +
  scale_shape_manual(values = c(16, 17), guide = guide_legend(order = 2)) +
  labs(
    x = "RMSE",
    y = NULL,
    color = "Method family",
    shape = "Condition"
  ) +
  theme_dissertation()
save_plot(robustness_plot, "robustness_method_rmse_comparison.png", width = 9.8, height = 5.4)

# --- Figure B1: pooled main effects ---
main_effects <- bind_rows(
  robustness_rds$init_effect %>% transmute(factor = "Initialization", level = recode(initialization, censored_mle = "Censored MLE", staircase = "Staircase"), RMSE, coverage),
  robustness_rds$strat_effect %>% transmute(factor = "Stratification", level = if_else(stratified, "BAS", "No BAS"), RMSE, coverage),
  robustness_rds$exp_effect %>% transmute(factor = "Exposure control", level = if_else(exposure_control, "Sympson-Hetter", "No control"), RMSE, coverage)
) %>%
  mutate(level = fct_inorder(level))

main_effects_long <- main_effects %>%
  pivot_longer(cols = c(RMSE, coverage), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, RMSE = "RMSE", coverage = "Coverage (%)"))

main_effects_plot <- ggplot(main_effects_long, aes(x = value, y = level, color = factor)) +
  geom_point(size = 3) +
  facet_wrap(~ metric, scales = "free_x") +
  scale_color_manual(values = c(
    "Initialization" = "#1F4E79",
    "Stratification" = "#6B8E23",
    "Exposure control" = "#B85C38"
  )) +
  labs(
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme_dissertation()
save_plot(main_effects_plot, "robustness_main_effects_rmse.png", width = 10.5, height = 4.8)

# --- Figure B3: coverage by calibration-error level ---
error_long <- error_summary %>%
  transmute(sigma_err, `Noisy scoring` = coverage_noisy, `Oracle scoring` = coverage_oracle) %>%
  pivot_longer(cols = c(`Noisy scoring`, `Oracle scoring`), names_to = "condition", values_to = "coverage")

error_plot <- ggplot(error_long, aes(x = factor(sigma_err), y = coverage, group = condition, color = condition)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Noisy scoring" = "#1F4E79", "Oracle scoring" = "#B85C38")) +
  labs(
    x = expression(sigma[err]),
    y = "Coverage (%)",
    color = NULL
  ) +
  theme_dissertation()
save_plot(error_plot, "robustness_coverage_by_error.png", width = 7.2, height = 4.8)

# --- Figure B2: method-level RMSE by bank sparsity ---
method_by_sparsity <- as_tibble(robustness_rds$sparsity_method)

sparsity_interaction_data <- method_by_sparsity %>%
  mutate(method = method_label(selection_method), family = method_family(method)) %>%
  filter(method %in% c("InfoMax", "CS-UCB Bernstein", "CS-UCB Theory", "CS-UCB Mixture", "Progressive"))

sparsity_interaction_plot <- ggplot(
  sparsity_interaction_data,
  aes(x = sparsity, y = RMSE, group = method, color = family, linetype = method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4) +
  scale_color_manual(values = family_palette) +
  scale_linetype_manual(values = c(
    "InfoMax" = "solid",
    "CS-UCB Bernstein" = "solid",
    "CS-UCB Theory" = "dashed",
    "CS-UCB Mixture" = "dotted",
    "Progressive" = "solid"
  )) +
  scale_x_continuous(breaks = c(0.03, 0.05, 0.10, 0.30, 0.50), labels = label_percent(accuracy = 1)) +
  labs(
    x = "Bank sparsity level",
    y = "RMSE",
    color = "Method family",
    linetype = "Method"
  ) +
  theme_dissertation() +
  # Not in 09: with the Progressive line the two legends no longer fit on one row, so they are stacked.
  theme(legend.box = "vertical")
save_plot(sparsity_interaction_plot, "robustness_sparsity_interaction.png", width = 10, height = 5.4)
