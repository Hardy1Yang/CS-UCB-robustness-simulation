# =============================================================================
# 03_analyze.R - Robustness Study Analysis
# =============================================================================
# Description: Analysis of Robustness Study results with 2×2×2 factorial design
#              (per committee recommendation)
#
# Factorial Design Analysis:
#   - Main effects: Initialization, Stratification, Exposure Control
#   - Two-way interactions
#   - Three-way interaction
#
# Additional Analyses:
#   - Sparsity effect
#   - Parameter error effect
#   - Method comparison across conditions
#
# Prerequisites: Run scripts/02_run_simulation.R first
# Usage: Rscript scripts/03_analyze.R   (from the repository root)
# Output: output/robustness/analysis/ (summaries read by 04_tables.R and
#         05_figures.R) and exploratory plots in output/visualization/robustness/
# =============================================================================

if (!exists("robustness_settings")) source("scripts/00_setup.R")

print_section("Robustness Study Analysis (2×2×2 Factorial)")

# --- Load Results ---

if (!file.exists("output/robustness/robustness_results_full.rds")) {
  stop("Please run scripts/02_run_simulation.R first")
}

robustness_results <- readRDS("output/robustness/robustness_results_full.rds")

factorial_df <- robustness_results$factorial_perfect
error_df <- robustness_results$factorial_error
config <- robustness_results$config

cat("Loaded factorial results:", nrow(factorial_df), "observations\n")
cat("Loaded error results:", nrow(error_df), "observations\n")
cat("\nFactorial design:\n")
print(config$factorial_design)

# --- Create Output Directories ---

output_dirs <- c(
  "output/robustness/analysis",
  "output/visualization/robustness/factorial",
  "output/visualization/robustness/sparsity",
  "output/visualization/robustness/error",
  "output/visualization/robustness/methods"
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# =============================================================================
# 1. Main Effects Analysis (Perfect Parameters)
# =============================================================================

print_section("1. Main Effects Analysis")

# 1.1 Main Effect: Initialization
init_effect <- factorial_df %>%
  group_by(initialization) %>%
  summarise(
    n = n(),
    mean_bias = mean(bias),
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    mean_width = mean(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Main Effect: Initialization ===\n")
print(init_effect)

# 1.2 Main Effect: Stratification
strat_effect <- factorial_df %>%
  group_by(stratified) %>%
  summarise(
    n = n(),
    mean_bias = mean(bias),
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    mean_width = mean(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Main Effect: Stratification ===\n")
print(strat_effect)

# 1.3 Main Effect: Exposure Control
exp_effect <- factorial_df %>%
  group_by(exposure_control) %>%
  summarise(
    n = n(),
    mean_bias = mean(bias),
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    mean_width = mean(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Main Effect: Exposure Control ===\n")
print(exp_effect)

# =============================================================================
# 2. Two-Way Interactions
# =============================================================================

print_section("2. Two-Way Interactions")

# 2.1 Initialization × Stratification
init_strat <- factorial_df %>%
  group_by(initialization, stratified) %>%
  summarise(
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )

cat("\n=== Initialization × Stratification ===\n")
print(init_strat %>% pivot_wider(names_from = stratified,
                                  values_from = c(RMSE, coverage)))

# 2.2 Initialization × Exposure Control
init_exp <- factorial_df %>%
  group_by(initialization, exposure_control) %>%
  summarise(
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )

cat("\n=== Initialization × Exposure Control ===\n")
print(init_exp %>% pivot_wider(names_from = exposure_control,
                                values_from = c(RMSE, coverage)))

# 2.3 Stratification × Exposure Control
strat_exp <- factorial_df %>%
  group_by(stratified, exposure_control) %>%
  summarise(
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )

cat("\n=== Stratification × Exposure Control ===\n")
print(strat_exp %>% pivot_wider(names_from = exposure_control,
                                 values_from = c(RMSE, coverage)))

# =============================================================================
# 3. Three-Way Interaction (Full Factorial)
# =============================================================================

print_section("3. Full 2×2×2 Factorial")

factorial_summary <- factorial_df %>%
  group_by(initialization, stratified, exposure_control) %>%
  summarise(
    n = n(),
    mean_bias = mean(bias),
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    mean_width = mean(ci_width, na.rm = TRUE),
    sparse_ratio = mean(sparse_ratio) * 100,
    .groups = "drop"
  ) %>%
  arrange(initialization, stratified, exposure_control)

cat("\n=== Full Factorial Summary ===\n")
print(factorial_summary)

# Best and worst combinations
best_rmse <- factorial_summary %>% filter(RMSE == min(RMSE))
best_coverage <- factorial_summary %>% filter(coverage == max(coverage))

cat("\n=== Best Combination (RMSE) ===\n")
print(best_rmse)

cat("\n=== Best Combination (Coverage) ===\n")
print(best_coverage)

# =============================================================================
# 4. ANOVA for Factorial Design
# =============================================================================

print_section("4. Factorial ANOVA")

# Convert to factors for ANOVA
factorial_df_anova <- factorial_df %>%
  mutate(
    initialization = factor(initialization),
    stratified = factor(stratified),
    exposure_control = factor(exposure_control)
  )

# 4.1 ANOVA for RMSE (using squared error as dependent variable)
cat("\n=== Three-way ANOVA: Squared Error ===\n")
aov_rmse <- aov(bias^2 ~ initialization * stratified * exposure_control,
                data = factorial_df_anova)
print(summary(aov_rmse))

# 4.2 ANOVA for Coverage
cat("\n=== Three-way ANOVA: Coverage ===\n")
aov_coverage <- aov(ci_coverage ~ initialization * stratified * exposure_control,
                    data = factorial_df_anova)
print(summary(aov_coverage))

# Effect sizes (eta-squared)
cat("\n=== Effect Sizes (Eta-squared) for RMSE ===\n")
ss_total <- sum(summary(aov_rmse)[[1]]$`Sum Sq`)
eta_sq <- summary(aov_rmse)[[1]]$`Sum Sq` / ss_total
names(eta_sq) <- rownames(summary(aov_rmse)[[1]])
print(round(eta_sq, 4))

# =============================================================================
# 5. Sparsity Effect Analysis
# =============================================================================

print_section("5. Sparsity Effect Analysis")

sparsity_summary <- factorial_df %>%
  group_by(sparsity) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    mean_sparse_used = mean(n_sparse_used),
    sparse_ratio = mean(sparse_ratio) * 100,
    .groups = "drop"
  )

cat("\n=== Overall Sparsity Effect ===\n")
print(sparsity_summary)

# Sparsity × Selection Method
sparsity_method <- factorial_df %>%
  group_by(sparsity, selection_method) %>%
  summarise(
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )

cat("\n=== Sparsity × Selection Method ===\n")
print(sparsity_method %>%
        pivot_wider(names_from = selection_method, values_from = c(RMSE, coverage)))

# Sparsity × Factorial
sparsity_factorial <- factorial_df %>%
  group_by(sparsity, initialization, stratified, exposure_control) %>%
  summarise(
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )

cat("\n=== Sparsity Effect by Factorial Condition ===\n")
print(sparsity_factorial %>% filter(sparsity %in% c(0.03, 0.50)))

# =============================================================================
# 6. Parameter Error Effect Analysis
# =============================================================================

print_section("6. Parameter Error Effect Analysis")

error_summary <- error_df %>%
  group_by(sigma_err) %>%
  summarise(
    n = n(),
    RMSE_noisy = sqrt(mean(bias_noisy^2)),
    RMSE_oracle = sqrt(mean(bias_oracle^2)),
    coverage_noisy = mean(ci_coverage_noisy) * 100,
    coverage_oracle = mean(ci_coverage_oracle) * 100,
    width_noisy = mean(ci_width_noisy, na.rm = TRUE),
    width_oracle = mean(ci_width_oracle, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Parameter Error Effect ===\n")
print(error_summary)

# Degradation from oracle to noisy
error_summary <- error_summary %>%
  mutate(
    RMSE_degradation = RMSE_noisy - RMSE_oracle,
    coverage_degradation = coverage_oracle - coverage_noisy
  )

cat("\n=== Performance Degradation due to Parameter Error ===\n")
print(error_summary %>% select(sigma_err, RMSE_degradation, coverage_degradation))

# Error × Selection Method
error_method <- error_df %>%
  group_by(sigma_err, selection_method) %>%
  summarise(
    RMSE_noisy = sqrt(mean(bias_noisy^2)),
    coverage_noisy = mean(ci_coverage_noisy) * 100,
    .groups = "drop"
  )

cat("\n=== Parameter Error × Selection Method ===\n")
print(error_method %>%
        pivot_wider(names_from = selection_method,
                    values_from = c(RMSE_noisy, coverage_noisy)))

# Error × Factorial
error_factorial <- error_df %>%
  group_by(sigma_err, initialization, stratified, exposure_control) %>%
  summarise(
    RMSE_noisy = sqrt(mean(bias_noisy^2)),
    coverage_noisy = mean(ci_coverage_noisy) * 100,
    .groups = "drop"
  )

cat("\n=== Parameter Error Effect by Factorial Condition ===\n")
print(error_factorial %>%
        group_by(initialization, stratified, exposure_control) %>%
        summarise(
          RMSE_err_05 = RMSE_noisy[sigma_err == 0.05],
          RMSE_err_20 = RMSE_noisy[sigma_err == 0.20],
          RMSE_increase = RMSE_err_20 - RMSE_err_05,
          .groups = "drop"
        ))

# =============================================================================
# 7. Selection Method Comparison
# =============================================================================

print_section("7. Selection Method Comparison")

method_summary <- factorial_df %>%
  group_by(selection_method) %>%
  summarise(
    n = n(),
    mean_bias = mean(bias),
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    mean_width = mean(ci_width, na.rm = TRUE),
    sparse_ratio = mean(sparse_ratio) * 100,
    mean_cum_regret = mean(cum_regret, na.rm = TRUE),
    median_cum_regret = median(cum_regret, na.rm = TRUE),
    mean_selection_time = mean(mean_selection_time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(RMSE)

cat("\n=== Selection Method Summary (Perfect Parameters) ===\n")
print(method_summary)

# Method × Factorial interaction
method_factorial <- factorial_df %>%
  group_by(selection_method, initialization, stratified, exposure_control) %>%
  summarise(
    RMSE = sqrt(mean(.data$bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )

cat("\n=== Best Factorial Condition for Each Method ===\n")
best_by_method <- method_factorial %>%
  group_by(selection_method) %>%
  filter(RMSE == min(RMSE)) %>%
  ungroup()
print(best_by_method)

# Method comparison under parameter error
method_error <- error_df %>%
  group_by(selection_method) %>%
  summarise(
    RMSE_noisy = sqrt(mean(bias_noisy^2)),
    RMSE_oracle = sqrt(mean(bias_oracle^2)),
    robustness = 1 - (RMSE_noisy - RMSE_oracle) / RMSE_oracle,
    coverage_noisy = mean(ci_coverage_noisy) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(robustness))

cat("\n=== Method Robustness to Parameter Error ===\n")
print(method_error)

# =============================================================================
# 8. Visualizations
# =============================================================================

print_section("8. Creating Visualizations")

# 8.1 Main Effects Plot
main_effects_data <- bind_rows(
  init_effect %>% mutate(factor = "Initialization", level = initialization) %>%
    select(factor, level, RMSE, coverage),
  strat_effect %>% mutate(factor = "Stratification",
                          level = ifelse(stratified, "Yes", "No")) %>%
    select(factor, level, RMSE, coverage),
  exp_effect %>% mutate(factor = "Exposure Control",
                        level = ifelse(exposure_control, "Yes", "No")) %>%
    select(factor, level, RMSE, coverage)
)

p_main_rmse <- ggplot(main_effects_data, aes(x = level, y = RMSE, fill = factor)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
  facet_wrap(~factor, scales = "free_x") +
  labs(
    title = "Main Effects on RMSE",
    subtitle = "2×2×2 Factorial Design",
    x = "",
    y = "RMSE",
    fill = "Factor"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

ggsave("output/visualization/robustness/factorial/main_effects_rmse.png",
       p_main_rmse, width = 12, height = 6, dpi = 300)
cat("Saved: main_effects_rmse.png\n")

# 8.2 Interaction Plots
p_init_strat <- ggplot(init_strat, aes(x = initialization, y = RMSE,
                                        color = factor(stratified),
                                        group = stratified)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  labs(
    title = "Initialization × Stratification Interaction",
    x = "Initialization",
    y = "RMSE",
    color = "Stratified"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("output/visualization/robustness/factorial/interaction_init_strat.png",
       p_init_strat, width = 8, height = 6, dpi = 300)
cat("Saved: interaction_init_strat.png\n")

# 8.3 Full Factorial Heatmap
factorial_summary <- factorial_summary %>%
  mutate(
    condition = paste0(
      substr(initialization, 1, 5), "_",
      ifelse(stratified, "Str", "NoS"), "_",
      ifelse(exposure_control, "Exp", "NoE")
    )
  )

p_factorial_heatmap <- ggplot(factorial_summary,
                               aes(x = condition, y = "RMSE", fill = RMSE)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(RMSE, 3)), color = "white", size = 4) +
  scale_fill_gradient(low = "green", high = "red") +
  labs(
    title = "RMSE by Factorial Condition",
    x = "Condition",
    y = "",
    fill = "RMSE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("output/visualization/robustness/factorial/factorial_heatmap.png",
       p_factorial_heatmap, width = 12, height = 4, dpi = 300)
cat("Saved: factorial_heatmap.png\n")

# 8.4 Sparsity Effect Plot
p_sparsity <- ggplot(sparsity_method,
                     aes(x = factor(sparsity * 100), y = RMSE,
                         color = selection_method, group = selection_method)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(
    title = "Effect of Sparsity on RMSE",
    subtitle = "By Selection Method",
    x = "Sparsity (%)",
    y = "RMSE",
    color = "Selection Method"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("output/visualization/robustness/sparsity/sparsity_rmse.png",
       p_sparsity, width = 10, height = 6, dpi = 300)
cat("Saved: sparsity_rmse.png\n")

p_sparsity_coverage <- ggplot(sparsity_method,
                              aes(x = factor(sparsity * 100), y = coverage,
                                  color = selection_method, group = selection_method)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray50") +
  labs(
    title = "Effect of Sparsity on Coverage",
    subtitle = "By Selection Method (dashed line = 95% target)",
    x = "Sparsity (%)",
    y = "Coverage (%)",
    color = "Selection Method"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("output/visualization/robustness/sparsity/sparsity_coverage.png",
       p_sparsity_coverage, width = 10, height = 6, dpi = 300)
cat("Saved: sparsity_coverage.png\n")

# 8.5 Parameter Error Effect Plot
error_long <- error_summary %>%
  select(sigma_err, RMSE_noisy, RMSE_oracle) %>%
  pivot_longer(cols = c(RMSE_noisy, RMSE_oracle),
               names_to = "estimate_type", values_to = "RMSE") %>%
  mutate(estimate_type = ifelse(estimate_type == "RMSE_noisy", "Noisy", "Oracle"))

p_error <- ggplot(error_long,
                  aes(x = factor(sigma_err), y = RMSE,
                      fill = estimate_type)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
  labs(
    title = "Effect of Parameter Error on RMSE",
    subtitle = "Comparing Noisy vs Oracle Estimates",
    x = expression(sigma[err]),
    y = "RMSE",
    fill = "Estimate Type"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("output/visualization/robustness/error/error_rmse.png",
       p_error, width = 10, height = 6, dpi = 300)
cat("Saved: error_rmse.png\n")

# 8.6 Method Comparison
p_methods <- ggplot(method_summary,
                    aes(x = reorder(selection_method, -RMSE), y = RMSE)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Selection Method Comparison",
    subtitle = "Perfect Parameter Condition",
    x = "Selection Method",
    y = "RMSE"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("output/visualization/robustness/methods/method_comparison.png",
       p_methods, width = 10, height = 6, dpi = 300)
cat("Saved: method_comparison.png\n")

# 8.7 Method Robustness Plot
p_robustness <- ggplot(method_error,
                       aes(x = reorder(selection_method, robustness),
                           y = robustness * 100)) +
  geom_bar(stat = "identity", fill = "forestgreen", alpha = 0.8) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Method Robustness to Parameter Error",
    subtitle = "Higher = more robust (100% = no degradation)",
    x = "Selection Method",
    y = "Robustness Score (%)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("output/visualization/robustness/methods/method_robustness.png",
       p_robustness, width = 10, height = 6, dpi = 300)
cat("Saved: method_robustness.png\n")

# =============================================================================
# 9. Save Analysis Results
# =============================================================================

print_section("9. Saving Analysis Results")

analysis_results <- list(
  # Main effects
  init_effect = init_effect,
  strat_effect = strat_effect,
  exp_effect = exp_effect,

  # Interactions
  init_strat = init_strat,
  init_exp = init_exp,
  strat_exp = strat_exp,

  # Full factorial
  factorial_summary = factorial_summary,

  # ANOVA
  anova_rmse = summary(aov_rmse),
  anova_coverage = summary(aov_coverage),
  eta_squared = eta_sq,

  # Sparsity
  sparsity_summary = sparsity_summary,
  sparsity_method = sparsity_method,

  # Parameter error
  error_summary = error_summary,
  error_method = error_method,

  # Methods
  method_summary = method_summary,
  method_error = method_error,

  # Config
  config = config
)

saveRDS(analysis_results, "output/robustness/analysis/robustness_analysis.rds")
cat("Saved analysis results\n")

# Save summary tables
write.csv(factorial_summary, "output/robustness/analysis/factorial_summary.csv",
          row.names = FALSE)
write.csv(sparsity_summary, "output/robustness/analysis/sparsity_summary.csv",
          row.names = FALSE)
write.csv(error_summary, "output/robustness/analysis/error_summary.csv",
          row.names = FALSE)
write.csv(method_summary, "output/robustness/analysis/method_summary.csv",
          row.names = FALSE)

cat("Saved CSV tables\n")

# =============================================================================
# 10. Key Findings Summary
# =============================================================================

print_section("10. Key Findings Summary")

cat("\n========================================\n")
cat("ROBUSTNESS STUDY KEY FINDINGS\n")
cat("========================================\n\n")

cat("1. MAIN EFFECTS:\n")
cat("   - Initialization:", init_effect$initialization[which.min(init_effect$RMSE)],
    "is better (RMSE:", round(min(init_effect$RMSE), 4), ")\n")
strat_rmse_yes <- strat_effect$RMSE[strat_effect$stratified == TRUE]
strat_rmse_no  <- strat_effect$RMSE[strat_effect$stratified == FALSE]
strat_verdict <- if (length(strat_rmse_yes) > 0 && length(strat_rmse_no) > 0) {
  ifelse(strat_rmse_yes < strat_rmse_no, "helps", "hurts")
} else { "insufficient data" }
cat("   - Stratification:", strat_verdict, "performance\n")

exp_rmse_yes <- exp_effect$RMSE[exp_effect$exposure_control == TRUE]
exp_rmse_no  <- exp_effect$RMSE[exp_effect$exposure_control == FALSE]
exp_verdict <- if (length(exp_rmse_yes) > 0 && length(exp_rmse_no) > 0) {
  ifelse(exp_rmse_yes < exp_rmse_no, "helps", "hurts")
} else { "insufficient data" }
cat("   - Exposure Control:", exp_verdict, "performance\n\n")

cat("2. BEST FACTORIAL COMBINATION:\n")
cat("   - Initialization:", best_rmse$initialization, "\n")
cat("   - Stratified:", best_rmse$stratified, "\n")
cat("   - Exposure Control:", best_rmse$exposure_control, "\n")
cat("   - RMSE:", round(best_rmse$RMSE, 4), "\n\n")

cat("3. SPARSITY EFFECT:\n")
cat("   - RMSE increases with sparsity\n")
cat("   - From", round(sparsity_summary$RMSE[1], 4), "(3%) to",
    round(sparsity_summary$RMSE[5], 4), "(50%)\n\n")

cat("4. PARAMETER ERROR EFFECT:\n")
cat("   - RMSE degradation from error:\n")
for (i in 1:nrow(error_summary)) {
  cat("     σ_err =", error_summary$sigma_err[i], ":",
      round(error_summary$RMSE_degradation[i], 4), "\n")
}

cat("\n5. METHOD RANKING (by RMSE):\n")
for (i in 1:min(5, nrow(method_summary))) {
  cat("   ", i, ".", method_summary$selection_method[i],
      "(RMSE:", round(method_summary$RMSE[i], 4), ")\n")
}

cat("\n6. MOST ROBUST METHOD:\n")
cat("   ", method_error$selection_method[1],
    "(Robustness:", round(method_error$robustness[1] * 100, 1), "%)\n")

best_regret <- method_summary %>% arrange(mean_cum_regret)
cat("\n7. LOWEST CUMULATIVE REGRET:\n")
for (i in 1:min(3, nrow(best_regret))) {
  cat("   ", i, ".", best_regret$selection_method[i],
      "(mean regret:", round(best_regret$mean_cum_regret[i], 4), ")\n")
}

fastest <- method_summary %>% arrange(mean_selection_time)
cat("\n8. FASTEST SELECTION TIME:\n")
for (i in 1:min(3, nrow(fastest))) {
  cat("   ", i, ".", fastest$selection_method[i],
      "(mean time:", round(fastest$mean_selection_time[i], 6), "s)\n")
}

# --- Complete ---

print_section("Robustness Study Analysis Complete")
cat("\nOutput files:\n")
cat("  - Analysis: output/robustness/analysis/\n")
cat("  - Visualizations: output/visualization/robustness/\n")
