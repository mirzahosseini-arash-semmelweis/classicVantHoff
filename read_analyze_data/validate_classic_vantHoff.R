# -----------------------------------------------------------------------------
# validate_classic_vantHoff.R
#
# Validation of the automated decisions of
# classic_vantHoff_analyzer.R (version 1.1.0):
#
#   V1  Monte Carlo: false exclusions and detection power of the outlier screen
#       (v1.0 rule vs v1.1 baseline-aware rule) and the resulting dH error
#   V2  Monte Carlo: false-positive rate and power of the elution-order
#       crossing test
#   V3  Coelution completion: fits with vs without completed Peak 2 points
#   V4  Hold-up time: analytical and numerical t0 sensitivity, column- and
#       temperature-specific t0 scenarios
#   V5  Comparative evaluation: naive all-point OLS ("manual" regression),
#       Theil-Sen regression, v1.0 rule and v1.1 automated result
#   V6  Shared (condition-wide) residual check
#
# Run from this folder after analyze_classic_vantHoff.R has been checked to
# work. Runtime about 10 min on one core. All tables are written as CSV and
# the validation figure as PDF/SVG.
#
# -----------------------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(gridExtra)
})

source("../classicVantHoff/classic_vantHoff_analyzer.R")

set.seed(0)

cond_cols   <- c("compoundname", "columnID", "eluentID", "modifierID", "flow_mL_min")
series_cols <- c(cond_cols, "peak")
Rgas        <- 8.314462618

# ---------- Data exactly as in analyze_classic_vantHoff.R --------------------

src <- readLines("analyze_classic_vantHoff.R", encoding = "UTF-8")
i1  <- grep("^retention <- readRDS", src)
i2  <- grep("^vh_filt <- vh_fit\\(vh_temp_filt\\)", src)
eval(parse(text = src[i1:i2], encoding = "UTF-8"))

sigma_ref <- vh_filt$sigma_ref
cat("sigma_ref (ln k) =", signif(sigma_ref, 4),
    "from", vh_filt$sigma_ref_n, "standardized residuals\n")

temps <- c(10, 20, 30, 40)
x4    <- 1000 / (temps + 273.15)

# ---------- V1  Outlier screen: Monte Carlo ----------------------------------

sim_outlier <- function(delta_sd, n_sim = 2000, slope = 2.0) {

  map_dfr(seq_len(n_sim), function(i) {

    y   <- 1 + slope * x4 + rnorm(4, sd = sigma_ref)
    pos <- sample.int(4, 1)
    sgn <- sample(c(-1, 1), 1)
    y[pos] <- y[pos] + sgn * delta_sd * sigma_ref

    d <- tibble(temp_C = temps, invT_1000 = x4, ln_k = y)

    old <- vh_flag_outliers_one(d, sigma_ref = NA_real_)
    new <- vh_flag_outliers_one(d, sigma_ref = sigma_ref)

    dh <- function(dd) -Rgas * unname(coef(lm(ln_k ~ invT_1000, data = dd))[2])

    tibble(
      delta_sd     = delta_sd,
      excl_old     = any(old$.vh_outlier),
      excl_new     = any(new$.vh_outlier),
      hit_old      = old$.vh_outlier[pos],
      hit_new      = new$.vh_outlier[pos],
      dH_err_ols   = dh(d) - (-Rgas * slope),
      dH_err_old   = dh(old[!old$.vh_outlier, ]) - (-Rgas * slope),
      dH_err_new   = dh(new[!new$.vh_outlier, ]) - (-Rgas * slope)
    )
  })
}

delta_grid <- c(0, 1, 2, 3, 4, 5, 6, 8, 10, 15)

sim1 <- map_dfr(delta_grid, sim_outlier)

v1 <- sim1 |>
  group_by(delta_sd) |>
  summarise(
    delta_ln_k          = first(delta_sd) * sigma_ref,
    any_exclusion_old   = mean(excl_old),
    any_exclusion_new   = mean(excl_new),
    correct_hit_old     = mean(hit_old),
    correct_hit_new     = mean(hit_new),
    rmse_dH_ols         = sqrt(mean(dH_err_ols^2)),
    rmse_dH_old_rule    = sqrt(mean(dH_err_old^2)),
    rmse_dH_new_rule    = sqrt(mean(dH_err_new^2)),
    .groups = "drop"
  )

print(v1, width = Inf)
write.csv(v1, "validation_V1_outlier_simulation.csv", row.names = FALSE)

# ---------- V2  Crossing test: Monte Carlo -----------------------------------

run_cross <- function(y1, y2) {

  # assign labels by elution order (lower k = Peak 1), as a peak picker would
  lo <- pmin(y1, y2)
  hi <- pmax(y1, y2)

  d <- tibble(
    compoundname = "c", columnID = "c", eluentID = "c", modifierID = "c",
    flow_mL_min = 0.7,
    temp_C    = rep(temps, 2),
    invT_1000 = rep(x4, 2),
    peak      = rep(c("Peak 1", "Peak 2"), each = 4),
    ln_k      = c(lo, hi),
    .vh_peak2_imputed = FALSE
  )

  vh_crossing_test_one(d)$diagnostic$crossing_detected
}

# (a) False positives: true, non-crossing lines of every real two-peak
#     condition with two complete four-temperature fits
pairs_real <- vh_filt$results |>
  filter(peak %in% c("Peak 1", "Peak 2"), n_temp_total == 4, is.finite(slope)) |>
  select(all_of(cond_cols), peak, intercept, slope) |>
  pivot_wider(names_from = peak, values_from = c(intercept, slope)) |>
  filter(if_all(c(`intercept_Peak 1`, `intercept_Peak 2`), is.finite))

n_rep <- 300

v2a <- pairs_real |>
  mutate(
    false_crossing_rate = pmap_dbl(
      list(`intercept_Peak 1`, `slope_Peak 1`, `intercept_Peak 2`, `slope_Peak 2`),
      function(a1, b1, a2, b2) {
        mean(replicate(n_rep, run_cross(
          a1 + b1 * x4 + rnorm(4, sd = sigma_ref),
          a2 + b2 * x4 + rnorm(4, sd = sigma_ref)
        )))
      }
    ),
    ln_alpha_25C = abs((`intercept_Peak 2` + `slope_Peak 2` * 1000 / 298.15) -
                       (`intercept_Peak 1` + `slope_Peak 1` * 1000 / 298.15))
  )

cat("\nCrossing test, false-positive rate over", nrow(v2a), "real non-crossing pairs:",
    "mean", round(mean(v2a$false_crossing_rate), 4),
    "max",  round(max(v2a$false_crossing_rate), 4), "\n")

# (b) Power: two lines crossing at 25 C, difference in apparent dH varied
x_cross <- 1000 / 298.15

v2b <- map_dfr(c(1, 2, 3, 4, 5, 6, 8, 10, 15, 20), function(ddH) {
  db <- ddH / Rgas # slope difference corresponding to ddH (kJ/mol)
  b1 <- 2.0
  b2 <- b1 + db
  a1 <- 1 - b1 * x_cross
  a2 <- 1 - b2 * x_cross
  tibble(
    ddH_kJ_mol     = ddH,
    detection_rate = mean(replicate(1000, run_cross(
      a1 + b1 * x4 + rnorm(4, sd = sigma_ref),
      a2 + b2 * x4 + rnorm(4, sd = sigma_ref)
    )))
  )
})

print(v2b)
write.csv(v2a, "validation_V2a_crossing_false_positive.csv", row.names = FALSE)
write.csv(v2b, "validation_V2b_crossing_power.csv",          row.names = FALSE)

# ---------- V3  Coelution completion: with vs without ------------------------

vh_nocomp <- vh_fit(vh_temp_filt, complete_missing_peak2 = FALSE)

affected <- vh_filt$results |>
  filter(peak == "Peak 2", peak2_imputed_n > 0) |>
  select(all_of(cond_cols))

v3 <- vh_compare_fits(vh_filt, vh_nocomp, labels = c("completed", "observed_only")) |>
  semi_join(affected, by = cond_cols) |>
  filter(peak == "Peak 2") |>
  left_join(
    vh_filt$results |>
      filter(peak == "Peak 1") |>
      select(all_of(cond_cols), deltaH_peak1 = deltaH_kJ_mol),
    by = cond_cols
  ) |>
  left_join(
    vh_filt$results |>
      filter(peak == "Peak 2") |>
      select(all_of(cond_cols), peak2_imputed_n),
    by = cond_cols
  ) |>
  arrange(compoundname, columnID, eluentID)

print(v3 |> select(compoundname, columnID, eluentID, peak2_imputed_n,
                   n_temp_completed, n_temp_observed_only,
                   deltaH_kJ_mol_completed, deltaH_kJ_mol_observed_only,
                   deltaH_diff_kJ_mol, deltaH_peak1), n = Inf, width = Inf)
write.csv(v3, "validation_V3_coelution_completion.csv", row.names = FALSE)

# ---------- V4  Hold-up time sensitivity -------------------------------------

prep_args <- list(run_id_col = "run_key", strict_qc = TRUE)

# Marker-based t0 tables (median over wavelengths of the "T0" injections)
t0_marker <- retention_t0_marker |>
  filter(compoundname == "T0") |>
  group_by(columnID, eluentID, temp_C) |>
  summarise(t0 = median(tR1, na.rm = TRUE), .groups = "drop")

t0_marker_col <- t0_marker |>
  group_by(columnID, eluentID) |>
  summarise(t0 = mean(t0), .groups = "drop")

# Thermal-expansion scenario: at constant pump volumetric flow the column
# flow scales with 1/rho(T); t0 = V0 * rho(T) / (F * rho_pump). A volumetric
# expansion coefficient of about 1.2e-3 1/K (order of magnitude for methanol
# and 2-propanol) gives a ~3.6 % decrease of t0 between 10 and 40 C.
t0_expansion <- tibble(temp_C = temps, t0 = 2.54 * (1 - 1.2e-3 * (temps - 25)))

# Marker used only as a run-specific RELATIVE time reference: the pycnometric
# t0 is multiplied by exp(residual) of a smooth ln(t_marker) vs 1000/T fit per
# column x eluent. This removes temperature-specific timing deviations shared
# by all injections of a condition without importing the marker's own
# (retention-related) temperature trend.
t0_marker_detr <- t0_marker |>
  group_by(columnID, eluentID) |>
  mutate(t0 = 2.54 * exp(resid(lm(log(t0) ~ I(1000 / (temp_C + 273.15)))))) |>
  ungroup()

t0_scen <- list(
  pycnometric_2.54    = 2.54,
  minus_2pct          = 2.54 * 0.98,
  plus_2pct           = 2.54 * 1.02,
  minus_5pct          = 2.54 * 0.95,
  plus_5pct           = 2.54 * 1.05,
  thermal_expansion   = t0_expansion,
  marker_column_const = t0_marker_col,
  marker_column_temp  = t0_marker,
  marker_timing_detrended = t0_marker_detr
)

sens <- vh_t0_sensitivity(retention_filt, t0_scen, prepare_args = prep_args)

base <- sens |>
  filter(t0_scenario == "pycnometric_2.54") |>
  select(all_of(series_cols), dH_base = deltaH_kJ_mol, se_base = deltaH_se_kJ_mol,
         sens_pct = deltaH_t0_sens_kJ_mol_per_pct)

v4 <- sens |>
  filter(t0_scenario != "pycnometric_2.54") |>
  inner_join(base, by = series_cols) |>
  mutate(diff = deltaH_kJ_mol - dH_base) |>
  group_by(t0_scenario) |>
  summarise(
    n_series          = sum(is.finite(diff)),
    median_abs_diff   = median(abs(diff), na.rm = TRUE),
    p95_abs_diff      = quantile(abs(diff), 0.95, na.rm = TRUE),
    max_abs_diff      = max(abs(diff), na.rm = TRUE),
    median_rel_to_SE  = median(abs(diff) / se_base, na.rm = TRUE),
    .groups = "drop"
  )

print(v4, width = Inf)
write.csv(v4, "validation_V4_t0_sensitivity.csv", row.names = FALSE)

# Analytical first-order sensitivity vs numerical +/- 2 % refits
chk <- sens |>
  filter(t0_scenario %in% c("minus_2pct", "plus_2pct")) |>
  select(t0_scenario, all_of(series_cols), deltaH_kJ_mol) |>
  pivot_wider(names_from = t0_scenario, values_from = deltaH_kJ_mol) |>
  inner_join(base, by = series_cols) |>
  mutate(numeric_per_pct = (plus_2pct - minus_2pct) / 4) |>
  filter(is.finite(numeric_per_pct), is.finite(sens_pct))

cat("\nAnalytical vs numerical dH sensitivity to t0 (kJ/mol per 1 %):",
    "median analytical", round(median(abs(chk$sens_pct)), 3),
    "| max |analytical - numerical|", signif(max(abs(chk$sens_pct - chk$numeric_per_pct)), 3),
    "| range", paste(round(range(chk$sens_pct), 3), collapse = " to "), "\n")

# Analyte injections eluting BEFORE the hold-up marker of the same
# column x eluent x temperature (k < 0 if the marker time were used as t0)
early <- retention_filt |>
  filter(peak_detection_ok %in% TRUE, is.finite(tR1)) |>
  distinct(run_key, .keep_all = TRUE) |>
  inner_join(t0_marker |> rename(t0_marker = t0), by = c("columnID", "eluentID", "temp_C")) |>
  summarise(n_injections = n(), n_before_marker = sum(tR1 < t0_marker))
print(early)

write.csv(t0_marker |> pivot_wider(names_from = temp_C, values_from = t0, names_prefix = "t0_"),
          "validation_V4_marker_t0_table.csv", row.names = FALSE)

# ---------- V5  Comparative evaluation ---------------------------------------

# (a) "manual" analysis: every observed temperature, ordinary least squares,
#     no coelution completion, no outlier screen, no crossing test
vh_naive <- vh_fit(vh_temp_filt, complete_missing_peak2 = FALSE,
                   detect_peak_crossing = FALSE, max_outliers = 0,
                   outlier_scale = "none")

# (b) v1.0 behaviour (studentized residual only)
vh_v10 <- vh_fit(vh_temp_filt, outlier_scale = "none")

# (c) Theil-Sen (median of pairwise slopes) on the same observed points
theil_sen <- function(x, y) {
  ij <- combn(length(x), 2)
  median((y[ij[2, ]] - y[ij[1, ]]) / (x[ij[2, ]] - x[ij[1, ]]))
}

ts <- vh_temp_filt |>
  group_by(across(all_of(series_cols))) |>
  filter(n_distinct(temp_C) >= 3) |>
  summarise(deltaH_theil_sen = -Rgas * theil_sen(invT_1000, ln_k), .groups = "drop")

cmp <- vh_filt$results |>
  select(all_of(series_cols), dH_auto = deltaH_kJ_mol, se_auto = deltaH_se_kJ_mol,
         n_outliers, status) |>
  left_join(vh_naive$results |> select(all_of(series_cols), dH_naive = deltaH_kJ_mol),
            by = series_cols) |>
  left_join(vh_v10$results   |> select(all_of(series_cols), dH_v10 = deltaH_kJ_mol),
            by = series_cols) |>
  left_join(ts, by = series_cols)

summ_cmp <- function(a, b, lab) {
  ok <- is.finite(a) & is.finite(b)
  d  <- b[ok] - a[ok]
  tibble(
    comparison          = lab,
    n_series            = sum(ok),
    identical_lt_1e_9   = sum(abs(d) < 1e-9),
    median_abs_diff     = median(abs(d)),
    p95_abs_diff        = unname(quantile(abs(d), 0.95)),
    max_abs_diff        = max(abs(d)),
    n_diff_gt_2kJ       = sum(abs(d) > 2)
  )
}

v5 <- bind_rows(
  summ_cmp(cmp$dH_auto, cmp$dH_naive,         "automated v1.1 vs all-point OLS (manual)"),
  summ_cmp(cmp$dH_auto, cmp$deltaH_theil_sen, "automated v1.1 vs Theil-Sen"),
  summ_cmp(cmp$dH_auto, cmp$dH_v10,           "automated v1.1 vs v1.0 rule"),
  summ_cmp(cmp$dH_naive, cmp$dH_v10,          "all-point OLS vs v1.0 rule")
)

print(v5, width = Inf)
write.csv(v5,  "validation_V5_comparison_summary.csv", row.names = FALSE)
write.csv(cmp, "validation_V5_comparison_by_series.csv", row.names = FALSE)

# Exported temperature-level table for an independent re-implementation
write.csv(
  vh_filt$data |>
    filter(!.vh_outlier) |>
    select(all_of(series_cols), temp_C, invT_1000, ln_k),
  "validation_V5_fit_input.csv", row.names = FALSE
)

# Points excluded by v1.1 and v1.0
excl_new <- vh_filt$data |> filter(.vh_outlier) |>
  select(all_of(series_cols), temp_C, ln_k, .vh_rstudent, .vh_z_ref)
excl_old <- vh_v10$data  |> filter(.vh_outlier) |>
  select(all_of(series_cols), temp_C, ln_k, .vh_rstudent)
print(excl_new, width = Inf)
cat("v1.0 rule excluded", nrow(excl_old), "points; by temperature:\n")
print(table(excl_old$eluentID, excl_old$temp_C))
write.csv(excl_new, "validation_V5_excluded_points_v1.1.csv", row.names = FALSE)
write.csv(excl_old, "validation_V5_excluded_points_v1.0.csv", row.names = FALSE)

# ---------- V6  Shared residual check, with pycnometric and marker t0 ---------

v6_base <- vh_shared_residual_check(vh_filt, by = "eluentID") |>
  mutate(t0 = "pycnometric constant")

run_m   <- vh_prepare(retention_filt, run_id_col = "run_key", t0_external = t0_marker)
vh_mark <- vh_fit(vh_temperature_summary(run_m))

v6_mark <- vh_shared_residual_check(vh_mark, by = "eluentID") |>
  mutate(t0 = "marker, column x temperature")

run_d   <- vh_prepare(retention_filt, run_id_col = "run_key", t0_external = t0_marker_detr)
vh_detr <- vh_fit(vh_temperature_summary(run_d))

v6_detr <- vh_shared_residual_check(vh_detr, by = "eluentID") |>
  mutate(t0 = "pycnometric x detrended marker timing")

v6 <- bind_rows(v6_base, v6_mark, v6_detr)
print(v6, n = Inf, width = Inf)
cat("sigma_ref pycnometric:", signif(vh_filt$sigma_ref, 3),
    "| marker-based:", signif(vh_mark$sigma_ref, 3),
    "| detrended marker timing:", signif(vh_detr$sigma_ref, 3), "\n")
cat("median SE of dH: pycnometric", round(median(vh_filt$results$deltaH_se_kJ_mol, na.rm = TRUE), 3),
    "| detrended marker timing", round(median(vh_detr$results$deltaH_se_kJ_mol, na.rm = TRUE), 3), "\n")
write.csv(v6, "validation_V6_shared_residuals.csv", row.names = FALSE)

# Marker-time residual pattern (log t0_marker vs 1000/T per column x eluent)
v6_marker <- t0_marker |>
  group_by(columnID, eluentID) |>
  mutate(res = resid(lm(log(t0) ~ I(1000 / (temp_C + 273.15))))) |>
  group_by(eluentID, temp_C) |>
  summarise(mean_resid_ln_t0 = mean(res), frac_pos = mean(res > 0), .groups = "drop")
print(v6_marker)
write.csv(v6_marker, "validation_V6_marker_residuals.csv", row.names = FALSE)

# ---------- Validation figure ------------------------------------------------

cols_rule <- c("v1.0 rule (|r*| > 3)" = "#B2182B",
               "v1.1 rule (|r*| > 3 and |z| > 3)" = "#2166AC",
               "no screening (OLS)" = "grey40")

pA <- v1 |>
  select(delta_sd, any_exclusion_old, any_exclusion_new,
         correct_hit_old, correct_hit_new) |>
  pivot_longer(-delta_sd) |>
  mutate(
    rule = if_else(grepl("old", name), names(cols_rule)[1], names(cols_rule)[2]),
    what = if_else(grepl("any", name), "any point excluded", "displaced point excluded")
  ) |>
  ggplot(aes(delta_sd, value, colour = rule, linetype = what)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = cols_rule) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title    = "(A) Outlier screen, 4 temperatures",
    x        = expression("Displacement of one point ("*sigma[ref]*" units)"),
    y        = "Fraction of simulated series",
    colour   = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.box = "vertical",
        panel.grid.minor = element_blank())

pB <- v1 |>
  select(delta_sd, rmse_dH_ols, rmse_dH_old_rule, rmse_dH_new_rule) |>
  pivot_longer(-delta_sd) |>
  mutate(rule = recode(name,
                       rmse_dH_ols      = names(cols_rule)[3],
                       rmse_dH_old_rule = names(cols_rule)[1],
                       rmse_dH_new_rule = names(cols_rule)[2])) |>
  ggplot(aes(delta_sd, value, colour = rule)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = cols_rule) +
  labs(
    title  = expression("(B) Error of "*Delta*italic(H)[vH*","*app]),
    x      = expression("Displacement of one point ("*sigma[ref]*" units)"),
    y      = expression("RMSE of "*Delta*italic(H)~"(kJ mol"^-1*")"),
    colour = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        panel.grid.minor = element_blank())

pC <- ggplot(v2b, aes(ddH_kJ_mol, detection_rate)) +
  geom_hline(yintercept = mean(v2a$false_crossing_rate), linetype = "dashed") +
  geom_line(linewidth = 0.7, colour = "#2166AC") +
  geom_point(size = 1.6, colour = "#2166AC") +
  annotate("text", x = max(v2b$ddH_kJ_mol), y = mean(v2a$false_crossing_rate),
           label = "false-positive rate,\nreal non-crossing pairs",
           hjust = 1, vjust = -0.3, size = 2.8) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "(C) Elution-order crossing test",
    x     = expression("|"*Delta*Delta*italic(H)[vH*","*app]*"| of crossing pair (kJ mol"^-1*")"),
    y     = "Crossing detected"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

fig_val <- arrangeGrob(pA, pB, pC, nrow = 1, widths = c(1.15, 1, 1))

ggsave("plot_validation.pdf", fig_val, width = 13, height = 4.8)
ggsave("plot_validation.svg", fig_val, width = 13, height = 4.8)

sessionInfo()
