# Clear environment and load functions, data ----------------------------------

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("../classicVantHoff/classic_vantHoff_analyzer.R")

retention <- readRDS("retention_data.rds")

retention <- retention |>
  mutate(run_key = sub("__[0-9]+$", "", logical_key))

retention_t0_marker <- retention |>
  filter(compoundname == "T0")

retention <- retention |> 
  filter(compoundname != "T0")

# Extracting retention data for classic van't Hoff analysis -------------------

vh_run <- vh_prepare(
  retention,
  run_id_col  = "run_key",
  strict_qc   = TRUE,
  t0_external = 2.54
)

vh_run |>
  count(
    compoundname,
    columnID,
    eluentID,
    modifierID,
    flow_mL_min,
    peak,
    temp_C
  )

vh_temp <- vh_temperature_summary(vh_run)

# Linear fit ------------------------------------------------------------------

vh <- vh_fit(vh_temp)

plot_vanthoff(
  vh,
  facet_vars = c("compoundname", "columnID"),
  facet_ncol = 6,
  colour_by  = "eluentID"
)

# Manual filtering and data cleaning ------------------------------------------

excluded_compounds <- c("Labetalol", "Propranolol", "Timolol")
retention_filt <- retention |>
  filter(!compoundname %in% excluded_compounds)

update_lookup <- tibble::tribble(
  ~logical_key,          ~tR1,      ~fwhm1,    ~tR2,      ~fwhm2,    ~peak_detection_ok, ~manual_override,
  "20260811__7354__210", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7354__230", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7354__254", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7402__210", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7402__230", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7402__254", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7546__210", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7546__230", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7546__254", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7841__210", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7841__230", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7841__254", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7833__210", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7833__230", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7833__254", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7278__210", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7278__230", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7278__254", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7326__210", 3.266664,  0.106674,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7326__230", 3.266664,  0.106674,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7326__254", 3.266664,  0.106674,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7374__210", 3.319997,  0.106674,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7374__230", 3.319997,  0.106674,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7374__254", 3.319997,  0.106674,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7837__210", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7837__230", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7837__254", NA_real_,  NA_real_,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__6954__210", 47.03909,  12.17852,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__6954__230", 47.03909,  12.17852,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__6954__254", 47.03909,  12.17852,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7582__210", 29.74539,  6.517592,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7582__230", 29.74539,  6.517592,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7582__254", 29.74539,  6.517592,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__8221__210", 19.16313,  3.433821,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__8221__230", 19.16313,  3.433821,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__8221__254", 19.16313,  3.433821,  NA_real_,  NA_real_,  TRUE,               TRUE,
  "20260811__7194__210", 5.643206,  0.926798,  10.87616,  1.496879,  TRUE,               TRUE,
  "20260811__7194__230", 5.643206,  0.926798,  10.87616,  1.496879,  TRUE,               TRUE,
  "20260811__7194__254", 5.643206,  0.926798,  10.87616,  1.496879,  TRUE,               TRUE,
  "20260811__7829__210", 5.046583,  0.730104,  8.249684,  0.930132,  TRUE,               TRUE,
  "20260811__7829__230", 5.046583,  0.730104,  8.249684,  0.930132,  TRUE,               TRUE,
  "20260811__7829__254", 5.046583,  0.730104,  8.249684,  0.930132,  TRUE,               TRUE,
  "20260811__8461__210", 4.566617,  0.550078,  6.506477,  0.646758,  TRUE,               TRUE,
  "20260811__8461__230", 4.566617,  0.550078,  6.506477,  0.646758,  TRUE,               TRUE,
  "20260811__8461__254", 4.566617,  0.550078,  6.506477,  0.646758,  TRUE,               TRUE
)

retention_filt <- retention_filt |>
  rows_update(update_lookup, by = "logical_key", unmatched = "ignore")

vh_run_filt <- vh_prepare(
  retention_filt,
  run_id_col  = "run_key",
  strict_qc   = TRUE,
  t0_external = 2.54
)

vh_temp_filt <- vh_temperature_summary(vh_run_filt)

vh_filt <- vh_fit(vh_temp_filt)

vh_filt$sigma_ref # dataset-level reference residual scale of ln k

fit_status <- vh_fit_status_summary(vh_filt)
fit_status

shared_resid <- vh_shared_residual_check(vh_filt, by = "eluentID")
shared_resid

write.csv(fit_status,   "fit_status_summary.csv",    row.names = FALSE)
write.csv(shared_resid, "shared_residual_check.csv", row.names = FALSE)
write.csv(vh_filt$results |> select(where(~ !is.list(.x))),
          "vh_results.csv", row.names = FALSE)

plot_fits <- plot_vanthoff(
  vh_filt,
  facet_vars = c("compoundname", "columnID"),
  facet_ncol = 6,
  colour_by  = "eluentID"
)

ggsave(
  filename = "plot_fits.svg",
  plot     =  plot_fits,
  width    = 14,
  height   = 13
)
ggsave(
  filename = "plot_fits.pdf",
  plot     =  plot_fits,
  width    = 14,
  height   = 13
)

# Temperature effects ---------------------------------------------------------

vh_temp_effect <- vh_k_temp_effect(vh_filt$results)

vh_temp_effect |>
  filter(peak == "Peak 1") |>
  select(
    compoundname,
    columnID,
    eluentID,
    modifierID,
    k_ref,
    delta_k_per_K,
    pct_k_per_K
  )

vh_alpha_effect <- vh_selectivity_temp_effect(vh_filt$results)

vh_alpha_effect |>
  select(
    compoundname,
    columnID,
    eluentID,
    modifierID,
    alpha_ref,
    delta_alpha_per_K,
    pct_alpha_per_K
  )

# Heatmap of results ----------------------------------------------------------

vh_results_status <- vh$results |>
  mutate(
    fit_ok = if_else(status == "ok", "ok", "flagged")
  )

hm_status <- prep_heatmap_data(
  data        = vh_results_status,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "fit_ok",
  filter_list = list(peak = "Peak 1"),
  agg_fun     = function(x, na.rm = TRUE) names(sort(table(x), decreasing = TRUE))[1]
)

ggplot(hm_status, aes(x = .col, y = .row, fill = value)) +
  geom_tile(color = "white") +
  facet_wrap(vars(), scales = "free") +
  scale_fill_manual(values = c(ok = "lightgreen", flagged = "tomato")) +
  theme_bw() +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

retention_filt <- retention |> 
  filter(!compoundname %in% excluded_compounds)

heatmap_tR <- make_heatmap(
  data        = retention_filt,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "tR1",
  facet_vars  = "temp_C",
  filter_list = list(excluded = FALSE),
  agg_fun     = median,
  fill_label  = expression(italic(t)[italic(R)*1]),
  x_label     = "Column | Eluent",
  y_label     = "Compound",
  text_digits = 1
)

ggsave(
  filename = "heatmap_tR.svg",
  plot     =  heatmap_tR,
  width    = 10,
  height   = 7
)
ggsave(
  filename = "heatmap_tR.pdf",
  plot     =  heatmap_tR,
  width    = 10,
  height   = 7
)

heatmap_Rs <- make_heatmap(
  data        = retention_filt,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "Rs12",
  facet_vars  = "temp_C",
  agg_fun     = median,
  fill_label  = expression(italic(R[s])),
  high        = "#66FF33",
  text_digits = 2
)

ggsave(
  filename = "heatmap_Rs.svg",
  plot     =  heatmap_Rs,
  width    = 10,
  height   = 7
)
ggsave(
  filename = "heatmap_Rs.pdf",
  plot     =  heatmap_Rs,
  width    = 10,
  height   = 7
)

heatmap_H <- make_heatmap(
  data        = vh_filt$results,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "deltaH_kJ_mol",
  filter_list = list(peak = "Peak 1"),
  agg_fun     = median,
  fill_label  = expression(Delta*H~"(kJ/mol)"),
  text_digits = 0,
  title       = expression("(A)"~Delta*H~"of peak 1")
)

vh_deltaH_diff <- vh_filt$results |>
  filter(peak %in% c("Peak 1", "Peak 2")) |>
  
  mutate(
    peak_short = recode(
      peak,
      "Peak 1" = "p1",
      "Peak 2" = "p2"
    ),
    peak_present = TRUE
  ) |>
  
  select(
    compoundname,
    columnID,
    eluentID,
    modifierID,
    flow_mL_min,
    peak_short,
    deltaH_kJ_mol,
    peak_present
  ) |>
  
  pivot_wider(
    names_from  = peak_short,
    values_from = c(deltaH_kJ_mol, peak_present),
    values_fill = list(
      peak_present = FALSE
    )
  ) |>
  
  mutate(
    deltaH_diff_abs = case_when(
      
      # Two successfully fitted peaks
      peak_present_p1 &
        peak_present_p2 &
        is.finite(deltaH_kJ_mol_p1) &
        is.finite(deltaH_kJ_mol_p2) ~
        abs(deltaH_kJ_mol_p1 - deltaH_kJ_mol_p2),
      
      # Genuine single-peak case
      peak_present_p1 &
        !peak_present_p2 &
        is.finite(deltaH_kJ_mol_p1) ~
        0,
      
      # Peak 1 missing/failed, or Peak 2 exists but its fit failed
      TRUE ~ NA_real_
    ),
    
    deltaH_diff_signed = case_when(
      peak_present_p1 &
        peak_present_p2 &
        is.finite(deltaH_kJ_mol_p1) &
        is.finite(deltaH_kJ_mol_p2) ~
        deltaH_kJ_mol_p2 - deltaH_kJ_mol_p1,
      
      peak_present_p1 &
        !peak_present_p2 &
        is.finite(deltaH_kJ_mol_p1) ~
        0,
      
      TRUE ~ NA_real_
    )
  )

heatmap_Hdiff <- make_heatmap(
  data        = vh_deltaH_diff,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "deltaH_diff_signed",
  fill_label  = expression(Delta*H[2] - Delta*H[1]~"(kJ/mol)"),
  text_digits = 0,
  title       = expression("(B) Signed difference of"~Delta*H)
)

heatmap_H_Hdiff <- grid.arrange(
  heatmap_H, heatmap_Hdiff,
  nrow = 1
)

ggsave(
  filename = "heatmap_H_Hdiff.svg",
  plot     =  heatmap_H_Hdiff,
  width    = 12,
  height   = 4
)
ggsave(
  filename = "heatmap_H_Hdiff.pdf",
  plot     =  heatmap_H_Hdiff,
  width    = 12,
  height   = 4
)

vh_temp_effect <- vh_temp_effect |> 
  mutate(
    delta_millik_per_K = 1000 * delta_k_per_K
  )
vh_alpha_effect <- vh_alpha_effect |> 
  mutate(
    delta_millialpha_per_K = 1000 * delta_alpha_per_K
  )

heatmap_k_temp <- make_heatmap(
  data        = vh_temp_effect,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "delta_millik_per_K",
  filter_list = list(peak = "Peak 1"),
  fill_label  = expression(10^3 %.% frac(Delta*italic(k), Delta*italic(T)) ~ "(K"^{-1}*")"),
  text_digits = 0,
  title       = "(A) Temperature dependence of retention factor at 25 °C (Peak 1)"
)

heatmap_alpha_temp <- make_heatmap(
  data        = vh_alpha_effect,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "delta_millialpha_per_K",
  fill_label  = expression(10^3 %.% frac(Delta*italic(alpha), Delta*italic(T))~"(K"^{-1}*")"),
  text_digits = 0,
  title       = "(B) Temperature dependence of selectivity at 25 °C"
)

heatmap_pctk_temp <- make_heatmap(
  data        = vh_temp_effect,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "pct_k_per_K",
  filter_list = list(peak = "Peak 1"),
  fill_label  = expression(100*frac(Delta*italic(k), italic(k)*Delta*italic(T))~"(% K"^{-1}*")"),
  text_digits = 1,
  title       = "(C) Relative temperature dependence of retention factor (Peak 1)"
)

heatmap_pctalpha_temp <- make_heatmap(
  data        = vh_alpha_effect,
  row_var     = "compoundname",
  col_vars    = c("columnID", "eluentID"),
  value_var   = "pct_alpha_per_K",
  fill_label  = expression(100*frac(Delta*italic(alpha), italic(alpha)*Delta*italic(T))~"(% K"^{-1}*")"),
  text_digits = 1,
  title       = "(D) Relative temperature dependence of selectivity"
)

heatmap_temp_effect <- grid.arrange(
  heatmap_k_temp,    heatmap_alpha_temp,
  heatmap_pctk_temp, heatmap_pctalpha_temp,
  nrow = 2
)

ggsave(
  filename = "heatmap_temp_effect.svg",
  plot     =  heatmap_temp_effect,
  width    = 14,
  height   = 10
)
ggsave(
  filename = "heatmap_temp_effect.pdf",
  plot     =  heatmap_temp_effect,
  width    = 14,
  height   = 10
)

# Synthetic example of elution order reversal ---------------------------------

# Construct two true molecular identities

toy_temps <- c(10, 20, 30, 40)

T_cross_C <- 25
x_cross   <- 1000 / (T_cross_C + 273.15)

# Two true van't Hoff lines:
#
# Identity A: ln(k) = aA + bA * 1000/T
# Identity B: ln(k) = aB + bB * 1000/T
#
# Both are forced to cross at ln(k) = 1 at 25 °C.

bA <- 2
bB <- 5

aA <- 1 - bA * x_cross
aB <- 1 - bB * x_cross

toy_truth <- crossing(
  temp_C        = toy_temps,
  true_identity = c("Identity A", "Identity B")
) |>
  mutate(
    T_K       = temp_C + 273.15,
    invT_1000 = 1000 / T_K,
    
    ln_k = case_when(
      true_identity == "Identity A" ~ aA + bA * invT_1000,
      true_identity == "Identity B" ~ aB + bB * invT_1000
    ),
    
    # Add small deterministic experimental deviations so that
    # the crossed model is much better, but not mathematically perfect.
    ln_k = ln_k + case_when(
      true_identity == "Identity A" & temp_C == 20 ~  0.020,
      true_identity == "Identity B" & temp_C == 20 ~ -0.010,
      true_identity == "Identity A" & temp_C == 30 ~  0.015,
      true_identity == "Identity B" & temp_C == 30 ~ -0.015,
      TRUE ~ 0
    ),
    
    k = exp(ln_k)
  )

toy_obs <- toy_truth |>
  group_by(temp_C) |>
  arrange(k, .by_group = TRUE) |>
  mutate(
    peak = paste("Peak", row_number())
  ) |>
  ungroup() |>
  mutate(
    compoundname = "Compound",
    columnID     = "Column",
    eluentID     = "Eluent",
    modifierID   = "DEA",
    flow_mL_min  = 0.7,
    n_runs       = 1L,
    phase_ratio  = NA_real_
  ) |>
  select(
    compoundname,
    columnID,
    eluentID,
    modifierID,
    flow_mL_min,
    temp_C,
    T_K,
    invT_1000,
    peak,
    k,
    ln_k,
    n_runs,
    phase_ratio
  ) |>
  group_by(
    compoundname,
    columnID,
    eluentID,
    modifierID,
    flow_mL_min,
    peak
  ) |>
  mutate(
    series_id = cur_group_id()
  ) |>
  ungroup()

toy_vh <- vh_fit(
  toy_obs,
  complete_missing_peak2    = FALSE,
  detect_peak_crossing      = TRUE,
  crossing_sse_ratio_cutoff = 0.5,
  max_outliers              = 0 # Demonstration focused purely on identity crossing
)

toy_vh$peak_identity_diagnostics

cols_peak <- c(
  "Peak 1" = "#F05837",
  "Peak 2" = "#595775"
)

p_before <- ggplot(data    = toy_obs,
                   mapping = aes(
                     x     = invT_1000,
                     y     = ln_k,
                     color = peak,
                     shape = peak,
                     group = peak
                   )) +
  geom_line(
    linewidth = 0.7,
    linetype  = "dashed"
  ) +
  geom_point(
    size = 3
  ) +
  scale_color_manual(
    values = cols_peak
  ) +
  scale_shape_manual(
    values = c(
      "Peak 1" = 16,
      "Peak 2" = 17
    )
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    title    = "(A) Peak identity based on elution order",
    subtitle = "Apparent Peak 1 and Peak 2 trajectories cannot cross",
    x        = expression(1000/italic(T)~(K^{-1})),
    y        = expression(ln(italic(k))),
    color    = "Observed peak",
    shape    = "Observed peak"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

crossing_C <- toy_vh$peak_identity_diagnostics |>
  filter(crossing_detected) |>
  pull(crossing_temp_C) |>
  first()

crossing_x <- 1000 / (crossing_C + 273.15)

p_after <- ggplot(data    = toy_vh$data,
                  mapping = aes(
                    x = invT_1000,
                    y = ln_k
                  )) +
  
  # Final fitted identity lines
  geom_line(
    data = toy_vh$predictions,
    aes(
      y     = ln_k_fit,
      color = peak,
      group = peak
    ),
    linewidth = 0.9
  ) +
  
  # Points colored by recovered identity,
  # but shaped by original observed elution order
  geom_point(
    aes(
      color = peak,
      shape = .vh_peak_elution_order
    ),
    size = 3
  ) +
  
  geom_vline(
    xintercept = crossing_x,
    linetype   = "dotted",
    linewidth  = 0.7
  ) +
  
  annotate(
    "text",
    x = crossing_x,
    y = Inf,
    label = sprintf(
      "Detected crossing\n%.1f °C",
      crossing_C
    ),
    vjust = 1.3,
    hjust = -0.08,
    size  = 3.2
  ) +
  
  scale_color_manual(
    values = cols_peak,
    labels = c(
      "Peak 1" = "Analyte A",
      "Peak 2" = "Analyte B"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Peak 1" = 16,
      "Peak 2" = 17
    ),
    labels = c(
      "Peak 1" = "First eluted peak",
      "Peak 2" = "Second eluted peak"
    )
  ) +
  
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  
  labs(
    title    = "(B) Recovered peak identity",
    subtitle = "Crossed assignment selected by SSE comparison",
    x        = expression(1000/italic(T)~(K^{-1})),
    y        = expression(ln(italic(k))),
    color    = "Recovered\nidentity",
    shape    = "Observed\nelution order"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.box       = "horizontal"
  )

plot_crossing <- grid.arrange(
  p_before, p_after,
  nrow = 1
)

ggsave(
  filename = "plot_crossing.svg",
  plot     =  plot_crossing,
  width    = 10,
  height   = 4.5
)
ggsave(
  filename = "plot_crossing.pdf",
  plot     =  plot_crossing,
  width    = 10,
  height   = 4.5
)
