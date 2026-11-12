## ---------------------------------------------------------------------------
## classic_vantHoff_analyzer.R -- Folder-level HPLC thermodynamic analysis.
##
## Repository: https://github.com/mirzahosseini-arash-semmelweis/classicVantHoff
## ---------------------------------------------------------------------------

VH_ANALYZER_VERSION <- "1.1.0"

suppressWarnings(
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(ggplot2)
    library(signs)
    library(rlang)
  })
)

# ---------- Prepare data -----------------------------------------------------

# Join a user-supplied lookup table (t0 or phase ratio).
# The table may be keyed by ANY subset of `condition_cols` (e.g. only
# columnID) and, if `allow_temp = TRUE`, additionally by temp_C, so that
# column-specific and temperature-dependent hold-up times can be supplied.
.vh_join_lookup <- function(
    x,
    tab,
    value_col,
    new_name,
    condition_cols,
    arg_name,
    allow_temp = TRUE
) {

  candidate_keys <- if (allow_temp) c(condition_cols, "temp_C") else condition_cols
  key_cols       <- intersect(candidate_keys, names(tab))

  if (!value_col %in% names(tab) || length(key_cols) == 0L) {
    stop(
      "`", arg_name, "` must contain '", value_col, "' and at least one key ",
      "column out of: ", paste(candidate_keys, collapse = ", ")
    )
  }

  tab <- tab |>
    select(all_of(c(key_cols, value_col))) |>
    distinct()

  n_dup <- tab |>
    count(across(all_of(key_cols))) |>
    filter(n > 1L) |>
    nrow()

  if (n_dup > 0L) {
    stop("`", arg_name, "` contains more than one value for ", n_dup,
         " key combination(s) of: ", paste(key_cols, collapse = ", "))
  }

  if (new_name != value_col) {
    tab <- tab |> rename(!!new_name := all_of(value_col))
  }

  x <- x |>
    select(-any_of(new_name)) |>
    left_join(tab, by = key_cols)

  n_miss <- sum(!is.finite(x[[new_name]]))
  if (n_miss > 0L) {
    warning(n_miss, " row(s) have no finite value in `", arg_name,
            "` after joining by ", paste(key_cols, collapse = ", "), ".")
  }

  x
}

vh_assign_peak_identity_by_area <- function(
    data,
    peak_area_ratio,
    group_cols,
    area_cols      = c("area1", "area2"),
    min_log_margin = 0.2 # rougly ±20% tolerance
) {
  
  if (
    !is.numeric(peak_area_ratio)  ||
    length(peak_area_ratio) != 1L ||
    !is.finite(peak_area_ratio)   ||
    peak_area_ratio <= 0
  ) {
    stop("`peak_area_ratio` must be a single finite positive number.")
  }
  
  if (length(area_cols) != 2L || !all(area_cols %in% names(data))) {
    stop(
      "`area_cols` must identify two existing peak-area columns. Missing: ",
      paste(setdiff(area_cols, names(data)), collapse = ", ")
    )
  }
  
  a1 <- area_cols[1]
  a2 <- area_cols[2]
  
  # ------------------------------------------------------------
  # Determine ONE identity decision per physical run.
  # Multiple wavelengths therefore cannot receive conflicting labels.
  # ------------------------------------------------------------
  
  identity_by_run <- data |>
    mutate(
      .vh_area_ratio_row = case_when(
        is.finite(.data[[a1]]) &
          is.finite(.data[[a2]]) &
          .data[[a1]] > 0 &
          .data[[a2]] > 0 ~
          .data[[a1]] / .data[[a2]],
        
        TRUE ~ NA_real_
      )
    ) |>
    group_by(across(all_of(group_cols))) |>
    summarise(
      .vh_area_ratio_observed =
        if (any(is.finite(.vh_area_ratio_row))) {
          # Median on logarithmic scale = robust geometric centre
          exp(
            median(
              log(.vh_area_ratio_row[is.finite(.vh_area_ratio_row)]),
              na.rm = TRUE
            )
          )
        } else {
          NA_real_
        },
      
      .vh_area_ratio_n = sum(is.finite(.vh_area_ratio_row)),
      
      .groups = "drop"
    ) |>
    mutate(
      .vh_peak_area_ratio_target = peak_area_ratio,
      
      .vh_area_identity_available =
        is.finite(.vh_area_ratio_observed) &
        .vh_area_ratio_observed > 0,
      
      # Error if current elution-order assignment is retained
      .vh_area_error_keep = if_else(
        .vh_area_identity_available,
        abs(log(.vh_area_ratio_observed) - log(peak_area_ratio)),
        NA_real_
      ),
      
      # Error if Peak 1 and Peak 2 identities are swapped
      .vh_area_error_swap = if_else(
        .vh_area_identity_available,
        abs(log(1 / .vh_area_ratio_observed) - log(peak_area_ratio)),
        NA_real_
      ),
      
      # How strongly one orientation is preferred
      .vh_area_identity_margin = abs(
        .vh_area_error_keep - .vh_area_error_swap
      ),
      
      .vh_area_identity_ambiguous =
        .vh_area_identity_available &
        .vh_area_identity_margin <= min_log_margin,
      
      .vh_area_identity_swapped =
        .vh_area_identity_available &
        !.vh_area_identity_ambiguous &
        .vh_area_error_swap < .vh_area_error_keep,
      
      .vh_area_identity_used =
        .vh_area_identity_available &
        !.vh_area_identity_ambiguous
    )
  
  data |>
    left_join(
      identity_by_run,
      by = group_cols
    )
}

vh_prepare <- function(
    data,
    condition_cols = c(
      "compoundname",
      "columnID",
      "eluentID",
      "modifierID",
      "flow_mL_min"
    ),
    tr_pattern                   = "^tR[0-9]+$",
    k_cols                       = NULL,
    run_id_col                   = "logical_key",
    collapse_wavelengths         = TRUE,
    strict_qc                    = TRUE,
    t0_external                  = NULL,
    phase_ratio                  = NULL,
    peak_area_ratio              = NULL,
    area_cols                    = c("area1", "area2"),
    area_identity_min_log_margin = 0.2
) {
  
  x <- data
  
  # ------------------------------------------------------------
  # Basic QC
  # ------------------------------------------------------------
  
  if ("excluded" %in% names(x)) {
    x <- x |> filter(is.na(excluded) | !excluded)
  }
  
  if (strict_qc) {
    if ("peak_detection_ok" %in% names(x)) {
      x <- x |> filter(peak_detection_ok %in% TRUE)
    }
    if ("deconvolution_failed_any" %in% names(x)) {
      x <- x |> filter(is.na(deconvolution_failed_any) | !deconvolution_failed_any)
    }
  }
  
  # ------------------------------------------------------------
  # Handle External t0 & Dynamically Calculate retention factor k
  # ------------------------------------------------------------
  
  tr_cols <- grep(tr_pattern, names(x), value = TRUE)
  
  if (!is.null(t0_external)) {
    if (is.numeric(t0_external) && length(t0_external) == 1L) {
      x <- x |> mutate(t0 = t0_external)
    } else if (is.data.frame(t0_external)) {
      x <- .vh_join_lookup(
        x              = x,
        tab            = t0_external,
        value_col      = "t0",
        new_name       = "t0",
        condition_cols = condition_cols,
        arg_name       = "t0_external",
        allow_temp     = TRUE
      )
    } else {
      stop("`t0_external` must be a single number or a data frame.")
    }
    
    # Calculate k dynamically for all tR columns found
    if (length(tr_cols) > 0 && "t0" %in% names(x)) {
      x <- x |> 
        mutate(across(all_of(tr_cols), ~ (.x - t0) / t0, .names = "k_{.col}"))
      
      # Map newly created k columns
      k_cols <- paste0("k_", tr_cols)
    }
  }
  
  # Fallback to default k_cols if no tR dynamic calculation took place
  if (is.null(k_cols)) {
    k_cols <- grep("^k[0-9]+$", names(x), value = TRUE)
  }
  
  if (length(k_cols) == 0) {
    stop("No retention time (tR) or retention factor (k) columns were detected in the input data.")
  }
  
  # ------------------------------------------------------------
  # Handle Phase Ratio (phi)
  # ------------------------------------------------------------
  
  if (!is.null(phase_ratio)) {
    if (is.numeric(phase_ratio) && length(phase_ratio) == 1L) {
      x <- x |> mutate(phase_ratio = phase_ratio)
    } else if (is.data.frame(phase_ratio)) {
      x <- .vh_join_lookup(
        x              = x,
        tab            = phase_ratio,
        value_col      = "phi",
        new_name       = "phase_ratio",
        condition_cols = condition_cols,
        arg_name       = "phase_ratio",
        allow_temp     = FALSE
      )
    }
  }
  
  # ------------------------------------------------------------
  # Optional area-ratio-based identity assignment
  # ------------------------------------------------------------
  
  run_group_cols <- unique(c(
    condition_cols,
    "temp_C",
    run_id_col
  ))
  
  if (!is.null(peak_area_ratio)) {
    
    if (length(k_cols) != 2L) {
      stop(
        "`peak_area_ratio` can only be used when exactly two ",
        "retention-factor columns are present."
      )
    }
    
    x <- vh_assign_peak_identity_by_area(
      data            = x,
      peak_area_ratio = peak_area_ratio,
      group_cols      = run_group_cols,
      area_cols       = area_cols,
      min_log_margin  = area_identity_min_log_margin
    )
    
  } else {
    
    x <- x |>
      mutate(
        .vh_peak_area_ratio_target  = NA_real_,
        .vh_area_ratio_observed     = NA_real_,
        .vh_area_ratio_n            = 0L,
        .vh_area_identity_available = FALSE,
        .vh_area_error_keep         = NA_real_,
        .vh_area_error_swap         = NA_real_,
        .vh_area_identity_margin    = NA_real_,
        .vh_area_identity_ambiguous = FALSE,
        .vh_area_identity_swapped   = FALSE,
        .vh_area_identity_used      = FALSE
      )
  }
  
  # ------------------------------------------------------------
  # Dynamic k -> long format conversion
  # ------------------------------------------------------------
  
  x <- x |>
    pivot_longer(
      cols      = all_of(k_cols),
      names_to  = "peak_raw",
      values_to = "k"
    ) |>
    mutate(
      .vh_peak_number = as.integer(gsub("[^0-9]", "", peak_raw)),
      
      # Identity label, not necessarily elution-order label
      peak = case_when(
        .vh_area_identity_swapped &
          .vh_peak_number == 1L ~ "Peak 2",
        
        .vh_area_identity_swapped &
          .vh_peak_number == 2L ~ "Peak 1",
        
        TRUE ~ paste("Peak", .vh_peak_number)
      ),
      
      T_K       = temp_C + 273.15,
      invT_1000 = 1000 / T_K
    ) |>
    filter(
      is.finite(temp_C),
      is.finite(T_K),
      is.finite(k),
      k > 0
    ) |>
    mutate(
      ln_k = log(k)
    )
  
  # ------------------------------------------------------------
  # Collapse wavelength channels from the same physical run
  # ------------------------------------------------------------
  
  if (collapse_wavelengths) {
    
    
    has_phase_ratio   <- "phase_ratio"   %in% names(x)
    has_wavelength_nm <- "wavelength_nm" %in% names(x)
    
    group_run <- unique(c(
      condition_cols,
      "temp_C",
      run_id_col,
      "peak"
    ))
    
    x <- x |>
      group_by(across(all_of(group_run))) |>
      summarise(
        T_K       = first(T_K),
        invT_1000 = first(invT_1000),
        
        # Retain phase_ratio column if it exists in data
        phase_ratio = if (.env$has_phase_ratio) {
          first(phase_ratio)
        } else {
          NA_real_
        },
        
        k         = median(k, na.rm = TRUE),
        ln_k      = median(ln_k, na.rm = TRUE),
        
        k_wavelength_sd = if (n() > 1) {
          sd(k, na.rm = TRUE)
        } else {
          NA_real_
        },
        
        n_wavelengths_used = if (.env$has_wavelength_nm) {
          n_distinct(wavelength_nm)
        } else {
          n()
        },
        
        peak_area_ratio_target = first(.vh_peak_area_ratio_target),
        
        area_ratio_observed = if (any(is.finite(.vh_area_ratio_observed))) {
            median(.vh_area_ratio_observed[is.finite(.vh_area_ratio_observed)], na.rm = TRUE)
          } else {
            NA_real_
          },
        
        area_identity_available = any(.vh_area_identity_available),
        area_identity_used      = any(.vh_area_identity_used),
        area_identity_swapped   = any(.vh_area_identity_swapped),
        area_identity_ambiguous = any(.vh_area_identity_ambiguous),
        
        .groups = "drop"
      )
  }
  
  # Unique van't Hoff series identifier
  x <- x |>
    group_by(
      across(all_of(c(condition_cols, "peak")))
    ) |>
    mutate(
      series_id = cur_group_id()
    ) |>
    ungroup()
  
  x
}

# ---------- Extract summaries ------------------------------------------------

vh_temperature_summary <- function(
    data,
    condition_cols = c(
      "compoundname",
      "columnID",
      "eluentID",
      "modifierID",
      "flow_mL_min"
    )
) {
  
  has_phase_ratio <- "phase_ratio" %in% names(data)
  
  group_cols <- c(
    condition_cols,
    "peak",
    "temp_C"
  )
  
  data |>
    group_by(across(all_of(group_cols))) |>
    summarise(
      T_K       = first(T_K),
      invT_1000 = first(invT_1000),
      
      phase_ratio = if (.env$has_phase_ratio) {
        first(phase_ratio)
      } else {
        NA_real_
      },
      
      sd_ln_k = if (n() > 1) {
        sd(ln_k, na.rm = TRUE)
      } else {
        NA_real_
      },
      ln_k    = mean(ln_k, na.rm = TRUE),

      k      = exp(ln_k),
      n_runs = n(),
      
      peak_area_ratio_target = first(peak_area_ratio_target),
      
      area_ratio_observed = if (any(is.finite(area_ratio_observed))) {
          median(area_ratio_observed, na.rm = TRUE)
        } else {
          NA_real_
        },
      
      area_identity_available = any(area_identity_available),
      area_identity_used      = any(area_identity_used),
      area_identity_swapped   = any(area_identity_swapped),
      area_identity_ambiguous = any(area_identity_ambiguous),
      
      .groups = "drop"
    ) |>
    group_by(
      across(all_of(c(condition_cols, "peak")))
    ) |>
    mutate(
      series_id = cur_group_id()
    ) |>
    ungroup()
}

# ---------- Fitter wrappers --------------------------------------------------

vh_complete_missing_peak2 <- function(
    data,
    condition_cols,
    enabled                  = TRUE,
    min_peak2_observed_temps = 1L
) {
  
  if (!".vh_peak2_imputed" %in% names(data)) {
    data$.vh_peak2_imputed <- FALSE
  }
  
  if (!enabled) {
    return(data)
  }
  
  data |>
    group_by(across(all_of(condition_cols))) |>
    group_modify(
      function(d, key) {
        
        # Temperatures with a finite observed Peak 1
        t_p1 <- d |>
          filter(
            peak == "Peak 1",
            is.finite(ln_k)
          ) |>
          pull(temp_C) |>
          unique()
        
        # Temperatures with a finite observed Peak 2
        t_p2 <- d |>
          filter(
            peak == "Peak 2",
            is.finite(ln_k),
            !.vh_peak2_imputed
          ) |>
          pull(temp_C) |>
          unique()
        
        # Only infer missing Peak 2 if the condition genuinely
        # contains Peak 2 at another temperature.
        if (length(t_p2) < min_peak2_observed_temps) {
          return(d)
        }
        
        missing_t <- setdiff(t_p1, t_p2)
        
        if (length(missing_t) == 0L) {
          return(d)
        }
        
        add <- d |>
          filter(
            peak == "Peak 1",
            temp_C %in% missing_t,
            is.finite(ln_k)
          ) |>
          mutate(
            peak = "Peak 2",
            .vh_peak2_imputed = TRUE
          )
        
        bind_rows(d, add)
      }
    ) |>
    ungroup()
}

vh_fit_one <- function(
    d,
    min_temps            = 3,
    min_span_C           = 10,
    r2_warn              = 0.95,
    sigma_ref            = NA_real_,
    excess_scatter_alpha = 0.01,
    t0_rel_uncertainty   = NULL,     # e.g. 0.02 = 2 % standard uncertainty
    conf_level           = 0.95
) {

  Rgas <- 8.314462618# J / (mol * K)

  n_temp    <- n_distinct(d$temp_C)
  temp_span <- diff(range(d$temp_C))

  if (n_temp < min_temps) {
    return(
      list(
        model = NULL,
        stats = tibble(
          n_temp                        = n_temp,
          n_runs                        = sum(d$n_runs),
          temp_min_C                    = min(d$temp_C),
          temp_max_C                    = max(d$temp_C),
          temp_span_C                   = temp_span,
          intercept                     = NA_real_,
          intercept_se                  = NA_real_,
          slope                         = NA_real_,
          slope_se                      = NA_real_,
          slope_p                       = NA_real_,
          deltaH_kJ_mol                 = NA_real_,
          deltaH_se_kJ_mol              = NA_real_,
          deltaH_ci_low_kJ_mol          = NA_real_,
          deltaH_ci_high_kJ_mol         = NA_real_,
          deltaH_t0_sens_kJ_mol_per_pct = NA_real_,
          deltaH_se_total_kJ_mol        = NA_real_,
          deltaS_J_mol_K                = NA_real_,
          deltaS_se_J_mol_K             = NA_real_,
          deltaG_kJ_mol                 = NA_real_,
          deltaG_se_kJ_mol              = NA_real_,
          r2                            = NA_real_,
          adj_r2                        = NA_real_,
          rmse_ln_k                     = NA_real_,
          scatter_ratio                 = NA_real_,
          scatter_p                     = NA_real_,
          status                        = "too_few_temperatures"
        )
      )
    )
  }
  
  # ------------------------------------------------------------
  # Classical linear van't Hoff model
  # ------------------------------------------------------------
  
  mod <- lm(
    ln_k ~ invT_1000,
    data = d
  )
  
  sm <- summary(mod)
  co <- coef(sm)
  
  intercept    <- co["(Intercept)", "Estimate"]
  intercept_se <- co["(Intercept)", "Std. Error"]
  slope        <- co["invT_1000", "Estimate"]
  slope_se     <- co["invT_1000", "Std. Error"]
  
  # invT_1000 = 1000/T => slope = -deltaH / (1000 * R)
  deltaH    <- -slope * 1000 * Rgas / 1000 # returns kJ/mol
  deltaH_se <- slope_se * 1000 * Rgas / 1000

  # Significance of the slope and t-based confidence interval of deltaH.
  df_res  <- df.residual(mod)
  slope_p <- co["invT_1000", "Pr(>|t|)"]
  t_q     <- if (df_res > 0) qt(1 - (1 - conf_level) / 2, df_res) else NA_real_
  deltaH_ci_low  <- deltaH - t_q * deltaH_se
  deltaH_ci_high <- deltaH + t_q * deltaH_se

  # First-order propagation of a hold-up-time error into deltaH.
  # For k = (tR - t0)/t0: d ln k / d ln t0 = -(1 + k)/k. A relative t0 error
  # that is common to all temperatures of a series shifts each ln k by a
  # k-dependent amount; regressing that shift on 1000/T gives the change of
  # the slope, hence of deltaH, per unit relative t0 error.
  g_t0 <- -(1 + d$k) / d$k
  dslope_dlnt0 <- if (all(is.finite(g_t0))) {
    unname(coef(lm(g_t0 ~ d$invT_1000))[2])
  } else {
    NA_real_
  }
  deltaH_t0_sens <- -Rgas * dslope_dlnt0 * 0.01 # kJ/mol per +1 % in t0

  deltaH_se_total <- if (!is.null(t0_rel_uncertainty) && is.finite(deltaH_t0_sens)) {
    sqrt(deltaH_se^2 + (deltaH_t0_sens * 100 * t0_rel_uncertainty)^2)
  } else {
    deltaH_se
  }

  deltaS       <- NA_real_
  deltaS_se    <- NA_real_
  deltaG       <- NA_real_
  deltaG_se    <- NA_real_
  
  # Check if phase_ratio is present in data subset d
  if ("phase_ratio" %in% names(d) && !all(is.na(d$phase_ratio))) {
    phi <- d$phase_ratio[1L]
    
    # Entropy: J / (mol * K)
    deltaS    <- (intercept - log(phi)) * Rgas
    deltaS_se <- intercept_se * Rgas
    
    # Gibbs Free Energy at std temp (kJ/mol)
    T_std_K <- 298.15
    deltaG  <- deltaH - (T_std_K * (deltaS / 1000))
    
    # Variance-Covariance matrix for exact SE propagation:
    cov_mat <- vcov(mod)
    var_b0  <- cov_mat["(Intercept)", "(Intercept)"]
    var_b1  <- cov_mat["invT_1000", "invT_1000"]
    cov_b01 <- cov_mat["(Intercept)", "invT_1000"]
    
    invT_std_1000 <- 1000 / T_std_K
    var_ln_k_std  <- var_b0 + (invT_std_1000^2) * var_b1 + (2 * invT_std_1000 * cov_b01)
    
    deltaG_se <- (Rgas * T_std_K * sqrt(var_ln_k_std)) / 1000
  }
  
  # ------------------------------------------------------------
  # Automated warning flags
  # ------------------------------------------------------------
  
  flags <- character()

  if (n_temp == 3)            {flags <- c(flags, "only_3_temperatures")}
  if (temp_span < min_span_C) {flags <- c(flags, "narrow_temperature_span")}
  if (sm$r.squared < r2_warn) {flags <- c(flags, "low_R2")}

  # Baseline-aware lack-of-fit check. R^2 depends on the slope
  # magnitude and is a poor linearity criterion for 3-4 points. Instead, the
  # residual variance of the series is compared with the dataset-level
  # reference variance: (n - 2) s^2 / sigma_ref^2 ~ chi^2(n - 2).
  scatter_ratio <- NA_real_
  scatter_p     <- NA_real_

  if (is.finite(sigma_ref) && sigma_ref > 0 && df_res > 0) {
    scatter_ratio <- sigma(mod) / sigma_ref
    scatter_p     <- pchisq(df_res * scatter_ratio^2, df = df_res, lower.tail = FALSE)

    if (scatter_p < excess_scatter_alpha) {
      flags <- c(flags, "excess_scatter")
    }
  }

  status <- if (length(flags) == 0) {
    "ok"
  }
  else {
    paste(flags, collapse = "; ")
  }
  
  list(
    model = mod,
    
    stats = tibble(
      n_temp = n_temp,
      n_runs = sum(d$n_runs),
      
      temp_min_C  = min(d$temp_C),
      temp_max_C  = max(d$temp_C),
      temp_span_C = temp_span,
      
      intercept    = intercept,
      intercept_se = intercept_se,
      slope        = slope,
      slope_se     = slope_se,
      slope_p      = slope_p,

      deltaH_kJ_mol    = deltaH,
      deltaH_se_kJ_mol = deltaH_se,

      deltaH_ci_low_kJ_mol  = deltaH_ci_low,
      deltaH_ci_high_kJ_mol = deltaH_ci_high,

      deltaH_t0_sens_kJ_mol_per_pct = deltaH_t0_sens,
      deltaH_se_total_kJ_mol        = deltaH_se_total,

      deltaS_J_mol_K    = deltaS,
      deltaS_se_J_mol_K = deltaS_se,
      
      deltaG_kJ_mol     = deltaG,
      deltaG_se_kJ_mol  = deltaG_se,
      
      r2        = sm$r.squared,
      adj_r2    = sm$adj.r.squared,
      rmse_ln_k = sigma(mod),

      scatter_ratio = scatter_ratio,
      scatter_p     = scatter_p,

      status = status
    )
  )
}

vh_flag_outliers_one <- function(
    d,
    min_temps       = 3,
    rstudent_cutoff = 3,
    max_outliers    = 1,
    sigma_ref       = NA_real_,
    ref_cutoff      = 3
) {

  # With 4 temperatures the externally studentized residual
  # has a single residual degree of freedom (t distribution with 1 df), so
  # |r*| > 3 alone is exceeded by the most deviating point of ~73 % of
  # perfectly linear, normally distributed series. When `sigma_ref` (a
  # dataset-level reference scale of ln k residuals) is finite, a point is
  # only excluded if it is ALSO improbable relative to that reference:
  #   |e_i| / (sigma_ref * sqrt(1 - h_ii)) > ref_cutoff.
  # sigma_ref = NA reproduces the version 1.0 behaviour.

  d <- d |>
    mutate(
      .vh_outlier           = FALSE,
      .vh_outlier_iteration = NA_integer_,
      .vh_rstudent          = NA_real_,
      .vh_z_ref             = NA_real_, # [NEW v1.1]
      .vh_cooks             = NA_real_
    )
  
  if (
    nrow(d) < 4              ||
    n_distinct(d$temp_C) < 4 ||
    max_outliers < 1
  ) {
    
    return(
      d |>
        mutate(.vh_used = TRUE)
    )
  }
  
  active    <- rep(TRUE, nrow(d))
  iteration <- 0L
  
  repeat {
    
    idx <- which(active)
    
    # Protect minimum final model size
    if (
      length(idx) <= min_temps               ||
      n_distinct(d$temp_C[idx]) <= min_temps ||
      iteration >= max_outliers
    ) {
      break
    }
    
    mod <- lm(
      ln_k ~ invT_1000,
      data = d[idx, , drop = FALSE]
    )
    
    rs <- suppressWarnings(rstudent(mod))
    cd <- cooks.distance(mod)

    d$.vh_rstudent[idx] <- rs
    d$.vh_cooks[idx]    <- cd

    # Standardized residual relative to the reference scale
    use_ref <- is.finite(sigma_ref) && sigma_ref > 0
    z_ref   <- if (use_ref) {
      resid(mod) / (sigma_ref * sqrt(1 - hatvalues(mod)))
    } else {
      rep(NA_real_, length(idx))
    }
    d$.vh_z_ref[idx] <- z_ref

    score                    <- abs(rs)
    score[!is.finite(score)] <- -Inf

    # Only points passing the reference gate are candidates
    if (use_ref) {
      score[!(abs(z_ref) > ref_cutoff)] <- -Inf
    }

    j_local <- which.max(score)

    if (
      length(j_local) == 0 ||
      score[j_local] <= rstudent_cutoff
    ) {
      break
    }
    
    j <- idx[j_local]
    
    if (sum(active) - 1 < min_temps) {
      break
    }
    
    iteration <- iteration + 1L
    
    d$.vh_outlier[j]           <- TRUE
    d$.vh_outlier_iteration[j] <- iteration
    
    active[j] <- FALSE
  }
  
  d |>
    mutate(
      .vh_used = !.vh_outlier
    )
}

vh_crossing_test_one <- function(
    d,
    sse_ratio_cutoff = 0.5,
    min_temps        = 4L,
    area_ratio_target_cols = c(
      "peak_area_ratio_target",
      ".vh_peak_area_ratio_target"
    )
) {
  
  # ------------------------------------------------------------
  # Was external area-based identity information supplied?
  #
  # If yes, crossing inference from elution order is unnecessary.
  # ------------------------------------------------------------
  
  area_cols_present <- intersect(
    area_ratio_target_cols,
    names(d)
  )
  
  area_identity_supplied <- FALSE
  
  if (length(area_cols_present) > 0) {
    area_identity_supplied <- any(
      vapply(
        area_cols_present,
        function(nm) {
          any(is.finite(d[[nm]]))
        },
        logical(1)
      )
    )
  }
  
  # Defaults
  diagnostic <- tibble(
    crossing_test_status = "not_evaluated",
    crossing_sse_order   = NA_real_,
    crossing_sse_switch  = NA_real_,
    crossing_sse_ratio   = NA_real_,
    crossing_temp_C      = NA_real_,
    crossing_detected    = FALSE,
    crossing_n_eval      = 0L,
    peak2_imputed_n      = sum(d$.vh_peak2_imputed, na.rm = TRUE)
  )
  
  d <- d |>
    mutate(
      .vh_peak_elution_order  = peak,
      .vh_crossing_relabelled = FALSE
    )
  
  if (area_identity_supplied) {
    diagnostic$crossing_test_status <- "skipped_external_area_identity"
    
    return(
      list(
        data       = d,
        diagnostic = diagnostic
      )
    )
  }
  
  # ------------------------------------------------------------
  # Need Peak 1 and Peak 2
  # ------------------------------------------------------------
  
  d2 <- d |>
    filter(
      peak %in% c("Peak 1", "Peak 2"),
      is.finite(temp_C),
      is.finite(invT_1000),
      is.finite(ln_k)
    ) |>
    mutate(
      peak_short = recode(
        peak,
        "Peak 1" = "p1",
        "Peak 2" = "p2"
      ),
      imputed = .vh_peak2_imputed
    )
  
  if (n_distinct(d2$temp_C) < min_temps) {
    diagnostic$crossing_test_status <- "not_evaluable_too_few_temperatures"
    
    return(
      list(
        data       = d,
        diagnostic = diagnostic
      )
    )
  }
  
  # ------------------------------------------------------------
  # One row per temperature
  # ------------------------------------------------------------
  
  w <- d2 |>
    select(
      temp_C,
      invT_1000,
      peak_short,
      ln_k,
      imputed
    ) |>
    pivot_wider(
      names_from  = peak_short,
      values_from = c(ln_k, imputed)
    ) |>
    arrange(temp_C)
  
  if (
    !all(c(
      "ln_k_p1",
      "ln_k_p2",
      "imputed_p2"
    ) %in% names(w))
  ) {
    diagnostic$crossing_test_status <- "not_evaluable_incomplete_pair"
    
    return(
      list(
        data       = d,
        diagnostic = diagnostic
      )
    )
  }
  
  i_lo <- 1L
  i_hi <- nrow(w)
  
  # ------------------------------------------------------------
  # The four corners must be OBSERVED.
  #
  # An imputed Peak 2 may contribute at intermediate temperatures,
  # but must not define the endpoint connection hypotheses.
  # ------------------------------------------------------------
  
  endpoints_ok <- all(
    is.finite(c(
      w$ln_k_p1[i_lo],
      w$ln_k_p2[i_lo],
      w$ln_k_p1[i_hi],
      w$ln_k_p2[i_hi]
    ))
  ) &&
    !isTRUE(w$imputed_p2[i_lo]) &&
    !isTRUE(w$imputed_p2[i_hi])
  
  if (!endpoints_ok) {
    diagnostic$crossing_test_status <- "not_evaluable_endpoint_peak2_missing"
    
    return(
      list(
        data       = d,
        diagnostic = diagnostic
      )
    )
  }
  
  # ------------------------------------------------------------
  # Helper: straight line through two endpoints
  # ------------------------------------------------------------
  
  line_between <- function(x, x1, y1, x2, y2) {
    y1 + (y2 - y1) * (x - x1) / (x2 - x1)
  }
  
  x_lo  <- w$invT_1000[i_lo]
  x_hi  <- w$invT_1000[i_hi]
  
  p1_lo <- w$ln_k_p1[i_lo]
  p2_lo <- w$ln_k_p2[i_lo]
  
  p1_hi <- w$ln_k_p1[i_hi]
  p2_hi <- w$ln_k_p2[i_hi]
  
  # ============================================================
  # Hypothesis I:
  # retain elution-order identities
  # ============================================================
  
  pred_order_1 <- line_between(
    w$invT_1000,
    x_lo, p1_lo,
    x_hi, p1_hi
  )
  
  pred_order_2 <- line_between(
    w$invT_1000,
    x_lo, p2_lo,
    x_hi, p2_hi
  )
  
  # ============================================================
  # Hypothesis II:
  # identities cross
  # ============================================================
  
  pred_cross_A <- line_between(
    w$invT_1000,
    x_lo, p1_lo,
    x_hi, p2_hi
  )
  
  pred_cross_B <- line_between(
    w$invT_1000,
    x_lo, p2_lo,
    x_hi, p1_hi
  )
  
  # Under the crossing hypothesis, observed Peak 1 is always the
  # earlier/lower-k peak, irrespective of identity.
  pred_cross_p1 <- pmin(
    pred_cross_A,
    pred_cross_B
  )
  
  pred_cross_p2 <- pmax(
    pred_cross_A,
    pred_cross_B
  )
  
  # ------------------------------------------------------------
  # Only intermediate temperatures discriminate the models
  # ------------------------------------------------------------
  
  mid <- seq_len(nrow(w))
  mid <- mid[mid != i_lo & mid != i_hi]
  
  if (length(mid) == 0L) {
    diagnostic$crossing_test_status <- "not_evaluable_no_intermediate_points"
    
    return(
      list(
        data       = d,
        diagnostic = diagnostic
      )
    )
  }
  
  sq <- function(obs, pred) {
    ok <- is.finite(obs) & is.finite(pred)
    
    list(
      sse = sum((obs[ok] - pred[ok])^2),
      n   = sum(ok)
    )
  }
  
  q01 <- sq(w$ln_k_p1[mid], pred_order_1[mid])
  q02 <- sq(w$ln_k_p2[mid], pred_order_2[mid])
  q11 <- sq(w$ln_k_p1[mid], pred_cross_p1[mid])
  q12 <- sq(w$ln_k_p2[mid], pred_cross_p2[mid])
  
  sse_order  <- q01$sse + q02$sse
  sse_switch <- q11$sse + q12$sse
  
  n_eval <- q01$n + q02$n
  
  if (n_eval == 0L) {
    diagnostic$crossing_test_status <- "not_evaluable_no_finite_intermediate_points"
    
    return(
      list(
        data       = d,
        diagnostic = diagnostic
      )
    )
  }
  
  # ------------------------------------------------------------
  # SSE ratio
  # ------------------------------------------------------------
  
  sse_ratio <- if (
    is.finite(sse_order) &&
    sse_order > .Machine$double.eps
  ) {
    sse_switch / sse_order
  } else {
    Inf
  }
  
  # ------------------------------------------------------------
  # Exact crossing of the two crossed corner-lines
  # ------------------------------------------------------------
  
  m_A <- (p2_hi - p1_lo) / (x_hi - x_lo)
  b_A <- p1_lo - m_A * x_lo
  m_B <- (p1_hi - p2_lo) / (x_hi - x_lo)
  b_B <- p2_lo - m_B * x_lo
  
  x_cross <- if (
    is.finite(m_A - m_B) &&
    abs(m_A - m_B) > .Machine$double.eps
  ) {
    (b_B - b_A) / (m_A - m_B)
  } else {
    NA_real_
  }
  
  T_cross_K <- if (
    is.finite(x_cross) &&
    x_cross > 0
  ) {
    1000 / x_cross
  } else {
    NA_real_
  }
  
  temp_cross_C <- T_cross_K - 273.15
  
  crossing_inside <- is.finite(temp_cross_C) && temp_cross_C > min(w$temp_C) && temp_cross_C < max(w$temp_C)
  
  # ------------------------------------------------------------
  # Conservative decision:
  #
  # crossed assignment must reduce SSE by ratio cutoff
  # ------------------------------------------------------------
  
  crossing_detected <- crossing_inside && is.finite(sse_ratio) && sse_ratio < sse_ratio_cutoff
  
  diagnostic <- tibble(
    crossing_test_status =
      if (crossing_detected) {
        "crossing_detected"
      } else {
        "tested_no_crossing"
      },
    
    crossing_sse_order  = sse_order,
    crossing_sse_switch = sse_switch,
    crossing_sse_ratio  = sse_ratio,
    crossing_temp_C     = temp_cross_C,
    crossing_detected   = crossing_detected,
    crossing_n_eval     = n_eval,
    peak2_imputed_n     = sum(d$.vh_peak2_imputed, na.rm = TRUE)
  )
  
  # ------------------------------------------------------------
  # If accepted, relabel according to IDENTITY rather than order
  #
  # Identity Peak 1 is defined as the species corresponding to
  # observed Peak 1 at the LOWEST temperature.
  # ------------------------------------------------------------
  
  if (crossing_detected) {
    
    pred_A_all <- line_between(
      d$invT_1000,
      x_lo, p1_lo,
      x_hi, p2_hi
    )
    
    pred_B_all <- line_between(
      d$invT_1000,
      x_lo, p2_lo,
      x_hi, p1_hi
    )
    
    old_peak <- d$peak
    
    A_is_earlier <- pred_A_all <= pred_B_all
    
    d$peak <- case_when(
      
      old_peak == "Peak 1" &  A_is_earlier ~ "Peak 1",
      old_peak == "Peak 2" &  A_is_earlier ~ "Peak 2",
      old_peak == "Peak 1" & !A_is_earlier ~ "Peak 2",
      old_peak == "Peak 2" & !A_is_earlier ~ "Peak 1",
      
      TRUE ~ old_peak
    )
    
    d$.vh_crossing_relabelled <- d$peak != old_peak
  }
  
  list(
    data       = d,
    diagnostic = diagnostic
  )
}

vh_reconcile_two_peak_identity <- function(
    data,
    condition_cols,
    complete_missing_peak2    = TRUE,
    min_peak2_observed_temps  = 1L,
    detect_crossing           = TRUE,
    crossing_sse_ratio_cutoff = 0.5,
    crossing_min_temps        = 4L
) {
  
  # ------------------------------------------------------------
  # First: infer coelution where Peak 2 disappears
  # ------------------------------------------------------------
  
  data2 <- vh_complete_missing_peak2(
    data                     = data,
    condition_cols           = condition_cols,
    enabled                  = complete_missing_peak2,
    min_peak2_observed_temps = min_peak2_observed_temps
  )
  
  if (!detect_crossing) {

    data2 <- data2 |>
      select(-any_of("series_id")) |>
      group_by(
        across(
          all_of(c(condition_cols, "peak"))
        )
      ) |>
      mutate(
        series_id = cur_group_id()
      ) |>
      ungroup()

    diagnostics <- data2 |>
      group_by(across(all_of(condition_cols))) |>
      summarise(
        crossing_test_status = "disabled",
        crossing_detected    = FALSE,
        peak2_imputed_n      = sum(.vh_peak2_imputed %in% TRUE),
        .groups              = "drop"
      )

    return(
      list(
        data        = data2,
        diagnostics = diagnostics
      )
    )
  }
  
  # ------------------------------------------------------------
  # One crossing test per fixed chromatographic condition
  # ------------------------------------------------------------
  
  nested <- data2 |>
    group_by(
      across(all_of(condition_cols))
    ) |>
    nest() |>
    ungroup()
  
  nested <- nested |>
    mutate(
      reconciliation = map(
        data,
        vh_crossing_test_one,
        sse_ratio_cutoff = crossing_sse_ratio_cutoff,
        min_temps        = crossing_min_temps
      ),
      
      data = map(
        reconciliation,
        "data"
      ),
      
      diagnostic = map(
        reconciliation,
        "diagnostic"
      )
    )
  
  corrected <- nested |>
    select(
      all_of(condition_cols),
      data
    ) |>
    unnest(data) |>
    select(-any_of("series_id")) |>
    group_by(
      across(
        all_of(c(condition_cols, "peak"))
      )
    ) |>
    mutate(
      series_id = cur_group_id()
    ) |>
    ungroup()
  
  diagnostics <- nested |>
    select(
      all_of(condition_cols),
      diagnostic
    ) |>
    unnest(diagnostic)
  
  list(
    data        = corrected,
    diagnostics = diagnostics
  )
}

# Dataset-level ("baseline") residual scale of ln k.
#
# Every series with >= `min_temps` observed temperatures is fitted by OLS and
# its residuals are standardized by sqrt(1 - h_ii), which gives them the
# common variance sigma^2 under the linear model. The pooled scale is the
# normal-consistent median absolute value (1.4826 * median |.|), which is
# insensitive to the few genuine artefacts present in the pool. Points that
# were imputed by coelution completion are not used.
vh_reference_scale <- function(
    data,
    series_cols,
    min_temps = 4L
) {

  d <- data
  if (".vh_peak2_imputed" %in% names(d)) {
    d <- d |> filter(!.vh_peak2_imputed %in% TRUE)
  }

  z <- d |>
    filter(is.finite(ln_k), is.finite(invT_1000)) |>
    group_by(across(all_of(series_cols))) |>
    group_modify(function(s, key) {
      if (n_distinct(s$temp_C) < min_temps) return(tibble(z = numeric()))
      m <- lm(ln_k ~ invT_1000, data = s)
      tibble(z = resid(m) / sqrt(1 - hatvalues(m)))
    }) |>
    ungroup() |>
    pull(z)

  z <- z[is.finite(z)]

  list(
    sigma_ref = if (length(z) >= 10L) 1.4826 * median(abs(z)) else NA_real_,
    n_points  = length(z)
  )
}

vh_fit <- function(
    data,
    condition_cols = c(
      "compoundname",
      "columnID",
      "eluentID",
      "modifierID",
      "flow_mL_min"
    ),
    min_temps                 = 3,
    min_span_C                = 10,
    r2_warn                   = 0.95,
    prediction_points         = 100,
    rstudent_cutoff           = 3,
    max_outliers              = 1,
    complete_missing_peak2    = TRUE,
    min_peak2_observed_temps  = 1L,
    detect_peak_crossing      = TRUE,
    crossing_sse_ratio_cutoff = 0.5,
    crossing_min_temps        = 4L,
    outlier_scale             = "pooled", # "pooled", a number (ln k units), or "none"
    ref_cutoff                = 3,
    excess_scatter_alpha      = 0.01,
    t0_rel_uncertainty        = NULL,
    conf_level                = 0.95
) {

  # ------------------------------------------------------------
  # Two-peak identity reconciliation
  # ------------------------------------------------------------
  
  identity_check <- vh_reconcile_two_peak_identity(
    data                      = data,
    condition_cols            = condition_cols,
    complete_missing_peak2    = complete_missing_peak2,
    min_peak2_observed_temps  = min_peak2_observed_temps,
    detect_crossing           = detect_peak_crossing,
    crossing_sse_ratio_cutoff = crossing_sse_ratio_cutoff,
    crossing_min_temps        = crossing_min_temps
  )
  
  data <- identity_check$data
  
  peak_identity_diagnostics <- identity_check$diagnostics
  
  # ------------------------------------------------------------
  # Data grouping and fit
  # ------------------------------------------------------------
  
  series_cols <- c(
    condition_cols,
    "peak"
  )

  # ------------------------------------------------------------
  # Reference (baseline) residual scale
  # ------------------------------------------------------------

  sigma_ref_n <- NA_integer_

  sigma_ref <- if (is.numeric(outlier_scale) && length(outlier_scale) == 1L) {
    outlier_scale
  } else if (identical(outlier_scale, "pooled")) {
    rs <- vh_reference_scale(data, series_cols = series_cols, min_temps = 4L)
    sigma_ref_n <- rs$n_points
    rs$sigma_ref
  } else if (identical(outlier_scale, "none")) {
    NA_real_
  } else {
    stop("`outlier_scale` must be \"pooled\", \"none\" or a single positive number.")
  }

  nested <- data |>
    group_by(across(all_of(series_cols))) |>
    nest() |>
    ungroup()

  fitted <- nested |>
    mutate(

      # ----------------------------------------------------------
      # Outlier diagnostics within each van't Hoff series
      # ----------------------------------------------------------

      data = map(
        data,
        vh_flag_outliers_one,
        min_temps       = min_temps,
        rstudent_cutoff = rstudent_cutoff,
        max_outliers    = max_outliers,
        sigma_ref       = sigma_ref,
        ref_cutoff      = ref_cutoff
      ),

      # ----------------------------------------------------------
      # Data used for the FINAL model
      # ----------------------------------------------------------

      fit_data = map(
        data,
        ~ .x |>
          filter(!.vh_outlier)
      ),

      # ----------------------------------------------------------
      # Final van't Hoff model
      # ----------------------------------------------------------

      fit_object = map(
        fit_data,
        vh_fit_one,
        min_temps            = min_temps,
        min_span_C           = min_span_C,
        r2_warn              = r2_warn,
        sigma_ref            = sigma_ref,
        excess_scatter_alpha = excess_scatter_alpha,
        t0_rel_uncertainty   = t0_rel_uncertainty,
        conf_level           = conf_level
      ),

      model = map(
        fit_object,
        "model"
      ),

      data = map2(
        data,
        model,
        function(d, m) {

          if (is.null(m)) {
            return(
              d |>
                mutate(
                  .vh_final_fitted = NA_real_,
                  .vh_final_resid = NA_real_
                )
            )
          }

          fitted_all <- predict(
            m,
            newdata = d
          )

          d |>
            mutate(
              .vh_final_fitted = fitted_all,
              .vh_final_resid = ln_k - fitted_all
            )
        }
      ),

      stats = map(
        fit_object,
        "stats"
      ),

      # QC counts
      n_temp_total = map_int(
        data,
        ~ n_distinct(.x$temp_C)
      ),

      n_outliers = map_int(
        data,
        ~ sum(.x$.vh_outlier, na.rm = TRUE)
      ),

      # Coelution-completed points used in the final model
      n_imputed = map_int(
        fit_data,
        ~ {if (".vh_peak2_imputed" %in% names(.x)) sum(.x$.vh_peak2_imputed %in% TRUE) else 0L}
      )
    ) |>
    select(-fit_object) |>
    unnest(stats)

  # ------------------------------------------------------------
  # Fitted lines for plotting
  # ------------------------------------------------------------

  predictions <- fitted |>
    filter(
      map_lgl(model, ~ !is.null(.x))
    ) |>
    mutate(
      prediction = map2(
        model,
        fit_data,
        function(m, d) {

          newdata <- tibble(
            invT_1000 = seq(
              min(d$invT_1000),
              max(d$invT_1000),
              length.out = prediction_points
            )
          )

          newdata$ln_k_fit <- predict(
            m,
            newdata = newdata
          )

          newdata
        }
      )
    ) |>
    select(
      all_of(series_cols),
      prediction
    ) |>
    unnest(prediction)

  # Recreate common series ID in both tables
  diagnostic_data <- fitted |>
    select(
      all_of(series_cols),
      data
    ) |>
    unnest(data)

  keys <- diagnostic_data |>
    distinct(
      across(all_of(series_cols)),
      series_id
    )

  predictions <- predictions |>
    left_join(
      keys,
      by = series_cols
    )

  results <- fitted |>
    select(
      -data,
      -fit_data,
      -model
    )
  
  if (nrow(peak_identity_diagnostics) > 0) {
    
    results <- results |>
      left_join(
        peak_identity_diagnostics,
        by = condition_cols
      )
  }
  
  results <- results |>
    rowwise() |>
    mutate(
      status = {
        
        fl <- if (
          is.na(status) ||
          status == "ok"
        ) {
          character()
        } else {
          strsplit(status, "; ", fixed = TRUE)[[1]]
        }
        
        if (isTRUE(get0("crossing_detected", ifnotfound = FALSE))) {
          fl <- c(fl, "peak_order_crossing_corrected")
        }

        if (n_imputed > 0) {
          fl <- c(fl, "peak2_coelution_imputed")
        }

        # Completed points are copies of Peak 1. When they are at
        # least as many as the independently observed Peak 2 points, the
        # Peak 2 slope is dominated by Peak 1 and Delta Delta H is biased
        # towards zero.
        if (n_imputed > 0 && is.finite(n_temp) && n_imputed >= n_temp - n_imputed) {
          fl <- c(fl, "peak2_mostly_imputed")
        }
        
        fl <- unique(fl)
        
        if (length(fl) == 0) {
          "ok"
        } else {
          paste(fl, collapse = "; ")
        }
      }
    ) |>
    ungroup()

  list(
    data                      = diagnostic_data,
    results                   = results,
    models                    = fitted,
    predictions               = predictions,
    peak_identity_diagnostics = peak_identity_diagnostics,
    condition_cols            = condition_cols,
    sigma_ref                 = sigma_ref,
    sigma_ref_n               = sigma_ref_n,
    version                   = VH_ANALYZER_VERSION
  )
}

# ---------- Diagnostics and summaries ----------------------------------------

# Fit-status summary: one row per category, counts of van't Hoff series.
vh_fit_status_summary <- function(vh) {

  r <- vh$results

  has_flag <- function(flag) grepl(flag, r$status, fixed = TRUE)
  fitted   <- is.finite(r$slope)

  imp_pts <- sum(r$n_imputed, na.rm = TRUE)

  tibble::tribble(
    ~category,                                               ~n_series,
    "Series (compound x column x eluent x peak)",            nrow(r),
    "  Peak 1 / Peak 2",                                     NA_integer_,
    "Observed temperatures: 4",                              sum(r$n_temp_total == 4),
    "Observed temperatures: 3",                              sum(r$n_temp_total == 3),
    "Observed temperatures: < 3 (not fitted)",               sum(r$n_temp_total < 3),
    "Fitted series",                                         sum(fitted),
    "  one point excluded as anomalous (4 -> 3)",            sum(r$n_outliers > 0),
    "  final fit on 4 temperatures",                         sum(fitted & r$n_temp == 4),
    "  final fit on 3 temperatures",                         sum(fitted & r$n_temp == 3),
    "  R^2 < r2_warn (low_R2)",                              sum(has_flag("low_R2")),
    "  excess residual scatter (excess_scatter)",            sum(has_flag("excess_scatter")),
    "  slope not significant (p > 0.05)",                    sum(fitted & r$slope_p > 0.05, na.rm = TRUE),
    "  Peak 2 series with coelution-completed points",       sum(has_flag("peak2_coelution_imputed")),
    "    of which completed points >= observed points",      sum(has_flag("peak2_mostly_imputed")),
    "  series relabelled after elution-order crossing",      sum(has_flag("peak_order_crossing_corrected")),
    "  status 'ok' (no flag)",                               sum(r$status == "ok"),
    "Coelution-completed Peak 2 points used in fits",        imp_pts
  ) |>
    mutate(
      n_series = if_else(
        category == "  Peak 1 / Peak 2",
        NA_integer_,
        as.integer(n_series)
      ),
      detail = if_else(
        category == "  Peak 1 / Peak 2",
        paste(sum(r$peak == "Peak 1"), "/", sum(r$peak == "Peak 2")),
        NA_character_
      )
    )
}

# Shared (condition-wide) residual check. Final-model residuals of retained
# observations of all series with >= `min_temps` temperatures are pooled by
# `by` x temperature. A median residual that is consistently signed across
# many independent series indicates a systematic, analyte-independent effect
# (e.g. a temperature-specific flow, temperature or hold-up deviation) rather
# than an analyte-specific anomaly. A two-sided sign test is reported.
vh_shared_residual_check <- function(
    vh,
    by        = "eluentID",
    min_temps = 4L
) {

  series_cols <- c(vh$condition_cols, "peak")

  keep <- vh$results |>
    filter(n_temp >= min_temps) |>
    select(all_of(series_cols))

  d <- vh$data |>
    semi_join(keep, by = series_cols) |>
    filter(!.vh_outlier, is.finite(.vh_final_resid))

  if (".vh_peak2_imputed" %in% names(d)) {
    d <- d |> filter(!.vh_peak2_imputed %in% TRUE)
  }

  d |>
    group_by(across(all_of(c(by, "temp_C")))) |>
    summarise(
      n_series      = n(),
      median_resid  = median(.vh_final_resid),
      mean_resid    = mean(.vh_final_resid),
      frac_positive = mean(.vh_final_resid > 0),
      sign_test_p   = if (n() > 0) {
        binom.test(sum(.vh_final_resid > 0), n(), 0.5)$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

# Series-by-series comparison of two vh_fit() results, e.g. with/without
# coelution completion, outlier screening or a different t0.
vh_compare_fits <- function(
    vh_ref,
    vh_alt,
    labels = c("ref", "alt")
) {

  series_cols <- c(vh_ref$condition_cols, "peak")
  keep_cols   <- c(series_cols, "n_temp", "n_outliers", "deltaH_kJ_mol",
                   "deltaH_se_kJ_mol", "r2", "status")

  a <- vh_ref$results |> select(any_of(keep_cols))
  b <- vh_alt$results |> select(any_of(keep_cols))

  full_join(
    a, b,
    by     = series_cols,
    suffix = paste0("_", labels)
  ) |>
    mutate(
      deltaH_diff_kJ_mol =
        .data[[paste0("deltaH_kJ_mol_", labels[2])]] -
        .data[[paste0("deltaH_kJ_mol_", labels[1])]]
    )
}

# Re-run the complete preparation and fitting chain for several hold-up time
# specifications (each element of `t0_scenarios` is anything accepted by
# vh_prepare(t0_external = ...)). Returns the stacked result tables.
vh_t0_sensitivity <- function(
    data,
    t0_scenarios,
    prepare_args = list(),
    fit_args     = list()
) {

  if (is.null(names(t0_scenarios)) || any(names(t0_scenarios) == "")) {
    stop("`t0_scenarios` must be a named list.")
  }

  purrr::imap_dfr(t0_scenarios, function(t0_spec, nm) {
    run <- do.call(vh_prepare, c(list(data = data, t0_external = t0_spec), prepare_args))
    tmp <- vh_temperature_summary(run)
    fit <- do.call(vh_fit, c(list(data = tmp), fit_args))
    fit$results |> mutate(t0_scenario = nm, .before = 1)
  })
}

# ---------- Fit plotter ------------------------------------------------------

plot_vanthoff <- function(
    vh,
    facet_vars  = "compoundname",
    facet_nrow  = NULL,
    facet_ncol  = NULL,
    colour_by   = "columnID"
) {
  
  peak_levels <- unique(vh$results$peak)
  
  d <- vh$data |>
    mutate(
      point_type = if_else(
        .vh_outlier,
        paste0(peak, " excluded"),
        peak
      ),
      point_type = factor(
        point_type,
        levels = c(
          peak_levels,
          paste0(peak_levels, " excluded")
        )
      )
    )
  
  pred <- vh$predictions
  
  facet_layer <- if (length(facet_vars) >= 2) {
    facet_grid(
      rows     = vars(!!sym(facet_vars[1])),
      cols     = vars(!!sym(facet_vars[2])),
      scales   = "free_y",
      switch   = "y",
      labeller = labeller(.multi_line = FALSE)
    )
  } else {
    facet_wrap(
      vars(!!!syms(facet_vars)),
      scales   = "free_y",
      labeller = labeller(.multi_line = FALSE),
      nrow     = facet_nrow,
      ncol     = facet_ncol
    )
  }
  
  ggplot(data    = d,
         mapping = aes(
           x     = invT_1000,
           y     = ln_k,
           group = series_id
         )) +
    
    geom_point(
      aes(
        color = .data[[colour_by]],
        shape = point_type
      ),
      size = 2
    ) +
    
    geom_line(
      data = pred,
      aes(
        y        = ln_k_fit,
        color    = .data[[colour_by]],
        linetype = peak
      ),
      linewidth   = 0.7,
      show.legend = FALSE
    ) +
    
    scale_shape_manual(
      values = c(
        "Peak 1"          = 16,
        "Peak 1 excluded" = 1,
        "Peak 2"          = 17,
        "Peak 2 excluded" = 2
      )
    ) +
    
    scale_y_continuous(labels = signs::signs_format()) +
    
    facet_layer +
    
    labs(
      x        = expression(1000/italic(T)~(K^-1)),
      y        = expression(ln(italic(k))),
      color    = colour_by,
      shape    = NULL,
      linetype = NULL
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      strip.placement  = "outside"
    )
}

# ---------- Heatmaps ---------------------------------------------------------

prep_heatmap_data <- function(
    data,
    row_var,
    col_vars,
    value_var,
    facet_vars    = NULL,
    filter_list   = list(),
    agg_fun       = median,
    na_rm         = TRUE,
    drop_na_value = FALSE,
    col_sep       = " | "
) {
  
  x <- data
  
  # ------------------------------------------------------------
  # Apply user-defined filters
  # ------------------------------------------------------------
  
  if (length(filter_list) > 0 && !is.null(names(filter_list))) {
    for (nm in names(filter_list)) {
      if (nm != "" && nm %in% names(x)) {
        val <- filter_list[[nm]]
        x <- x |> filter(.data[[nm]] %in% !!val)
      }
    }
  }
  
  all_rows <- unique(x[[row_var]])
  
  # ------------------------------------------------------------
  # Optionally drop missing values in the target variable
  # ------------------------------------------------------------
  
  if (drop_na_value) {
    x <- x |> filter(!is.na(.data[[value_var]]))
  }
  
  # ------------------------------------------------------------
  # Construct composite column label
  # ------------------------------------------------------------
  
  x <- x |>
    mutate(
      .row = .data[[row_var]],
      .col = interaction(
        !!!syms(col_vars),
        sep       = col_sep,
        drop      = TRUE,
        lex.order = TRUE
      ) |> as.character()
    )
  
  group_vars <- c(".row", ".col", facet_vars)
  
  # ------------------------------------------------------------
  # Aggregate within each heatmap cell
  # ------------------------------------------------------------
  
  safe_agg <- function(v) {
    if (all(is.na(v))) return(NA_real_)
    suppressWarnings(agg_fun(v, na.rm = na_rm))
  }
  
  x_sum <- x |>
    group_by(across(all_of(group_vars))) |>
    summarise(
      value   = safe_agg(.data[[value_var]]),
      n_cell  = n(),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------
  # Complete grid so missing/all-NA compounds still get grid cells
  # ------------------------------------------------------------
  if (!is.null(facet_vars) && length(facet_vars) > 0) {
    x_sum <- x_sum |>
      complete(
        .row = all_rows[!is.na(all_rows)],
        nesting(!!!syms(c(".col", facet_vars)))
      )
  } else {
    x_sum <- x_sum |>
      complete(
        .row = all_rows[!is.na(all_rows)],
        .col
      )
  }
  
  x_sum
}

plot_heatmap <- function(
    data,
    fill_label  = "Value",
    facet_vars  = NULL,
    x_label     = NULL,
    y_label     = NULL,
    title       = NULL,
    text_digits = NULL,
    text_size   = 3,
    rotate_x    = TRUE,
    low         = "#2166AC",
    mid         = "white",
    high        = "#B2182B",
    midpoint    = 0,
    limits      = NULL,
    na_fill     = "grey90"
) {
  
  # Build dynamic format string
  d_count <- text_digits %||% 2
  fmt_str <- sprintf("%%0.%df", d_count)
  
  # Formatting function with mathematical minus sign (U+2212)
  fmt_val <- function(val) {
    
    val_round <- round(val, d_count)
    
    # Avoid negative zero after rounding
    val_round[val_round == 0] <- 0
    
    out <- sprintf(fmt_str, abs(val_round))
    
    neg <- !is.na(val_round) & val_round < 0
    out[neg] <- paste0("\u2212", out[neg])
    
    out[is.na(val_round)] <- ""
    
    out
  }
  
  p <- ggplot(data    = data,
              mapping = aes(
                x    = .col,
                y    = .row,
                fill = value
              )) +
    geom_tile(
      color     = "white",
      linewidth = 0.3
    ) +
    scale_fill_gradient2(
      low      = low,
      mid      = mid,
      high     = high,
      midpoint = midpoint,
      limits   = limits,
      na.value = na_fill,
      name     = fill_label,
      labels   = fmt_val
    ) +
    labs(
      x = x_label %||% "Condition",
      y = y_label %||% "Compound",
      title = title
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(hjust = 0.5),
      panel.grid       = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      axis.text.x      = element_text(
        angle = if (rotate_x) 45 else 0,
        hjust = 1,
        vjust = 1
      )
    )
  
  if (!is.null(text_digits)) {
    p <- p +
      geom_text(
        aes(label = fmt_val(value)),
        size = text_size
      )
  }
  
  if (!is.null(facet_vars) && length(facet_vars) > 0) {
    p <- p +
      facet_wrap(
        vars(!!!syms(facet_vars)),
        scales   = "free",
        labeller = label_both
      )
  }
  
  p
}

make_heatmap <- function(
    data,
    row_var,
    col_vars,
    value_var,
    facet_vars    = NULL,
    filter_list   = list(),
    agg_fun       = median,
    na_rm         = TRUE,
    drop_na_value = FALSE,
    col_sep       = " | ",
    fill_label    = NULL,
    x_label       = NULL,
    y_label       = NULL,
    title         = NULL,
    text_digits   = NULL,
    text_size     = 3,
    rotate_x      = TRUE,
    low           = "#2166AC",
    mid           = "white",
    high          = "#B2182B",
    midpoint      = 0,
    limits        = NULL,
    na_fill       = "grey90"
) {
  
  hm <- prep_heatmap_data(
    data          = data,
    row_var       = row_var,
    col_vars      = col_vars,
    value_var     = value_var,
    facet_vars    = facet_vars,
    filter_list   = filter_list,
    agg_fun       = agg_fun,
    na_rm         = na_rm,
    drop_na_value = drop_na_value,
    col_sep       = col_sep
  )
  
  plot_heatmap(
    data        = hm,
    fill_label  = fill_label %||% value_var,
    facet_vars  = facet_vars,
    x_label     = x_label,
    y_label     = y_label,
    title       = title,
    text_digits = text_digits,
    text_size   = text_size,
    rotate_x    = rotate_x,
    low         = low,
    mid         = mid,
    high        = high,
    midpoint    = midpoint,
    limits      = limits,
    na_fill     = na_fill
  )
}

order_heatmap_levels <- function(
    data,
    by         = "value",
    decreasing = TRUE
  ) {
  
  row_levels <- data |>
    group_by(.row) |>
    summarise(score = mean(.data[[by]], na.rm = TRUE), .groups = "drop") |>
    arrange(if (decreasing) desc(score) else score) |>
    pull(.row)
  
  col_levels <- data |>
    group_by(.col) |>
    summarise(score = mean(.data[[by]], na.rm = TRUE), .groups = "drop") |>
    arrange(if (decreasing) desc(score) else score) |>
    pull(.col)
  
  data |>
    mutate(
      .row = factor(.row, levels = row_levels),
      .col = factor(.col, levels = col_levels)
    )
}

# ---------- Temperature dependence -------------------------------------------

vh_k_temp_effect <- function(
    vh_results,
    ref_temp_C = 25,
    delta_T    = 1
) {
  
  T1 <- ref_temp_C + 273.15
  T2 <- T1 + delta_T
  
  vh_results |>
    mutate(
      
      # fitted k at reference temperature
      k_ref = exp(intercept + slope * 1000 / T1),
      
      # fitted k after +delta_T
      k_plus = exp(intercept + slope * 1000 / T2),
      
      # absolute change
      delta_k_per_K = (k_plus - k_ref) / delta_T,
      
      # relative change
      pct_k_per_K = 100 * ((k_plus / k_ref) - 1) / delta_T)
}

vh_selectivity_temp_effect <- function(
    vh_results,
    ref_temp_C = 25,
    delta_T    = 1
) {
  
  T1 <- ref_temp_C + 273.15
  T2 <- T1 + delta_T
  
  id_cols <- c(
    "compoundname",
    "columnID",
    "eluentID",
    "modifierID",
    "flow_mL_min"
  )
  
  x <- vh_results |>
    filter(peak %in% c("Peak 1", "Peak 2")) |>
    mutate(
      peak = recode(
        peak,
        "Peak 1" = "p1",
        "Peak 2" = "p2"
      )
    ) |>
    select(
      all_of(id_cols),
      peak,
      intercept,
      slope
    ) |>
    pivot_wider(
      names_from  = peak,
      values_from = c(intercept, slope)
    )
  
  x |>
    mutate(
      
      # k values at reference T
      k1_ref = exp(intercept_p1 + slope_p1 * 1000 / T1),
      
      k2_ref = exp(intercept_p2 + slope_p2 * 1000 / T1),
      
      # k values at T + delta_T
      k1_plus = exp(intercept_p1 + slope_p1 * 1000 / T2),
      
      k2_plus = exp(intercept_p2 + slope_p2 * 1000 / T2),
      
      # selectivity
      # Identity-specific ratio; can cross 1
      ratio21_ref  = k2_ref  / k1_ref,
      ratio21_plus = k2_plus / k1_plus,
      
      # Conventional order-independent selectivity; always >= 1
      alpha_ref  = pmax(k1_ref, k2_ref)   / pmin(k1_ref, k2_ref),
      alpha_plus = pmax(k1_plus, k2_plus) / pmin(k1_plus, k2_plus),
      
      # absolute change in alpha per K
      delta_alpha_per_K = (alpha_plus - alpha_ref) / delta_T,
      
      # relative % change in alpha per K
      pct_alpha_per_K = 100 * ((alpha_plus / alpha_ref) - 1) / delta_T)
}
