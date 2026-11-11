## ---------------------------------------------------------------------------
## classic_vantHoff_analyzer.R -- Folder-level HPLC thermodynamic analysis.
## ---------------------------------------------------------------------------

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
      if (!all(c(condition_cols, "t0") %in% names(t0_external))) {
        stop("`t0_external` must contain 't0' and all specified `condition_cols`: ", 
             paste(setdiff(c(condition_cols, "t0"), names(t0_external)), collapse = ", "))
      }
      join_by_cols <- intersect(condition_cols, names(t0_external))
      x <- x |> 
        left_join(
          t0_external |> select(all_of(c(join_by_cols, "t0"))), 
          by = join_by_cols
        )
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
      if (!all(c(condition_cols, "phi") %in% names(phase_ratio))) {
        stop("`phase_ratio` must contain 'phi' and all specified `condition_cols`: ", 
             paste(setdiff(c(condition_cols, "phi"), names(phase_ratio)), collapse = ", "))
      }
      join_by_cols <- intersect(condition_cols, names(phase_ratio))
      x <- x |> 
        left_join(
          phase_ratio |> 
            select(all_of(c(join_by_cols, "phi"))) |> 
            rename(phase_ratio = phi), 
          by = join_by_cols
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
      
      ln_k    = mean(ln_k, na.rm = TRUE),
      sd_ln_k = if (n() > 1) {
        sd(ln_k, na.rm = TRUE)
      } else {
        NA_real_
      },
      
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
    min_temps  = 3,
    min_span_C = 10,
    r2_warn    = 0.95
) {
  
  Rgas <- 8.314462618# J / (mol * K)
  
  n_temp    <- n_distinct(d$temp_C)
  temp_span <- diff(range(d$temp_C))
  
  if (n_temp < min_temps) {
    return(
      list(
        model = NULL,
        stats = tibble(
          n_temp             = n_temp,
          n_runs             = sum(d$n_runs),
          temp_min_C         = min(d$temp_C),
          temp_max_C         = max(d$temp_C),
          temp_span_C        = temp_span,
          intercept          = NA_real_,
          intercept_se       = NA_real_,
          slope              = NA_real_,
          slope_se           = NA_real_,
          deltaH_kJ_mol      = NA_real_,
          deltaH_se_kJ_mol   = NA_real_,
          deltaS_J_mol_K     = NA_real_,
          deltaS_se_J_mol_K  = NA_real_,
          deltaG_kJ_mol      = NA_real_,
          deltaG_se_kJ_mol   = NA_real_,
          r2                 = NA_real_,
          adj_r2             = NA_real_,
          rmse_ln_k          = NA_real_,
          status             = "too_few_temperatures"
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
      
      deltaH_kJ_mol    = deltaH,
      deltaH_se_kJ_mol = deltaH_se,
      
      deltaS_J_mol_K    = deltaS,
      deltaS_se_J_mol_K = deltaS_se,
      
      deltaG_kJ_mol     = deltaG,
      deltaG_se_kJ_mol  = deltaG_se,
      
      r2        = sm$r.squared,
      adj_r2    = sm$adj.r.squared,
      rmse_ln_k = sigma(mod),
      
      status = status
    )
  )
}

vh_flag_outliers_one <- function(
    d,
    min_temps       = 3,
    rstudent_cutoff = 3,
    max_outliers    = 1
) {
  
  d <- d |>
    mutate(
      .vh_outlier           = FALSE,
      .vh_outlier_iteration = NA_integer_,
      .vh_rstudent          = NA_real_,
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
    
    score                    <- abs(rs)
    score[!is.finite(score)] <- -Inf
    
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
    
    return(
      list(
        data        = data2,
        diagnostics = tibble()
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
    crossing_min_temps        = 4L
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
        max_outliers    = max_outliers
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
        min_temps  = min_temps,
        min_span_C = min_span_C,
        r2_warn    = r2_warn
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
        
        if (isTRUE(crossing_detected)) {
          fl <- c(fl, "peak_order_crossing_corrected")
        }
        
        if (
          is.finite(peak2_imputed_n) &&
          peak2_imputed_n > 0
        ) {
          fl <- c(fl, "peak2_coelution_imputed")
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
    condition_cols            = condition_cols
  )
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