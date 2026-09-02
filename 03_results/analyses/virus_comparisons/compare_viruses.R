################################################################################
# Compare viral kinetic times across respiratory viruses
#
# Estimates and compares:
#   - time from first detection to viral peak
#   - time from viral peak to viral clearance
#   - total detectable shedding duration
#
# Uncertainty is propagated across Monte Carlo simulation draws.
# Kinetic outcomes are computed within each draw before summarising across draws.
################################################################################

rm(list = ls())


# ==============================================================================
# 1. Packages
# ==============================================================================

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)


# ==============================================================================
# 2. Settings
# ==============================================================================

Ct_LOD <- 40

simulation_dir <- "/home/laura.mulas/Monolix/SH/results"
output_dir <- "/home/laura.mulas/Monolix/SH/results"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

virus_levels <- c(
  "SARS-CoV-2",
  "IAV",
  "IBV",
  "RSV"
)


# ==============================================================================
# 3. Simulation files
# ==============================================================================

simulation_files <- list(
  "SARS-CoV-2" = file.path(
    simulation_dir,
    paste0(
      "simulation_TV_covid_",
      c("1_250", "251_500", "501_750", "751_1000"),
      "_IC_traj_vfinal.rds"
    )
  ),

  "IAV" = file.path(
    simulation_dir,
    paste0(
      "simulation_TV_flu_A_",
      c("1_250", "251_500", "501_750", "751_1000"),
      "_IC_traj_vfinal.rds"
    )
  ),

  "IBV" = file.path(
    simulation_dir,
    "simulation_TV_flu_B_1000_traj_IC_vfinal.rds"
  ),

  "RSV" = file.path(
    simulation_dir,
    "simulation_TV_VRS_1000_IC_traj_vfinal.rds"
  )
)


# ==============================================================================
# 4. Helper functions
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Linear interpolation of a Ct threshold crossing
# ------------------------------------------------------------------------------

interpolate_crossing <- function(
    t1,
    t2,
    ct1,
    ct2,
    threshold = Ct_LOD) {

  if (
    !is.finite(t1) ||
      !is.finite(t2) ||
      !is.finite(ct1) ||
      !is.finite(ct2)
  ) {
    return(NA_real_)
  }

  if (ct1 == ct2) {
    return(mean(c(t1, t2)))
  }

  t1 +
    (threshold - ct1) *
    (t2 - t1) /
    (ct2 - ct1)
}


# ------------------------------------------------------------------------------
# 4.2 Extract kinetic times from one simulated trajectory
# ------------------------------------------------------------------------------

extract_kinetic_times <- function(
    df_trajectory,
    threshold = Ct_LOD) {

  df_trajectory <- df_trajectory %>%
    transmute(
      time = time_since_symptoms_onset,
      Ct = Ct_value
    ) %>%
    filter(
      is.finite(time),
      is.finite(Ct)
    ) %>%
    arrange(time)

  if (nrow(df_trajectory) < 3) {
    return(
      tibble(
        t_detection = NA_real_,
        t_peak = NA_real_,
        t_clearance = NA_real_,
        detection_to_peak = NA_real_,
        peak_to_clearance = NA_real_,
        detectable_duration = NA_real_
      )
    )
  }

  # Peak viral load corresponds to the minimum Ct.
  idx_peak <- which.min(df_trajectory$Ct)[1]
  t_peak <- df_trajectory$time[idx_peak]

  # First detection before peak: Ct >= threshold -> Ct < threshold.
  if (idx_peak > 1) {
    idx_detection <- which(
      head(df_trajectory$Ct[seq_len(idx_peak)], -1) >= threshold &
        tail(df_trajectory$Ct[seq_len(idx_peak)], -1) < threshold
    )
  } else {
    idx_detection <- integer(0)
  }

  if (length(idx_detection) > 0) {
    i_det <- tail(idx_detection, 1)

    t_detection <- interpolate_crossing(
      t1 = df_trajectory$time[i_det],
      t2 = df_trajectory$time[i_det + 1],
      ct1 = df_trajectory$Ct[i_det],
      ct2 = df_trajectory$Ct[i_det + 1],
      threshold = threshold
    )
  } else {
    t_detection <- NA_real_
  }

  # Viral clearance after peak: Ct < threshold -> Ct >= threshold.
  if (idx_peak < nrow(df_trajectory)) {
    post_peak_index <- idx_peak:nrow(df_trajectory)

    idx_clearance_local <- which(
      head(df_trajectory$Ct[post_peak_index], -1) < threshold &
        tail(df_trajectory$Ct[post_peak_index], -1) >= threshold
    )
  } else {
    idx_clearance_local <- integer(0)
  }

  if (length(idx_clearance_local) > 0) {
    i_clear <- idx_peak + idx_clearance_local[1] - 1

    t_clearance <- interpolate_crossing(
      t1 = df_trajectory$time[i_clear],
      t2 = df_trajectory$time[i_clear + 1],
      ct1 = df_trajectory$Ct[i_clear],
      ct2 = df_trajectory$Ct[i_clear + 1],
      threshold = threshold
    )
  } else {
    t_clearance <- NA_real_
  }

  tibble(
    t_detection = t_detection,
    t_peak = t_peak,
    t_clearance = t_clearance,
    detection_to_peak = ifelse(
      is.finite(t_detection),
      t_peak - t_detection,
      NA_real_
    ),
    peak_to_clearance = ifelse(
      is.finite(t_clearance),
      t_clearance - t_peak,
      NA_real_
    ),
    detectable_duration = ifelse(
      is.finite(t_detection) & is.finite(t_clearance),
      t_clearance - t_detection,
      NA_real_
    )
  )
}


# ------------------------------------------------------------------------------
# 4.3 Reduce one simulation file
# ------------------------------------------------------------------------------

reduce_simulation_file <- function(file) {

  message("Reading: ", basename(file))

  df <- readRDS(file)

  required_columns <- c(
    "sim",
    "time_since_symptoms_onset",
    "Ct"
  )

  missing_columns <- setdiff(
    required_columns,
    names(df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing column(s) in ",
      basename(file),
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  # Population-median Ct at each time point within each Monte Carlo draw.
  df_small <- df %>%
    group_by(
      sim,
      time_since_symptoms_onset
    ) %>%
    summarise(
      Ct_value = median(Ct, na.rm = TRUE),
      .groups = "drop"
    )

  rm(df)
  gc()

  df_small
}


# ------------------------------------------------------------------------------
# 4.4 Process one virus
# ------------------------------------------------------------------------------

calculate_kinetics_virus <- function(
    files,
    virus_name,
    threshold = Ct_LOD) {

  message("\nProcessing ", virus_name)

  reduced_list <- vector(
    "list",
    length(files)
  )

  for (i in seq_along(files)) {

    if (!file.exists(files[i])) {
      stop(
        "Simulation file not found: ",
        files[i]
      )
    }

    reduced_list[[i]] <- reduce_simulation_file(
      files[i]
    )
  }

  df_reduced <- bind_rows(reduced_list)

  rm(reduced_list)
  gc()

  result <- df_reduced %>%
    group_by(sim) %>%
    group_modify(
      ~ extract_kinetic_times(
        .x,
        threshold = threshold
      )
    ) %>%
    ungroup() %>%
    mutate(
      virus = virus_name,
      .before = 1
    )

  rm(df_reduced)
  gc()

  result
}


# ------------------------------------------------------------------------------
# 4.5 Summarise uncertainty
# ------------------------------------------------------------------------------

summarise_metric <- function(x) {

  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(
      tibble(
        median = NA_real_,
        low95 = NA_real_,
        high95 = NA_real_,
        n_valid = 0,
        n_unique = 0
      )
    )
  }

  tibble(
    median = median(x),
    low95 = unname(
      quantile(
        x,
        probs = 0.025,
        type = 8
      )
    ),
    high95 = unname(
      quantile(
        x,
        probs = 0.975,
        type = 8
      )
    ),
    n_valid = length(x),
    n_unique = n_distinct(x)
  )
}


# ==============================================================================
# 5. Calculate kinetic times for all viruses
# ==============================================================================

kinetics_by_draw <- map_dfr(
  names(simulation_files),
  function(virus_name) {
    calculate_kinetics_virus(
      files = simulation_files[[virus_name]],
      virus_name = virus_name,
      threshold = Ct_LOD
    )
  }
) %>%
  mutate(
    virus = factor(
      virus,
      levels = virus_levels
    )
  )


# ==============================================================================
# 6. Summarise uncertainty across Monte Carlo draws
# ==============================================================================

kinetics_long <- kinetics_by_draw %>%
  select(
    virus,
    sim,
    detection_to_peak,
    peak_to_clearance,
    detectable_duration
  ) %>%
  pivot_longer(
    cols = c(
      detection_to_peak,
      peak_to_clearance,
      detectable_duration
    ),
    names_to = "metric",
    values_to = "value"
  )

kinetics_summary <- kinetics_long %>%
  group_by(
    virus,
    metric
  ) %>%
  group_modify(
    ~ summarise_metric(.x$value)
  ) %>%
  ungroup() %>%
  mutate(
    metric = recode(
      metric,
      "detection_to_peak" = "Detection to peak",
      "peak_to_clearance" = "Peak to clearance",
      "detectable_duration" = "Detectable shedding duration"
    )
  )

print(kinetics_summary)


# ==============================================================================
# 7. Diagnostic check
# ==============================================================================

# If n_unique = 1, identical median/low/high values are expected because
# the Monte Carlo draws contain no variability for that metric.

diagnostic_table <- kinetics_summary %>%
  select(
    virus,
    metric,
    median,
    low95,
    high95,
    n_valid,
    n_unique
  )

print(
  diagnostic_table,
  n = Inf
)


# ==============================================================================
# 8. Plot: time from viral peak to clearance
# ==============================================================================

clearance_summary <- kinetics_summary %>%
  filter(
    metric == "Peak to clearance"
  )

clearance_plot <- ggplot(
  clearance_summary,
  aes(
    x = virus,
    y = median
  )
) +
  geom_errorbar(
    aes(
      ymin = low95,
      ymax = high95
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    x = NULL,
    y = "Time from viral peak to clearance (days)"
  ) +
  theme_classic(
    base_size = 14
  )

clearance_plot


# ==============================================================================
# 9. Plot: all kinetic times
# ==============================================================================

kinetics_plot <- ggplot(
  kinetics_summary,
  aes(
    x = virus,
    y = median
  )
) +
  geom_errorbar(
    aes(
      ymin = low95,
      ymax = high95
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "Time (days)"
  ) +
  theme_classic(
    base_size = 14
  )

kinetics_plot


# ==============================================================================
# 10. Save results
# ==============================================================================

write.csv(
  kinetics_by_draw,
  file.path(
    output_dir,
    "virus_kinetics_by_draw.csv"
  ),
  row.names = FALSE
)

write.csv(
  kinetics_summary,
  file.path(
    output_dir,
    "virus_kinetics_summary.csv"
  ),
  row.names = FALSE
)

ggsave(
  filename = file.path(
    output_dir,
    "time_to_clearance_comparison.svg"
  ),
  plot = clearance_plot,
  width = 7,
  height = 5
)

ggsave(
  filename = file.path(
    output_dir,
    "virus_kinetic_times_comparison.svg"
  ),
  plot = kinetics_plot,
  width = 10,
  height = 5
)

message("Virus comparison analysis completed.")
