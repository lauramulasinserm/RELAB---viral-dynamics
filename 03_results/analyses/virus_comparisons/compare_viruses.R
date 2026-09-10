################################################################################
# Compare viral kinetic times across respiratory viruses
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

simulation_dir <- ".../results"
output_dir <- ".../results"

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

extract_clearance_time <- function(
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
        clearance_time = NA_real_
      )
    )
  }

  # Viral peak = minimum Ct
  idx_peak <- which.min(df_trajectory$Ct)[1]
  t_peak <- df_trajectory$time[idx_peak]

  # Search for clearance after the viral peak:
  # Ct < threshold -> Ct >= threshold
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

    clearance_time <- t_clearance - t_peak

  } else {

    clearance_time <- NA_real_

  }

  tibble(
    clearance_time = clearance_time
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

calculate_clearance_virus <- function(
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
      ~ extract_clearance_time(
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
    calculate_clearance_virus(
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

clearance_summary <- clearance_by_draw %>%
  group_by(virus) %>%
  group_modify(
    ~ summarise_metric(.x$clearance_time)
  ) %>%
  ungroup()

print(clearance_summary)


# ==============================================================================
# 7. Plot: time from viral peak to clearance
# ==============================================================================

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
# 8. Save results
# ==============================================================================

write.csv(
  clearance_by_draw,
  file.path(
    output_dir,
    "virus_clearance_by_draw.csv"
  ),
  row.names = FALSE
)

write.csv(
  clearance_summary,
  file.path(
    output_dir,
    "virus_clearance_summary.csv"
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

message("Virus clearance comparison completed.")
