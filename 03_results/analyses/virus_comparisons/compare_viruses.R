################################################################################################
# Comparison of time to viral peak across respiratory viruses
# Author: Laura MULAS
# Purpose:
#   Estimate and compare the time from first viral detection (Ct < 40)
#   to peak viral load across SARS-CoV-2, IAV, IBV, and RSV.
#
#   Uncertainty is propagated from the Monte Carlo trajectory simulations.
################################################################################################

rm(list = ls())

# ---------------------------------------------------------------------------
# 0. Libraries
# ---------------------------------------------------------------------------

library(dplyr)
library(ggplot2)


# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

Ct_LOD <- 40

simulation_dir <- "03_results_analysis/results/simulations"

output_dir <- "03_results_analysis/results/time_to_peak"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ---------------------------------------------------------------------------
# 2. Simulation files
# ---------------------------------------------------------------------------

simulation_files <- list(

  "SARS-CoV-2" = file.path(
    simulation_dir,
    "COVID",
    paste0(
      "COVID_trajectories_",
      c("1_250", "251_500", "501_750", "751_1000"),
      ".rds"
    )
  ),

  "IAV" = file.path(
    simulation_dir,
    "IAV",
    paste0(
      "IAV_trajectories_",
      c("1_250", "251_500", "501_750", "751_1000"),
      ".rds"
    )
  ),

  "IBV" = file.path(
    simulation_dir,
    "IBV",
    paste0(
      "IBV_trajectories_",
      c("1_250", "251_500", "501_750", "751_1000"),
      ".rds"
    )
  ),

  "RSV" = file.path(
    simulation_dir,
    "RSV",
    paste0(
      "RSV_trajectories_",
      c("1_250", "251_500", "501_750", "751_1000"),
      ".rds"
    )
  )
)


# ---------------------------------------------------------------------------
# 3. Helper functions
# ---------------------------------------------------------------------------

# Linear interpolation of the time at which Ct crosses the detection threshold.
interpolate_crossing <- function(
    t1,
    t2,
    ct1,
    ct2,
    threshold = 40) {

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


# Extract detection-to-peak time from one simulated trajectory.
extract_time_to_peak <- function(
    df_trajectory,
    threshold = 40) {

  df_trajectory <- df_trajectory %>%
    select(
      time_since_symptoms_onset,
      Ct_value
    ) %>%
    filter(
      is.finite(time_since_symptoms_onset),
      is.finite(Ct_value)
    ) %>%
    arrange(
      time_since_symptoms_onset
    )

  if (nrow(df_trajectory) < 2) {
    return(NA_real_)
  }


  # Viral peak corresponds to the minimum Ct.
  idx_peak <- which.min(
    df_trajectory$Ct_value
  )

  peak_time <- df_trajectory$
    time_since_symptoms_onset[idx_peak]


  # Restrict to the ascending phase before the viral peak.
  df_before_peak <- df_trajectory[
    seq_len(idx_peak),
  ]


  if (nrow(df_before_peak) < 2) {
    return(NA_real_)
  }


  # Identify threshold crossing from undetectable to detectable:
  # Ct >= 40 -> Ct < 40.
  idx_detection <- which(
    head(
      df_before_peak$Ct_value,
      -1
    ) >= threshold &
      tail(
        df_before_peak$Ct_value,
        -1
      ) < threshold
  )


  if (length(idx_detection) == 0) {
    return(NA_real_)
  }


  # Use the last crossing before the peak.
  i <- tail(
    idx_detection,
    1
  )

  detection_time <- interpolate_crossing(
    t1 = df_before_peak$
      time_since_symptoms_onset[i],

    t2 = df_before_peak$
      time_since_symptoms_onset[i + 1],

    ct1 = df_before_peak$
      Ct_value[i],

    ct2 = df_before_peak$
      Ct_value[i + 1],

    threshold = threshold
  )


  if (
    !is.finite(detection_time) ||
    !is.finite(peak_time)
  ) {
    return(NA_real_)
  }


  peak_time - detection_time
}


# Reduce a large trajectory file immediately to:
# simulation draw x time x median Ct.
reduce_simulation_file <- function(file) {

  message(
    "Reading: ",
    basename(file)
  )

  df <- readRDS(file)

  df_small <- df %>%
    group_by(
      sim,
      time_since_symptoms_onset
    ) %>%
    summarise(
      Ct_value = median(
        Ct,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  rm(df)
  gc()

  df_small
}


# Process all trajectory files from one virus.
calculate_time_to_peak_virus <- function(
    files,
    virus_name,
    threshold = 40) {

  message(
    "\nProcessing ",
    virus_name
  )


  # Read and reduce each large file sequentially.
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

    reduced_list[[i]] <-
      reduce_simulation_file(
        files[i]
      )
  }


  df_reduced <- bind_rows(
    reduced_list
  )

  rm(reduced_list)
  gc()


  # Estimate time to peak separately for each Monte Carlo draw.
  sim_ids <- sort(
    unique(df_reduced$sim)
  )

  time_to_peak <- vapply(
    sim_ids,
    function(sim_id) {

      extract_time_to_peak(
        df_reduced %>%
          filter(
            sim == sim_id
          ),
        threshold = threshold
      )
    },
    numeric(1)
  )


  result <- tibble(
    sim = sim_ids,
    virus = virus_name,
    time_to_peak = time_to_peak
  )


  rm(df_reduced)
  gc()

  result
}


# ---------------------------------------------------------------------------
# 4. Calculate time to peak for each virus
# ---------------------------------------------------------------------------

time_to_peak_by_draw <- bind_rows(

  lapply(
    names(simulation_files),
    function(virus_name) {

      calculate_time_to_peak_virus(
        files = simulation_files[[virus_name]],
        virus_name = virus_name,
        threshold = Ct_LOD
      )
    }
  )

) %>%
  mutate(
    virus = factor(
      virus,
      levels = c(
        "SARS-CoV-2",
        "IAV",
        "IBV",
        "RSV"
      )
    )
  )


# ---------------------------------------------------------------------------
# 5. Summarize uncertainty
# ---------------------------------------------------------------------------

time_to_peak_summary <- time_to_peak_by_draw %>%
  group_by(
    virus
  ) %>%
  summarise(

    time_to_peak = median(
      time_to_peak,
      na.rm = TRUE
    ),

    low95 = quantile(
      time_to_peak,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),

    high95 = quantile(
      time_to_peak,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),

    n_valid = sum(
      is.finite(time_to_peak)
    ),

    .groups = "drop"
  )


print(
  time_to_peak_summary
)


# ---------------------------------------------------------------------------
# 6. Plot comparison
# ---------------------------------------------------------------------------

time_to_peak_plot <- ggplot(
  time_to_peak_summary,
  aes(
    x = virus,
    y = time_to_peak
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
    y = "Time from first detection to viral peak (days)"
  ) +
  theme_classic(
    base_size = 14
  )


time_to_peak_plot


# ---------------------------------------------------------------------------
# 7. Save results
# ---------------------------------------------------------------------------

write.csv(
  time_to_peak_by_draw,
  file.path(
    output_dir,
    "time_to_peak_by_draw.csv"
  ),
  row.names = FALSE
)

write.csv(
  time_to_peak_summary,
  file.path(
    output_dir,
    "time_to_peak_summary.csv"
  ),
  row.names = FALSE
)

ggsave(
  filename = file.path(
    output_dir,
    "time_to_peak_comparison.svg"
  ),
  plot = time_to_peak_plot,
  width = 7,
  height = 5
)


message(
  "Time-to-peak analysis completed."
)
