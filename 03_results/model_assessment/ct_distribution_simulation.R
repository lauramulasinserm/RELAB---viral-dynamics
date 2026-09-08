# ============================================================
# Prepare lightweight RELAB data for Ct distribution plots
# ============================================================
#
# Purpose
# -------
# This script is intended to be run on the computing cluster.
# It combines observed Ct data with simulated viral-kinetic
# trajectories and exports only the lightweight summaries needed
# to reproduce the final Ct distribution figure locally.
#
# Outputs
# -------
#   - RELAB_plot_observed.csv
#   - RELAB_plot_violin.csv
#   - RELAB_plot_sim_stats.csv
#   - RELAB_plot_N.csv
#
# Required package
# ----------------
#   - dplyr
#
# Notes
# -----
# The heavy simulation files remain on the cluster. Only observed
# values, violin-density coordinates, simulated summary statistics,
# and sample sizes are exported.
# ============================================================

rm(list = ls())

library(dplyr)

set.seed(123)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

# Adapt these paths if needed.
data_dir <- "/home/laura.mulas/Monolix/data_monolix"
results_dir <- "/home/laura.mulas/Monolix/SH/results"

# ------------------------------------------------------------
# 2. Global settings
# ------------------------------------------------------------

tss_levels <- c("[0, 1]", "[2, 5]", "[6, +∞[")
ct_min <- 10
ct_max <- 40
violin_width <- 0.25

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

# Add time-since-symptom-onset categories.
add_tss_category <- function(df, time_var = "time_since_symptoms_onset") {
  df %>%
    mutate(
      TSS_category = case_when(
        .data[[time_var]] >= 0 & .data[[time_var]] <= 1 ~ "[0, 1]",
        .data[[time_var]] > 1 & .data[[time_var]] <= 5  ~ "[2, 5]",
        .data[[time_var]] > 5                          ~ "[6, +∞[",
        TRUE                                           ~ NA_character_
      ),
      TSS_category = factor(
        TSS_category,
        levels = tss_levels
      )
    )
}

# Convert Ct_diff back to Ct and restrict values to the plotting range.
prepare_observed_data <- function(data_obs, virus_name) {
  data_obs %>%
    mutate(
      time_since_symptoms_onset = time,
      Ct_obs = pmin(
        pmax(50 - Ct_diff, ct_min),
        ct_max
      )
    ) %>%
    add_tss_category() %>%
    filter(
      !is.na(TSS_category),
      !is.na(Ct_obs)
    ) %>%
    transmute(
      virus = virus_name,
      TSS_category,
      Ct = Ct_obs
    )
}

# Read one or several simulation files while retaining only trajectories
# corresponding to observed infection/time combinations.
load_trajectory_files <- function(files, data_obs) {
  if (length(files) == 0) {
    stop("No simulation trajectory files were found.")
  }

  obs_keys <- data_obs %>%
    select(
      ID_infection,
      time_since_symptoms_onset
    ) %>%
    distinct()

  pieces <- vector("list", length(files))

  for (i in seq_along(files)) {
    message(
      "  Reading simulation file ",
      i,
      "/",
      length(files),
      ": ",
      basename(files[i])
    )

    pieces[[i]] <- readRDS(files[i]) %>%
      select(any_of(c(
        "ID_infection",
        "time_since_symptoms_onset",
        "Ct",
        "sim"
      ))) %>%
      semi_join(
        obs_keys,
        by = c(
          "ID_infection",
          "time_since_symptoms_onset"
        )
      )

    gc()
  }

  bind_rows(pieces)
}

# Keep simulated Ct values within the plotting range.
filter_simulated_ct <- function(df_sim) {
  df_sim %>%
    filter(
      !is.na(TSS_category),
      !is.na(Ct),
      Ct >= ct_min,
      Ct <= ct_max
    )
}

# Calculate simulated median and interquartile range by TSS category.
summarise_simulated_ct <- function(df_sim, virus_name) {
  df_sim %>%
    filter_simulated_ct() %>%
    group_by(TSS_category) %>%
    summarise(
      med_sim = median(Ct, na.rm = TRUE),
      q25_sim = quantile(Ct, 0.25, na.rm = TRUE),
      q75_sim = quantile(Ct, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      virus = virus_name,
      .before = 1
    )
}

# Create the polygon coordinates used to reconstruct simulated violins.
# This avoids exporting millions of simulated Ct values.
make_violin_density <- function(
    df_sim,
    virus_name,
    width = violin_width
) {
  df_sim <- filter_simulated_ct(df_sim)

  violin_data <- lapply(seq_along(tss_levels), function(i) {
    tss_level <- tss_levels[i]
    ct_values <- df_sim$Ct[df_sim$TSS_category == tss_level]

    if (length(ct_values) < 2) {
      return(NULL)
    }

    dens <- density(
      ct_values,
      from = ct_min,
      to = ct_max,
      n = 512,
      na.rm = TRUE
    )

    dens_scaled <- dens$y / max(dens$y, na.rm = TRUE) * width

    bind_rows(
      data.frame(
        virus = virus_name,
        TSS_category = tss_level,
        side = "left",
        x = i - dens_scaled,
        Ct = dens$x
      ),
      data.frame(
        virus = virus_name,
        TSS_category = tss_level,
        side = "right",
        x = rev(i + dens_scaled),
        Ct = rev(dens$x)
      )
    )
  })

  bind_rows(violin_data) %>%
    mutate(
      TSS_category = factor(
        TSS_category,
        levels = tss_levels
      )
    )
}

# Process one virus from observed data through lightweight plot outputs.
process_virus <- function(
    virus_name,
    observed_file,
    simulation_pattern = NULL,
    simulation_file = NULL
) {
  message("Processing ", virus_name, "...")

  # Observed data
  data_obs <- read.csv(
    file.path(data_dir, observed_file),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      time_since_symptoms_onset = time
    )

  observed <- prepare_observed_data(
    data_obs = data_obs,
    virus_name = virus_name
  )

  n_obs <- observed %>%
    count(
      virus,
      TSS_category,
      name = "N"
    )

  # Simulation files
  if (!is.null(simulation_file)) {
    simulation_files <- file.path(
      results_dir,
      simulation_file
    )
  } else {
    simulation_files <- list.files(
      path = results_dir,
      pattern = simulation_pattern,
      full.names = TRUE
    )
  }

  if (length(simulation_files) == 0 || !all(file.exists(simulation_files))) {
    stop("Simulation file(s) not found for ", virus_name, ".")
  }

  simulated <- load_trajectory_files(
    files = simulation_files,
    data_obs = data_obs
  ) %>%
    add_tss_category()

  sim_stats <- summarise_simulated_ct(
    df_sim = simulated,
    virus_name = virus_name
  )

  violin <- make_violin_density(
    df_sim = simulated,
    virus_name = virus_name
  )

  rm(simulated)
  gc()

  list(
    observed = observed,
    violin = violin,
    sim_stats = sim_stats,
    n = n_obs
  )
}

# ------------------------------------------------------------
# 4. Virus-specific configuration
# ------------------------------------------------------------

virus_config <- list(
  IBV = list(
    observed_file = "data_for_monolix_grippe_B_24_26_vfinal.csv",
    simulation_file = "simulation_TV_flu_B_1000_traj_error_vfinal.rds"
  ),
  IAV = list(
    observed_file = "data_for_monolix_grippe_A_24_26_vfinal.csv",
    simulation_pattern = "^simulation_TV_flu_A_[0-9]+_[0-9]+_error_traj_vfinal\\.rds$"
  ),
  `SARS-CoV-2` = list(
    observed_file = "data_for_monolix_covid_24_26_vfinal.csv",
    simulation_pattern = "^simulation_TV_covid_[0-9]+_[0-9]+_error_traj_vfinal\\.rds$"
  ),
  RSV = list(
    observed_file = "data_for_monolix_VRS_24_26_vfinal.csv",
    simulation_pattern = "^simulation_TV_VRS_[0-9]+_[0-9]+_error_traj_vfinal\\.rds$"
  )
)

# ------------------------------------------------------------
# 5. Process all viruses
# ------------------------------------------------------------

results <- lapply(names(virus_config), function(virus_name) {
  cfg <- virus_config[[virus_name]]

  process_virus(
    virus_name = virus_name,
    observed_file = cfg$observed_file,
    simulation_pattern = cfg$simulation_pattern,
    simulation_file = cfg$simulation_file
  )
})

names(results) <- names(virus_config)

# ------------------------------------------------------------
# 6. Combine outputs
# ------------------------------------------------------------

plot_observed <- bind_rows(lapply(results, `[[`, "observed"))
plot_violin <- bind_rows(lapply(results, `[[`, "violin"))
plot_sim_stats <- bind_rows(lapply(results, `[[`, "sim_stats"))
plot_n <- bind_rows(lapply(results, `[[`, "n"))

# ------------------------------------------------------------
# 7. Export lightweight files
# ------------------------------------------------------------

output_files <- c(
  observed = "RELAB_plot_observed.csv",
  violin = "RELAB_plot_violin.csv",
  sim_stats = "RELAB_plot_sim_stats.csv",
  n = "RELAB_plot_N.csv"
)

write.csv(
  plot_observed,
  file.path(results_dir, output_files["observed"]),
  row.names = FALSE
)

write.csv(
  plot_violin,
  file.path(results_dir, output_files["violin"]),
  row.names = FALSE
)

write.csv(
  plot_sim_stats,
  file.path(results_dir, output_files["sim_stats"]),
  row.names = FALSE
)

write.csv(
  plot_n,
  file.path(results_dir, output_files["n"]),
  row.names = FALSE
)

message("")
message("============================================================")
message("DONE")
message("Lightweight plotting files written to:")
for (file_name in unname(output_files)) {
  message("  ", file.path(results_dir, file_name))
}
message("============================================================")
