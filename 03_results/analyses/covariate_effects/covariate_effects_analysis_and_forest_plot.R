#################################################################################################
# COVARIATE EFFECTS ON VIRAL KINETICS
#
# Purpose
#   1. Reduce large simulation files in a memory-efficient way
#   2. Estimate covariate effects on peak Ct and time to viral clearance
#   3. Save virus-specific summary tables
#   4. Build the combined forest plot across viruses
#
# Viruses
#   - SARS-CoV-2
#   - Influenza A virus (IAV)
#   - Influenza B virus (IBV)
#   - Respiratory syncytial virus (RSV)
#
# Notes for GitHub
#   - Edit only the paths in SECTION 1.
#   - Raw simulation files are not included in the repository.
#   - Large *.rds files should be excluded with .gitignore.
#################################################################################################

rm(list = ls())
set.seed(123)

# ================================================================================================
# 0) PACKAGES
# ================================================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(forcats)
library(rlang)

library(ggplot2)
library(ggtext)
library(patchwork)
library(cowplot)
library(grid)
library(svglite)

# ================================================================================================
# 1) USER SETTINGS AND PATHS
# ================================================================================================

# Detection threshold used to define viral clearance
Ct_LOD <- 40

# Directory in which summary CSV files are written.
# Change this path for your own machine/server.
results_dir <- "/home/laura.mulas/Monolix/R_simulation/results"

# Final figure path.
# Change this path for your own machine.
figure_file <- file.path(results_dir, "RELAB_cov_effects_FINAL.svg")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Simulation files are defined separately in each virus section below.
# Keeping them explicit makes it easy to see which files are used in each analysis.

# Model-derived kinetic outcomes
# Generic Monte Carlo propagation for covariate effects
###########################################################################################################################

#################################################################################################
# MEMORY-EFFICIENT ANALYSIS OF COVARIATE EFFECTS
# ALL RESPIRATORY VIRUSES
#################################################################################################


# ================================================================================================
# 1. Extract peak and clearance from one median Ct trajectory
# ================================================================================================

get_peak_and_clearance <- function(df, threshold = Ct_LOD) {
  
  df <- df %>%
    arrange(time_since_symptoms_onset)
  
  
  # ---------------------------------------------------------------------------
  # Peak
  # ---------------------------------------------------------------------------
  
  i_peak <- which.min(df$Ct)[1]
  
  t_peak  <- df$time_since_symptoms_onset[i_peak]
  Ct_peak <- df$Ct[i_peak]
  
  
  # ---------------------------------------------------------------------------
  # Internal function: Ct threshold crossing
  # ---------------------------------------------------------------------------
  
  get_t_cross <- function(x, y, threshold = Ct_LOD, t_min = -Inf) {
    
    keep <- x >= t_min
    
    x <- x[keep]
    y <- y[keep]
    
    
    if (length(x) < 2) {
      return(NA_real_)
    }
    
    
    idx <- which(
      head(y, -1) < threshold &
        tail(y, -1) >= threshold
    )[1]
    
    
    if (is.na(idx)) {
      return(NA_real_)
    }
    
    
    x1 <- x[idx]
    x2 <- x[idx + 1]
    
    y1 <- y[idx]
    y2 <- y[idx + 1]
    
    
    if (y2 == y1) {
      return(x1)
    }
    
    
    x1 +
      (threshold - y1) *
      (x2 - x1) /
      (y2 - y1)
  }
  
  
  # ---------------------------------------------------------------------------
  # Clearance
  # ---------------------------------------------------------------------------
  
  t_clearance <- get_t_cross(
    x = df$time_since_symptoms_onset,
    y = df$Ct,
    threshold = threshold,
    t_min = t_peak
  )
  
  
  time_to_clearance <- t_clearance - t_peak
  
  
  tibble(
    t_peak = t_peak,
    Ct_peak = Ct_peak,
    t_clearance = t_clearance,
    time_to_clearance = time_to_clearance
  )
}



# ================================================================================================
# 2. Read large simulation files ONE BY ONE
#    and immediately reduce them to median trajectories
#
# IMPORTANT:
# Only the selected covariate is retained.
#
# Therefore, for age:
# sim × age × time
#
# For season:
# sim × season × time
#
# etc.
# ================================================================================================

read_and_reduce_simulations <- function(
    files,
    covariate
) {
  
  covariate <- ensym(covariate)
  
  reduced_list <- vector(
    "list",
    length(files)
  )
  
  
  for (i in seq_along(files)) {
    
    cat(
      "\n--------------------------------------------------\n",
      "Processing file ",
      i,
      "/",
      length(files),
      "\n",
      basename(files[i]),
      "\n--------------------------------------------------\n",
      sep = ""
    )
    
    
    # -------------------------------------------------------------------------
    # Read ONLY one large file
    # -------------------------------------------------------------------------
    
    tmp <- readRDS(files[i])
    
    
    cat(
      "Raw rows:",
      nrow(tmp),
      "\n"
    )
    
    
    # -------------------------------------------------------------------------
    # Age categories
    # -------------------------------------------------------------------------
    
    if ("age_cat" %in% names(tmp)) {
      
      tmp <- tmp %>%
        mutate(
          age_cat = factor(
            age_cat,
            levels = c(
              "<5",
              "5-18",
              "18-65",
              ">65"
            )
          )
        )
    }
    
    
    # -------------------------------------------------------------------------
    # Immediately reduce size
    # -------------------------------------------------------------------------
    
    reduced_list[[i]] <- tmp %>%
      group_by(
        sim,
        !!covariate,
        time_since_symptoms_onset
      ) %>%
      summarise(
        Ct = median(
          Ct,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    
    cat(
      "Reduced rows:",
      nrow(reduced_list[[i]]),
      "\n"
    )
    
    
    # -------------------------------------------------------------------------
    # Delete large raw dataset
    # -------------------------------------------------------------------------
    
    rm(tmp)
    gc()
  }
  
  
  # ---------------------------------------------------------------------------
  # bind_rows is now performed ONLY on small datasets
  # ---------------------------------------------------------------------------
  
  result <- bind_rows(
    reduced_list
  )
  
  
  rm(reduced_list)
  gc()
  
  
  result
}



# ================================================================================================
# 3. Compute effects FROM ALREADY REDUCED trajectories
# ================================================================================================

compute_covariate_effects_reduced <- function(
    df_group_sim,
    covariate,
    reference,
    threshold = Ct_LOD
) {
  
  covariate <- ensym(covariate)
  
  
  # ==============================================================================================
  # 3.1 Outcomes for each simulation and subgroup
  # ==============================================================================================
  
  df_outcomes <- df_group_sim %>%
    group_by(
      sim,
      !!covariate
    ) %>%
    group_modify(
      ~ get_peak_and_clearance(
        .x,
        threshold = threshold
      )
    ) %>%
    ungroup()
  
  
  # ==============================================================================================
  # 3.2 Reference category
  # ==============================================================================================
  
  df_reference <- df_outcomes %>%
    filter(
      !!covariate == reference
    ) %>%
    transmute(
      sim,
      Ct_peak_ref = Ct_peak,
      time_to_clearance_ref = time_to_clearance
    )
  
  
  if (nrow(df_reference) == 0) {
    
    stop(
      paste0(
        "Reference category '",
        reference,
        "' was not found for covariate ",
        as_string(covariate),
        "."
      )
    )
  }
  
  
  # ==============================================================================================
  # 3.3 Difference relative to reference WITHIN simulation
  # ==============================================================================================
  
  df_differences <- df_outcomes %>%
    left_join(
      df_reference,
      by = "sim"
    ) %>%
    mutate(
      
      delta_Ct_peak =
        Ct_peak -
        Ct_peak_ref,
      
      delta_clearance =
        time_to_clearance -
        time_to_clearance_ref
    )
  
  
  # ==============================================================================================
  # 3.4 Absolute outcome summary
  # ==============================================================================================
  
  df_outcome_summary <- df_outcomes %>%
    group_by(
      !!covariate
    ) %>%
    summarise(
      
      Ct_peak_median =
        median(
          Ct_peak,
          na.rm = TRUE
        ),
      
      Ct_peak_low =
        quantile(
          Ct_peak,
          probs = 0.025,
          na.rm = TRUE,
          names = FALSE
        ),
      
      Ct_peak_high =
        quantile(
          Ct_peak,
          probs = 0.975,
          na.rm = TRUE,
          names = FALSE
        ),
      
      
      clearance_median =
        median(
          time_to_clearance,
          na.rm = TRUE
        ),
      
      clearance_low =
        quantile(
          time_to_clearance,
          probs = 0.025,
          na.rm = TRUE,
          names = FALSE
        ),
      
      clearance_high =
        quantile(
          time_to_clearance,
          probs = 0.975,
          na.rm = TRUE,
          names = FALSE
        ),
      
      
      n_valid_peak =
        sum(
          is.finite(Ct_peak)
        ),
      
      n_valid_clearance =
        sum(
          is.finite(time_to_clearance)
        ),
      
      .groups = "drop"
    )
  
  
  # ==============================================================================================
  # 3.5 Effect relative to reference
  # ==============================================================================================
  
  df_effect_summary <- df_differences %>%
    group_by(
      !!covariate
    ) %>%
    summarise(
      
      delta_Ct_peak_median =
        median(
          delta_Ct_peak,
          na.rm = TRUE
        ),
      
      delta_Ct_peak_low =
        quantile(
          delta_Ct_peak,
          probs = 0.025,
          na.rm = TRUE,
          names = FALSE
        ),
      
      delta_Ct_peak_high =
        quantile(
          delta_Ct_peak,
          probs = 0.975,
          na.rm = TRUE,
          names = FALSE
        ),
      
      
      delta_clearance_median =
        median(
          delta_clearance,
          na.rm = TRUE
        ),
      
      delta_clearance_low =
        quantile(
          delta_clearance,
          probs = 0.025,
          na.rm = TRUE,
          names = FALSE
        ),
      
      delta_clearance_high =
        quantile(
          delta_clearance,
          probs = 0.975,
          na.rm = TRUE,
          names = FALSE
        ),
      
      
      n_valid_peak =
        sum(
          is.finite(delta_Ct_peak)
        ),
      
      n_valid_clearance =
        sum(
          is.finite(delta_clearance)
        ),
      
      .groups = "drop"
    )
  
  
  # ==============================================================================================
  # 3.6 Median trajectory + uncertainty interval
  # ==============================================================================================
  
  df_trajectory_summary <- df_group_sim %>%
    group_by(
      !!covariate,
      time_since_symptoms_onset
    ) %>%
    summarise(
      
      Ct_median =
        median(
          Ct,
          na.rm = TRUE
        ),
      
      Ct_low =
        quantile(
          Ct,
          probs = 0.025,
          na.rm = TRUE,
          names = FALSE
        ),
      
      Ct_high =
        quantile(
          Ct,
          probs = 0.975,
          na.rm = TRUE,
          names = FALSE
        ),
      
      .groups = "drop"
    )
  
  
  # ==============================================================================================
  # RETURN
  # ==============================================================================================
  
  list(
    
    trajectories_by_sim =
      df_group_sim,
    
    outcomes_by_sim =
      df_outcomes,
    
    differences_by_sim =
      df_differences,
    
    outcome_summary =
      df_outcome_summary,
    
    effect_summary =
      df_effect_summary,
    
    trajectory_summary =
      df_trajectory_summary
  )
}



#################################################################################################
#################################################################################################
#
#                                   SARS-CoV-2
#
#################################################################################################
#################################################################################################


# ================================================================================================
# 4. COVID simulation files
# ================================================================================================

sim_files_cov <- c(
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_covid_1_250_IC_traj_vfinal.rds",
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_covid_251_500_IC_traj_vfinal.rds",
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_covid_501_750_IC_traj_vfinal.rds",
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_covid_751_1000_IC_traj_vfinal.rds"
)



# ================================================================================================
# 4.1 AGE
# ================================================================================================

cat("\n\nCOVID - AGE\n")


df_age_cov <- read_and_reduce_simulations(
  files = sim_files_cov,
  covariate = age_cat
)


results_age_cov <- compute_covariate_effects_reduced(
  df_group_sim = df_age_cov,
  covariate = age_cat,
  reference = "18-65",
  threshold = Ct_LOD
)


# Free memory
rm(df_age_cov)
gc()



# ================================================================================================
# 4.2 EPIDEMIC SEASON
# ================================================================================================

cat("\n\nCOVID - SEASON\n")


df_season_cov <- read_and_reduce_simulations(
  files = sim_files_cov,
  covariate = season
)


results_season_cov <- compute_covariate_effects_reduced(
  df_group_sim = df_season_cov,
  covariate = season,
  reference = "24_25",
  threshold = Ct_LOD
)


rm(df_season_cov)
gc()



# ================================================================================================
# 4.3 EPIDEMIC PERIOD
# ================================================================================================

cat("\n\nCOVID - EPIDEMIC PERIOD\n")


df_period_cov <- read_and_reduce_simulations(
  files = sim_files_cov,
  covariate = periode_epi
)


results_period_cov <- compute_covariate_effects_reduced(
  df_group_sim = df_period_cov,
  covariate = periode_epi,
  reference = "Phase 2",
  threshold = Ct_LOD
)


rm(df_period_cov)
gc()



# ================================================================================================
# 4.4 COVID RESULTS TABLE
# ================================================================================================

results_age_table_cov <- results_age_cov$effect_summary %>%
  rename(
    category = age_cat
  ) %>%
  mutate(
    virus = "SARS-CoV-2",
    covariate = "Age",
    .before = 1
  )


results_season_table_cov <- results_season_cov$effect_summary %>%
  rename(
    category = season
  ) %>%
  mutate(
    virus = "SARS-CoV-2",
    covariate = "Epidemic season",
    .before = 1
  )


results_period_table_cov <- results_period_cov$effect_summary %>%
  rename(
    category = periode_epi
  ) %>%
  mutate(
    virus = "SARS-CoV-2",
    covariate = "Epidemic period",
    .before = 1
  )


results_effects_table_cov <- bind_rows(
  results_age_table_cov,
  results_season_table_cov,
  results_period_table_cov
)

write.csv(
  results_effects_table_cov,
  file.path(results_dir, "COV_covariate_effects.csv"),
  row.names = FALSE
)

#################################################################################################
#################################################################################################
#
#                                   INFLUENZA A
#
#################################################################################################
#################################################################################################


# ================================================================================================
# IBV simulation files
# ================================================================================================

sim_files_A <- c(
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_flu_A_1_250_IC_traj_vfinal.rds",
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_flu_A_251_500_IC_traj_vfinal.rds",
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_flu_A_501_750_IC_traj_vfinal.rds",
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_flu_A_751_1000_IC_traj_vfinal.rds"
)



# ================================================================================================
# 5.1 AGE
# ================================================================================================

cat("\n\nIAV - AGE\n")


df_age_A <- read_and_reduce_simulations(
  files = sim_files_A,
  covariate = age_cat
)


results_age_A <- compute_covariate_effects_reduced(
  df_group_sim = df_age_A,
  covariate = age_cat,
  reference = "18-65",
  threshold = Ct_LOD
)


rm(df_age_A)
gc()



# ================================================================================================
# 5.2 EPIDEMIC SEASON
# ================================================================================================

cat("\n\nIAV - SEASON\n")


df_season_A <- read_and_reduce_simulations(
  files = sim_files_A,
  covariate = season
)


results_season_A <- compute_covariate_effects_reduced(
  df_group_sim = df_season_A,
  covariate = season,
  reference = "24_25",
  threshold = Ct_LOD
)


rm(df_season_A)
gc()



# ================================================================================================
# 5.3 EPIDEMIC PERIOD
# ================================================================================================

cat("\n\nIAV - EPIDEMIC PERIOD\n")


df_period_A <- read_and_reduce_simulations(
  files = sim_files_A,
  covariate = periode_epi
)


results_period_A <- compute_covariate_effects_reduced(
  df_group_sim = df_period_A,
  covariate = periode_epi,
  reference = "Phase 2",
  threshold = Ct_LOD
)


rm(df_period_A)
gc()



# ================================================================================================
# IBV RESULTS TABLE
# ================================================================================================

results_age_table_A <- results_age_A$effect_summary %>%
  rename(
    category = age_cat
  ) %>%
  mutate(
    virus = "Influenza A",
    covariate = "Age",
    .before = 1
  )


results_season_table_A <- results_season_A$effect_summary %>%
  rename(
    category = season
  ) %>%
  mutate(
    virus = "Influenza A",
    covariate = "Epidemic season",
    .before = 1
  )


results_period_table_A <- results_period_A$effect_summary %>%
  rename(
    category = periode_epi
  ) %>%
  mutate(
    virus = "Influenza A",
    covariate = "Epidemic period",
    .before = 1
  )


results_effects_table_A <- bind_rows(
  results_age_table_A,
  results_season_table_A,
  results_period_table_A
)


write.csv(
  results_effects_table_A,
  file.path(results_dir, "IAV_covariate_effects.csv"),
  row.names = FALSE
)

#################################################################################################
#################################################################################################
#
#                                   INFLUENZA B
#
#################################################################################################
#################################################################################################


# ================================================================================================
# 5. IAV simulation files
# ================================================================================================

sim_files_B <- c(
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_flu_B_1000_traj_IC_vfinal.rds"
)



# ================================================================================================
# 5.1 Vaccin
# ================================================================================================

cat("\n\nIBV - VACCIN\n")


df_vac_B <- read_and_reduce_simulations(
  files = sim_files_B,
  covariate = vaccin
)


results_vac_B <- compute_covariate_effects_reduced(
  df_group_sim = df_vac_B,
  covariate = vaccin,
  reference = 0,
  threshold = Ct_LOD
)


rm(df_vac_B)
gc()



# ================================================================================================
# 5.2 SEXE
# ================================================================================================

cat("\n\nIBV - SEX\n")


df_sexe_B <- read_and_reduce_simulations(
  files = sim_files_B,
  covariate = sexe
)


results_sexe_B <- compute_covariate_effects_reduced(
  df_group_sim = df_sexe_B,
  covariate = sexe,
  reference = "Female",
  threshold = Ct_LOD
)


rm(df_sexe_B)
gc()


# ================================================================================================
# 5.4 IBV RESULTS TABLE
# ================================================================================================

results_vac_table_B <- results_vac_B$effect_summary %>%
  rename(
    category = vaccin
  ) %>%
  mutate(
    virus = "Influenza B",
    covariate = "Vaccination status",
    .before = 1
  )


results_sexe_table_B <- results_sexe_B$effect_summary %>%
  rename(
    category = sexe
  ) %>%
  mutate(
    virus = "Influenza B",
    covariate = "Sex",
    .before = 1
  )


results_effects_table_B <- bind_rows(
  results_vac_table_B,
  results_sexe_table_B
)

write.csv(
  results_effects_table_B,
  file.path(results_dir, "IBV_covariate_effects.csv"),
  row.names = FALSE
)

#################################################################################################
#################################################################################################
#
#                                   RSV
#
#################################################################################################
#################################################################################################


# ================================================================================================
# RSV simulation files
# ================================================================================================

sim_files_RSV <- c(
  
  "/home/laura.mulas/Monolix/SH/results/simulation_TV_VRS_1000_IC_traj_vfinal.rds"
)



# ================================================================================================
# 5.1 AGE
# ================================================================================================

cat("\n\nRSV - AGE\n")


df_age_RSV <- read_and_reduce_simulations(
  files = sim_files_RSV,
  covariate = age_cat
)


results_age_RSV <- compute_covariate_effects_reduced(
  df_group_sim = df_age_RSV,
  covariate = age_cat,
  reference = "18-65",
  threshold = Ct_LOD
)


rm(df_age_RSV)
gc()



# ================================================================================================
# 5.2 EPIDEMIC SEASON
# ================================================================================================

cat("\n\nRSV - SEASON\n")


df_season_RSV <- read_and_reduce_simulations(
  files = sim_files_RSV,
  covariate = season
)


results_season_RSV <- compute_covariate_effects_reduced(
  df_group_sim = df_season_RSV,
  covariate = season,
  reference = "24_25",
  threshold = Ct_LOD
)


rm(df_season_RSV)
gc()



# ================================================================================================
# 5.3 EPIDEMIC PERIOD
# ================================================================================================

cat("\n\nRSV - EPIDEMIC PERIOD\n")


df_period_RSV <- read_and_reduce_simulations(
  files = sim_files_RSV,
  covariate = periode_epi
)


results_period_RSV <- compute_covariate_effects_reduced(
  df_group_sim = df_period_RSV,
  covariate = periode_epi,
  reference = "Phase 2",
  threshold = Ct_LOD
)


rm(df_period_RSV)
gc()



# ================================================================================================
# RSV RESULTS TABLE
# ================================================================================================

results_age_table_RSV <- results_age_RSV$effect_summary %>%
  rename(
    category = age_cat
  ) %>%
  mutate(
    virus = "RSV",
    covariate = "Age",
    .before = 1
  )


results_season_table_RSV <- results_season_RSV$effect_summary %>%
  rename(
    category = season
  ) %>%
  mutate(
    virus = "RSV",
    covariate = "Epidemic season",
    .before = 1
  )


results_period_table_RSV <- results_period_RSV$effect_summary %>%
  rename(
    category = periode_epi
  ) %>%
  mutate(
    virus = "RSV",
    covariate = "Epidemic period",
    .before = 1
  )


results_effects_table_RSV <- bind_rows(
  results_age_table_RSV,
  results_season_table_RSV,
  results_period_table_RSV
)



write.csv(
  results_effects_table_RSV,
  file.path(results_dir, "RSV_covariate_effects.csv"),
  row.names = FALSE
)


# ================================================================================================
# COMBINE VIRUS-SPECIFIC EFFECT TABLES FOR THE FINAL FIGURE
# ================================================================================================

results_effects_table_all <- bind_rows(
  results_effects_table_cov,
  results_effects_table_A,
  results_effects_table_B,
  results_effects_table_RSV
)

# ================================================================================================
# PLOT 1) BASIC CLEANING
# ================================================================================================
df_effects <- results_effects_table_all %>%
  mutate(
    
    # Empty character strings -> NA
    across(
      where(is.character),
      ~ na_if(trimws(.x), "")
    ),
    
    # Ensure numeric variables are numeric
    across(
      c(
        delta_Ct_peak_median,
        delta_Ct_peak_low,
        delta_Ct_peak_high,
        delta_clearance_median,
        delta_clearance_low,
        delta_clearance_high
      ),
      ~ as.numeric(gsub(",", ".", .))
    )
  )


# ================================================================================================
# 3) VIRUS COLOURS
# ================================================================================================

virus_cols <- c(
  "SARS-CoV-2" = "#F76548",
  "IAV"        = "#34C7BE",
  "IBV"        = "#4A80D9",
  "RSV"        = "#F5BE3B"
)


virus_levels <- c(
  "SARS-CoV-2",
  "IAV",
  "IBV",
  "RSV"
)


# ================================================================================================
# 4) CLEAN LABELS
# ================================================================================================

df_effects <- df_effects %>%
  mutate(
    
    virus = case_when(
      virus == "Influenza A" ~ "IAV",
      virus == "Influenza B" ~ "IBV",
      TRUE ~ virus
    ),
    
    # --------------------------------------------------
    # Covariate groups
    # --------------------------------------------------
    
    variable_group = recode(
      covariate,
      
      "Age"    = "Age category",
      "Sex"   = "Sex",
      "Epidemic period" = "Epidemic period",
      "Epidemic season" = "Epidemic season",
      "Vaccination" = "Vaccination status",
      "Vaccination status" = "Vaccination status",
      
      .default = covariate
    ),
    
    
    # --------------------------------------------------
    # Compared category
    # --------------------------------------------------
    
    term_label = recode(
      category,
      
      # Age
      "<5"    = "<5 years",
      "5-18"  = "5-18 years",
      "18-65" = "19-65 years",
      ">65"   = ">65 years",
      
      # Sex
      "male"   = "Male",
      "Male"   = "Male",
      "female" = "Female",
      "Female" = "Female",
      
      # Epidemic period
      "Phase 1" = "Early",
      "Phase 2"   = "Mid",
      "Phase 3"  = "Late",
      
      # Season
      "24_25" = "2024/2025",
      "24-25" = "2024/2025",
      "2024/25" = "2024/2025",
      
      "25_25" = "2025/2026",
      "25_26" = "2025/2026",
      "25-26" = "2025/2026",
      "2025/26" = "2025/2026",
      
      # Vaccination
      "0" = "Unvaccinated",
      "1" = "Vaccinated",
      
      .default = category
    ),
    
    
    # --------------------------------------------------
    # Reference category
    # --------------------------------------------------
    
    ref_label = recode(
      covariate,
      
      "Age"    = "19-65 years",
      "Sex"   = "Female",
      "Epidemic period" = "Mid",
      "Epidemic season" = "2024/2025",
      "Vaccination" = "Unvaccinated",
      "Vaccination status" = "Unvaccinated",
      
      .default = covariate
    ),
  )


# ================================================================================================
# 5) ORDER OF GROUPS / VIRUSES
# ================================================================================================

facet_order <- c(
  "Age category",
  "Vaccination status",
  "Sex",
  "Epidemic period",
  "Epidemic season"
)


df_effects <- df_effects %>%
  mutate(
    
    variable_group = factor(
      variable_group,
      levels = facet_order
    ),
    
    virus = factor(
      virus,
      levels = virus_levels
    )
  )


# ================================================================================================
# 6) LEFT-HAND LABELS
# ================================================================================================

# One label for each category.
# We take the first reference label, because the reference should normally
# be identical across viruses for a given comparison.

label_df <- df_effects %>%
  group_by(
    variable_group,
    term_label
  ) %>%
  summarise(
    
    ref_label = first(
      ref_label[!is.na(ref_label)]
    ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    label_md = paste0(
      "<b>",
      term_label,
      "</b><br>",
      "<span style='color:black;'>(ref ",
      ref_label,
      ")</span>"
    )
  )


# ================================================================================================
# 7) Y POSITIONS
# ================================================================================================

# Age now has THREE comparisons:
#
# <5      vs 18-65
# 5-18    vs 18-65
# >65     vs 18-65
#
# We leave some vertical space between them.

y_map <- bind_rows(
  
  # --------------------------------------------------
  # AGE
  # --------------------------------------------------
  
  tibble(
    variable_group = "Age category",
    term_label = c(
      "<5 years",
      "5-18 years",
      ">65 years"
    ),
    y = c(
      7,
      4,
      1
    )
  ),
  
  
  # --------------------------------------------------
  # VACCINATION
  # --------------------------------------------------
  
  tibble(
    variable_group = "Vaccination status",
    term_label = "Vaccinated",
    y = 1
  ),
  
  
  # --------------------------------------------------
  # SEX
  # --------------------------------------------------
  
  tibble(
    variable_group = "Sex",
    term_label = "Male",
    y = 1
  ),
  
  
  # --------------------------------------------------
  # EPIDEMIC PERIOD
  # --------------------------------------------------
  
  tibble(
    variable_group = "Epidemic period",
    term_label = c(
      "Early",
      "Late"
    ),
    y = c(
      4,
      1
    )
  ),
  
  
  # --------------------------------------------------
  # SEASON
  # --------------------------------------------------
  
  tibble(
    variable_group = "Epidemic season",
    term_label = "2025/2026",
    y = 1
  )
)


# Convert group to character temporarily for joining

label_df <- label_df %>%
  mutate(
    variable_group = as.character(variable_group)
  )


df_effects <- df_effects %>%
  mutate(
    variable_group = as.character(variable_group)
  )


# Join Y positions

label_df_plot <- label_df %>%
  left_join(
    y_map,
    by = c(
      "variable_group",
      "term_label"
    )
  ) %>%
  filter(!is.na(y))


effects_plot <- df_effects %>%
  left_join(
    y_map,
    by = c(
      "variable_group",
      "term_label"
    )
  )


# ================================================================================================
# 8) OFFSET VIRUSES VERTICALLY
# ================================================================================================

effects_plot <- effects_plot %>%
  mutate(
    
    y_offset = case_when(
      
      virus == "SARS-CoV-2" ~  0.24,
      virus == "IAV"        ~  0.08,
      virus == "IBV"        ~ -0.08,
      virus == "RSV"        ~ -0.24,
      
      TRUE ~ 0
    ),
    
    y_plot = y + y_offset
  ) %>%
  filter(!is.na(y_plot))


# ================================================================================================
# 9) CHECK THAT ALL CATEGORIES HAVE A Y POSITION
# ================================================================================================

if (any(is.na(label_df_plot$y))) {
  
  cat(
    "\nWARNING: Some labels have no y position:\n\n"
  )
  
  print(
    label_df_plot %>%
      filter(is.na(y)) %>%
      select(
        variable_group,
        term_label,
        ref_label
      )
  )
}


if (any(is.na(effects_plot$y_plot))) {
  
  cat(
    "\nWARNING: Some effects have no y position:\n\n"
  )
  
  print(
    effects_plot %>%
      filter(is.na(y_plot)) %>%
      distinct(
        variable_group,
        term_label
      )
  )
}


# ================================================================================================
# 10) CAP VALUES AT GRAPH LIMITS
# ================================================================================================

# Ct axis: -2 to +2
# Clearance axis: -12 to +12

effects_plot <- effects_plot %>%
  mutate(
    
    # --------------------------------------------------
    # Peak Ct
    # --------------------------------------------------
    
    delta_ct_peak_disp = pmax(
      pmin(
        delta_Ct_peak_median,
        2
      ),
      -2
    ),
    
    delta_ct_peak_low_disp = pmax(
      delta_Ct_peak_low,
      -2
    ),
    
    delta_ct_peak_high_disp = pmin(
      delta_Ct_peak_high,
      2
    ),
    
    
    # --------------------------------------------------
    # Time to clearance
    # --------------------------------------------------
    
    delta_time_clearance_disp = pmax(
      pmin(
        delta_clearance_median,
        12
      ),
      -12
    ),
    
    delta_time_clearance_low_disp = pmax(
      delta_clearance_low,
      -12
    ),
    
    delta_time_clearance_high_disp = pmin(
      delta_clearance_high,
      12
    )
  )


# ================================================================================================
# 11) COMMON STRIP
# ================================================================================================

make_strip <- function(label) {
  
  ggplot() +
    
    annotate(
      "rect",
      
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      
      fill = "grey95",
      colour = "grey20",
      linewidth = 0.6
    ) +
    
    annotate(
      "text",
      
      x = 0,
      y = 0,
      
      label = label,
      
      size = 6,
      fontface = "bold"
    ) +
    
    xlim(
      -1,
      1
    ) +
    
    ylim(
      -1,
      1
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(
        0,
        0,
        0,
        0
      )
    )
}


# ================================================================================================
# 12) COLUMN TITLES
# ================================================================================================

make_title <- function(label) {
  
  ggplot() +
    
    annotate(
      "text",
      
      x = 0.5,
      y = 0.5,
      
      label = label,
      
      size = 7,
      fontface = "bold"
    ) +
    
    xlim(
      0,
      1
    ) +
    
    ylim(
      0,
      1
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(
        0,
        0,
        2,
        0
      )
    )
}


title_ct <- make_title(
  "Peak Ct"
)


title_clear <- make_title(
  "Time to viral clearance (days)"
)


# ================================================================================================
# 13) FUNCTION: Y LIMITS FOR EACH GROUP
# ================================================================================================

get_y_limits <- function(group_name) {
  
  if (group_name == "Age category") {
    
    return(
      c(
        0.2,
        7.8
      )
    )
    
  }
  
  
  if (group_name == "Epidemic period") {
    
    return(
      c(
        0.2,
        4.8
      )
    )
    
  }
  
  
  return(
    c(
      0.2,
      1.8
    )
  )
}


# ================================================================================================
# 14) LEFT LABEL PLOT
# ================================================================================================

make_label_plot <- function(group_name) {
  
  dd <- label_df_plot %>%
    filter(
      variable_group == group_name
    )
  
  
  ylim_use <- get_y_limits(
    group_name
  )
  
  
  ggplot(
    dd,
    aes(
      x = 1,
      y = y
    )
  ) +
    
    geom_richtext(
      aes(
        label = label_md
      ),
      
      hjust = 1,
      vjust = 0.5,
      
      fill = NA,
      label.color = NA,
      
      size = 6,
      lineheight = 1.05
    ) +
    
    scale_x_continuous(
      limits = c(
        0,
        1
      ),
      
      expand = c(
        0,
        0
      )
    ) +
    
    scale_y_continuous(
      limits = ylim_use,
      
      breaks = NULL,
      
      expand = c(
        0,
        0
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(
        2,
        0,
        2,
        30
      )
    )
}


# ================================================================================================
# 15) FOREST PLOT: DIFFERENCE IN PEAK Ct
# ================================================================================================

make_ct_plot <- function(group_name) {
  
  dd <- effects_plot %>%
    filter(
      variable_group == group_name
    )
  
  
  ylim_use <- get_y_limits(
    group_name
  )
  
  
  p <- ggplot(
    dd,
    aes(
      x = delta_ct_peak_disp,
      y = y_plot,
      colour = virus
    )
  ) +
    
    # Zero
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    # Confidence / uncertainty intervals
    geom_errorbarh(
      aes(
        xmin = delta_ct_peak_low_disp,
        xmax = delta_ct_peak_high_disp
      ),
      
      height = 0,
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    
    # Point estimates
    geom_point(
      size = 4,
      na.rm = TRUE
    ) +
    
    # Colours
    scale_color_manual(
      values = virus_cols,
      breaks = virus_levels,
      limits = virus_levels,
      name = "Virus",
      drop = FALSE
    ) +
    
    # Reversed Ct axis:
    # positive delta Ct = lower viral load -> shown on LEFT
    scale_x_reverse(
      
      name = NULL,
      
      limits = c(
        2,
        -2
      ),
      
      breaks = seq(
        -2,
        2,
        by = 0.5
      ),
      
      expand = c(
        0,
        0
      )
    ) +
    
    scale_y_continuous(
      
      limits = ylim_use,
      
      breaks = NULL,
      
      expand = c(
        0,
        0
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_classic(
      base_size = 13
    ) +
    
    theme(
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      
      legend.position = "none",
      
      plot.margin = margin(
        2,
        15,
        2,
        0
      )
    )
  
  
  # ==============================================================================================
  # Only show X axis for last row = Season
  # ==============================================================================================
  
  if (group_name != "Epidemic season") {
    
    p <- p +
      
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )
    
  } else {
    
    p <- p +
      
      theme(
        
        axis.text.x = element_text(
          size = 16
        ),
        
        axis.ticks.x = element_line(),
        
        axis.line.x = element_line(),
        
        plot.margin = margin(
          2,
          15,
          40,
          0
        )
      ) +
      
      
      # ------------------------------------------------------------------------------------------
    # Left arrow = lower viral load
    # ------------------------------------------------------------------------------------------
    
    annotation_custom(
      
      grob = segmentsGrob(
        
        x0 = unit(
          0.49,
          "npc"
        ),
        
        x1 = unit(
          0.25,
          "npc"
        ),
        
        y0 = unit(
          -0.90,
          "npc"
        ),
        
        y1 = unit(
          -0.90,
          "npc"
        ),
        
        arrow = arrow(
          length = unit(
            0.18,
            "cm"
          )
        ),
        
        gp = gpar(
          col = "black",
          lwd = 1
        )
      )
    ) +
      
      
      annotation_custom(
        
        grob = textGrob(
          
          "Lower viral load",
          
          x = unit(
            0.37,
            "npc"
          ),
          
          y = unit(
            -0.70,
            "npc"
          ),
          
          gp = gpar(
            fontsize = 17
          )
        )
      ) +
      
      
      # ------------------------------------------------------------------------------------------
    # Right arrow = higher viral load
    # ------------------------------------------------------------------------------------------
    
    annotation_custom(
      
      grob = segmentsGrob(
        
        x0 = unit(
          0.51,
          "npc"
        ),
        
        x1 = unit(
          0.75,
          "npc"
        ),
        
        y0 = unit(
          -0.90,
          "npc"
        ),
        
        y1 = unit(
          -0.90,
          "npc"
        ),
        
        arrow = arrow(
          length = unit(
            0.18,
            "cm"
          )
        ),
        
        gp = gpar(
          col = "black",
          lwd = 1
        )
      )
    ) +
      
      
      annotation_custom(
        
        grob = textGrob(
          
          "Higher viral load",
          
          x = unit(
            0.63,
            "npc"
          ),
          
          y = unit(
            -0.70,
            "npc"
          ),
          
          gp = gpar(
            fontsize = 17
          )
        )
      )
  }
  
  
  p
}


# ================================================================================================
# 16) FOREST PLOT: DIFFERENCE IN TIME TO CLEARANCE
# ================================================================================================

make_clear_plot <- function(group_name) {
  
  dd <- effects_plot %>%
    filter(
      variable_group == group_name
    )
  
  
  ylim_use <- get_y_limits(
    group_name
  )
  
  
  p <- ggplot(
    dd,
    aes(
      x = delta_time_clearance_disp,
      y = y_plot,
      colour = virus
    )
  ) +
    
    # Zero
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    # Confidence / uncertainty intervals
    geom_errorbarh(
      aes(
        xmin = delta_time_clearance_low_disp,
        xmax = delta_time_clearance_high_disp
      ),
      
      height = 0,
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    
    # Point estimates
    geom_point(
      size = 4,
      na.rm = TRUE
    ) +
    
    # Colours
    scale_color_manual(
      values = virus_cols,
      breaks = virus_levels,
      limits = virus_levels,
      name = "Virus",
      drop = FALSE
    ) +
    
    # Same orientation as your previous figure:
    # negative = shorter -> LEFT
    # positive = longer -> RIGHT
    #
    # limits c(12, -12) reverses the physical axis,
    # exactly as in your previous code.
    scale_x_continuous(
      
      name = NULL,
      
      limits = c(
        12,
        -12
      ),
      
      breaks = seq(
        -12,
        12,
        by = 2
      ),
      
      expand = c(
        0,
        0
      )
    ) +
    
    scale_y_continuous(
      
      limits = ylim_use,
      
      breaks = NULL,
      
      expand = c(
        0,
        0
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_classic(
      base_size = 13
    ) +
    
    theme(
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      
      legend.position = "none",
      
      plot.margin = margin(
        2,
        20,
        2,
        0
      )
    )
  
  
  # ==============================================================================================
  # Only show X axis for Season
  # ==============================================================================================
  
  if (group_name != "Epidemic season") {
    
    p <- p +
      
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )
    
  } else {
    
    p <- p +
      
      theme(
        
        axis.text.x = element_text(
          size = 16
        ),
        
        axis.ticks.x = element_line(),
        
        axis.line.x = element_line(),
        
        plot.margin = margin(
          2,
          20,
          40,
          0
        )
      ) +
      
      
      # ------------------------------------------------------------------------------------------
    # Left arrow = shorter clearance
    # ------------------------------------------------------------------------------------------
    
    annotation_custom(
      
      grob = segmentsGrob(
        
        x0 = unit(
          0.49,
          "npc"
        ),
        
        x1 = unit(
          0.15,
          "npc"
        ),
        
        y0 = unit(
          -0.90,
          "npc"
        ),
        
        y1 = unit(
          -0.90,
          "npc"
        ),
        
        arrow = arrow(
          length = unit(
            0.18,
            "cm"
          )
        ),
        
        gp = gpar(
          col = "black",
          lwd = 1
        )
      )
    ) +
      
      
      annotation_custom(
        
        grob = textGrob(
          
          "Shorter time to clearance",
          
          x = unit(
            0.32,
            "npc"
          ),
          
          y = unit(
            -0.70,
            "npc"
          ),
          
          gp = gpar(
            fontsize = 17
          )
        )
      ) +
      
      
      # ------------------------------------------------------------------------------------------
    # Right arrow = longer clearance
    # ------------------------------------------------------------------------------------------
    
    annotation_custom(
      
      grob = segmentsGrob(
        
        x0 = unit(
          0.52,
          "npc"
        ),
        
        x1 = unit(
          0.85,
          "npc"
        ),
        
        y0 = unit(
          -0.90,
          "npc"
        ),
        
        y1 = unit(
          -0.90,
          "npc"
        ),
        
        arrow = arrow(
          length = unit(
            0.18,
            "cm"
          )
        ),
        
        gp = gpar(
          col = "black",
          lwd = 1
        )
      )
    ) +
      
      
      annotation_custom(
        
        grob = textGrob(
          
          "Longer time to clearance",
          
          x = unit(
            0.68,
            "npc"
          ),
          
          y = unit(
            -0.70,
            "npc"
          ),
          
          gp = gpar(
            fontsize = 17
          )
        )
      )
  }
  
  
  p
}

# ================================================================================================
# 16bis) LÉGENDE UNIQUE (4 virus)
# ================================================================================================

make_legend_grob <- function() {
  
  dummy <- tibble(
    virus = factor(virus_levels, levels = virus_levels),
    x = 1,
    y = 1
  )
  
  p <- ggplot(dummy, aes(x = x, y = y, colour = virus)) +
    geom_point(size = 4) +
    scale_color_manual(
      values = virus_cols,
      breaks = virus_levels,
      limits = virus_levels,
      name   = "Virus",
      drop   = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text  = element_text(size = 15)
    )
  
  cowplot::get_legend(p)
}

legend_grob <- make_legend_grob()

# ================================================================================================
# 17) BUILD A COMPLETE BLOCK FOR ONE COVARIATE
#
# Design:
#
#             ┌────────────────────────────────────────────────────┐
#             │                        AGE                         │
#             └────────────────────────────────────────────────────┘
#
# labels          Peak Ct column            Clearance column
#
#
# IMPORTANT:
# "ABB"
# "CDE"
#
# means that plot B (the grey strip) SPANS columns 2 AND 3.
# ================================================================================================
make_vertical_separator <- function() {
  
  ggplot() +
    
    annotate(
      "segment",
      x = 0.5,
      xend = 0.5,
      y = 0,
      yend = 1,
      colour = "grey80",
      linewidth = 0.6
    ) +
    
    xlim(0, 1) +
    ylim(0, 1) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
}

p_sep <- make_vertical_separator()

make_group_block <- function(
    group_name,
    body_height = 1
) {
  
  p_strip <- make_strip(group_name)
  p_label <- make_label_plot(group_name)
  p_ct <- make_ct_plot(group_name)
  p_clear <- make_clear_plot(group_name)
  
  
  # A = espace à gauche
  # B = strip commun sur Ct + espace + clearance
  # C = labels
  # D = Peak Ct
  # E = espace entre les deux graphiques
  # F = Clearance
  
  design_group <- "
ABBB
CDEF
"
  
  wrap_plots(
    
    plot_spacer(),
    p_strip,
    
    p_label,
    p_ct,
    p_sep,
    p_clear,
    
    design = "
ABBB
CDEF
",
    
    widths = c(
      1,
      5,
      0.5,
      5
    ),
    
    heights = c(
      0.50,
      body_height
    )
  )
}


# ================================================================================================
# 18) BUILD THE FIVE COVARIATE BLOCKS
# ================================================================================================

# Age has 3 rows -> larger

block_age <- make_group_block(
  "Age category",
  body_height = 3
)


# Vaccination has 1 row

block_vaccination <- make_group_block(
  "Vaccination status",
  body_height = 1
)


# Sex has 1 row

block_sex <- make_group_block(
  "Sex",
  body_height = 1
)


# Period has 2 rows

block_period <- make_group_block(
  "Epidemic period",
  body_height = 2
)


# Season has 1 row

block_season <- make_group_block(
  "Epidemic season",
  body_height = 1
)


# ================================================================================================
# 19) TOP TITLES
# ================================================================================================

title_row <- wrap_plots(
  
  plot_spacer(),
  title_ct,
  plot_spacer(),
  title_clear,
  
  nrow = 1,
  
  widths = c(
    1,
    5,
    0.5,
    5
  )
)

# ================================================================================================
# 20) FINAL FIGURE
# ================================================================================================

final_combined <- (
  
  title_row /
    block_age /
    block_vaccination /
    block_sex /
    block_period /
    block_season /
    wrap_elements(full = legend_grob)
  
) +
  
  plot_layout(
    heights = c(
      0.65,   # Titres
      3.30,   # Age
      1.30,   # Vaccination
      1.30,   # Sex
      2.30,   # Period
      1.30,   # Season
      0.30    # Légende
    )
  )

# Display
final_combined


# ================================================================================================
# 21) SAVE SVG
# ================================================================================================

ggsave(
  
  filename = figure_file,
  
  plot = final_combined,
  
  device = svglite::svglite,
  
  width = 1700 / 96,
  
  height = 950 / 96
)

