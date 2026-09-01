################################################################################
# Simulation of SARS-CoV-2 Ct trajectories with residual error
# Author: Laura Mulas
# Purpose: Simulate individual Ct trajectories from final Monolix estimates,
#          including inter-individual variability and additive residual error.
################################################################################
# ==================== Libraries ====================
library(readr)
library(dplyr)
library(deSolve)

# Input data and parameter file
data_csv <- "02_monolix/data/data_for_monolix_SARS_CoV_2.csv"
popparam_txt <- "03_results_analysis/model_outputs/COVID/populationParameters.txt"

# Output directory
outdir <- "03_results_analysis/results/simulations_error/COVID"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ==================== Load and prepare data ====================
data_monolix <- read.csv(data_csv, sep=",", stringsAsFactors = FALSE)

# Rename and harmonize variable names
data <- data_monolix %>%
  mutate(
    vaccin = statut_vaccin,
    periode_epi = epidemic_period,
    season = season,
    time_since_symptoms_onset = time,
    Ct_norm = obs,
    ID = ID_infection
  )

# Convert categorical variables to factors and set reference levels
data$sexe <- as.factor(data$sexe)
data$age_cat <- as.factor(data$age_cat)
data$vaccin <- as.factor(data$vaccin)
data$periode_epi <- as.factor(data$periode_epi)
data$season <- as.factor(data$season)
data$age_cat <- relevel(data$age_cat, ref = "18-65")
data$periode_epi <- relevel(data$periode_epi, ref = "Phase 2")
data$season <- relevel(data$season, ref = "24_25")

# Create a simplified dataset for simulation
df_sim <- data %>%
  distinct(ID_infection, age_cat, season, periode_epi, vaccin) %>%
  mutate(ID = row_number())

n_ind <- nrow(df_sim)
n_sim <- 1000
random_seed <- 123
set.seed(random_seed)

# ==================== Population parameters ====================
# Read population-level parameters (fixed effects + omega)
df_params <- readr::read_delim(popparam_txt, show_col_types = FALSE)

# ==================== Population parameters ====================
sigma_a     <- df_params$value[df_params$parameter == "a"]
omega_beta  <- df_params$value[df_params$parameter == "omega_beta"]
omega_delta <- df_params$value[df_params$parameter == "omega_delta"]

# Fixed effects on the original scale
beta_pop  <- df_params$value[df_params$parameter == "beta_pop"]
c_pop     <- df_params$value[df_params$parameter == "c_pop"]
p_pop     <- df_params$value[df_params$parameter == "p_pop"]
eta_pop   <- df_params$value[df_params$parameter == "eta_pop"]
delta_pop <- df_params$value[df_params$parameter == "delta_pop"]
tinc_pop  <- df_params$value[df_params$parameter == "tinc_pop"]

beta_beta_age_cat__5              <- df_params$value[df_params$parameter == "beta_beta_age_cat__5"]
beta_beta_age_cat_5_18	          <- df_params$value[df_params$parameter == "beta_beta_age_cat_5_18"]
beta_beta_age_cat__65             <- df_params$value[df_params$parameter == "beta_beta_age_cat__65"]
beta_beta_season_25_26            <- df_params$value[df_params$parameter == "beta_beta_season_25_26"]
beta_delta_age_cat__5             <- df_params$value[df_params$parameter == "beta_delta_age_cat__5"]
beta_delta_age_cat_5_18           <- df_params$value[df_params$parameter == "beta_delta_age_cat_5_18"]
beta_delta_age_cat__65            <- df_params$value[df_params$parameter == "beta_delta_age_cat__65"]
beta_delta_season_25_26           <- df_params$value[df_params$parameter == "beta_delta_season_25_26"]
beta_delta_period_fus_G_Phase_3   <- df_params$value[df_params$parameter == "beta_delta_period_fus_G_Phase_3"]


# ==================== Random effects ====================
# Monolix omega parameter used for random-effect dispersion on the log scale
ranef_beta <- matrix(
  rnorm(n_ind * n_sim, mean = 0, sd = sqrt(omega_beta)),
  nrow = n_ind, ncol = n_sim
)
ranef_delta <- matrix(
  rnorm(n_ind * n_sim, mean = 0, sd = sqrt(omega_delta)),
  nrow = n_ind, ncol = n_sim
)

# ==================== Fixed constants ====================
T0    <- 8e7
TIMES <- seq(-5, 30, by = 0.1)

# ==================== ODE model ====================
model_virus <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    dT   <- -beta * VI * T
    dVI  <-  eta * (p / c) * beta * VI * T - delta * VI
    dVNI <- (1 - eta) * (p / c) * beta * VI * T - delta * VNI
    V    <- VI + VNI
    list(c(dT, dVI, dVNI), c(V = V))
  })
}

# =============================================================================
# MAIN SIMULATION LOOP BY CHUNKS
# =============================================================================

chunk_size <- 250

chunk_starts <- seq(1, n_sim, by = chunk_size)

for (chunk_start in chunk_starts) {
  
  chunk_end <- min(chunk_start + chunk_size - 1, n_sim)
  
  cat("\n")
  cat("============================================================\n")
  cat("SIMULATIONS", chunk_start, "TO", chunk_end, "\n")
  cat("============================================================\n")
  
  # ===========================================================================
  # INITIALISE LISTS FOR CURRENT CHUNK ONLY
  # ===========================================================================
  
  n_sim_chunk <- chunk_end - chunk_start + 1
  
  sim_list_chunk     <- vector("list", n_sim_chunk * n_ind)
  results_list_chunk <- vector("list", n_sim_chunk * n_ind)
  
  idx <- 1
  
  # ===========================================================================
  # RUN SIMULATIONS FOR CURRENT CHUNK
  # ===========================================================================
  
  for (i in chunk_start:chunk_end) {
    for (j in seq_len(n_ind)) {
      
      # --- Individual parameters according to the Monolix model ---
      log_beta_i <- log(beta_pop) +
        beta_beta_age_cat__5 * (df_sim$age_cat[j] == "<5") +
        beta_beta_age_cat_5_18 * (df_sim$age_cat[j] == "5-18") +
        beta_beta_age_cat__65 * (df_sim$age_cat[j] == ">65") +
        beta_beta_season_25_26 * (df_sim$season[j] == "25_26") +
        ranef_beta[j, i]
      
      beta_i <- exp(log_beta_i)
      
      # log(c) = log(c_pop): no inter-individual variability
      c_i     <- c_pop
      
      # log(p) = log(p_pop): no inter-individual variability
      p_i     <- p_pop
      
      # log(eta) = log(eta_pop): no inter-individual variability
      eta_i   <- eta_pop
      
      # log(tinc) = log(tinc_pop): no inter-individual variability
      tinc_i  <- tinc_pop
      
      log_delta_i <- log(delta_pop) +
        beta_delta_age_cat__5 * (df_sim$age_cat[j] == "<5") +
        beta_delta_age_cat_5_18 * (df_sim$age_cat[j] == "5-18") +
        beta_delta_age_cat__65 * (df_sim$age_cat[j] == ">65") +
        beta_delta_season_25_26 * (df_sim$season[j] == "25_26") +
        beta_delta_period_fus_G_Phase_3 * (df_sim$periode_epi[j] == "Phase 3") +
        ranef_delta[j, i]
      
      delta_i <- exp(log_delta_i)
      
      # --- Initial conditions ---
      VI0       <- eta_i * p_i / c_i
      state_ini <- c(T = T0, VI = VI0, VNI = 0)
      pars_i    <- c(beta  = beta_i,
                     c     = c_i,
                     p     = p_i,
                     eta   = eta_i,
                     delta = delta_i,
                     tinc  = tinc_i)
      
      # --- R0 ---
      R0 <- eta_i * T0 * p_i * beta_i / (delta_i * c_i)
      
      # --- ODE solution ---
      if (R0 <= 1) {
        
        sim_i <- NULL
        
      } else {
        
        sim_i <- tryCatch(
          as.data.frame(ode(y      = state_ini,
                            times  = TIMES,
                            func   = model_virus,
                            parms  = pars_i,
                            method = "lsoda")),
          error = function(e) NULL
        )
        
        if (!is.null(sim_i)) {
          
          # Structural prediction (ipred)
          V_total  <- pmax(sim_i$VI + sim_i$VNI, 1e-2)
          Ct_ipred <- 49 - 3 * log10(V_total / 30)
          
          # Additive residual error: y = ipred + a * epsilon
          epsilon  <- rnorm(nrow(sim_i), mean = 0, sd = 1)
          sim_i$Ct <- Ct_ipred + sigma_a * epsilon
          
          # Censoring: Ct > 40 is considered undetectable
          sim_i$Ct_obs <- ifelse(sim_i$Ct > 40, NA, sim_i$Ct)
          
          # Metadata
          sim_i$sim          <- i
          sim_i$ind          <- j
          sim_i$ID           <- df_sim$ID[j]
          sim_i$ID_infection <- df_sim$ID_infection[j]
          sim_i$age_cat      <- df_sim$age_cat[j]
          sim_i$periode_epi  <- df_sim$periode_epi[j]
          sim_i$season       <- df_sim$season[j]
          sim_i$vaccin       <- df_sim$vaccin[j]
        }
      }
      
      # =========================================================================
      # STORE TRAJECTORY
      # =========================================================================
      
      sim_list_chunk[[idx]] <- sim_i
      
      # =========================================================================
      # STORE PEAK RESULTS
      # =========================================================================
      
      results_list_chunk[[idx]] <- data.frame(
        sim          = i,
        ind          = j,
        ID           = df_sim$ID[j],
        ID_infection = df_sim$ID_infection[j],
        age_cat      = df_sim$age_cat[j],
        periode_epi  = df_sim$periode_epi[j],
        season       = df_sim$season[j],
        vaccin       = df_sim$vaccin[j],
        R0           = R0
      )
      
      idx <- idx + 1
    }
    
    cat("Simulation", i, "completed\n")
  }
  
  # ===========================================================================
  # REMOVE NULL TRAJECTORIES
  # ===========================================================================
  
  sim_list_chunk <- Filter(Negate(is.null), sim_list_chunk)
  
  # ===========================================================================
  # COMBINE TRAJECTORIES FOR CURRENT CHUNK
  # ===========================================================================
  
  df_traj_chunk <- dplyr::bind_rows(sim_list_chunk)
  
  df_long_sim_chunk <- df_traj_chunk %>%
    rename(time_since_symptoms_onset = time) %>%
    dplyr::select(sim, ID, ID_infection, age_cat, periode_epi, season, vaccin,
                  time_since_symptoms_onset, Ct, Ct_obs, T, VI, VNI, V)
  
  # ===========================================================================
  # COMBINE PEAK RESULTS
  # ===========================================================================
  
  df_peak_chunk <- dplyr::bind_rows(results_list_chunk)
  
  # ===========================================================================
  # SAVE TRAJECTORIES
  # ===========================================================================
  
  outfile_traj <- file.path(
    outdir,
    paste0("simulation_TV_covid_", chunk_start, "_", chunk_end, "_error_traj_vfinal.rds")
  )
  
  saveRDS(df_long_sim_chunk, file = outfile_traj)
  
  cat("Trajectories saved ->", outfile_traj, "\n")
  
  # ===========================================================================
  # SAVE PEAK RESULTS
  # ===========================================================================
  
  outfile_peak <- file.path(
    outdir,
    paste0("simulation_TV_covid_", chunk_start, "_", chunk_end, "_error_peak_vfinal.rds")
  )
  
  saveRDS(df_peak_chunk, file = outfile_peak)
  
  cat("Peak results saved ->", outfile_peak, "\n")
  
  # ===========================================================================
  # FREE MEMORY BEFORE NEXT CHUNK
  # ===========================================================================
  
  rm(sim_list_chunk, results_list_chunk, df_traj_chunk, df_long_sim_chunk, df_peak_chunk)
  
  gc()
  
  cat("\nChunk", chunk_start, "-", chunk_end, "completed and memory cleared.\n")
}
