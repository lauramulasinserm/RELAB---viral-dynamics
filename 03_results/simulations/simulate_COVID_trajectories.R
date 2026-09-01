####################################################################################################
# Simulation of SARS-CoV-2 viral trajectories
# Author: Laura MULAS
# Purpose: Propagate parameter uncertainty and inter-individual variability
#          to simulate individual Ct trajectories and peak viral kinetics.
####################################################################################################


# ==================== Libraries ====================
library(readr)
library(dplyr)
library(tidyr)
library(MASS)
library(deSolve)

# Input data and parameter files
data_csv <- "02_monolix/data/data_for_monolix_SARS_CoV_2.csv"
popparam_txt <- "03_results_analysis/model_outputs/COVID/populationParameters.txt"
cov_txt <- "03_results_analysis/model_outputs/COVID/FisherInformation/covarianceEstimatesLin.txt"

# Output directory
outdir <- "03_results_analysis/results/simulations/COVID"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ==================== Load and Prepare Data ====================
data_monolix <- read.csv(data_csv, sep=",", stringsAsFactors = FALSE)

# Rename and harmonize variable names
data <- data_monolix %>%
  mutate(
    vaccin = statut_vaccin,
    periode_epi = epidemic_period,
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
chunk_size <- 250
random_seed <- 123

set.seed(random_seed)

# ==================== Population Parameters ====================
# Read population-level parameters (fixed effects + omega)
df_params <- readr::read_delim(popparam_txt)

# Extract fixed effects and omega
mu_full <- df_params$value[1:17]
names(mu_full) <- df_params$parameter[1:17]

# Fixed parameters are excluded from the MVN draw because they are not in the covariance matrix.
fixed_params_names <- c("eta_pop", "c_pop", "tinc_pop")

# Apply log-transform for parameters defined on log-scale
mu_full_log <- mu_full
log_params <- c("beta_pop", "p_pop", "delta_pop")
mu_full_log[log_params] <- log(mu_full[log_params])

mu_est <- mu_full_log[setdiff(names(mu_full_log), fixed_params_names)]

# Read covariance matrix (Fisher information estimates)
Sigma_full <- as.matrix(readr::read_delim(cov_txt, col_names = FALSE)[1:14,2:15])
rownames(Sigma_full) <- names(mu_est)
colnames(Sigma_full) <- names(mu_est)

# ==================== Simulate Population Parameters ====================
# Draw random samples of parameters from the multivariate normal distribution
sim_est <- MASS::mvrnorm(n = n_sim, mu = mu_est, Sigma = Sigma_full)
sim_est <- as.data.frame(sim_est)

# Reconstruct the full parameter set: sampled estimated parameters + fixed constants.
sim_full_df <- matrix(rep(mu_full_log, each = n_sim), nrow = n_sim, byrow = TRUE)
colnames(sim_full_df) <- names(mu_full_log)
sim_df <- as.data.frame(sim_full_df)

sim_df[, names(mu_est)] <- sim_est

# Extract omega and ensure positivity
omega_beta  <- sim_df$omega_beta
omega_delta <- sim_df$omega_delta

omega_beta[omega_beta <= 0]   <- 1e-6
omega_delta[omega_delta <= 0] <- 1e-6

sd_beta  <- sqrt(omega_beta)
sd_delta <- sqrt(omega_delta)

# Simulate individual random effects for beta and delta for each subject and parameter draw.
ranef_beta <- matrix(
  rnorm(n_ind * n_sim, mean = 0, sd = rep(sd_beta,  each = n_ind)),
  nrow = n_ind, ncol = n_sim
)

ranef_delta <- matrix(
  rnorm(n_ind * n_sim, mean = 0, sd = rep(sd_delta, each = n_ind)),
  nrow = n_ind, ncol = n_sim
)

# Keep only fixed effects
omega_names <- c("omega_beta", "omega_delta")


theta_names <- setdiff(colnames(sim_df), omega_names)
theta_sim <- sim_df[, theta_names, drop = FALSE]



# ==================== Main Simulation Loop ====================
beta_sim  <- matrix(NA, nrow = n_ind, ncol = n_sim)
p_sim     <- matrix(NA, nrow = n_ind, ncol = n_sim)
delta_sim <- matrix(NA, nrow = n_ind, ncol = n_sim)

for (i in 1:n_sim) {
  theta_i <- theta_sim[i, ]  # Fixed effects for simulation i
  
  for (j in 1:nrow(df_sim)) {
    ID_j <- df_sim$ID[j]
    
    log_beta <- theta_i["beta_pop"] +
      (if("beta_beta_age_cat__5"   %in% names(theta_i) && df_sim$age_cat[j] == "<5")   theta_i["beta_beta_age_cat__5"]   else 0) +
      (if("beta_beta_age_cat_5_18"   %in% names(theta_i) && df_sim$age_cat[j] == "5-18")   theta_i["beta_beta_age_cat_5_18"]   else 0) +
      (if("beta_beta_age_cat__65"   %in% names(theta_i) && df_sim$age_cat[j] == ">65")   theta_i["beta_beta_age_cat__65"]   else 0) +
      (if("beta_beta_season_25_26"   %in% names(theta_i) && df_sim$season[j] == "25_26")   theta_i["beta_beta_season_25_26"]   else 0) +
      ranef_beta[j, i]
    
    beta <- exp(as.numeric(log_beta))
    
    
    log_p <- theta_i["p_pop"]
    
    p <- exp(as.numeric(log_p))
    
    log_delta <- theta_i["delta_pop"] +
      (if("beta_delta_age_cat__5"   %in% names(theta_i) && df_sim$age_cat[j] == "<5")   theta_i["beta_delta_age_cat__5"]   else 0) +
      (if("beta_delta_age_cat_5_18"   %in% names(theta_i) && df_sim$age_cat[j] == "5-18")   theta_i["beta_delta_age_cat_5_18"]   else 0) +
      (if("beta_delta_age_cat__65"   %in% names(theta_i) && df_sim$age_cat[j] == ">65")   theta_i["beta_delta_age_cat__65"]   else 0) +
      (if("beta_delta_period_fus_G_Phase_3"   %in% names(theta_i) && df_sim$periode_epi[j] == "Phase 3")   theta_i["beta_delta_period_fus_G_Phase_3"]   else 0) +
      (if("beta_delta_season_25_26"   %in% names(theta_i) && df_sim$season[j] == "25_26")   theta_i["beta_delta_season_25_26"]   else 0) +
      ranef_delta[j, i]
    
    delta <- exp(as.numeric(log_delta))
    
    beta_sim[j, i]  <- beta
    p_sim[j, i]     <- p
    delta_sim[j, i] <- delta
    
  }
  
  cat("Simulation", i, "completed\n")
}



# Fixed constants (not sampled).
ETA  <- as.numeric(mu_full["eta_pop"])
C    <- as.numeric(mu_full["c_pop"])
TINC <- as.numeric(mu_full["tinc_pop"])
T0   <- 8e7
TIMES <- seq(-5, 30, by = 0.1)

# ODE model
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
# MAIN SIMULATION LOOP
# =============================================================================

chunk_starts <- seq(
  1,
  n_sim,
  by = chunk_size
)


for (chunk_start in chunk_starts) {
  
  chunk_end <- min(
    chunk_start + chunk_size - 1,
    n_sim
  )
  
  cat("\n")
  cat("============================================================\n")
  cat(
    "SIMULATIONS",
    chunk_start,
    "TO",
    chunk_end,
    "\n"
  )
  cat("============================================================\n")
  
  
  # ===========================================================================
  # INITIALISE LISTS FOR CURRENT CHUNK ONLY
  # ===========================================================================
  
  n_sim_chunk <- chunk_end - chunk_start + 1
  
  sim_list_chunk <- vector(
    "list",
    n_sim_chunk * n_ind
  )
  
  results_list_chunk <- vector(
    "list",
    n_sim_chunk * n_ind
  )
  
  idx <- 1
  
  
  # ===========================================================================
  # RUN SIMULATIONS FOR CURRENT CHUNK
  # ===========================================================================
  
  for (i in chunk_start:chunk_end) {
    
    for (j in 1:n_ind) {
      
      beta  <- beta_sim[j, i]
      p     <- p_sim[j, i]
      delta <- delta_sim[j, i]
      
      
      # -----------------------------------------------------------------------
      # Initial infectious viral load and R0
      # -----------------------------------------------------------------------
      
      VI0 <- ETA * p / C
      
      R0 <- ETA * T0 * p * beta /
        (delta * C)
      
      
      # -----------------------------------------------------------------------
      # R0 <= 1
      # -----------------------------------------------------------------------
      
      if (R0 <= 1) {
        
        Vpeak_ana    <- NA
        Vpeak_ode_VI <- NA
        Vpeak_ode_V  <- NA
        Ct_peak_ode  <- NA
        
        sim_i <- NULL
        
        
        # -----------------------------------------------------------------------
        # R0 > 1
        # -----------------------------------------------------------------------
        
      } else {
        
        Vpeak_ana <-
          VI0 +
          (ETA * p / C) * T0 -
          (delta / beta) * (log(R0) + 1)
        
        
        # Initial conditions
        
        state_ini <- c(
          T   = T0,
          VI  = VI0,
          VNI = 0
        )
        
        
        # Parameters
        
        pars_i <- c(
          beta  = beta,
          c     = C,
          p     = p,
          eta   = ETA,
          delta = delta,
          tinc  = TINC
        )
        
        
        # ---------------------------------------------------------------------
        # Solve ODE
        # ---------------------------------------------------------------------
        
        sim_i <- tryCatch(
          
          as.data.frame(
            
            ode(
              y      = state_ini,
              times  = TIMES,
              func   = model_virus,
              parms  = pars_i,
              method = "lsoda"
            )
          ),
          
          error = function(e) NULL
        )
        
        
        # ---------------------------------------------------------------------
        # ODE failed
        # ---------------------------------------------------------------------
        
        if (is.null(sim_i)) {
          
          Vpeak_ode_VI <- NA
          Vpeak_ode_V  <- NA
          Ct_peak_ode  <- NA
          
          
          # ---------------------------------------------------------------------
          # ODE successful
          # ---------------------------------------------------------------------
          
        } else {
          
          Vpeak_ode_VI <- max(
            sim_i$VI,
            na.rm = TRUE
          )
          
          
          Vpeak_ode_V <- max(
            sim_i$VI + sim_i$VNI,
            na.rm = TRUE
          )
          
          
          # -------------------------------------------------------------------
          # Convert viral load to Ct
          # -------------------------------------------------------------------
          
          sim_i$Ct <- 49 -
            3 * log10(
              pmax(
                (sim_i$VI + sim_i$VNI) / 30,
                1e-2
              )
            )
          
          
          Ct_peak_ode <- min(
            sim_i$Ct,
            na.rm = TRUE
          )
          
          
          # -------------------------------------------------------------------
          # Add metadata
          # -------------------------------------------------------------------
          
          sim_i$sim <- i
          
          sim_i$ind <- j
          
          sim_i$ID <- df_sim$ID[j]
          
          sim_i$ID_infection <-
            df_sim$ID_infection[j]
          
          sim_i$age_cat <-
            df_sim$age_cat[j]
          
          sim_i$periode_epi <-
            df_sim$periode_epi[j]
          
          sim_i$season <-
            df_sim$season[j]
          
          sim_i$vaccin <-
            df_sim$vaccin[j]
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
        
        sim = i,
        
        ind = j,
        
        ID = df_sim$ID[j],
        
        ID_infection =
          df_sim$ID_infection[j],
        
        age_cat =
          df_sim$age_cat[j],
        
        periode_epi =
          df_sim$periode_epi[j],
        
        season =
          df_sim$season[j],
        
        vaccin =
          df_sim$vaccin[j],
        
        R0 = R0,
        
        Vpeak_ana =
          Vpeak_ana,
        
        Vpeak_ode_VI =
          Vpeak_ode_VI,
        
        Vpeak_ode_V =
          Vpeak_ode_V,
        
        Ct_peak_ode =
          Ct_peak_ode
      )
      
      
      idx <- idx + 1
    }
    
    
    cat(
      "Simulation",
      i,
      "completed\n"
    )
  }
  
  
  # ===========================================================================
  # REMOVE NULL TRAJECTORIES
  # ===========================================================================
  
  sim_list_chunk <- Filter(
    Negate(is.null),
    sim_list_chunk
  )
  
  
  # ===========================================================================
  # COMBINE TRAJECTORIES FOR CURRENT CHUNK
  # ===========================================================================
  
  df_traj_chunk <- dplyr::bind_rows(
    sim_list_chunk
  )
  
  
  df_long_sim_chunk <- df_traj_chunk %>%
    
    mutate(
      time_since_symptoms_onset = time,
      y_sim = Ct
    ) %>%
    
    dplyr::select(
      sim,
      ID,
      ID_infection,
      age_cat,
      periode_epi,
      season,
      vaccin,
      time_since_symptoms_onset,
      Ct,
      T,
      VI,
      VNI,
      V
    )
  
  
  # ===========================================================================
  # COMBINE PEAK RESULTS
  # ===========================================================================
  
  df_peak_chunk <- dplyr::bind_rows(
    results_list_chunk
  )
  
  
  # ===========================================================================
  # SAVE TRAJECTORIES
  # ===========================================================================
  
  outfile_traj <- file.path(
    
    outdir,
    
    paste0(
      "COVID_trajectories_",
      chunk_start,
      "_",
      chunk_end,
      ".rds"
    )
  )
  
  
  saveRDS(
    df_long_sim_chunk,
    file = outfile_traj
  )
  
  
  cat(
    "Trajectories saved ->",
    outfile_traj,
    "\n"
  )
  
  
  # ===========================================================================
  # SAVE PEAK RESULTS
  # ===========================================================================
  
  outfile_peak <- file.path(
    
    outdir,
    
    paste0(
      "COVID_peaks_",
      chunk_start,
      "_",
      chunk_end,
      ".rds"
    )
  )
  
  
  saveRDS(
    df_peak_chunk,
    file = outfile_peak
  )
  
  
  cat(
    "Peak results saved ->",
    outfile_peak,
    "\n"
  )
  
  
  # ===========================================================================
  # FREE MEMORY BEFORE NEXT CHUNK
  # ===========================================================================
  
  rm(
    sim_list_chunk,
    results_list_chunk,
    df_traj_chunk,
    df_long_sim_chunk,
    df_peak_chunk
  )
  
  
  gc()
  
  
  cat(
    "\nChunk",
    chunk_start,
    "-",
    chunk_end,
    "completed and memory cleared.\n"
  )
}