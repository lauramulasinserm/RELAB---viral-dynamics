################################################################################################
# Simulation of Influenza B virus (IBV) trajectories
# Author: Laura MULAS
# Purpose:
#   Propagate parameter uncertainty and inter-individual variability to simulate
#   individual Ct trajectories and peak viral kinetics.
################################################################################################

rm(list = ls())

# ---------------------------------------------------------------------------
# 0. Libraries
# ---------------------------------------------------------------------------

library(readr)
library(dplyr)
library(MASS)
library(deSolve)


# ---------------------------------------------------------------------------
# 1. File paths and simulation settings
# ---------------------------------------------------------------------------

data_csv <- "02_monolix/data/data_for_monolix_IBV.csv"

popparam_txt <- paste0(
  "03_results_analysis/model_outputs/IBV/",
  "populationParameters.txt"
)

cov_txt <- paste0(
  "03_results_analysis/model_outputs/IBV/",
  "FisherInformation/covarianceEstimatesLin.txt"
)

outdir <- "03_results_analysis/results/simulations/IBV"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

n_sim <- 1000
chunk_size <- 250
random_seed <- 123

set.seed(random_seed)


# ---------------------------------------------------------------------------
# 2. Load and prepare Monolix data
# ---------------------------------------------------------------------------

data <- read_csv(
  data_csv,
  show_col_types = FALSE
) %>%
  mutate(
    vaccin = factor(statut_vaccin),
    sexe = factor(sexe)
  )

df_sim <- data %>%
  distinct(
    ID_infection,
    sexe,
    vaccin
  ) %>%
  mutate(
    ID = row_number()
  )

n_ind <- nrow(df_sim)


# ---------------------------------------------------------------------------
# 3. Population parameters and uncertainty
# ---------------------------------------------------------------------------

df_params <- read_delim(
  popparam_txt,
  show_col_types = FALSE
)

mu_full <- df_params$value[1:10]
names(mu_full) <- df_params$parameter[1:10]

fixed_params_names <- c(
  "eta_pop",
  "c_pop",
  "tinc_pop"
)

mu_full_log <- mu_full

log_params <- c(
  "beta_pop",
  "p_pop",
  "delta_pop"
)

mu_full_log[log_params] <- log(
  mu_full[log_params]
)

mu_est <- mu_full_log[
  setdiff(
    names(mu_full_log),
    fixed_params_names
  )
]


Sigma_full <- as.matrix(
  read_delim(
    cov_txt,
    col_names = FALSE,
    show_col_types = FALSE
  )[1:7, 2:8]
)

rownames(Sigma_full) <- names(mu_est)
colnames(Sigma_full) <- names(mu_est)


sim_est <- MASS::mvrnorm(
  n = n_sim,
  mu = mu_est,
  Sigma = Sigma_full
) %>%
  as.data.frame()


sim_df <- matrix(
  rep(
    mu_full_log,
    each = n_sim
  ),
  nrow = n_sim,
  byrow = TRUE
) %>%
  as.data.frame()

colnames(sim_df) <- names(mu_full_log)

sim_df[, names(mu_est)] <- sim_est


# ---------------------------------------------------------------------------
# 4. Inter-individual variability
# ---------------------------------------------------------------------------

omega_beta <- pmax(
  sim_df$omega_beta,
  1e-6
)

omega_delta <- pmax(
  sim_df$omega_delta,
  1e-6
)

ranef_beta <- matrix(
  rnorm(
    n_ind * n_sim,
    mean = 0,
    sd = rep(
      sqrt(omega_beta),
      each = n_ind
    )
  ),
  nrow = n_ind,
  ncol = n_sim
)

ranef_delta <- matrix(
  rnorm(
    n_ind * n_sim,
    mean = 0,
    sd = rep(
      sqrt(omega_delta),
      each = n_ind
    )
  ),
  nrow = n_ind,
  ncol = n_sim
)


theta_sim <- sim_df[
  ,
  setdiff(
    colnames(sim_df),
    c(
      "omega_beta",
      "omega_delta"
    )
  ),
  drop = FALSE
]


# ---------------------------------------------------------------------------
# 5. Individual parameter values
# ---------------------------------------------------------------------------

beta_sim <- matrix(
  NA_real_,
  nrow = n_ind,
  ncol = n_sim
)

p_sim <- matrix(
  NA_real_,
  nrow = n_ind,
  ncol = n_sim
)

delta_sim <- matrix(
  NA_real_,
  nrow = n_ind,
  ncol = n_sim
)


for (i in seq_len(n_sim)) {

  theta_i <- theta_sim[i, ]

  for (j in seq_len(n_ind)) {

    log_beta <-
      theta_i["beta_pop"] +
      ranef_beta[j, i]

    beta_sim[j, i] <- exp(
      as.numeric(log_beta)
    )


    p_sim[j, i] <- exp(
      as.numeric(
        theta_i["p_pop"]
      )
    )


    log_delta <-
      theta_i["delta_pop"] +

      if (
        "beta_delta_sexe_Male" %in% names(theta_i) &&
          df_sim$sexe[j] == "Male"
      ) {
        theta_i["beta_delta_sexe_Male"]
      } else {
        0
      } +

      if (
        "beta_delta_statut_vaccin_1" %in% names(theta_i) &&
          df_sim$vaccin[j] == "1"
      ) {
        theta_i["beta_delta_statut_vaccin_1"]
      } else {
        0
      } +

      ranef_delta[j, i]

    delta_sim[j, i] <- exp(
      as.numeric(log_delta)
    )
  }
}


# ---------------------------------------------------------------------------
# 6. Viral dynamics model
# ---------------------------------------------------------------------------

ETA <- as.numeric(
  mu_full["eta_pop"]
)

C <- as.numeric(
  mu_full["c_pop"]
)

TINC <- as.numeric(
  mu_full["tinc_pop"]
)

T0 <- 4e6

TIMES <- seq(
  -2,
  30,
  by = 0.1
)


model_virus <- function(t, state, pars) {

  with(
    as.list(
      c(
        state,
        pars
      )
    ),
    {

      dT <- -beta * VI * T

      dVI <-
        eta *
        (p / c) *
        beta *
        VI *
        T -
        delta * VI

      dVNI <-
        (1 - eta) *
        (p / c) *
        beta *
        VI *
        T -
        delta * VNI

      V <- VI + VNI

      list(
        c(
          dT,
          dVI,
          dVNI
        ),
        c(V = V)
      )
    }
  )
}


# ---------------------------------------------------------------------------
# 7. Simulate trajectories in chunks
# ---------------------------------------------------------------------------

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

  message(
    "IBV simulations ",
    chunk_start,
    "-",
    chunk_end
  )


  n_sim_chunk <-
    chunk_end -
    chunk_start +
    1

  sim_list_chunk <- vector(
    "list",
    n_sim_chunk * n_ind
  )

  results_list_chunk <- vector(
    "list",
    n_sim_chunk * n_ind
  )

  idx <- 1


  for (i in chunk_start:chunk_end) {

    for (j in seq_len(n_ind)) {

      beta <- beta_sim[j, i]
      p <- p_sim[j, i]
      delta <- delta_sim[j, i]

      VI0 <- ETA * p / C

      R0 <-
        ETA *
        T0 *
        p *
        beta /
        (delta * C)


      if (R0 <= 1) {

        Vpeak_ana <- NA_real_
        Vpeak_ode_VI <- NA_real_
        Vpeak_ode_V <- NA_real_
        Ct_peak_ode <- NA_real_
        sim_i <- NULL

      } else {

        Vpeak_ana <-
          VI0 +
          (ETA * p / C) * T0 -
          (delta / beta) *
          (log(R0) + 1)


        state_ini <- c(
          T = T0,
          VI = VI0,
          VNI = 0
        )

        pars_i <- c(
          beta = beta,
          c = C,
          p = p,
          eta = ETA,
          delta = delta,
          tinc = TINC
        )


        sim_i <- tryCatch(
          as.data.frame(
            ode(
              y = state_ini,
              times = TIMES,
              func = model_virus,
              parms = pars_i,
              method = "lsoda"
            )
          ),
          error = function(e) NULL
        )


        if (is.null(sim_i)) {

          Vpeak_ode_VI <- NA_real_
          Vpeak_ode_V <- NA_real_
          Ct_peak_ode <- NA_real_

        } else {

          Vpeak_ode_VI <- max(
            sim_i$VI,
            na.rm = TRUE
          )

          Vpeak_ode_V <- max(
            sim_i$VI + sim_i$VNI,
            na.rm = TRUE
          )

          sim_i$Ct <-
            49 -
            3 *
            log10(
              pmax(
                (
                  sim_i$VI +
                  sim_i$VNI
                ) / 30,
                1e-2
              )
            )

          Ct_peak_ode <- min(
            sim_i$Ct,
            na.rm = TRUE
          )


          sim_i <- sim_i %>%
            mutate(
              sim = i,
              ind = j,
              ID = df_sim$ID[j],
              ID_infection =
                df_sim$ID_infection[j],
              vaccin =
                df_sim$vaccin[j],
              sexe =
                df_sim$sexe[j]
            )
        }
      }


      sim_list_chunk[[idx]] <- sim_i

      results_list_chunk[[idx]] <- data.frame(
        sim = i,
        ind = j,
        ID = df_sim$ID[j],
        ID_infection =
          df_sim$ID_infection[j],
        vaccin =
          df_sim$vaccin[j],
        sexe =
          df_sim$sexe[j],
        R0 = R0,
        Vpeak_ana = Vpeak_ana,
        Vpeak_ode_VI = Vpeak_ode_VI,
        Vpeak_ode_V = Vpeak_ode_V,
        Ct_peak_ode = Ct_peak_ode
      )

      idx <- idx + 1
    }
  }


  sim_list_chunk <- Filter(
    Negate(is.null),
    sim_list_chunk
  )

  df_traj_chunk <- bind_rows(
    sim_list_chunk
  )

  df_long_sim_chunk <- df_traj_chunk %>%
    mutate(
      time_since_symptoms_onset = time
    ) %>%
    select(
      sim,
      ID,
      ID_infection,
      vaccin,
      sexe,
      time_since_symptoms_onset,
      Ct,
      T,
      VI,
      VNI,
      V
    )

  df_peak_chunk <- bind_rows(
    results_list_chunk
  )


  outfile_traj <- file.path(
    outdir,
    paste0(
      "IBV_trajectories_",
      chunk_start,
      "_",
      chunk_end,
      ".rds"
    )
  )

  outfile_peak <- file.path(
    outdir,
    paste0(
      "IBV_peaks_",
      chunk_start,
      "_",
      chunk_end,
      ".rds"
    )
  )


  saveRDS(
    df_long_sim_chunk,
    outfile_traj
  )

  saveRDS(
    df_peak_chunk,
    outfile_peak
  )


  rm(
    sim_list_chunk,
    results_list_chunk,
    df_traj_chunk,
    df_long_sim_chunk,
    df_peak_chunk
  )

  gc()
}

message("IBV simulations completed.")
