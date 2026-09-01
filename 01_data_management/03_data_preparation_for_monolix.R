################################################################################################
# Prepare RELAB data for Monolix
# Purpose:
#   Create virus-specific Monolix datasets for SARS-CoV-2, influenza A,
#   influenza B, and RSV from the normalized RELAB dataset.
################################################################################################

rm(list = ls())

# ---------------------------------------------------------------------------
# 0. Libraries
# ---------------------------------------------------------------------------

library(tidyverse)


# ---------------------------------------------------------------------------
# 1. File paths
# ---------------------------------------------------------------------------

input_file <- "data/processed/RELAB_24_26_clean_start_norm.csv"

output_dir <- "02_monolix/data"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ---------------------------------------------------------------------------
# 2. Load normalized RELAB dataset
# ---------------------------------------------------------------------------

data <- read_csv(
  input_file,
  show_col_types = FALSE
)


# ---------------------------------------------------------------------------
# 3. Epidemic-period functions
# ---------------------------------------------------------------------------

get_epidemic_period_covid <- function(season, week) {
  case_when(
    season == "24_25" & week %in% 35:42 ~ "Phase 1",
    season == "24_25" & week %in% 43:48 ~ "Phase 2",
    season == "24_25" & week %in% c(49:52, 1:22) ~ "Phase 3",

    season == "25_26" & week %in% 35:39 ~ "Phase 1",
    season == "25_26" & week %in% 40:43 ~ "Phase 2",
    season == "25_26" & week %in% c(44:52, 1:22) ~ "Phase 3",

    TRUE ~ NA_character_
  )
}


get_epidemic_period_iav <- function(season, week) {
  case_when(
    season == "24_25" & week %in% c(35:52, 1) ~ "Phase 1",
    season == "24_25" & week %in% 2:4 ~ "Phase 2",
    season == "24_25" & week %in% 5:25 ~ "Phase 3",

    season == "25_26" & week %in% 35:51 ~ "Phase 1",
    season == "25_26" & week %in% c(52, 1:3) ~ "Phase 2",
    season == "25_26" & week %in% 4:25 ~ "Phase 3",

    TRUE ~ NA_character_
  )
}


get_epidemic_period_ibv <- function(season, week) {
  case_when(
    season == "24_25" & week %in% c(35:52, 1:3) ~ "Phase 1",
    season == "24_25" & week %in% 4:6 ~ "Phase 2",
    season == "24_25" & week %in% 7:25 ~ "Phase 3",

    season == "25_26" & week %in% c(35:52, 1:4) ~ "Phase 1",
    season == "25_26" & week %in% 5:6 ~ "Phase 2",
    season == "25_26" & week %in% 7:25 ~ "Phase 3",

    TRUE ~ NA_character_
  )
}


get_epidemic_period_rsv <- function(season, week) {
  case_when(
    season == "24_25" & week %in% 35:51 ~ "Phase 1",
    season == "24_25" & week %in% c(52, 1:2) ~ "Phase 2",
    season == "24_25" & week %in% 3:25 ~ "Phase 3",

    season == "25_26" & week %in% 35:50 ~ "Phase 1",
    season == "25_26" & week %in% c(51:52, 1:2) ~ "Phase 2",
    season == "25_26" & week %in% 3:25 ~ "Phase 3",

    TRUE ~ NA_character_
  )
}


# ---------------------------------------------------------------------------
# 4. Generic Monolix preparation function
# ---------------------------------------------------------------------------

prepare_monolix_data <- function(
    data,
    virus,
    infection_pattern,
    pcr_result,
    ct_normalized,
    epidemic_period_fun,
    vaccination_variable = NULL,
    influenza_type = NULL,
    obs_scale = c("Ct", "log10_proxy"),
    include_symptoms = FALSE) {

  obs_scale <- match.arg(obs_scale)

  # -------------------------------------------------------------------------
  # 4.1 Select infection episodes
  # -------------------------------------------------------------------------

  data_virus <- data %>%
    filter(
      str_detect(
        coalesce(ID_infection, ""),
        regex(infection_pattern, ignore_case = TRUE)
      )
    )

  # For influenza, retain episodes containing the requested influenza type.
  if (!is.null(influenza_type)) {
    data_virus <- data_virus %>%
      group_by(ID_infection) %>%
      filter(
        any(
          PCR_grippe_typage_MultiName == influenza_type,
          na.rm = TRUE
        )
      ) %>%
      ungroup()
  }


  # -------------------------------------------------------------------------
  # 4.2 Define normalized Ct
  # -------------------------------------------------------------------------

  data_virus <- data_virus %>%
    mutate(
      Ct_norm = if_else(
        .data[[pcr_result]] == "Pos",
        .data[[ct_normalized]],
        40
      )
    )


  # -------------------------------------------------------------------------
  # 4.3 Symptoms
  # -------------------------------------------------------------------------

  if (include_symptoms) {
    data_virus <- data_virus %>%
      group_by(ID_infection) %>%
      mutate(
        symptoms = case_when(
          any(
            fievre == "Yes" &
              signes_respiratoires == "Yes",
            na.rm = TRUE
          ) ~ "Both",

          any(
            fievre == "Yes",
            na.rm = TRUE
          ) ~ "Fever only",

          any(
            signes_respiratoires == "Yes",
            na.rm = TRUE
          ) ~ "Respiratory signs only",

          any(
            fievre == "No" &
              signes_respiratoires == "No",
            na.rm = TRUE
          ) ~ "None",

          TRUE ~ "Unknown"
        )
      ) %>%
      ungroup()
  }


  # -------------------------------------------------------------------------
  # 4.4 Prepare Monolix variables
  # -------------------------------------------------------------------------

  data_monolix <- data_virus %>%
    mutate(
      time = time_since_symptoms_onset,

      # Keep the observation definitions used in the original virus-specific
      # scripts.
      obs = case_when(
        obs_scale == "Ct" ~ Ct_norm,
        obs_scale == "log10_proxy" ~ (49 - Ct_norm) / 3
      ),

      censor = case_when(
        .data[[pcr_result]] == "Neg" ~ 1,
        .data[[pcr_result]] == "Pos" ~ 0,
        TRUE ~ NA_real_
      ),

      Ct_diff = 50 - Ct_norm,

      tag_inf = 1,

      tSS = 0,

      epidemic_period = epidemic_period_fun(
        epidemic_season,
        Periode.semaine
      ),

      season = epidemic_season
    )


  # -------------------------------------------------------------------------
  # 4.5 Vaccination status
  # -------------------------------------------------------------------------

  if (!is.null(vaccination_variable)) {
    data_monolix <- data_monolix %>%
      mutate(
        statut_vaccin = .data[[vaccination_variable]]
      ) %>%
      filter(
        !is.na(statut_vaccin)
      )
  }


  # -------------------------------------------------------------------------
  # 4.6 Keep valid observations
  # -------------------------------------------------------------------------

  data_monolix <- data_monolix %>%
    filter(
      !is.na(ID_infection),
      !is.na(Ct_diff),
      !is.na(obs)
    )


  # -------------------------------------------------------------------------
  # 4.7 Select output columns
  # -------------------------------------------------------------------------

  common_columns <- c(
    "ID_infection",
    "time",
    "obs",
    "censor",
    "Ct_diff",
    "tag_inf",
    "sexe",
    "age_cat",
    "season",
    "epidemic_period"
  )

  if (!is.null(vaccination_variable)) {
    common_columns <- c(
      common_columns,
      "statut_vaccin"
    )
  }

  if (include_symptoms) {
    common_columns <- c(
      common_columns,
      "symptoms"
    )
  }

  common_columns <- c(
    common_columns,
    "tSS"
  )

  data_monolix <- data_monolix %>%
    select(
      all_of(common_columns)
    )


  # -------------------------------------------------------------------------
  # 4.8 Summary
  # -------------------------------------------------------------------------

  message(
    virus,
    ": ",
    n_distinct(data_monolix$ID_infection),
    " infection episodes; ",
    nrow(data_monolix),
    " observations."
  )

  return(data_monolix)
}


# ---------------------------------------------------------------------------
# 5. SARS-CoV-2
# ---------------------------------------------------------------------------

data_monolix_covid <- prepare_monolix_data(
  data = data,
  virus = "SARS-CoV-2",
  infection_pattern = "Cov",
  pcr_result = "PCR_covid",
  ct_normalized = "Ct_covid_norm",
  epidemic_period_fun = get_epidemic_period_covid,
  vaccination_variable = "vaccin_covid_statut",
  obs_scale = "Ct",
  include_symptoms = TRUE
)


# ---------------------------------------------------------------------------
# 6. Influenza A
# ---------------------------------------------------------------------------

data_monolix_iav <- prepare_monolix_data(
  data = data,
  virus = "Influenza A",
  infection_pattern = "Inf",
  pcr_result = "PCR_grippe",
  ct_normalized = "Ct_grippe_norm",
  epidemic_period_fun = get_epidemic_period_iav,
  vaccination_variable = "vaccin_grippe_statut",
  influenza_type = "Grippe A",
  obs_scale = "log10_proxy"
)


# ---------------------------------------------------------------------------
# 7. Influenza B
# ---------------------------------------------------------------------------

data_monolix_ibv <- prepare_monolix_data(
  data = data,
  virus = "Influenza B",
  infection_pattern = "Inf",
  pcr_result = "PCR_grippe",
  ct_normalized = "Ct_grippe_norm",
  epidemic_period_fun = get_epidemic_period_ibv,
  vaccination_variable = "vaccin_grippe_statut",
  influenza_type = "Grippe B",
  obs_scale = "log10_proxy"
)


# ---------------------------------------------------------------------------
# 8. RSV
# ---------------------------------------------------------------------------

data_monolix_rsv <- prepare_monolix_data(
  data = data,
  virus = "RSV",
  infection_pattern = "RSV",
  pcr_result = "PCR_VRS",
  ct_normalized = "Ct_VRS_norm",
  epidemic_period_fun = get_epidemic_period_rsv,
  vaccination_variable = NULL,
  obs_scale = "Ct"
)


# ---------------------------------------------------------------------------
# 9. Save Monolix datasets
# ---------------------------------------------------------------------------

write_csv(
  data_monolix_covid,
  file.path(
    output_dir,
    "data_for_monolix_SARS_CoV_2.csv"
  )
)

write_csv(
  data_monolix_iav,
  file.path(
    output_dir,
    "data_for_monolix_IAV.csv"
  )
)

write_csv(
  data_monolix_ibv,
  file.path(
    output_dir,
    "data_for_monolix_IBV.csv"
  )
)

write_csv(
  data_monolix_rsv,
  file.path(
    output_dir,
    "data_for_monolix_RSV.csv"
  )
)

message("Monolix datasets saved to: ", output_dir)
