################################################################################################
# RELAB data management
# Author: Laura MULAS
# Purpose:
#   Clean and harmonize RELAB data for the 2024/25 and 2025/26 epidemic seasons.
#
#   The same pipeline is applied to both seasons, with season-specific settings
#   for the study period and symptom-onset information.
################################################################################################

rm(list = ls())

# ---------------------------------------------------------------------------
# 0. Libraries
# ---------------------------------------------------------------------------

library(tidyverse)
library(lubridate)


# ---------------------------------------------------------------------------
# 1. Season configuration
# ---------------------------------------------------------------------------

season_config <- tibble(
  season = c("24_25", "25_26"),
  input_file = c(
    "data/raw/Export_data_RELAB_24_25.csv",
    "data/raw/Export_data_RELAB_25_26.csv"
  ),
  start_date = as.Date(c(
    "2024-09-01",
    "2025-09-01"
  )),
  end_date = as.Date(c(
    "2025-05-31",
    "2026-05-31"
  )),
  use_fever_tss = c(
    FALSE,
    TRUE
  ),
  output_file = c(
    "data/processed/RELAB_24_25_clean_start.csv",
    "data/processed/RELAB_25_26_clean_start.csv"
  )
)


# ---------------------------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------------------------

# Find a representative symptom-onset date from repeated interval-censored
# observations within the same infection episode.
get_tss_episode <- function(lower, upper) {

  ok <- !is.na(lower) & !is.na(upper)

  lower <- lower[ok]
  upper <- upper[ok]

  if (length(lower) == 0) {
    return(as.Date(NA))
  }

  # Intersection shared by all intervals.
  common_lower <- max(lower)
  common_upper <- min(upper)

  if (common_lower <= common_upper) {
    return(
      as.Date(
        (as.numeric(common_lower) + as.numeric(common_upper)) / 2,
        origin = "1970-01-01"
      )
    )
  }

  # If no global intersection exists, find the largest subset of compatible
  # intervals.
  n <- length(lower)

  if (n >= 2) {
    for (k in n:2) {

      combinations <- combn(
        seq_len(n),
        k,
        simplify = FALSE
      )

      valid_combinations <- lapply(
        combinations,
        function(idx) {

          lower_intersection <- max(lower[idx])
          upper_intersection <- min(upper[idx])

          if (lower_intersection <= upper_intersection) {
            tibble(
              lower = lower_intersection,
              upper = upper_intersection
            )
          } else {
            NULL
          }
        }
      ) %>%
        bind_rows()

      if (nrow(valid_combinations) > 0) {

        return(
          valid_combinations %>%
            mutate(
              midpoint = as.Date(
                (as.numeric(lower) + as.numeric(upper)) / 2,
                origin = "1970-01-01"
              )
            ) %>%
            arrange(midpoint) %>%
            pull(midpoint) %>%
            first()
        )
      }
    }
  }

  # If none of the intervals overlap, retain the earliest possible
  # symptom-onset date.
  min(lower)
}


# Propagate the last known vaccination response after vaccination has first
# been reported for an individual.
propagate_vaccination <- function(vaccine, sampling_date) {

  valid <- vaccine %in% 3:6

  if (!any(valid, na.rm = TRUE)) {
    return(vaccine)
  }

  first_vacc_date <- min(
    sampling_date[valid],
    na.rm = TRUE
  )

  valid_indices <- which(valid)

  last_index <- valid_indices[
    which.max(sampling_date[valid_indices])
  ]

  last_valid_vaccine <- vaccine[last_index]

  replace <- (
    is.na(vaccine) |
      vaccine %in% c(1, 2)
  ) &
    sampling_date > first_vacc_date

  vaccine[replace] <- last_valid_vaccine

  vaccine
}


# ---------------------------------------------------------------------------
# 3. Main cleaning function
# ---------------------------------------------------------------------------

prepare_season <- function(
    input_file,
    season,
    start_date,
    end_date,
    use_fever_tss = FALSE,
    output_file = NULL) {

  message("Preparing RELAB season ", season, "...")

  # -------------------------------------------------------------------------
  # 3.1 Import and basic cleaning
  # -------------------------------------------------------------------------

  df <- read_delim(
    input_file,
    delim = ";",
    show_col_types = FALSE,
    trim_ws = TRUE
  ) %>%
    distinct() %>%
    mutate(
      ID = match(
        patient_ID,
        unique(patient_ID)
      ),
      date_prelevement = dmy(date_prelevement)
    ) %>%
    filter(
      date_prelevement >= start_date,
      date_prelevement <= end_date
    ) %>%
    mutate(
      year = year(date_prelevement),
      week = isoweek(date_prelevement),
      year = if_else(
        week == 53,
        year + 1L,
        year
      ),
      week = if_else(
        week == 53,
        1L,
        week
      )
    ) %>%
    relocate(
      ID,
      date_prelevement,
      .before = Labo_nom
    ) %>%
    arrange(
      ID,
      date_prelevement
    )


  # -------------------------------------------------------------------------
  # 3.2 Age and sex
  # -------------------------------------------------------------------------

  df <- df %>%
    mutate(
      annee_naissance = case_when(
        annee_naissance %in% 1900:2026 ~ annee_naissance,
        TRUE ~ mois_naissance
      )
    ) %>%
    filter(
      !is.na(annee_naissance),
      annee_naissance > 1900
    ) %>%
    mutate(
      age = as.numeric(Periode.annee) -
        as.numeric(annee_naissance),

      sexe = case_when(
        sexe == "F" ~ "Female",
        sexe == "M" ~ "Male",
        TRUE ~ NA_character_
      ),

      age_cat = case_when(
        age <= 5 ~ "<5",
        age > 5 & age <= 18 ~ "5-18",
        age > 18 & age <= 65 ~ "18-65",
        age > 65 ~ ">65",
        TRUE ~ NA_character_
      ),

      age_cat = factor(
        age_cat,
        levels = c(
          "<5",
          "5-18",
          "18-65",
          ">65"
        )
      )
    ) %>%
    filter(
      !is.na(age),
      !is.na(sexe)
    ) %>%
    relocate(
      age,
      sexe,
      .after = date_prelevement
    )


  # -------------------------------------------------------------------------
  # 3.3 PCR results and Ct values
  # -------------------------------------------------------------------------

  valid_pcr_values <- c("Pos", "Neg")

  df <- df %>%
    filter(
      PCR_covid %in% valid_pcr_values |
        PCR_grippe %in% valid_pcr_values |
        PCR_VRS %in% valid_pcr_values
    ) %>%
    mutate(
      across(
        c(
          PCR_covid_Ct,
          PCR_grippe_Ct,
          PCR_VRS_Ct
        ),
        ~ parse_number(
          as.character(.x),
          locale = locale(decimal_mark = ",")
        )
      ),

      PCR_covid_Ct = if_else(
        PCR_covid_Ct <= 10 |
          PCR_covid_Ct > 40,
        NA_real_,
        PCR_covid_Ct
      ),

      PCR_grippe_Ct = if_else(
        PCR_grippe_Ct <= 10 |
          PCR_grippe_Ct > 40,
        NA_real_,
        PCR_grippe_Ct
      ),

      PCR_VRS_Ct = if_else(
        PCR_VRS_Ct <= 10 |
          PCR_VRS_Ct > 40,
        NA_real_,
        PCR_VRS_Ct
      )
    ) %>%
    filter(
      (
        PCR_covid == "Neg" &
          PCR_grippe == "Neg" &
          PCR_VRS == "Neg"
      ) |
        (
          PCR_covid == "Pos" &
            !is.na(PCR_covid_Ct)
        ) |
        (
          PCR_grippe == "Pos" &
            !is.na(PCR_grippe_Ct)
        ) |
        (
          PCR_VRS == "Pos" &
            !is.na(PCR_VRS_Ct)
        )
    )


  # -------------------------------------------------------------------------
  # 3.4 Symptoms
  # -------------------------------------------------------------------------

  df <- df %>%
    mutate(
      fievre = if_else(
        fievre %in% c("O", "N"),
        fievre,
        NA_character_
      ),

      signes_respiratoires = if_else(
        signes_respiratoires %in% c("O", "N"),
        signes_respiratoires,
        NA_character_
      ),

      symptom_respi_onset = case_when(
        debut_signe_respi == 3 ~ "0-1",
        debut_signe_respi == 4 ~ "2-5",
        debut_signe_respi == 5 ~ "6-10",
        TRUE ~ NA_character_
      ),

      fievre = case_when(
        fievre == "O" ~ "Yes",
        fievre == "N" ~ "No",
        TRUE ~ NA_character_
      ),

      signes_respiratoires = case_when(
        signes_respiratoires == "O" ~ "Yes",
        signes_respiratoires == "N" ~ "No",
        TRUE ~ NA_character_
      )
    )


  # Fever-onset categories were available/used in 2025/26.
  if (use_fever_tss) {

    df <- df %>%
      mutate(
        symptom_fievre_onset = case_when(
          debut_signe_respi_fievre == 4 ~ "0-1",
          debut_signe_respi_fievre == 5 ~ "2-5",
          debut_signe_respi_fievre == 6 ~ "6-10",
          TRUE ~ NA_character_
        ),

        fievre = case_when(
          debut_signe_respi_fievre == 2 ~ "No",
          debut_signe_respi_fievre %in% c(4, 5, 6) ~ "Yes",
          TRUE ~ fievre
        )
      )
  }


  # -------------------------------------------------------------------------
  # 3.5 Healthcare establishment
  # -------------------------------------------------------------------------

  df <- df %>%
    mutate(
      etablissement_soins = if_else(
        etablissement_soins %in% c("O", "E"),
        "Yes",
        "No"
      )
    )


  # -------------------------------------------------------------------------
  # 3.6 Vaccination
  # -------------------------------------------------------------------------

  df <- df %>%
    group_by(ID) %>%
    arrange(
      date_prelevement,
      .by_group = TRUE
    ) %>%
    mutate(
      vaccin_covid = propagate_vaccination(
        vaccin_covid,
        date_prelevement
      ),

      vaccin_grippe = propagate_vaccination(
        vaccin_grippe,
        date_prelevement
      )
    ) %>%
    ungroup() %>%
    mutate(
      vaccin_covid_statut = case_when(
        vaccin_covid %in% c(4, 5) ~ 1,
        vaccin_covid %in% c(3, 6, 7) ~ 0,
        TRUE ~ NA_real_
      ),

      vaccin_grippe_statut = case_when(
        vaccin_grippe %in% c(4, 5) ~ 1,
        vaccin_grippe %in% c(3, 6, 7) ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(
      (
        PCR_covid == "Pos" &
          !is.na(vaccin_covid_statut)
      ) |
        (
          PCR_grippe == "Pos" &
            !is.na(vaccin_grippe_statut)
        ) |
        PCR_VRS == "Pos" |
        (
          PCR_covid != "Pos" &
            PCR_grippe != "Pos" &
            PCR_VRS != "Pos" &
            (
              !is.na(vaccin_grippe_statut) |
                !is.na(vaccin_covid_statut)
            )
        )
    )


  # -------------------------------------------------------------------------
  # 3.7 Define infection episodes
  # -------------------------------------------------------------------------

  df <- df %>%
    arrange(
      ID,
      date_prelevement
    ) %>%
    group_by(ID) %>%
    mutate(
      days_since_prev = as.numeric(
        date_prelevement -
          lag(
            date_prelevement,
            default = first(date_prelevement)
          )
      ),

      new_episode = if_else(
        row_number() == 1 |
          days_since_prev > 21,
        1L,
        0L
      ),

      ID_paire = cumsum(new_episode)
    ) %>%
    ungroup() %>%
    select(
      -days_since_prev,
      -new_episode
    )


  df <- df %>%
    arrange(
      ID,
      ID_paire,
      date_prelevement
    ) %>%
    group_by(
      ID,
      ID_paire
    ) %>%
    mutate(
      viruses_detected = pmap_chr(
        list(
          PCR_covid,
          PCR_grippe,
          PCR_VRS
        ),
        function(covid, influenza, rsv) {

          viruses <- c(
            if (covid == "Pos") "Covid" else NA_character_,
            if (influenza == "Pos") "Influenza" else NA_character_,
            if (rsv == "Pos") "RSV" else NA_character_
          )

          viruses <- viruses[
            !is.na(viruses)
          ]

          if (length(viruses) == 0) {
            "None"
          } else {
            paste(
              viruses,
              collapse = "_"
            )
          }
        }
      ),

      is_pos = viruses_detected != "None",

      new_infection =
        is_pos &
        (
          is.na(
            lag(date_prelevement)
          ) |
            as.numeric(
              date_prelevement -
                lag(date_prelevement)
            ) >= 21
        ),

      infection_episode = cumsum(
        replace_na(
          new_infection,
          FALSE
        )
      ),

      days_since_prev = as.numeric(
        date_prelevement -
          lag(date_prelevement)
      ),

      days_to_next = as.numeric(
        lead(date_prelevement) -
          date_prelevement
      ),

      attach_negative =
        !is_pos &
        (
          (
            !is.na(days_since_prev) &
              days_since_prev < 21
          ) |
            (
              !is.na(days_to_next) &
                days_to_next <= 5
            )
        ),

      episode_temp = if_else(
        is_pos |
          attach_negative,
        infection_episode,
        NA_integer_
      ),

      virus_for_episode = if_else(
        is_pos,
        viruses_detected,
        NA_character_
      )
    ) %>%
    fill(
      episode_temp,
      virus_for_episode,
      .direction = "down"
    ) %>%
    mutate(
      ID_infection = if_else(
        is.na(episode_temp) |
          is.na(virus_for_episode),

        paste(
          season,
          ID,
          ID_paire,
          "none",
          sep = "_"
        ),

        paste(
          season,
          ID,
          ID_paire,
          episode_temp,
          virus_for_episode,
          sep = "_"
        )
      )
    ) %>%
    ungroup() %>%
    select(
      -is_pos,
      -new_infection,
      -days_since_prev,
      -days_to_next,
      -attach_negative,
      -episode_temp,
      -virus_for_episode
    )


  # -------------------------------------------------------------------------
  # 3.8 Separate influenza A/B episodes when necessary
  # -------------------------------------------------------------------------

  df <- df %>%
    group_by(ID_infection) %>%
    mutate(
      n_typages_grippe = n_distinct(
        PCR_grippe_typage_MultiName[
          !is.na(PCR_grippe_typage_MultiName)
        ]
      ),

      split_infection_grippe =
        str_detect(
          ID_infection,
          "Influenza"
        ) &
        n_typages_grippe >= 2
    ) %>%
    ungroup() %>%
    mutate(
      ID_infection = if_else(
        split_infection_grippe &
          !is.na(PCR_grippe_typage_MultiName),

        paste(
          sub(
            "_Influenza$",
            "",
            ID_infection
          ),
          PCR_grippe_typage_MultiName,
          "Influenza",
          sep = "_"
        ),

        ID_infection
      )
    ) %>%
    select(
      -n_typages_grippe,
      -split_infection_grippe
    )


  # Remove duplicate records corresponding to the same medical record/date.
  df <- df %>%
    group_by(
      ID_infection,
      date_dossier
    ) %>%
    arrange(
      date_prelevement,
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup()


  # -------------------------------------------------------------------------
  # 3.9 Time since symptom onset
  # -------------------------------------------------------------------------

  df <- df %>%
    mutate(
      lower_TSS_respi = case_when(
        debut_signe_respi == 3 ~ date_prelevement - 1,
        debut_signe_respi == 4 ~ date_prelevement - 5,
        debut_signe_respi == 5 ~ date_prelevement - 10,
        TRUE ~ as.Date(NA)
      ),

      upper_TSS_respi = case_when(
        debut_signe_respi == 3 ~ date_prelevement,
        debut_signe_respi == 4 ~ date_prelevement - 2,
        debut_signe_respi == 5 ~ date_prelevement - 6,
        TRUE ~ as.Date(NA)
      )
    ) %>%
    group_by(
      ID,
      ID_infection
    ) %>%
    mutate(
      time_of_symptoms_respi = get_tss_episode(
        lower_TSS_respi,
        upper_TSS_respi
      )
    ) %>%
    ungroup()


  if (use_fever_tss) {

    df <- df %>%
      mutate(
        lower_TSS_fievre = case_when(
          debut_signe_respi_fievre == 4 ~ date_prelevement - 1,
          debut_signe_respi_fievre == 5 ~ date_prelevement - 5,
          debut_signe_respi_fievre == 6 ~ date_prelevement - 10,
          TRUE ~ as.Date(NA)
        ),

        upper_TSS_fievre = case_when(
          debut_signe_respi_fievre == 4 ~ date_prelevement,
          debut_signe_respi_fievre == 5 ~ date_prelevement - 2,
          debut_signe_respi_fievre == 6 ~ date_prelevement - 6,
          TRUE ~ as.Date(NA)
        )
      ) %>%
      group_by(
        ID,
        ID_infection
      ) %>%
      mutate(
        time_of_symptoms_fievre = get_tss_episode(
          lower_TSS_fievre,
          upper_TSS_fievre
        )
      ) %>%
      ungroup() %>%
      mutate(
        time_since_symptoms_onset_respi =
          as.numeric(
            date_prelevement -
              time_of_symptoms_respi
          ),

        time_since_symptoms_onset_fievre =
          as.numeric(
            date_prelevement -
              time_of_symptoms_fievre
          ),

        # Respiratory symptom onset is used preferentially; fever onset is
        # used when respiratory symptom onset is unavailable.
        time_since_symptoms_onset = coalesce(
          time_since_symptoms_onset_respi,
          time_since_symptoms_onset_fievre
        )
      ) %>%
      filter(
        !is.na(
          time_since_symptoms_onset
        )
      ) %>%
      select(
        -lower_TSS_respi,
        -upper_TSS_respi,
        -lower_TSS_fievre,
        -upper_TSS_fievre
      )

  } else {

    df <- df %>%
      mutate(
        time_since_symptoms_onset_respi =
          as.numeric(
            date_prelevement -
              time_of_symptoms_respi
          ),

        time_since_symptoms_onset =
          time_since_symptoms_onset_respi
      ) %>%
      filter(
        !is.na(
          time_since_symptoms_onset
        )
      ) %>%
      select(
        -lower_TSS_respi,
        -upper_TSS_respi
      )
  }


  # Number of observations available for each individual.
  df <- df %>%
    group_by(ID) %>%
    mutate(
      n_tests = n()
    ) %>%
    ungroup()


  # -------------------------------------------------------------------------
  # 3.10 Harmonize influenza typing
  # -------------------------------------------------------------------------

  df <- df %>%
    mutate(
      variant_grippe = case_when(
        !is.na(FLU_SeqClade_Lyon) &
          FLU_SeqClade_Lyon != "" ~ FLU_SeqClade_Lyon,
        TRUE ~ FLU_clade
      ),

      grippe_subtype = case_when(
        str_detect(
          coalesce(variant_grippe, ""),
          "5a\\.2a"
        ) ~ "H1N1",

        str_detect(
          coalesce(variant_grippe, ""),
          "2a\\.3a\\.1"
        ) ~ "H3N2",

        str_detect(
          coalesce(variant_grippe, ""),
          "V1A\\.3a\\.2"
        ) ~ "B/Victoria",

        TRUE ~ "Others"
      ),

      PCR_grippe_typage_MultiName = case_when(
        PCR_grippe_typage_MultiName %in%
          c(
            "A",
            "GRA",
            "Grippe A"
          ) ~ "Grippe A",

        PCR_grippe_typage_MultiName %in%
          c(
            "B",
            "GRB",
            "Grippe B"
          ) ~ "Grippe B",

        grippe_subtype %in%
          c(
            "H1N1",
            "H3N2"
          ) ~ "Grippe A",

        grippe_subtype ==
          "B/Victoria" ~ "Grippe B",

        TRUE ~ NA_character_
      )
    ) %>%
    select(
      -variant_grippe
    ) %>%
    group_by(ID_infection) %>%
    mutate(
      PCR_grippe_typage_MultiName = if_else(
        str_detect(
          ID_infection,
          regex(
            "influenza",
            ignore_case = TRUE
          )
        ) &
          is.na(
            PCR_grippe_typage_MultiName
          ) &
          any(
            !is.na(
              PCR_grippe_typage_MultiName
            )
          ),

        first(
          PCR_grippe_typage_MultiName[
            !is.na(
              PCR_grippe_typage_MultiName
            )
          ]
        ),

        PCR_grippe_typage_MultiName
      )
    ) %>%
    ungroup() %>%
    filter(
      !(
        PCR_grippe == "Pos" &
          is.na(
            PCR_grippe_typage_MultiName
          ) &
          PCR_covid == "Neg" &
          PCR_VRS == "Neg"
      )
    )


  # -------------------------------------------------------------------------
  # 3.11 Harmonize PCR technique labels
  # -------------------------------------------------------------------------

  eurobio_names <- c(
    "Eurobio - Réactif Allplex (Starlet/CFX) ",
    "Eurobio - Réactif Allplex (Starlet/CFX)",
    " Eurobio EBX 42"
  )

  df <- df %>%
    mutate(
      PCR_covid_technique_update = if_else(
        PCR_covid_technique %in%
          eurobio_names,
        "Eurobio",
        PCR_covid_technique
      ),

      PCR_grippe_technique_update = if_else(
        PCR_grippe_technique %in%
          eurobio_names,
        "Eurobio",
        PCR_grippe_technique
      ),

      PCR_VRS_technique_update = if_else(
        PCR_VRS_technique %in%
          eurobio_names,
        "Eurobio",
        PCR_VRS_technique
      )
    )


  # -------------------------------------------------------------------------
  # 3.12 Add epidemic-season label
  # -------------------------------------------------------------------------

  df <- df %>%
    mutate(
      epidemic_season = season
    )


  # -------------------------------------------------------------------------
  # 3.13 Summary and save
  # -------------------------------------------------------------------------

  message(
    season,
    ": ",
    n_distinct(df$patient_ID),
    " individuals; ",
    n_distinct(df$ID_infection),
    " infection episodes; ",
    nrow(df),
    " observations."
  )

  if (!is.null(output_file)) {

    dir.create(
      dirname(output_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    write_csv(
      df,
      output_file
    )

    message(
      "Saved to: ",
      output_file
    )
  }

  return(df)
}


# ---------------------------------------------------------------------------
# 4. Run both epidemic seasons
# ---------------------------------------------------------------------------

relab_24_25 <- prepare_season(
  input_file = season_config$input_file[
    season_config$season == "24_25"
  ],
  season = "24_25",
  start_date = season_config$start_date[
    season_config$season == "24_25"
  ],
  end_date = season_config$end_date[
    season_config$season == "24_25"
  ],
  use_fever_tss = FALSE,
  output_file = season_config$output_file[
    season_config$season == "24_25"
  ]
)


relab_25_26 <- prepare_season(
  input_file = season_config$input_file[
    season_config$season == "25_26"
  ],
  season = "25_26",
  start_date = season_config$start_date[
    season_config$season == "25_26"
  ],
  end_date = season_config$end_date[
    season_config$season == "25_26"
  ],
  use_fever_tss = TRUE,
  output_file = season_config$output_file[
    season_config$season == "25_26"
  ]
)


# ---------------------------------------------------------------------------
# 5. Combine seasons
# ---------------------------------------------------------------------------

relab_24_26 <- bind_rows(
  relab_24_25,
  relab_25_26
)

combined_output <- "data/processed/RELAB_24_26_clean_start.csv"

write_csv(
  relab_24_26,
  combined_output
)

message(
  "Combined dataset saved to: ",
  combined_output
)
