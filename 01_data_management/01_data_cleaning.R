# ==============================================================================
# RELAB data preparation
# ==============================================================================
# Purpose:
#   Clean and prepare RELAB community RT-PCR data for downstream analyses.
#
# Reproducibility:
#   - Raw individual-level data are NOT distributed with this repository.
#   - Paths are relative to the repository root.
#   - Place the authorised raw export in `data/raw/` before running this script.
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)

# ------------------------------------------------------------------------------
# File paths and analysis settings
# ------------------------------------------------------------------------------
input_file <- file.path(
  "data", "raw", "Export_data_RELAB_2025_26_week_09.csv"
)

output_file <- file.path(
  "data", "processed", "RELAB_25_26_clean_start.csv"
)

start_date <- as.Date("2025-09-01")
end_date   <- as.Date("2026-05-31")

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------
if (!file.exists(input_file)) {
  stop(
    "Raw RELAB file not found: ", input_file,
    "\nPlace the authorised input file in `data/raw/` or update `input_file`."
  )
}

df <- read.csv(
  input_file,
  sep = ";",
  stringsAsFactors = FALSE,
  row.names = NULL
)


df %>%
  summarise(nb_individus = n_distinct(patient_ID))


# ---------------------------------------------------------------------------
# Remove duplicates
# ---------------------------------------------------------------------------

df_no_dup <- df %>%
  distinct()

df <- df_no_dup

rm(df_no_dup)

# Number of individuals
df %>%
  summarise(nb_individus = n_distinct(patient_ID))

# ---------------------------------------------------------------------------
# ID
# ---------------------------------------------------------------------------

# Create an anonymised sequential analysis ID.
# `patient_ID` remains in memory here because it is required for linkage, but raw
# identifiers should never be committed to the public repository.
df <- df %>%
  mutate(ID = match(patient_ID, unique(patient_ID))) %>%
  relocate(ID)

# ---------------------------------------------------------------------------
# Sampling date
# ---------------------------------------------------------------------------

df_date <- df %>%
  mutate(date_prelevement = as.Date(date_prelevement, format = "%d/%m/%Y"))

summary(df_date$date_prelevement)

# Filter the dataset to the predefined study period.
df_date <- df_date %>%
  filter(date_prelevement >= start_date & date_prelevement <= end_date) |>
  
  # Adjust weeks that are 53 -> not an entire week
  mutate(
    year = year(date_prelevement),
    week = isoweek(date_prelevement),  
    year = if_else(week == 53, year + 1L, year),
    week = if_else(week == 53, 1L, week)
  )

# Reorganize columns and sort rows
df <- df_date %>%
  relocate(c(ID, date_prelevement), .before = "Labo_nom") %>%
  arrange(ID, date_prelevement)

# Clean up intermediate objects
rm(df_date)

df %>%
  summarise(nb_individus = n_distinct(patient_ID))

# ---------------------------------------------------------------------------
# Age and sex
# ---------------------------------------------------------------------------

df_age <- df |>
  mutate(
    annee_naissance = case_when(
      annee_naissance %in% c(1900:2026) ~ annee_naissance,
      TRUE ~ mois_naissance
    )
  ) |>
  filter(!is.na(annee_naissance), annee_naissance > 1900, Periode.annee %in% c(2025,2026))

df_age <- df_age |>
  mutate(age = as.numeric(`Periode.annee`) - as.numeric(annee_naissance),
         sexe = ifelse(sexe %in% c("F", "M"), sexe, NA),
         sexe = as.factor(ifelse(sexe == "F", "Female","Male"))) |>
  filter(!(is.na(age)), !(is.na(sexe)))

summary(as.factor(df_age$sexe))

df <- df_age |>
  relocate(c(age, sexe), .after = date_prelevement) |>
  as.data.frame()

rm(df_age)

df %>%
  summarise(nb_individus = n_distinct(patient_ID))

# Creation of age categories
# NOTE: The cut-points below reproduce the original analysis exactly.
#       In particular, age == 5 is assigned to "<5" and age == 18 to "5-18".
#       Confirm these boundary conventions before publication if the manuscript
#       describes the groups as <5, 5-18, 18-65, and >65.
df_age_cat <- df |> 
  mutate(age_cat = case_when(
    age <= 5 ~ "<5",
    age > 5 & age <= 18 ~ "5-18",
    age > 18 & age <= 65 ~ "18-65",
    age > 65 ~ ">65"
  )) |>
  as.data.frame()

df_age_cat$age_cat <- factor(df_age_cat$age_cat, levels = c("<5", "5-18", "18-65", ">65"))

table(df_age_cat$age_cat, useNA = "always")

df <- df_age_cat 

rm(df_age_cat)

# ---------------------------------------------------------------------------
# Remove rows with uninterpretable PCR results
# ---------------------------------------------------------------------------

# Define invalid values
valid_values <- c("Pos", "Neg")

df_valid_result <- df %>%
  # Keep rows where at least one PCR column is valid
  filter(
    (PCR_covid %in% valid_values |
        PCR_grippe %in% valid_values |
        PCR_VRS %in% valid_values)
  )

df <- df_valid_result

df %>%
  summarise(nb_individus = n_distinct(patient_ID))

rm(df_valid_result, valid_values)

# ---------------------------------------------------------------------------
# Ct values
# ---------------------------------------------------------------------------
# Monitor Ct values to keep only values between 0 and 40
df_Ct <- df %>%
  mutate(
    PCR_covid_Ct = parse_number(PCR_covid_Ct, locale = locale(decimal_mark = ",")),
    PCR_grippe_Ct = parse_number(PCR_grippe_Ct, locale = locale(decimal_mark = ",")),
    PCR_VRS_Ct = parse_number(PCR_VRS_Ct, locale = locale(decimal_mark = ","))
  )

df_Ct <- df_Ct |>
  mutate(PCR_covid_Ct = ifelse(PCR_covid_Ct <= 10 | PCR_covid_Ct > 40, NA, PCR_covid_Ct),
         PCR_grippe_Ct = ifelse(PCR_grippe_Ct <= 10 | PCR_grippe_Ct > 40, NA, PCR_grippe_Ct),
         PCR_VRS_Ct = ifelse(PCR_VRS_Ct <= 10 | PCR_VRS_Ct > 40, NA, PCR_VRS_Ct))


# Monitor Ct values to keep only positive tests with an associated Ct value

df_Ct <- df_Ct |>
  filter(
    # Case 1: All PCR tests are "neg"
    (PCR_covid  == "Neg" & PCR_grippe == "Neg" & PCR_VRS == "Neg") |
      # Case 2: At least one PCR is "Pos" and the corresponding Ct is available
      ( (PCR_covid == "Pos"  & !is.na(PCR_covid_Ct))  |
          (PCR_grippe == "Pos" & !is.na(PCR_grippe_Ct)) |
          (PCR_VRS == "Pos"     & !is.na(PCR_VRS_Ct)) )
  )

df <- df_Ct

df %>%
  summarise(nb_individus = n_distinct(patient_ID))

rm(df_Ct)

# ---------------------------------------------------------------------------
# Symptoms management
# ---------------------------------------------------------------------------

df_symp <- df |>
  mutate(fievre = ifelse(df$fievre %in% c("O", "N"), df$fievre, NA),
         signes_respiratoires = 
           ifelse(df$signes_respiratoires %in% c("O", "N"), df$signes_respiratoires, NA),
         symptom_respi_onset = case_when(
           debut_signe_respi == 1 | (debut_signe_respi == 2 & is.na(signes_respiratoires)) ~ NA_character_,
           debut_signe_respi %in% c(2,7) ~ NA_character_,
           debut_signe_respi == 3 ~ "0-1",
           debut_signe_respi == 4 ~ "2-5",
           debut_signe_respi == 5 ~ "6-10"),
         symptom_fievre_onset = case_when(
           debut_signe_respi_fievre == 1 | (debut_signe_respi_fievre == 2 & is.na(fievre)) ~ NA_character_,
           debut_signe_respi_fievre %in% c(2,3) ~ NA_character_,
           debut_signe_respi_fievre == 4 ~ "0-1",
           debut_signe_respi_fievre == 5 ~ "2-5",
           debut_signe_respi_fievre == 6 ~ "6-10"))

table(df_symp$fievre, useNA = "always")
table(df_symp$signes_respiratoires, useNA = "always")
table(df_symp$symptom_respi_onset, useNA = "always")
table(df_symp$symptom_fievre_onset, useNA = "always")
table(df$debut_signe_respi, useNA = "always")
table(df$debut_signe_respi_fievre, useNA = "always")

df_symp <- df_symp |>
  mutate(fievre = case_when(
           fievre == "N" | debut_signe_respi_fievre == 2 ~ "No",
           fievre == "O" | debut_signe_respi_fievre %in% c(4,5,6) ~ "Yes"),
         signes_respiratoires = case_when(
           signes_respiratoires == "N" ~ "No",
           signes_respiratoires == "O" ~ "Yes"))


df <- df_symp

df %>%
  summarise(nb_individus = n_distinct(patient_ID))

rm(df_symp)

# ---------------------------------------------------------------------------
# Healthcare establishment
# ---------------------------------------------------------------------------

table(df$etablissement_soins, useNA = "always")

df_s <- df |>
  mutate(etablissement_soins = ifelse(etablissement_soins %in% c("O", "E"), "Yes", "No"))

table(df_s$etablissement_soins,  useNA = "always")

df <- df_s

rm(df_s)


# ---------------------------------------------------------------------------
# Vaccine information
# ---------------------------------------------------------------------------

# --- Step 1: Update missing vaccination info after first known vaccine for COVID-19
# NOTE: This section preserves the original imputation rule. Review carefully if
#       vaccination status can change more than once within a patient's records.
df_vax <- df %>%
  group_by(ID) %>%
  mutate(
    first_vacc_date_covid = if (any(vaccin_covid %in% 3:6))
      min(date_prelevement[vaccin_covid %in% 3:6], na.rm = TRUE)
    else as.Date(NA),
    
    last_valid_vacc_covid = if (any(vaccin_covid %in% 3:6))
      vaccin_covid[which.max(date_prelevement[vaccin_covid %in%3:6])]
    else NA_integer_,
    
    vaccin_covid = if_else(
      (is.na(vaccin_covid) | vaccin_covid %in% c(1,2)) &
        !is.na(first_vacc_date_covid) &
        date_prelevement > first_vacc_date_covid,
      last_valid_vacc_covid,
      vaccin_covid
    )
  ) %>%
  ungroup() %>%
  select(-first_vacc_date_covid, -last_valid_vacc_covid)

# --- Step 2: Update missing vaccination info after first known vaccine for influenza
df_vax <- df_vax %>%
  group_by(ID) %>%
  mutate(
    first_vacc_date_grippe = if (any(vaccin_grippe %in% 3:6))
      min(date_prelevement[vaccin_grippe %in% 3:6], na.rm = TRUE)
    else as.Date(NA),
    
    last_valid_vacc_grippe = if (any(vaccin_grippe %in%3:6))
      vaccin_grippe[which.max(date_prelevement[vaccin_grippe %in% 3:6])]
    else NA_integer_,
    
    vaccin_grippe = if_else(
      (is.na(vaccin_grippe) | vaccin_grippe %in% c(1,2)) &
        !is.na(first_vacc_date_grippe) &
        date_prelevement > first_vacc_date_grippe,
      last_valid_vacc_grippe,
      vaccin_grippe
    )
  ) %>%
  ungroup() %>%
  select(-first_vacc_date_grippe, -last_valid_vacc_grippe)

# --- Step 3: Label vaccine status
df_vax <- df_vax %>%
  mutate(
    vaccin_covid_statut = case_when(
      vaccin_covid == 1 ~ NA,
      vaccin_covid == 2 ~ NA,
      vaccin_covid == 3 ~ 0,
      vaccin_covid == 4 ~ 1,
      vaccin_covid == 5 ~ 1,
      vaccin_covid == 6 ~ 0,
      vaccin_covid == 7 ~ 0,
      TRUE ~ NA
    ),
    vaccin_grippe_statut = case_when(
      vaccin_grippe == 1 ~ NA,
      vaccin_grippe == 2 ~ NA,
      vaccin_grippe == 3 ~ 0,
      vaccin_grippe == 4 ~ 1,
      vaccin_grippe == 5 ~ 1,
      vaccin_grippe == 6 ~ 0,
      vaccin_grippe == 7 ~ 0,
      TRUE ~ NA
    )
  )

# --- Step 4: Remove individuals with no vaccination information at all
df_vax <- df_vax %>%
  filter(
    # COVID positif : garder seulement si info vaccin COVID
    (PCR_covid == "Pos" & !is.na(vaccin_covid_statut)) |
      
      # Grippe positive : garder seulement si info vaccin grippe
      (PCR_grippe == "Pos" & !is.na(vaccin_grippe_statut)) |
      
      # RSV positif : garder, à adapter si tu as une variable vaccin RSV
      (PCR_VRS == "Pos") |
      
      # Négatifs : garder si tu veux les conserver comme témoins
      ((PCR_covid != "Pos" & PCR_grippe != "Pos" & PCR_VRS != "Pos") & (!is.na(vaccin_grippe_statut) | !is.na(vaccin_covid_statut)))
  )
df <- df_vax

df %>%
  summarise(nb_individus = n_distinct(patient_ID))

rm(df_vax)

# ---------------------------------------------------------------------------
# Time Since Symptoms Onset (TSS) Management
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 1. Create infection pair ID (ID_paire)
#    -> Groups tests from the same patient occurring within 21 days
# ---------------------------------------------------------------------------
df_tss <- df %>%
  arrange(ID, date_prelevement) %>%
  group_by(ID) %>%
  mutate(
    # Compute days since the previous test
    days_since_prev = as.numeric(date_prelevement - lag(date_prelevement, default = first(date_prelevement))),
    
    # Start a new episode if >21 days between two tests
    new_episode = if_else(row_number() == 1 | days_since_prev > 21, 1, 0),
    
    # Cumulative count of new episodes → episode identifier
    ID_paire = cumsum(new_episode)
  ) %>%
  ungroup() %>%
  select(-days_since_prev, -new_episode)  # clean up temporary columns

# ---------------------------------------------------------------------------
# 2. Create infection ID (ID_infection)
#    -> Each infection episode can have multiple viruses (Covid, Influenza, RSV)
# ---------------------------------------------------------------------------
df_tss$date_prelevement <- as.Date(df_tss$date_prelevement, format = "%Y-%m-%d")

df_tss <- df_tss %>%
  arrange(ID, ID_paire, date_prelevement) %>%
  group_by(ID, ID_paire) %>%
  mutate(
    # Virus détecté
    viruses_detected = pmap_chr(
      list(PCR_covid, PCR_grippe, PCR_VRS),
      ~ {
        viruses <- c(
          ifelse(..1 == "Pos", "Covid", NA),
          ifelse(..2 == "Pos", "Influenza", NA),
          ifelse(..3 == "Pos", "RSV", NA)
        )
        viruses <- viruses[!is.na(viruses)]
        if (length(viruses) == 0) "None" else paste(viruses, collapse = "_")
      }
    ),
    
    is_pos = viruses_detected != "None",
    
    # Nouveau épisode si positif et >21 jours depuis le précédent
    new_infection = is_pos & 
      (is.na(lag(date_prelevement)) | as.numeric(date_prelevement - lag(date_prelevement)) >= 21),
    
    infection_episode = cumsum(replace_na(new_infection, FALSE)),
    
    # Intervalles pour les négatifs
    days_since_prev = as.numeric(date_prelevement - lag(date_prelevement)),
    days_to_next    = as.numeric(lead(date_prelevement) - date_prelevement),
    
    # Déterminer quels négatifs sont attachés à un épisode
    attach_negative =
      !is_pos &
      (
        (!is.na(days_since_prev) & days_since_prev < 21) |
          (!is.na(days_to_next) & days_to_next <= 5)
      ),
    
    # Créer une colonne temporaire pour l'épisode
    episode_temp = if_else(is_pos | attach_negative, infection_episode, NA_integer_),
    
    # Virus associé au premier positif
    virus_for_episode = if_else(is_pos, viruses_detected, NA_character_)
  ) %>%
  # Propager l'épisode et le virus vers les tests négatifs attachés
  fill(episode_temp, virus_for_episode, .direction = "down") %>%
  mutate(
    # Construire l'ID infection final
    ID_infection = if_else(
      is.na(episode_temp) | is.na(virus_for_episode),
      paste("25_26_", ID, ID_paire, "none", sep = "_"),
      paste("25_26_", ID, ID_paire, episode_temp, virus_for_episode, sep = "_")
    )
  ) %>%
  ungroup() %>%
  select(-is_pos, -new_infection, -days_since_prev, -days_to_next, -attach_negative, -episode_temp, -virus_for_episode)

# Traiter le cas particulier où dans une même infection on a un cas de grippe B et un de grippe A

df_tss <- df_tss %>%
  group_by(ID_infection) %>%
  mutate(
    n_typages_grippe = n_distinct(
      PCR_grippe_typage_MultiName[!is.na(PCR_grippe_typage_MultiName)]
    ),
    split_infection_grippe =
      grepl("Influenza", ID_infection) &
      n_typages_grippe >= 2
  ) %>%
  ungroup() %>%
  mutate(
    ID_infection = if_else(
      split_infection_grippe & !is.na(PCR_grippe_typage_MultiName),
      paste(
        sub("_Influenza$", "", ID_infection),
        PCR_grippe_typage_MultiName,
        "Influenza",
        sep = "_"
      ),
      ID_infection
    )
  ) %>%
  select(-n_typages_grippe, -split_infection_grippe)


### Enlever les duplicats de date_dossier

df_no_dup <- df_tss %>%
  group_by(ID_infection, date_dossier) %>%
  arrange(date_prelevement, .by_group = TRUE) %>%   
  slice(1) %>%   # garde la ligne la plus ancienne pour chaque ID_infection/date_dossier
  ungroup()



# ---------------------------------------------------------------------------
# 3. Handle TSS (Time Since Symptoms Onset)
# ---------------------------------------------------------------------------

# --- Step 3.1: Assign initial time_of_symptoms for NA infections (no ID_infection)
df_tss <- df_no_dup %>%
  mutate(
    time_of_symptoms_respi = case_when(
      is.na(ID_infection) & debut_signe_respi == 3 ~ 1,
      is.na(ID_infection) & debut_signe_respi == 4 ~ 3,
      is.na(ID_infection) & debut_signe_respi == 5 ~ 8,
      is.na(ID_infection) ~ NA_real_,
      TRUE ~ NA_real_
    ),
    time_of_symptoms_fievre = case_when(
      is.na(ID_infection) & debut_signe_respi_fievre == 4 ~ 1,
      is.na(ID_infection) & debut_signe_respi_fievre == 5 ~ 3,
      is.na(ID_infection) & debut_signe_respi_fievre == 6 ~ 8,
      is.na(ID_infection) ~ NA_real_,
      TRUE ~ NA_real_
    )
  )

# --- Step 3.2: Define lower and upper TSS bounds based on respiratory symptom category
df_tss_update <- df_tss %>%
  mutate(date_prelevement = as.Date(date_prelevement)) %>%
  mutate(
    lower_TSS_respi = case_when(
      debut_signe_respi == 3 ~ date_prelevement - 1,
      debut_signe_respi == 4 ~ date_prelevement - 5,
      debut_signe_respi == 5 ~ date_prelevement - 10,
      TRUE ~ as.Date(NA)
    ),
    upper_TSS_respi = case_when(
      debut_signe_respi == 3 ~ date_prelevement - 0,
      debut_signe_respi == 4 ~ date_prelevement - 2,
      debut_signe_respi == 5 ~ date_prelevement - 6,
      TRUE ~ as.Date(NA)
    ),
    lower_TSS_fievre = case_when(
      debut_signe_respi_fievre == 4 ~ date_prelevement - 1,
      debut_signe_respi_fievre == 5 ~ date_prelevement - 5,
      debut_signe_respi_fievre == 6 ~ date_prelevement - 10,
      TRUE ~ as.Date(NA)
    ),
    upper_TSS_fievre = case_when(
      debut_signe_respi_fievre == 4 ~ date_prelevement - 0,
      debut_signe_respi_fievre == 5 ~ date_prelevement - 2,
      debut_signe_respi_fievre == 6 ~ date_prelevement - 6,
      TRUE ~ as.Date(NA)
    )
  )

# --- Step 3.3: Compute the final time_of_symptoms per infection

get_tss_episode <- function(lower, upper) {
  
  ok <- !is.na(lower) & !is.na(upper)
  lower <- lower[ok]
  upper <- upper[ok]
  
  if (length(lower) == 0) {
    return(as.Date(NA))
  }
  
  n <- length(lower)
  
  # 1) Intersection commune à tous les intervalles
  common_lower <- max(lower)
  common_upper <- min(upper)
  
  if (common_lower <= common_upper) {
    return(as.Date(
      (as.numeric(common_lower) + as.numeric(common_upper)) / 2,
      origin = "1970-01-01"
    ))
  }
  
  # 2) Sinon, chercher le plus grand sous-groupe d'intervalles
  # qui possède une intersection commune
  if (n >= 2) {
    for (k in n:2) {
      
      combs <- combn(seq_len(n), k, simplify = FALSE)
      
      valid_combs <- lapply(combs, function(idx) {
        
        l <- max(lower[idx])
        u <- min(upper[idx])
        
        if (l <= u) {
          data.frame(
            k = k,
            common_lower = l,
            common_upper = u
          )
        } else {
          NULL
        }
      })
      
      valid_combs <- dplyr::bind_rows(valid_combs)
      
      if (nrow(valid_combs) > 0) {
        
        valid_combs <- valid_combs %>%
          mutate(
            mid = as.Date(
              (as.numeric(common_lower) + as.numeric(common_upper)) / 2,
              origin = "1970-01-01"
            )
          ) %>%
          arrange(mid)
        
        return(valid_combs$mid[1])
      }
    }
  }
  
  # 3) Si aucun intervalle ne se recoupe :
  # prendre la date de symptômes la plus ancienne
  as.Date(min(lower), origin = "1970-01-01")
}


df_tss_update <- df_tss_update %>%
  group_by(ID, ID_infection) %>%
  mutate(
    n_obs = n(),
    
    time_of_symptoms_respi = if_else(
      !is.na(ID_infection),
      get_tss_episode(lower_TSS_respi, upper_TSS_respi),
      as.Date(NA)
    ),
    
    time_of_symptoms_fievre = if_else(
      !is.na(ID_infection),
      get_tss_episode(lower_TSS_fievre, upper_TSS_fievre),
      as.Date(NA)
    )
  ) %>%
  ungroup() %>%
  
  mutate(
    time_of_symptoms_respi = if_else(
      is.na(time_of_symptoms_respi) &
        !is.na(lower_TSS_respi) &
        !is.na(upper_TSS_respi),
      as.Date(
        (as.numeric(lower_TSS_respi) + as.numeric(upper_TSS_respi)) / 2,
        origin = "1970-01-01"
      ),
      time_of_symptoms_respi
    ),
    
    time_of_symptoms_fievre = if_else(
      is.na(time_of_symptoms_fievre) &
        !is.na(lower_TSS_fievre) &
        !is.na(upper_TSS_fievre),
      as.Date(
        (as.numeric(lower_TSS_fievre) + as.numeric(upper_TSS_fievre)) / 2,
        origin = "1970-01-01"
      ),
      time_of_symptoms_fievre
    )
  ) %>%
  select(-n_obs)

# --- Step 3.4: Compute TSS (days between sampling and symptom onset)
df_tss_final <- df_tss_update %>%
  mutate(
    TSS_respi = as.numeric(date_prelevement - time_of_symptoms_respi),
    TSS_fievre = as.numeric(date_prelevement - time_of_symptoms_fievre)
  ) 

#### Ajouter filtre sur TSS
df_tss_final <- df_tss_final |>
  filter(!is.na(TSS_respi) | !is.na(TSS_fievre))

df <- df_tss_final |>
  mutate(
    time_since_symptoms_onset_respi = TSS_respi,
    time_since_symptoms_onset_fievre = TSS_fievre,
    time_since_symptoms_onset = ifelse(is.na(TSS_respi), TSS_fievre, TSS_respi)
  ) |>
  dplyr::select(-TSS_respi, -TSS_fievre, -upper_TSS_respi, -lower_TSS_respi, -upper_TSS_fievre, -lower_TSS_fievre)

# Adding the number of tests by person
df <- df %>%
  group_by(ID) %>%
  mutate(n_tests = n()) %>%
  ungroup()

df %>%
  summarise(nb_individus = n_distinct(patient_ID))

df %>%
  summarise(nb_infections = n_distinct(ID_infection))

rm(df_tss, df_tss_final, df_tss_update, df_no_dup)

# ---------------------------------------------------------------------------
# Uniformize Influenza Typage
# ---------------------------------------------------------------------------
df_grippe <- df %>%
  mutate(
    variant_grippe = case_when(
      !is.na(FLU_SeqClade_Lyon) &
        FLU_SeqClade_Lyon != "" ~ FLU_SeqClade_Lyon,
      TRUE ~ FLU_clade
    )
  ) 

# Update variant names to standardized groups
df_grippe <- df_grippe %>%
  mutate(grippe_subtype = case_when(
    grepl("5a.2a", variant_grippe) ~ "H1N1",
    grepl("2a.3a.1", variant_grippe) ~ "H3N2",
    grepl("V1A.3a.2", variant_grippe) ~ "B/Victoria",
    TRUE ~ "Others"
  ))

# Uniformize Influenza types
df_grippe <- df_grippe %>%
  mutate(PCR_grippe_typage_MultiName = case_when(
    PCR_grippe_typage_MultiName %in% c("A", "GRA", "Grippe A") ~ "Grippe A",
    PCR_grippe_typage_MultiName %in% c("B", "GRB", "Grippe B") ~ "Grippe B",
    grippe_subtype %in% c("H1N1", "H3N2") ~ "Grippe A",
    grippe_subtype %in% c("B/Victoria") ~ "Grippe B",
    TRUE ~ NA 
  )) %>%
  select(-variant_grippe)

df <- df_grippe %>%
  group_by(ID_infection) %>%
  mutate(
    PCR_grippe_typage_MultiName = if_else(
      grepl("flu", ID_infection) &
        is.na(PCR_grippe_typage_MultiName) &
        any(!is.na(PCR_grippe_typage_MultiName)),
      PCR_grippe_typage_MultiName[!is.na(PCR_grippe_typage_MultiName)][1],
      PCR_grippe_typage_MultiName
    )
  ) %>%
  ungroup() %>%
  filter(!(PCR_grippe == "Pos" & is.na(PCR_grippe_typage_MultiName) & PCR_covid == "Neg" & PCR_VRS == "Neg"))

df %>%
  summarise(nb_individus = n_distinct(patient_ID))

df %>%
  summarise(nb_infections = n_distinct(ID_infection))

rm(df_grippe)



# Uniformize test techniques
df <- df %>%
  mutate(
    PCR_covid_technique_update = case_when(PCR_covid_technique %in% c("Eurobio - R\xe9actif Allplex (Starlet/CFX) ", "Eurobio - Réactif Allplex (Starlet/CFX)"," Eurobio EBX 42") ~ "Eurobio",TRUE ~ PCR_covid_technique),
    PCR_grippe_technique_update = case_when(PCR_grippe_technique %in% c("Eurobio - R\xe9actif Allplex (Starlet/CFX) ", "Eurobio - Réactif Allplex (Starlet/CFX)"," Eurobio EBX 42") ~ "Eurobio",TRUE ~ PCR_grippe_technique),
    PCR_VRS_technique_update = case_when(PCR_VRS_technique %in% c("Eurobio - R\xe9actif Allplex (Starlet/CFX) ", "Eurobio - Réactif Allplex (Starlet/CFX)"," Eurobio EBX 42") ~ "Eurobio",TRUE ~ PCR_VRS_technique)
    )




# ------------------------------------------------------------------------------
# Save processed data
# ------------------------------------------------------------------------------
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

write.csv(
  df,
  file = output_file,
  row.names = FALSE
)

message("Processed dataset written to: ", output_file)
