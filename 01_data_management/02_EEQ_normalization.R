################################################################################################
# EEQ-based Ct normalization for RELAB data
# Purpose:
#   1. Clean and harmonize External Quality Evaluation (EEQ) Ct data
#   2. Estimate technique- and virus-specific Ct normalization rules
#   3. Assign an EEQ PCR technique to each RELAB test
#   4. Normalize RELAB Ct values
################################################################################################

rm(list = ls())

# ---------------------------------------------------------------------------
# 0. Libraries
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(readr)


# ---------------------------------------------------------------------------
# 1. File paths
# ---------------------------------------------------------------------------

relab_file <- "data/processed/RELAB_24_25_clean_start.csv"
eeq_file   <- "data/raw/EEQ_2025_negatif.csv"
cerba_file <- "data/raw/cartographie_technique_relab_v3.csv"

rules_output <- "data/processed/EEQ_Ct_normalization_rules.csv"
relab_output <- "data/processed/RELAB_24_25_clean_start_norm.csv"


# ---------------------------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------------------------

parse_ct <- function(x) {
  x <- str_replace_all(as.character(x), ",", ".")
  x[x %in% c("", "0", "*", "NaN")] <- NA_character_
  suppressWarnings(as.numeric(x))
}


standardize_technique <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    str_detect(x, regex("pkamp", ignore_case = TRUE)) ~ "PKamp",
    str_detect(x, regex("alinity", ignore_case = TRUE)) &
      str_detect(x, regex("anatolia", ignore_case = TRUE)) ~ x,
    str_detect(x, regex("alinity", ignore_case = TRUE)) ~ "Alinity",
    str_detect(x, regex("anatolia", ignore_case = TRUE)) ~ "Anatolia",
    str_detect(x, regex("elitech", ignore_case = TRUE)) ~ "Elitech",
    str_detect(x, regex("eurobio|id solution|idsolution", ignore_case = TRUE)) ~ "Eurobio",
    str_detect(x, regex("allplex|seegene", ignore_case = TRUE)) ~ "Seegene",
    str_detect(x, regex("biosynex", ignore_case = TRUE)) ~ "Biosynex",
    str_detect(x, regex("appolon", ignore_case = TRUE)) ~ "Appolon",
    str_detect(x, regex("panther", ignore_case = TRUE)) ~ "Panther",
    TRUE ~ x
  )
}


# ---------------------------------------------------------------------------
# 3. Load input data
# ---------------------------------------------------------------------------

relab <- read_csv(relab_file, show_col_types = FALSE)

eeq <- read_delim(
  eeq_file,
  delim = ";",
  show_col_types = FALSE,
  trim_ws = TRUE
)

cerba_info <- read_delim(
  cerba_file,
  delim = ";",
  show_col_types = FALSE,
  trim_ws = TRUE
)


# ---------------------------------------------------------------------------
# 4. Clean and harmonize EEQ data
# ---------------------------------------------------------------------------

eeq_clean <- eeq %>%
  rename(test_center = centre_test) %>%
  select(-starts_with("X")) %>%
  mutate(
    Virus_test = case_when(
      str_detect(ID_EEQ_test, "^SC")  ~ "SARS-CoV-2",
      str_detect(ID_EEQ_test, "^GRA") ~ "Influenza A",
      str_detect(ID_EEQ_test, "^GRB") ~ "Influenza B",
      str_detect(ID_EEQ_test, "^VB")  ~ "RSV",
      TRUE ~ NA_character_
    ),

    Ct_range = case_when(
      str_detect(ID_EEQ_test, "^SC\\s?102|^GRA\\s?102|^GRB\\s?102|^VB[^0-9]*103") ~ "High",
      str_detect(ID_EEQ_test, "^SC\\s?104|^GRA\\s?104|^GRB\\s?104|^VB[^0-9]*105") ~ "Medium",
      str_detect(ID_EEQ_test, "^SC\\s?105|^SC\\s?106|^GRA\\s?106|^GRB\\s?106|^VB[^0-9]*107") ~ "Low",
      TRUE ~ NA_character_
    ),

    across(
      c(CT_SARS.CoV2, CT_Grippe, CT_Grippe_A, CT_Grippe_B, CT_VRS),
      parse_ct
    ),

    PCR_technique = case_when(
      test_center %in% c(
        "Cerballiance HDF",
        "Cerballiance NAQ",
        "Cerballiance Provence Marseille",
        "Cerballiance ARA"
      ) ~ "Alinity",

      test_center %in% c(
        "Alphabio",
        "Biogroup Bioesterel",
        "Biogroup Lorraine"
      ) ~ "Elitech",

      test_center %in% c(
        "Biogroup Unilians",
        "BIOLAM LCD Biogroup",
        "BIOLAM LCD Biogroup (fev-2026)",
        "Biogroup Unilians (fev-2026)"
      ) ~ "Eurobio",

      test_center %in% c(
        "Cerballiance LISSES",
        "Cerballiance IDF",
        "Cerba Spé (Frépillon)"
      ) ~ "PKamp",

      test_center == "Cerballiance Provence-Azur" ~ "Anatolia",
      test_center %in% c("Cerballiance Occitanie", "Inovie Genbio") ~ "Seegene",
      test_center == "BioLBS" ~ "Biosynex",
      test_center == "HCL" ~ "Panther",
      test_center == "Cerballiance St Martin" ~ "Appolon",
      test_center == "Cerballiance Normandie" ~ "Vitro respi",
      TRUE ~ NA_character_
    )
  ) %>%

  # Use the generic influenza Ct when the subtype-specific Ct is not the
  # appropriate output for these techniques/centers.
  mutate(
    CT_Grippe_A = if_else(
      test_center %in% c(
        "Cerballiance Provence-Azur",
        "Biogroup Unilians",
        "BIOLAM LCD Biogroup"
      ) &
        str_detect(ID_EEQ_test, "^GRA"),
      CT_Grippe,
      CT_Grippe_A
    ),

    CT_Grippe_B = if_else(
      test_center %in% c(
        "Cerballiance Provence-Azur",
        "Biogroup Unilians",
        "BIOLAM LCD Biogroup"
      ) &
        str_detect(ID_EEQ_test, "^GRB"),
      CT_Grippe,
      CT_Grippe_B
    )
  ) %>%

  # Exclude known problematic measurements.
  filter(
    !(test_center == "Cerballiance Normandie" & Extracteur == "N/A"),
    !test_center %in% c("Biogroup Unilians", "BIOLAM LCD Biogroup")
  ) %>%

  mutate(
    test_center = recode(
      test_center,
      "BIOLAM LCD Biogroup (fev-2026)" = "BIOLAM LCD Biogroup",
      "Biogroup Unilians (fev-2026)" = "Biogroup Unilians"
    ),

    CT_Grippe_B = if_else(
      test_center %in% c(
        "Cerballiance NAQ",
        "Cerballiance Provence Marseille"
      ) &
        ID_EEQ_test %in% c(
          "GRA106-LY230625",
          "GRA102-LY230625",
          "GRA104-LY230625"
        ),
      NA_real_,
      CT_Grippe_B
    ),

    CT_Grippe_A = if_else(
      test_center == "Cerballiance HDF" &
        ID_EEQ_test == "GRB104-LY230625",
      NA_real_,
      CT_Grippe_A
    )
  ) %>%

  filter(
    !is.na(Virus_test),
    !is.na(Ct_range),
    !is.na(PCR_technique)
  )


# ---------------------------------------------------------------------------
# 5. Reshape EEQ data and calculate normalization rules
# ---------------------------------------------------------------------------

eeq_long <- eeq_clean %>%
  select(
    ID_EEQ_test,
    test_center,
    PCR_technique,
    Virus_test,
    Ct_range,
    CT_SARS.CoV2,
    CT_Grippe_A,
    CT_Grippe_B,
    CT_VRS
  ) %>%
  pivot_longer(
    cols = c(CT_SARS.CoV2, CT_Grippe_A, CT_Grippe_B, CT_VRS),
    names_to = "Ct_variable",
    values_to = "Ct_value"
  ) %>%
  mutate(
    Virus_measurement = recode(
      Ct_variable,
      "CT_SARS.CoV2" = "SARS-CoV-2",
      "CT_Grippe_A"  = "Influenza A",
      "CT_Grippe_B"  = "Influenza B",
      "CT_VRS"       = "RSV"
    )
  ) %>%
  filter(
    Virus_test == Virus_measurement,
    !is.na(Ct_value)
  )


# Mean Ct for each technique, virus and viral-load range.
technique_means <- eeq_long %>%
  group_by(Virus_test, Ct_range, PCR_technique) %>%
  summarise(
    Ct_mean_technique = mean(Ct_value, na.rm = TRUE),
    .groups = "drop"
  )


# Mean Ct across techniques for each virus and viral-load range.
range_means <- eeq_long %>%
  group_by(Virus_test, Ct_range) %>%
  summarise(
    Ct_mean_range = mean(Ct_value, na.rm = TRUE),
    .groups = "drop"
  )


# Technique-specific bias:
#   bias = mean Ct for a technique - overall mean Ct
#
# The final correction is the mean bias across the available
# High / Medium / Low viral-load ranges.
eeq_rules <- technique_means %>%
  left_join(
    range_means,
    by = c("Virus_test", "Ct_range")
  ) %>%
  mutate(
    Ct_bias = Ct_mean_technique - Ct_mean_range
  ) %>%
  group_by(Virus_test, PCR_technique) %>%
  summarise(
    bias_mean_3ranges = mean(Ct_bias, na.rm = TRUE),
    n_ranges = sum(!is.na(Ct_bias)),
    .groups = "drop"
  ) %>%
  arrange(Virus_test, PCR_technique)


write_csv(eeq_rules, rules_output)


# ---------------------------------------------------------------------------
# 6. Assign analysis center and EEQ technique to RELAB observations
# ---------------------------------------------------------------------------

cerba_map <- cerba_info %>%
  transmute(
    id_group_cel = as.character(id_groupe_cel),
    cerba_center = selas,
    cerba_technique = technique
  ) %>%
  distinct(id_group_cel, .keep_all = TRUE)


relab_technique <- relab %>%
  mutate(
    id_group_cel = if_else(
      Labo_nom == "CERBA",
      str_remove(
        str_extract(as.character(patient_ID), "^\\d+"),
        "^0+"
      ),
      NA_character_
    )
  ) %>%
  left_join(cerba_map, by = "id_group_cel") %>%
  mutate(
    centre_analyse = case_when(
      Labo_nom == "BIOGROUP" & departement_centre_preleveur == 6  ~ "Biogroup Bioesterel",
      Labo_nom == "BIOGROUP" & departement_centre_preleveur == 57 ~ "Biogroup Lorraine",
      Labo_nom == "BIOGROUP" & departement_centre_preleveur == 69 ~ "Biogroup Unilians",
      Labo_nom == "BIOGROUP" & departement_centre_preleveur == 93 ~ "BIOLAM LCD Biogroup",
      Labo_nom == "BIOGROUP" & departement_centre_preleveur == 13 ~ "Alphabio",

      Labo_nom == "CERBA" & departement_centre_preleveur == 59 ~ "Cerballiance HDF",
      Labo_nom == "CERBA" & departement_centre_preleveur == 91 ~ "Cerballiance LISSES",
      Labo_nom == "CERBA" & departement_centre_preleveur == 33 ~ "Cerballiance NAQ",
      Labo_nom == "CERBA" & departement_centre_preleveur == 14 ~ "Cerballiance St Martin",
      Labo_nom == "CERBA" & departement_centre_preleveur == 31 ~ "Cerballiance Occitanie",
      Labo_nom == "CERBA" & departement_centre_preleveur == 6  ~ "Cerballiance Provence-Azur",
      Labo_nom == "CERBA" & departement_centre_preleveur == 13 ~ "Cerballiance Provence Marseille",
      Labo_nom == "CERBA" & departement_centre_preleveur == 69 ~ "Cerballiance ARA",
      Labo_nom == "CERBA" & departement_centre_preleveur == 93 ~ "Cerballiance IDF",
      Labo_nom == "CERBA" & departement_centre_preleveur == 95 ~ "Cerba Spé (Frépillon)",

      Labo_nom == "BioLBS" ~ "BioLBS",
      Labo_nom == "INOVIE" ~ "Inovie Genbio",

      str_detect(coalesce(PCR_covid_technique, ""), regex("Elitech", ignore_case = TRUE)) |
        str_detect(coalesce(PCR_grippe_technique, ""), regex("Elitech", ignore_case = TRUE)) |
        str_detect(coalesce(PCR_VRS_technique, ""), regex("Elitech", ignore_case = TRUE)) ~ "Elitech",

      TRUE ~ NA_character_
    ),

    technique_EEQ = case_when(
      !is.na(cerba_technique) ~ cerba_technique,

      centre_analyse %in% c(
        "Elitech",
        "Alphabio",
        "Biogroup Bioesterel",
        "Biogroup Lorraine"
      ) ~ "Elitech",

      centre_analyse %in% c(
        "Biogroup Unilians",
        "BIOLAM LCD Biogroup"
      ) ~ "Eurobio",

      centre_analyse == "Inovie Genbio" ~ "Seegene",
      centre_analyse == "BioLBS" ~ "Biosynex",
      centre_analyse == "Cerballiance Occitanie" ~ "Seegene",

      departement_centre_preleveur %in% c(34, 9) ~ "Seegene",

      str_detect(coalesce(PCR_covid_technique, ""), regex("Eurobio|Id Solution|IDSolution", ignore_case = TRUE)) ~ "Eurobio",
      str_detect(coalesce(PCR_covid_technique, ""), regex("Elitech", ignore_case = TRUE)) ~ "Elitech",

      TRUE ~ NA_character_
    ),

    technique_EEQ = standardize_technique(technique_EEQ),

    # Resolve mixed Cerballiance "Alinity/MGI-Anatolia" labels.
    technique_EEQ = case_when(
      str_detect(
        coalesce(technique_EEQ, ""),
        regex("Alinity/MGI-Anatolia", ignore_case = TRUE)
      ) &
        (
          centre_analyse == "Cerballiance Provence Marseille" |
            departement_centre_preleveur == 83 |
            is.na(departement_centre_preleveur)
        ) ~ "Alinity",

      str_detect(
        coalesce(technique_EEQ, ""),
        regex("Alinity/MGI-Anatolia", ignore_case = TRUE)
      ) &
        centre_analyse == "Cerballiance Provence-Azur" ~ "Anatolia",

      TRUE ~ technique_EEQ
    )
  )


# ---------------------------------------------------------------------------
# 7. Apply EEQ normalization rules to RELAB Ct values
# ---------------------------------------------------------------------------

rules_wide <- eeq_rules %>%
  select(Virus_test, PCR_technique, bias_mean_3ranges) %>%
  pivot_wider(
    names_from = Virus_test,
    values_from = bias_mean_3ranges
  ) %>%
  rename(
    bias_covid = `SARS-CoV-2`,
    bias_iav   = `Influenza A`,
    bias_ibv   = `Influenza B`,
    bias_rsv   = RSV
  )


relab_norm <- relab_technique %>%
  left_join(
    rules_wide,
    by = c("technique_EEQ" = "PCR_technique")
  ) %>%
  mutate(
    Ct_covid_norm = if_else(
      PCR_covid == "Pos" & !is.na(bias_covid),
      PCR_covid_Ct - bias_covid,
      PCR_covid_Ct
    ),

    Ct_grippe_norm = case_when(
      PCR_grippe == "Pos" &
        PCR_grippe_typage_MultiName == "Grippe A" &
        !is.na(bias_iav) ~ PCR_grippe_Ct - bias_iav,

      PCR_grippe == "Pos" &
        PCR_grippe_typage_MultiName == "Grippe B" &
        !is.na(bias_ibv) ~ PCR_grippe_Ct - bias_ibv,

      TRUE ~ PCR_grippe_Ct
    ),

    Ct_VRS_norm = if_else(
      PCR_VRS == "Pos" & !is.na(bias_rsv),
      PCR_VRS_Ct - bias_rsv,
      PCR_VRS_Ct
    )
  ) %>%

  # For negative follow-up tests belonging to a known infection episode,
  # set Ct to the assay limit of detection (Ct = 40).
  mutate(
    Ct_covid_norm = case_when(
      PCR_covid == "Neg" &
        str_detect(coalesce(ID_infection, ""), regex("covid|cov", ignore_case = TRUE)) ~ 40,
      Ct_covid_norm > 40 ~ 40,
      TRUE ~ Ct_covid_norm
    ),

    Ct_grippe_norm = case_when(
      PCR_grippe == "Neg" &
        str_detect(coalesce(ID_infection, ""), regex("influenza|inf", ignore_case = TRUE)) ~ 40,
      Ct_grippe_norm > 40 ~ 40,
      TRUE ~ Ct_grippe_norm
    ),

    Ct_VRS_norm = case_when(
      PCR_VRS == "Neg" &
        str_detect(coalesce(ID_infection, ""), regex("rsv", ignore_case = TRUE)) ~ 40,
      Ct_VRS_norm > 40 ~ 40,
      TRUE ~ Ct_VRS_norm
    )
  ) %>%

  filter(
    !is.na(Ct_grippe_norm) |
      !is.na(Ct_covid_norm) |
      !is.na(Ct_VRS_norm) |
      str_detect(coalesce(ID_infection, ""), regex("none", ignore_case = TRUE))
  ) %>%

  select(
    -id_group_cel,
    -cerba_center,
    -cerba_technique,
    -bias_covid,
    -bias_iav,
    -bias_ibv,
    -bias_rsv
  )


# ---------------------------------------------------------------------------
# 8. Save normalized RELAB dataset
# ---------------------------------------------------------------------------

write_csv(relab_norm, relab_output)


# ---------------------------------------------------------------------------
# 9. Summary
# ---------------------------------------------------------------------------

message("EEQ normalization rules saved to: ", rules_output)
message("Normalized RELAB dataset saved to: ", relab_output)
message("Number of EEQ normalization rules: ", nrow(eeq_rules))
message("Number of RELAB observations after normalization: ", nrow(relab_norm))
