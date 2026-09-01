################################################################################################
# RELAB study population characteristics
# Author: Laura MULAS
# Purpose:
#   Prepare infection-level data and generate the descriptive characteristics
#   table for SARS-CoV-2, influenza A, influenza B, and RSV infections.
################################################################################################

rm(list = ls())

# ---------------------------------------------------------------------------
# 0. Libraries
# ---------------------------------------------------------------------------

library(tidyverse)
library(gtsummary)
library(gt)


# ---------------------------------------------------------------------------
# 1. File paths
# ---------------------------------------------------------------------------

covid_file <- "data/processed/data_for_monolix_covid_24_26_vfinal_table_carac.csv"
iav_file   <- "data/processed/data_for_monolix_grippe_A_24_26_vfinal_table_carac.csv"
ibv_file   <- "data/processed/data_for_monolix_grippe_B_24_26_vfinal_table_carac.csv"
rsv_file   <- "data/processed/data_for_monolix_VRS_24_26_vfinal_table_carac.csv"

table_output <- "outputs/tables/characteristics_table.html"
vaccination_output <- "outputs/tables/vaccination_by_season.csv"


# ---------------------------------------------------------------------------
# 2. Load and combine virus-specific datasets
# ---------------------------------------------------------------------------

data_cov <- read_csv(covid_file, show_col_types = FALSE)
data_iav <- read_csv(iav_file, show_col_types = FALSE)
data_ibv <- read_csv(ibv_file, show_col_types = FALSE)
data_rsv <- read_csv(rsv_file, show_col_types = FALSE)

data <- bind_rows(
  data_cov %>% mutate(groupe_virus = "SARS-CoV-2"),
  data_iav %>% mutate(groupe_virus = "IAV"),
  data_ibv %>% mutate(groupe_virus = "IBV"),
  data_rsv %>% mutate(groupe_virus = "RSV")
) %>%
  mutate(
    n_tests = as.numeric(n_tests)
  )


# ---------------------------------------------------------------------------
# 3. Harmonize vaccination status
# ---------------------------------------------------------------------------

data <- data %>%
  mutate(
    vaccination_statut = case_when(
      statut_vaccin == 1 ~ "Vaccinated",
      statut_vaccin == 0 ~ "Not vaccinated",
      TRUE ~ "Unknown"
    ),
    vaccination_statut = factor(
      vaccination_statut,
      levels = c("Vaccinated", "Not vaccinated", "Unknown")
    )
  )


# ---------------------------------------------------------------------------
# 4. Identify co-infections
# ---------------------------------------------------------------------------

coinfection_info <- data %>%
  group_by(ID_infection) %>%
  summarise(
    n_virus = n_distinct(groupe_virus),
    coinfection = case_when(
      n_virus == 1 ~ "No co-infection",
      n_virus == 2 ~ "Co-infection with 1 other virus",
      n_virus >= 3 ~ "Co-infection with 2 other viruses"
    ),
    .groups = "drop"
  )


# ---------------------------------------------------------------------------
# 5. Create infection-virus level dataset
# ---------------------------------------------------------------------------
# One row per infection episode and virus.
# Co-infections therefore contribute one row for each detected virus.

data_infection_level <- data %>%
  group_by(ID_infection, groupe_virus) %>%
  summarise(
    epidemic_season = case_when(
      first(season) == "24_25" ~ "24/25",
      first(season) == "25_26" ~ "25/26",
      TRUE ~ as.character(first(season))
    ),
    sexe = first(sexe),
    age_cat = first(age_cat),
    symptoms = first(symptoms),
    n_tests = first(n_tests),
    Labo_nom = first(Labo_nom),

    # TRUE if at least one observation in the infection episode is censored/negative.
    n_tests_negatifs = any(censor %in% TRUE, na.rm = TRUE),

    vaccination_statut = case_when(
      any(vaccination_statut == "Vaccinated", na.rm = TRUE) ~ "Vaccinated",
      any(vaccination_statut == "Not vaccinated", na.rm = TRUE) ~ "Not vaccinated",
      TRUE ~ "Unknown"
    ),
    .groups = "drop"
  ) %>%

  mutate(
    groupe_virus = factor(
      groupe_virus,
      levels = c("SARS-CoV-2", "IAV", "IBV", "RSV")
    ),
    vaccination_statut = factor(
      vaccination_statut,
      levels = c("Vaccinated", "Not vaccinated", "Unknown")
    ),
    age_cat = factor(
      age_cat,
      levels = c("<5", "5-18", "18-65", ">65")
    ),
    epidemic_season = factor(
      epidemic_season,
      levels = c("24/25", "25/26")
    )
  ) %>%

  # Vaccination status is required for SARS-CoV-2, IAV and IBV analyses.
  # RSV vaccination status was not analysed.
  filter(
    groupe_virus == "RSV" |
      vaccination_statut != "Unknown"
  ) %>%

  left_join(
    coinfection_info,
    by = "ID_infection"
  )


# ---------------------------------------------------------------------------
# 6. Check sample sizes
# ---------------------------------------------------------------------------

virus_counts <- data_infection_level %>%
  count(groupe_virus, name = "n_infections")

print(virus_counts)


# ---------------------------------------------------------------------------
# 7. Generate characteristics table
# ---------------------------------------------------------------------------

table_generale <- data_infection_level %>%
  select(
    epidemic_season,
    sexe,
    age_cat,
    vaccination_statut,
    symptoms,
    coinfection,
    Labo_nom,
    n_tests,
    n_tests_negatifs,
    groupe_virus
  ) %>%

  tbl_summary(
    by = groupe_virus,

    label = list(
      epidemic_season ~ "Epidemic season",
      sexe ~ "Sex",
      age_cat ~ "Age category",
      vaccination_statut ~ "Vaccination against infection virus",
      symptoms ~ "Symptoms",
      coinfection ~ "Co-infection",
      Labo_nom ~ "Laboratory",
      n_tests ~ "Number of tests per infection",
      n_tests_negatifs ~ "Infections with a negative test"
    ),

    statistic = list(
      all_categorical() ~ "{n} ({p}%)",
      all_continuous() ~ "{mean} ({sd})"
    ),

    missing = "no"
  ) %>%

  add_overall(last = TRUE) %>%
  bold_labels() %>%
  modify_header(label = "**Characteristics**") %>%

  # Display only the Female level for sex.
  # The complementary Male percentage can be inferred directly.
  modify_table_body(
    ~ .x %>%
      filter(
        !(variable == "sexe" &
            row_type == "level" &
            label != "Female")
      ) %>%

      # An overall vaccination percentage is not meaningful because
      # vaccination status was not analysed for RSV.
      mutate(
        stat_0 = if_else(
          variable == "vaccination_statut" &
            row_type == "level",
          "–",
          stat_0
        )
      )
  )


# Display table in RStudio
table_generale


# ---------------------------------------------------------------------------
# 8. Vaccination coverage by virus and epidemic season
# ---------------------------------------------------------------------------

table_vaccination_season <- data_infection_level %>%
  filter(
    groupe_virus %in% c("SARS-CoV-2", "IAV", "IBV")
  ) %>%
  group_by(
    groupe_virus,
    epidemic_season
  ) %>%
  summarise(
    n_total = n(),
    n_vaccinated = sum(
      vaccination_statut == "Vaccinated",
      na.rm = TRUE
    ),
    pct_vaccinated = 100 * n_vaccinated / n_total,
    .groups = "drop"
  ) %>%
  mutate(
    vaccinated = sprintf(
      "%d (%.1f%%)",
      n_vaccinated,
      pct_vaccinated
    )
  ) %>%
  select(
    groupe_virus,
    epidemic_season,
    n_total,
    vaccinated
  )

print(table_vaccination_season)


# ---------------------------------------------------------------------------
# 9. Save outputs
# ---------------------------------------------------------------------------

dir.create(
  dirname(table_output),
  recursive = TRUE,
  showWarnings = FALSE
)

gtsave(
  data = as_gt(table_generale),
  filename = table_output
)

write_csv(
  table_vaccination_season,
  vaccination_output
)

message("Characteristics table saved to: ", table_output)
message("Vaccination-by-season table saved to: ", vaccination_output)
