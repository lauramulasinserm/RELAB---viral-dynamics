# ============================================================
# Vaccine effectiveness against infection
# SARS-CoV-2, Influenza A and Influenza B
# ============================================================

rm(list = ls())

# -------------------------
# 0. Packages
# -------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(broom)


# -------------------------
# 1. Load data
# -------------------------

data_VE <- read.csv(
  "C:/Users/laura.mulas/OneDrive - INSERM/Documents/PhD/Projet_RELAB/RELAB_DATA/RELAB_24_26_clean_start_norm_vfinal_25_06.csv",
  sep = ",",
  stringsAsFactors = FALSE
)


# ============================================================
# 2. Function to estimate VE by epidemic season
# ============================================================

estimate_VE_by_season <- function(data) {

  data %>%
    filter(
      !is.na(infection),
      !is.na(vaccin)
    ) %>%

    group_by(epidemic_season) %>%

    nest() %>%

    mutate(

      model = map(
        data,
        ~ glm(
          infection ~
            vaccin +
            age_cat +
            sexe +
            factor(week) +
            technique_EEQ,
          family = binomial(),
          data = .x
        )
      ),

      tidy = map(
        model,
        ~ broom::tidy(
          .x,
          conf.int = TRUE
        )
      )
    ) %>%

    unnest(tidy) %>%

    filter(
      term == "vaccinVaccinated"
    ) %>%

    mutate(

      OR = exp(estimate),

      OR_low = exp(conf.low),

      OR_high = exp(conf.high),

      VE = 1 - OR,

      VE_low = 1 - OR_high,

      VE_high = 1 - OR_low

    ) %>%

    select(
      epidemic_season,
      VE,
      VE_low,
      VE_high,
      p.value
    )
}


# ============================================================
# 3. Influenza B
# ============================================================

data_VE_B <- data_VE %>%

  mutate(

    vaccin = case_when(
      vaccin_grippe_statut == 0 ~ "Not vaccinated",
      vaccin_grippe_statut == 1 ~ "Vaccinated",
      TRUE ~ NA_character_
    ),

    vaccin = factor(
      vaccin,
      levels = c(
        "Not vaccinated",
        "Vaccinated"
      )
    ),

    infection = case_when(

      PCR_grippe == "Pos" &
        PCR_grippe_typage_MultiName == "Grippe B" ~ 1,

      PCR_grippe == "Neg" ~ 0,

      TRUE ~ NA_real_
    ),

    grippeB_reported =
      !is.na(PCR_grippe_typage_MultiName) &
      PCR_grippe_typage_MultiName == "Grippe B"
  ) %>%

  # Keep one record per infection/test episode,
  # prioritising the row explicitly reported as Influenza B.
  arrange(
    ID_infection,
    desc(grippeB_reported)
  ) %>%

  distinct(
    ID_infection,
    .keep_all = TRUE
  ) %>%

  select(
    -grippeB_reported
  )


res_VE_B <- estimate_VE_by_season(
  data_VE_B
)

res_VE_B


# ============================================================
# 4. Influenza A
# ============================================================

data_VE_A <- data_VE %>%

  mutate(

    vaccin = case_when(
      vaccin_grippe_statut == 0 ~ "Not vaccinated",
      vaccin_grippe_statut == 1 ~ "Vaccinated",
      TRUE ~ NA_character_
    ),

    vaccin = factor(
      vaccin,
      levels = c(
        "Not vaccinated",
        "Vaccinated"
      )
    ),

    infection = case_when(

      PCR_grippe == "Pos" &
        PCR_grippe_typage_MultiName == "Grippe A" ~ 1,

      PCR_grippe == "Neg" ~ 0,

      TRUE ~ NA_real_
    )
  ) %>%

  filter(
    !is.na(infection)
  )


res_VE_A <- estimate_VE_by_season(
  data_VE_A
)

res_VE_A


# ============================================================
# 5. SARS-CoV-2
# ============================================================

data_VE_cov <- data_VE %>%

  mutate(

    vaccin = case_when(
      vaccin_covid_statut == 0 ~ "Not vaccinated",
      vaccin_covid_statut == 1 ~ "Vaccinated",
      TRUE ~ NA_character_
    ),

    vaccin = factor(
      vaccin,
      levels = c(
        "Not vaccinated",
        "Vaccinated"
      )
    ),

    infection = case_when(
      PCR_covid == "Pos" ~ 1,
      PCR_covid == "Neg" ~ 0,
      TRUE ~ NA_real_
    ),

    covid_reported = grepl(
      "Cov",
      ID_infection
    )
  ) %>%

  # Keep one record per infection/test episode,
  # prioritising the COVID-associated row when duplicated.
  arrange(
    ID_infection,
    desc(covid_reported)
  ) %>%

  distinct(
    ID_infection,
    .keep_all = TRUE
  ) %>%

  select(
    -covid_reported
  )


res_VE_cov <- estimate_VE_by_season(
  data_VE_cov
)

res_VE_cov


# ============================================================
# 6. Combine results
# ============================================================

VE_results <- bind_rows(

  res_VE_cov %>%
    mutate(
      virus = "SARS-CoV-2"
    ),

  res_VE_A %>%
    mutate(
      virus = "IAV"
    ),

  res_VE_B %>%
    mutate(
      virus = "IBV"
    )

) %>%

  select(
    virus,
    epidemic_season,
    VE,
    VE_low,
    VE_high,
    p.value
  ) %>%

  arrange(
    virus,
    epidemic_season
  )


print(
  VE_results
)
