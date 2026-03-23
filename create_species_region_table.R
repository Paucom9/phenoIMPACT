library(dplyr)
library(stringr)
library(tidyr)
library(here)

#### Data ----
phenology_estimates <- read.csv(here::here("output", "pheno_estimates_allspp.csv"))
CBMS_traits         <- read.csv(here::here("data", "Traits_CBMS.csv"), sep = ";")
species_traits <- read.csv(here::here("data", "species_trait_table.csv"), sep = ";", dec = ".")
vol_se <- read.csv(here::here("data","voltinism", "species_country_voltinism_SE.csv"), sep = ";", dec = ".")
vol_nl <- read.csv(here::here("data","voltinism", "species_country_voltinism_NL.csv"), sep = ";", dec = ".")
vol_fi <- read.csv(here::here("data","voltinism", "species_country_voltinism_FI.csv"), sep = ";", dec = ".")


#### Species × site → country ----
df_species_site <- phenology_estimates %>%
  select(SPECIES, SITE_ID) %>%
  distinct() %>%
  mutate(country = str_extract(SITE_ID, "^[^.]+"))

#### Counts table ----
species_country_table <- df_species_site %>%
  count(SPECIES, country) %>%
  pivot_wider(names_from = country, values_from = n, values_fill = 0)

#### Voltinism (traits) ----
traits_volt <- species_traits %>%
  transmute(
    SPECIES = gsub("_", " ", Taxon),
    vol_max = Vol_max
  )

#### Apply univoltine logic ----
species_country_volt <- species_country_table %>%
  left_join(traits_volt, by = "SPECIES") %>%
  mutate(
    across(
      -c(SPECIES, vol_max),
      ~ case_when(
        . == 0 ~ "not present",
        vol_max <= 1 ~ "univoltine",
        TRUE ~ as.character(.)
      )
    )
  )

#### CBMS voltinism (overwrite ES-CTBMS) ----
traits_cbms <- CBMS_traits %>%
  transmute(
    SPECIES = gsub("_", " ", Butterfly.species),
    voltinism_2 = tolower(voltinism_2)
  )

species_country_volt <- species_country_volt %>%
  left_join(traits_cbms, by = "SPECIES") %>%
  mutate(
    `ES-CTBMS` = ifelse(
      `ES-CTBMS` == "not present",
      "not present",
      voltinism_2
    )
  ) %>%
  mutate(
    across(
      -c(SPECIES),
      ~ ifelse(grepl("^[0-9]+$", .), "", .)
    )
  ) %>%
  select(-voltinism_2, -vol_max)

# Merge with SE, FI and NL

vol_se <- vol_se %>%
  select(SPECIES, SEBMS)
vol_fi <- vol_fi %>%
  select(SPECIES, FIBMS)
vol_nl <- vol_nl %>%
  select(SPECIES,  NLBMS)

str(species_country_volt)
str(vol_se)
str(vol_fi)
str(vol_nl)

# Complete Sweeden
species_country_volt <- species_country_volt %>%
  left_join(vol_se, by = "SPECIES", suffix = c("", "_new")) %>%
  mutate(
    SEBMS = ifelse(SEBMS == "" | is.na(SEBMS), SEBMS_new, SEBMS)
  ) %>%
  select(-SEBMS_new)

# Complete Finland
species_country_volt <- species_country_volt %>%
  left_join(vol_fi, by = "SPECIES", suffix = c("", "_new")) %>%
  mutate(
    FIBMS = ifelse(FIBMS == "" | is.na(FIBMS), FIBMS_new, FIBMS)
  ) %>%
  select(-FIBMS_new)

# Complete Netherlands
species_country_volt <- species_country_volt %>%
  left_join(vol_nl, by = "SPECIES", suffix = c("", "_new")) %>%
  mutate(
    NLBMS = ifelse(NLBMS == "" | is.na(NLBMS), NLBMS_new, NLBMS)
  ) %>%
  select(-NLBMS_new)

# Complete Luxembourgh
species_country_volt <- species_country_volt %>%
  mutate(
    LUBMS = ifelse(LUBMS == "" | is.na(LUBMS), "multivoltine", LUBMS)
  )

#### Save ----
write.csv(
  species_country_volt,
  here::here("output", "species_country_voltinism.csv"),
  row.names = FALSE
)
