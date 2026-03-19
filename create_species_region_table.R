library(dplyr)
library(stringr)

#### Data Import and Preparation ####
# ---
here::here() # Check the current working directory
# ---
phenology_estimates  <- read.csv(here::here("output", "pheno_estimates_allspp.csv"), sep = ",", dec = ".")

df_species_site <- phenology_estimates %>%
  select(SPECIES, SITE_ID) %>%
  distinct() %>%
  left_join(
    species_traits %>%
      transmute(
        SPECIES = gsub("_", " ", Taxon),
        vol_min = Vol_min,
        vol_max = Vol_max
      ) %>%
      filter(SPECIES %in% unique(phenology_estimates$SPECIES)),
    by = "SPECIES"
  )

str(df_species_site)


df_species_site <- df_species_site %>%
  mutate(
    country = str_extract(SITE_ID, "^[^.]+")
  )

unique(df_species_site$country)

df_counts <- df_species_site %>%
  count(SPECIES, country)

library(tidyr)

species_country_table <- df_counts %>%
  pivot_wider(
    names_from = country,
    values_from = n,
    values_fill = 0
  )
