############################################################ #
# Script name: "4_3_ecological_predictors.R"
# Objective: 
# Input: Data/filtered_interaction_data_all_cols.rds
#
# Output: Data/ 
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####
library(dplyr)

#### Load data ####
filtered_interaction_data = readRDS("Data/filtered_interaction_data_all_cols.rds")

####### Extract taxonomic information per species ############
species_data = filtered_interaction_data %>%
  dplyr::select(Pollinator_order, Pollinator_yes, Pollinator_family, Pollinator_genus, pollinator) %>%
  distinct()
species_data = species_data %>%
  mutate(Pollinator_order2 = if_else(Pollinator_order %in% c("Diptera", "Lepidoptera", "Hymenoptera"), Pollinator_order, "Other"))

######### Compute plant richness per network #############
plants_per_network = filtered_interaction_data %>%
  filter(!is.na(plant)) %>%
  group_by(id) %>%
  summarise(n_plants = n_distinct(plant)) %>%
  arrange(desc(n_plants))
plants_per_network_df<-as.data.frame(plants_per_network)

######### Compute intercations per network #############
degree_interactions = filtered_interaction_data %>%
  dplyr::mutate(interaction = as.numeric(interaction)) %>%
  dplyr::group_by(id, pollinator) %>%
  dplyr::summarise(
    sp_interactions = sum(interaction, na.rm = TRUE),
    sp_degree       = dplyr::n_distinct(plant[interaction > 0]),
    .groups = "drop")

############## Modify habitat categories ###########
habitat_net = filtered_interaction_data %>%
  group_by(id) %>%
  summarise(habitat = first(EuPPollNet_habitat))
#Check
unique(habitat_net$habitat)

habitat_net = habitat_net %>%
  mutate(habitat_type = case_when(
    habitat %in% c(
      "Forest/woodland",
      "Riparian vegetation",
      "Semi-natural grasslands",
      "Moors and heathland",
      "Sclerophyllous vegetation",
      "Montane to alpine grasslands",
      "Beaches, dunes, sands",
      "Riparian"
    ) ~ "Natural",
    habitat %in% c(
      "Ruderal vegetation",
      "Intensive grasslands",
      "Agricultural land",
      "Agricultural margins",
      "Green urban areas"
    ) ~ "Anthropic",
    TRUE ~ NA_character_  # in case there's any unclassified habitat
  ))

#### Combine all predictors into one list and save ####
ecological_predictors <- list(
  species_data = species_data,
  plants_per_network = plants_per_network_df,
  degree_interactions = degree_interactions,
  habitat_net = habitat_net)

#### Save data ####
saveRDS(ecological_predictors, "Data/ecological_predictors.rds")
