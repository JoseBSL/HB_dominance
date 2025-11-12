#Manuscript checks
library(dplyr)

filtered_interaction_data = readRDS("Data/filtered_interaction_data_all_cols.rds")

#Quick check of number of species
pol_check = filtered_interaction_data %>% 
  select(pollinator, Pollinator_yes) %>% 
  mutate(Pol_name1 = word(pollinator, 1)) %>% 
  mutate(Pol_name2 = word(pollinator, 2)) %>% 
  filter(!is.na(Pol_name2)) %>%
  filter(!Pollinator_yes == "no") %>% 
  summarise(n_pols = n_distinct(pollinator))

#Quick check of number of sampling events
filtered_interaction_data %>% 
  summarise(n_distinct(id))

#Independent networks
filtered_interaction_data %>% 
  mutate(study_network_id = paste0(Study_id, Network_id)) %>% 
  summarise(n_distinct(study_network_id))



