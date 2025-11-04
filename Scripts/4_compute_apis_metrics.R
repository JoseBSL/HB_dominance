############################################################ #
# Script name: "4_compute_apis_metrics.R"
# Objective: Compute Apis mellifera network-level metrics
#           Total number of interactions, Degree (number of plant partners)
#           Interaction evenness (Shannon-based), Simpson dominance index
#
# Input:  Data/dt_simple.rds
# Output: apis_summary (data frame in environment)
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####
library(dplyr)
library(purrr)

#### Load data ####
dt_simple = readRDS("Data/dt_simple.rds")

#### Compute total Apis interactions and degree per network ####
apis_interactions = dt_simple %>%
  filter(pollinator == "Apis mellifera") %>%
  group_by(id) %>%
  summarise(
    Apis_interactions = sum(interaction, na.rm = TRUE),
    Apis_degree = n_distinct(plant))

#### Compute Apis interaction evenness per network ####
apis_interaction_evenness <- dt_simple %>%
  filter(pollinator == "Apis mellifera") %>%
  group_by(id, pollinator, plant) %>%
  summarise(w_ij = sum(interaction), .groups = "drop") %>%
  group_by(id, pollinator) %>%
  mutate(
    total_w = sum(w_ij),
    p_ij = w_ij / total_w,
    log_k = log(n_distinct(plant))  # log of partner count
  ) %>%
  summarise(
    H = -sum(p_ij * log(p_ij)),  # Shannon entropy
    log_k = first(log_k),
    Apis_interaction_evenness = ifelse(log_k > 0, H / log_k, NA),
    .groups = "drop")


#### Compute Apis Simpson dominance per network ####
apis_simpson = dt_simple %>%
  filter(pollinator == "Apis mellifera") %>%
  group_by(id, plant) %>%
  summarise(n = sum(interaction), 
            .groups = "drop") %>%
  group_by(id) %>%
  mutate(p = n / sum(n)) %>%
  summarise(Apis_simpson = sum(p^2))

#### Compute total number of all interactions per network and total number of pollinator species ####
all_interactions = dt_simple %>%
  group_by(id) %>%
  summarise(
    All_interactions = sum(interaction, 
                           na.rm = TRUE),
    pollinator_n = n_distinct(pollinator),
    plant_n = n_distinct(plant))

#### Join Apis interactions, Apis interaction evenness, and all interactions ####
apis_summary = list(
  apis_interactions,
  all_interactions,
  apis_simpson,
  apis_interaction_evenness) %>%
  reduce(left_join, by = "id") %>%
  mutate(
    Dominance = log(Apis_interactions / (All_interactions - Apis_interactions)),
    Apis_ndegree = Apis_degree / plant_n,
    Apis_1_D = 1 - Apis_simpson)

#### Save data ####
saveRDS(apis_summary, "Data/apis_.rds")
