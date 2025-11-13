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

#Percentage of plants visited by the HB
all_plants = filtered_interaction_data %>% 
  summarise(n_distinct(plant))

apis_plants = filtered_interaction_data %>% 
  filter(pollinator == "Apis mellifera") %>% 
  summarise(n_distinct(plant))

apis_plants/all_plants * 100


#Percentage of plants visited by the HB

colnames(filtered_interaction_data)

hb_percent = filtered_interaction_data %>% 
  mutate(Study_network_id = paste0(Study_id, Network_id)) %>% 
  mutate(group = if_else(grepl("Apis", pollinator, ignore.case = TRUE), "Honeybee", "Other")) %>% 
  group_by(Study_network_id, group) %>% 
  summarise(
    n_interactions = sum(interaction, na.rm = TRUE),   # sum of interactions per group
    .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = n_interactions,
    values_fill = 0) %>%
  mutate(
    total = Honeybee + Other,
    perc_honeybee = 100 * Honeybee / total)
  
max(hb_percent$perc_honeybee)
min(hb_percent$perc_honeybee)
mean(hb_percent$perc_honeybee)

ggplot(hb_percent) +
  geom_histogram(aes(perc_honeybee))

#Percentage of plants visited by the HB
# Calculate the mean
mean_hb <- mean(hb_percent$perc_honeybee, na.rm = TRUE)

# Plot histogram + mean line
ggplot(hb_percent, aes(x = perc_honeybee)) +
  geom_histogram(
    bins = 30,                     # adjust bin number if needed
    fill = "skyblue", 
    color = "white"
  ) +
  geom_vline(
    xintercept = mean_hb,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = mean_hb,
    y = 1,
    label = paste0("Mean = ", round(mean_hb, 1), "%"),
    vjust = -0.5,
    hjust = -0.1,
    color = "red",
    size = 4
  ) +
  labs(
    x = "Honeybee interactions (%) per network",
    y = "Number of networks",
    title = "Distribution of honeybee interaction dominance across networks"
  ) +
  theme_minimal(base_size = 12)


hb_exclusive_global <- filtered_interaction_data %>%
  mutate(group = if_else(grepl("Apis", pollinator, ignore.case = TRUE),
                         "Honeybee", "Other")) %>%
  # aggregate across the WHOLE dataset (no network grouping)
  group_by(plant, group) %>%
  summarise(n_interactions = sum(interaction, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n_interactions, values_fill = 0) %>%
  mutate(exclusive_hb  = (Honeybee > 0 & Other == 0))

# Summary numbers
hb_exclusive_summary <- hb_exclusive_global %>%
  summarise(n_plants_hb_exclusive = sum(exclusive_hb),
    perc_exclusive_hb     = 100 * n_plants_hb_exclusive / n_plants_total)

hb_exclusive_summary
