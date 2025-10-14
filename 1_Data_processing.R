#### Load libraries ####
library(data.table)
library(dplyr)

#### Load data ####
# Version: 1.3
# DOI: 10.5281/zenodo.15183272
lanuza = readRDS("Data/Interaction_data.rds")
# Summarize the total number of interactions
sum(lanuza$Interaction)


#### Data preparation ####
# Remove all network IDs that do not meet the inclusion criteria
# Ensure that is in datatable format
setDT(lanuza)
#Create a study x network x date ID
lanuza$id = as.factor(
  paste(lanuza$Study_id, 
        lanuza$Network_id, 
        lanuza$Date, 
        sep = "//")) # id used for further analyses
#Keep only network IDs that include Apis mellifera
ids_with_apis = lanuza[na.omit(Pollinator_accepted_name) == "Apis mellifera", unique(id)]
#Subset the full dataset to only those networks
df_withApis = lanuza[id %chin% ids_with_apis]


#### Checks ####
#Total number of networks (all)
uniqueN(lanuza$id)
#Number of networks that contain Apis mellifera
uniqueN(df_withApis$id)
#Total interactions within networks that contain Apis mellifera
df_withApis[, sum(Interaction, na.rm = TRUE)]
#Define the bee families
bee_families = c("Apidae", "Halictidae", "Andrenidae", "Colletidae", "Megachilidae", "Melittidae")
#Create a Pollinator_yes column based on criteria of poll. functional groups
df_withApis = df_withApis %>%
    mutate(
    Pollinator_yes = case_when(
      Pollinator_order == "Coleoptera" ~ "Coleo.",
      Pollinator_order == "Diptera" & Pollinator_family == "Syrphidae" ~ "Syrph.",
      Pollinator_order == "Diptera" & Pollinator_family != "Syrphidae" ~ "Non-syr.",
      Pollinator_order == "Hymenoptera" & Pollinator_family %in% bee_families ~ "Bees",
      Pollinator_order == "Hymenoptera" & !(Pollinator_family %in% bee_families) ~ "Hymen.",
      Pollinator_order == "Lepidoptera" ~ "Lepid.",
      TRUE ~ "no"))
#Sum the number of interactions for each plant–pollinator pair within each site and date
#Only for networks containing Apis
dfs = df_withApis %>%
  group_by(Bioregion, 
           Country, 
           Study_id, 
           EuPPollNet_habitat, 
           Network_id, 
           Date, 
           id,
           Pollinator_yes,
           Pollinator_order, 
           Pollinator_family, 
           Pollinator_genus, 
           Plant_accepted_name, 
           Pollinator_accepted_name) %>%
  summarise(Interactions = sum(Interaction, na.rm = TRUE), .groups = "drop")
#Rename variables
dfs = dfs %>%
  rename(pollinator = Pollinator_accepted_name) %>% 
  rename(plant = Plant_accepted_name) %>% 
  rename(interaction = Interactions)
#Convert to datatable
DT = as.data.table(dfs)

#### Filter data ####
#Step 1: High-interaction networks
#Step 2: Single-plant networks
#Step 3: Single-pollinator networks

#Step 1: Remove ids with interaction > 1000
ids_high_interaction = DT[interaction > 1000, 
                          unique(id)]
#Step 2: Remove ids with only one plant
ids_one_plant = DT[, .(n_plant = uniqueN(plant)), 
                   by = id][n_plant == 1, 
                            id]
#Step 3: Remove ids with only one pollinator
ids_one_pollinator = DT[, .(n_pollinator = uniqueN(pollinator)), 
                        by = id][n_pollinator == 1, 
                                 id]
#Combine all sets of problematic ids
ids_to_remove = unique(c(ids_high_interaction, 
                         ids_one_plant, 
                         ids_one_pollinator))
#Filter them out
DT_filtered = DT[!id %in% ids_to_remove] ########### WORKING DATASET 
sum(DT_filtered$interaction)
sum(DT$interaction)
#Keep only required columns
dt_simple = unique(DT_filtered [, .(plant, pollinator, interaction, id)])
dt_simple[, `:=`(plant = as.character(plant),
                 pollinator = as.character(pollinator),
                 id = as.character(id),
                 interaction = as.numeric(interaction))]

#### Save data ####
saveRDS(dt_simple, "Data/dt_simple.rds")
