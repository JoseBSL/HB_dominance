
#### Load libraries ####
library(data.table)

#### Load data ####
# Version: 1.3
# DOI: 10.5281/zenodo.15183272
lanuza = readRDS("Interaction_data.rds")
# Summarize the total number of interactions
sum(lanuza$Interaction)


#### Data preparation ####
# Remove all network IDs that do not meet the inclusion criteria
# Ensure that is a datatable format
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

#Total number of networks (all)
uniqueN(lanuza$id)
#Number of networks that contain Apis mellifera
uniqueN(df_withApis$id)

#Total interactions within networks that contain Apis mellifera
df_withApis[, sum(Interaction, na.rm = TRUE)]

#Define the bee families
bee_families = c("Apidae", "Halictidae", "Andrenidae", "Colletidae", "Megachilidae", "Melittidae")

