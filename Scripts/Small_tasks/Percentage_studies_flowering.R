# Task 1: Methods — justify exclusion of floral abundance
#   - Explain in Methods why floral abundance was not used.
#   - Compute and report the % of studies in the final dataset 
#     that include floral abundance data.


#### Load libraries ####
library(dplyr)
library(tidyr)
#### Load data ####
#Load filtered data
df = readRDS("Data/filtered_interaction_data.rds")
#Load raw data
lanuza = readRDS("Data/raw_interaction_data.rds")

#Extract study id's
id = df %>%
  separate(id, into = c("Study_id", "Location", "Date"), sep = "//") %>% 
  distinct(Study_id) %>% 
  pull()

#Filter now studies of raw data and check % with flowering info
lanuza %>%
  filter(Study_id %in% id) %>%
  distinct(Study_id, Flower_data) %>%   # one record per study
  count(Flower_data) %>%
  mutate(percentage = 100 * n / sum(n))
    
#74% of the selected studies 3/4 contain flowering information
#or other way around for the text, we will need to exclude 1/4 of the studies
#From 46 studies