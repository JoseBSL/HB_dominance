############################################################ #
# Script name: "4_2_extract_climatic_data.R"
# Objective: Extract climatic data from https://www.chelsa-climate.org/ 
#
# Input1:  Data/CHELSA_EUR11_tas_norm_1981-2005_V1.1.nc
#         Note: file name renames intact to match the original source
# Input2: Data/data_to_extract_climatic.rds
#Output: Data/temperature.rds 
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####
library(terra)
library(dplyr)
library(sf)
library(exactextractr)

#### Load data ####
data_to_extract_climatic = readRDS("Data/data_to_extract_climatic.rds")

r = rast("Data/CHELSA_EUR11_tas_norm_1981-2005_V1.1.nc")
Temp_summer = (r[[5]]+r[[6]]+r[[7]]+r[[8]])/4 ## average summer temperature (May-Aug)
summary(Temp_summer)

#### Extract climatic data from coordinates ####
#Obtain coordinates per network
network_coords = data_to_extract_climatic %>%
  select(Latitude, Longitude, id) %>%
  distinct()
#Convert to sf with WGS84
points_sf = st_as_sf(network_coords, coords = c("Longitude", "Latitude"), crs = 4326)
#Create 5 km buffer (in meters)
points_buf = st_buffer(points_sf, dist = 5000)
#Extract for all months using a stack
extracted = as.data.frame(exact_extract(Temp_summer, points_buf, 'mean'))
names(extracted)[1] = "Temperature"
#Add site IDs
extracted$id = network_coords$id
#Summarise to have one value per site
temperature_summary = extracted %>%
  group_by(id) %>%
  dplyr::summarise(
    Temperature = mean(Temperature, na.rm = TRUE))

#### Save data ####
saveRDS(temperature_summary, "Data/temperature.rds")
