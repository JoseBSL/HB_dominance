# Install if needed:
# install.packages(c("dplyr", "tidyr", "sf", "ggplot2", "giscoR", "scatterpie"))

library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(giscoR)
library(scatterpie)

#-----------------------------
# 1. One coordinate per Study_id (in 3035)
#-----------------------------
#Load filtered data
df = readRDS("Data/filtered_interaction_data.rds")
#Extract study id's
id = df %>%
  separate(id, into = c("Study_id", "Location", "Date"), sep = "//") %>% 
  distinct(Study_id) %>% 
  pull()
#Load raw data
lanuza = readRDS("Data/raw_interaction_data.rds")

coordinates_sf <- lanuza %>%
  filter(Study_id %in% id) %>%
  select(Study_id, Latitude, Longitude) %>%
  distinct() %>%
  group_by(Study_id) %>%
  slice(1) %>%                # keep one coordinate per study
  ungroup() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%  # WGS84
  st_transform(3035)          # match map CRS

#-----------------------------
# 2. Apis vs Other counts per Study_id (wide format)
#-----------------------------

#apis_wide <- df %>%
#  #filter(Study_id %in% id) %>%
#  mutate(
#    is_apis = if_else(
#      grepl("Apis", pollinator, ignore.case = TRUE),
#      "Apis",
#      "Other"
#    )
#  ) %>%
#  count(Study_id, is_apis, name = "n") %>%   # one row per Study_id x group
#  pivot_wider(
#    names_from  = is_apis,
#    values_from = n,
#    values_fill = 0
#  )
#  


apis_wide <- df %>%
  separate(id, into = c("Study_id", "Location", "Date"), sep = "//") %>% 
  # filter(Study_id %in% id) %>%   # keep if needed
  mutate(
    is_apis = if_else(
      grepl("Apis", pollinator, ignore.case = TRUE),
      "Apis",
      "Other"
    )
  ) %>%
  group_by(Study_id, is_apis) %>%
  summarise(
    n_interactions = sum(interaction, na.rm = TRUE),   # sum interactions
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = is_apis,
    values_from = n_interactions,
    values_fill = 0
  )
# apis_wide now has columns: Study_id, Apis, Other

#-----------------------------
# 3. Join coordinates + Apis/Other, extract x/y for pies
#-----------------------------

pie_data <- coordinates_sf %>%
  left_join(apis_wide, by = "Study_id") %>%
  mutate(
    x = st_coordinates(geometry)[, 1],
    y = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry()

# Optional sanity check:
# head(pie_data)

#-----------------------------
# 4. Countries layer (Europe, EPSG:3035)
#-----------------------------

cntries <- gisco_get_countries(region = "Europe", year = "2016") %>%
  st_transform(3035)

#-----------------------------
# 5. Plot: map + pie charts at study locations
#-----------------------------
ggplot() +
  geom_sf(
    data = cntries,
    fill = "grey85",
    color = "grey60",
    size = 0.1
  ) +
  geom_scatterpie(
    data = pie_data,
    aes(x = x, y = y),
    cols = c("Apis", "Other"),
    color = NA,
    pie_scale = 2   # adjust size if needed
  ) +
  coord_sf(
    xlim = c(2200000, 7150000),
    ylim = c(1380000, 5500000),
    expand = FALSE
  ) +
  theme(
    panel.grid.major = element_line(
      color = gray(0.5),
      linetype = "dashed",
      size = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "bottom"
  ) +
  labs(fill = "Pollinator group")


