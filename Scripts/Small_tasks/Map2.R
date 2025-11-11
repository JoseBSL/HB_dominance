library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(giscoR)

#--------------------------------------------------
# 0. Assumptions from your previous code
# df: filtered_interaction_data
# lanuza: raw_interaction_data
# pie_data: has Study_id, Apis, Other, x, y in EPSG:3035
# cntries: Europe polygons in EPSG:3035
#--------------------------------------------------

# If you don't have pie_data yet, quick recap:
# (keep your existing version if it already works)

df <- readRDS("Data/filtered_interaction_data.rds")
lanuza <- readRDS("Data/raw_interaction_data.rds")

id <- df %>%
  separate(id, into = c("Study_id", "Location", "Date"), sep = "//") %>%
  distinct(Study_id) %>%
  pull()

coordinates_sf <- lanuza %>%
  filter(Study_id %in% id) %>%
  select(Study_id, Latitude, Longitude) %>%
  distinct() %>%
  group_by(Study_id) %>%
  slice(1) %>%
  ungroup() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(3035)

apis_wide <- df %>%
  separate(id, into = c("Study_id", "Location", "Date"), sep = "//") %>%
  mutate(
    is_apis = if_else(
      grepl("Apis", pollinator, ignore.case = TRUE),
      "Apis",
      "Other"
    )
  ) %>%
  group_by(Study_id, is_apis) %>%
  summarise(
    n_interactions = sum(interaction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = is_apis,
    values_from = n_interactions,
    values_fill = 0
  )

pie_data <- coordinates_sf %>%
  left_join(apis_wide, by = "Study_id") %>%
  mutate(
    x = st_coordinates(geometry)[, 1],
    y = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry()

# Convert to sf with geometry again for spatial ops
pie_sf <- st_as_sf(pie_data,
                   coords = c("x", "y"),
                   crs = 3035,
                   remove = FALSE)

# Countries layer
cntries <- gisco_get_countries(region = "Europe", year = "2016") %>%
  st_transform(3035)

#--------------------------------------------------
# 1. Build FULL rectangular grid over Europe extent
#    (NO clipping by political borders)
#--------------------------------------------------

# Use your plotting extent
bbox_grid <- st_bbox(c(
  xmin = 2200000,
  xmax = 7150000,
  ymin = 1380000,
  ymax = 5500000
), crs = st_crs(3035))

grid_poly <- st_make_grid(
  st_as_sfc(bbox_grid),
  cellsize = 100000,   # 200 km grid; change as needed
  square = TRUE
)

grid_sf <- st_sf(
  grid_id = seq_along(grid_poly),
  geometry = grid_poly
)

#--------------------------------------------------
# 2. Assign studies (pie points) to grid cells & aggregate
#--------------------------------------------------

# Join each point to the grid cell containing it
join <- st_join(pie_sf, grid_sf, join = st_within)

# Summarise Apis/Other interactions per grid cell
grid_agg <- join %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    Apis  = sum(Apis,  na.rm = TRUE),
    Other = sum(Other, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total = Apis + Other,
    prop_apis = if_else(total > 0, Apis / total, NA_real_)
  )

# Attach aggregated values back to ALL grid cells
grid_full <- grid_sf %>%
  left_join(grid_agg, by = "grid_id")

#--------------------------------------------------
# 3. Plot: GRIDDed map (not cut by borders)
#--------------------------------------------------

ggplot() +

  # country borders on top, just as a reference
  geom_sf(
    data = cntries,
    fill = "grey85",
    color = "black",
    size = 0.1
  ) +
  # full grid, colored where data exist
  geom_sf(
    data = grid_full,
    aes(fill = prop_apis),
    color = NA,
    size = 0.01
  ) +
  coord_sf(
    xlim = c(2200000, 7150000),
    ylim = c(1380000, 5500000),
    expand = FALSE
  ) +
  scale_fill_viridis_c(
    option = "C",
    na.value = NA,
    name = "Apis interactions (%)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "right"
  )

