############################################################ #
# Script name: "6_3_visualization_dominance_by_plant.R"
# Objective: 
#
# Inputs:
#   - Data/fit_SES_nd.rds
#   - Data/fit_SES_mh.rds
#   - Data/fit_SES_pdi.rds
#
# Outputs:
#   - Combined figure showing Dominance × plant effects
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####
library(GGally)
library(ggplot2)

#### Load data ####
apis_network_metrics = readRDS("Data/apis_network_metrics.rds")


######## FIG. 1 COVARIATION BETWEEN DOMINANCE, ND, N POLL AND N PLANT ####

apis_summary2_no1 = apis_network_metrics%>%
  filter(Apis_interactions>1)

vars <- c("plant_n","pollinator_n","Apis_ndegree", "Dominance")
df_subset <- apis_summary2_no1[, vars, drop = FALSE]

# Custom lower panel: GAM + linear
lower_gam_lm <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "gray40", size = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed", ...) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "red", ...)
}

ggpairs(
  df_subset,
  columns = c("plant_n","pollinator_n","Apis_ndegree",  "Dominance"),
  columnLabels = c("Plant richness", "Pollinator richness", "Hb normalized degree",  "Hb dominance"),
  lower = list(continuous = wrap(lower_gam_lm)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5, fill = "lightblue")),
  upper = list(continuous = wrap("cor", size = 3))
  #diag = list(continuous = wrap(custom_density_diag))
)

#####################  FIGURE RANK-ABUNDANCE  ##########

library(ggrepel)

df = readRDS("Data/filtered_interaction_data.rds")

# Frequency table
pollinator_freq <- df %>%
  distinct(id, pollinator) %>%
  count(pollinator, name = "networks_present_in") %>%
  arrange(desc(networks_present_in)) %>%
  mutate(rank = row_number())

pollinator_freq <- pollinator_freq %>%
  mutate(
    label = str_trim(pollinator),                         # remove leading/trailing spaces
    label = str_split(label, "\\s+"),                     # split into words
    label = sapply(label, function(words) {
      paste0(substr(words, 1, 3), collapse = " ")         # take first 3 letters of each word
    })
  )

# Add labels only for top 20
pollinator_freq <- pollinator_freq %>%
  mutate(label = ifelse(rank <= 50, label, NA))

# Plot
ggplot(pollinator_freq, aes(x = rank, y = log(networks_present_in))) +
  geom_point(col="grey55") +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf,         # Allow all labels to be shown
    box.padding = 0.5,          # Space around labels
    point.padding = 0.5,        # Space between label and point
    nudge_x = 50,             # Push labels upward
    segment.size = 0.2,         # Line thickness from point to label
    segment.size = 0.2,
    force_pull = 2,
    segment.color = "grey80",     # light grey color
    segment.alpha = 0.5           # Line transparency
  ) +
  scale_y_continuous(
    name = "Species count in the networks",
    breaks = log(c(1, 2, 10, 100, 500, 1500)),  # choose meaningful original values
    labels = c(1, 2, 10, 100, 500, 1500)        # show original scale on axis
  ) +
  labs(
    x = "Pollinator species rank",
    y = "Log(Frequency in networks)",
    title = ""
  ) +
  theme_minimal()





