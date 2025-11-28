############################################################ #
# Script name: "6_4_visualization_rank_abundance.R"
# Objective: Visualize rank-abundance distribution of pollinator species
#            Only names of top 50 are shown
# Inputs: Data/filtered_interaction_data.rds
#
# Outputs: Rank-abundance plot with top 50 species labeled
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####
library(ggrepel)
library(ggplot2)
library(dplyr)
library(stringr)

#### Load data ####
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





