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

vars = c("plant_n","pollinator_n","Apis_ndegree", "Dominance")
df_subset = apis_summary2_no1[, vars, drop = FALSE]

# Custom lower panel: GAM + linear
lower_gam_lm = function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(shape = 16, color = "lightblue", size = 0.8, alpha = 0.9) +  # shape=16 removes border
    geom_smooth(method = "lm", se = FALSE, color = "orange", linetype = "dashed", ...) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "black", ...) +
    scale_x_continuous(n.breaks = 4) +   # limit x-axis ticks
    scale_y_continuous(n.breaks = 4) 
}

p = ggpairs(
  df_subset,
  columns = c("plant_n","pollinator_n","Apis_ndegree",  "Dominance"),
  columnLabels = c("Plant richness", "Pollinator richness", "Hb normalized degree",  "Hb dominance"),
  lower = list(continuous = wrap(lower_gam_lm)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5, fill = "lightblue")),
  upper = list(continuous = wrap("cor", size = 3))
  #diag = list(continuous = wrap(custom_density_diag))
) 

custom_theme = theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

p +
  theme_set(custom_theme)


