############################################################ #
# Script name: "6_2_visualization_dominance_by_plant.R"
# Objective: Visualize the effect of honeybee dominance on 
#             SES metrics across different plant species
# Inputs:
#   - Data/fit_SES_nd_over1.rds
#   - Data/fit_SES_mh_over1.rds
#   - Data/fit_SES_pdi_over1.rds
#
# Outputs:
#   - Combined figure showing Dominance × plant effects
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####
library(ggeffects)
library(ggplot2)

#### Load data ####
fit_SES_pdi = readRDS("Data/fit_SES_pdi_over1.rds")
fit_SES_nd  = readRDS("Data/fit_SES_nd_over1.rds")
fit_SES_mh  = readRDS("Data/fit_SES_mh_over1.rds")

################## DOMINANCE X PLANTS
graphics.off()
palette("default")

# get marginal effects
gpp_mh <- ggpredict(fit_SES_mh, terms = c("Dominance[all]", "Plant"), type = "fixed")
# plot with custom order
mhp<-plot(gpp_mh) +
  labs(x = "", 
       y = "SES Morisita–Horn",
       color = "Plant",
       main="") +
  ggtitle("")+
  theme_classic()


gpp_pdi <- ggpredict(fit_SES_pdi, terms = c("Dominance[all]", "Plant"), type = "fixed")

# plot with custom order
pdip<-plot(gpp_pdi) +
  labs(x = "Honeybee dominance (ln-ratio)", 
       y = "SES Paired difference index",
       color = "Plant",
       main="") +
  ggtitle("")+
  theme_classic()

gpp_nd <- ggpredict(fit_SES_nd, terms = c("Dominance[all]", "Plant"), type = "fixed")

# plot with custom order
ndp<-plot(gpp_nd) +
  labs(x = "", 
       y = "SES normalised degree",
       color = "Plant",
       main="") +
  ggtitle("")+
  theme_classic()

library(patchwork)




# One row, equal panel sizes, single legend on the right, automatic A)–D) tags
p2 <- (mhp | pdip | ndp) +
  plot_layout(ncol = 3, nrow=1) 

p2 + plot_annotation(tag_levels = "A", tag_suffix = ")")


p3 <- (mhp |pdip | ndp) +
  plot_layout(ncol = 3, nrow=1, guides = "collect") &
  theme(legend.position = "right")

p3 + plot_annotation(tag_levels = "A", tag_suffix = ")")

