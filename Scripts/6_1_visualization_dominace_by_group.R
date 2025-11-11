############################################################ #
# Script name: "6_1_visualization_dominance_by_group.R"
# Objective: Visualize the effect of honeybee dominance on
#            SES network metrics (ND, MH, PDI) by group.
#
# Inputs:
#   - Data/fit_SES_nd.rds
#   - Data/fit_SES_mh.rds
#   - Data/fit_SES_pdi.rds
#
# Outputs:
#   - Combined figure showing Dominance × Group effects
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####

library(ggeffects)
library(ggplot2)
library(patchwork)

#### Load data ####
fit_SES_pdi = readRDS("Data/fit_SES_pdi.rds")
fit_SES_nd  = readRDS("Data/fit_SES_nd.rds")
fit_SES_mh  = readRDS("Data/fit_SES_mh.rds")

###################### DOMINANCE X GROUP ################ 
#graphics.off()
#palette("default")

# get marginal effects
gp_mh = ggpredict(fit_SES_mh, terms = c("Dominance[all]", "Group"), type = "fixed")
# plot with custom order
mh<-plot(gp_mh) +
  labs(x = "", 
       y = "SES Morisita–Horn",
       color = "Group",
       main="") +
  ggtitle("")+
  theme_classic()

gp_pdi = ggpredict(fit_SES_pdi, terms = c("Dominance[all]", "Group"), type = "fixed")

# plot with custom order
pdi = plot(gp_pdi) +
  labs(x = "Honeybee dominance (ln-ratio)", 
       y = "SES Paired difference index",
       color = "Group",
       main="") +
  ggtitle("")+
  theme_classic()

gp_nd <- ggpredict(fit_SES_nd, terms = c("Dominance[all]", "Group"), type = "fixed")

# plot with custom order
nd = plot(gp_nd) +
  labs(x = "", 
       y = "SES normalised degree",
       color = "Group",
       main="") +
  ggtitle("")+
  theme_classic()

# One row, equal panel sizes, single legend on the right, automatic A)–D) tags
p = (mh |pdi | nd) +
  plot_layout(ncol = 3, nrow=1, guides = "collect") &
  theme(legend.position = "right")

p + plot_annotation(tag_levels = "A", tag_suffix = ")")

