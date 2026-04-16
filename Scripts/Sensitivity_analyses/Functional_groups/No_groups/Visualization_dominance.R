############################################################ #
# Script name: "6_1_visualization_dominance.R"
# Objective: Visualize the effect of honeybee dominance on
#            SES network metrics (ND, MH, PDI) - all pollinators pooled.
#
# Inputs:
#   - Data/fit_SES_nd.rds
#   - Data/fit_SES_mh.rds
#   - Data/fit_SES_pdi.rds
#
# Outputs:
#   - Combined figure showing Dominance effects (all groups pooled)
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
fit_SES_pdi = readRDS("Data/fit_SES_pdi_nogroup.rds")
fit_SES_nd  = readRDS("Data/fit_SES_nd_nogroup.rds")
fit_SES_mh  = readRDS("Data/fit_SES_mh_nogroup.rds")

###################### DOMINANCE (POOLED) ######################

# Marginal effects - single term, no Group
gp_mh  <- ggpredict(fit_SES_mh,  terms = "Dominance[all]", type = "fixed")
gp_pdi <- ggpredict(fit_SES_pdi, terms = "Dominance[all]", type = "fixed")
gp_nd  <- ggpredict(fit_SES_nd,  terms = "Dominance[all]", type = "fixed")

# Plots
mh <- plot(gp_mh) +
  labs(x = "",
       y = "SES Morisita-Horn") +
  ggtitle("") +
  theme_classic()

pdi <- plot(gp_pdi) +
  labs(x = "Honeybee dominance (ln-ratio)",
       y = "SES Paired difference index") +
  ggtitle("") +
  theme_classic()

nd <- plot(gp_nd) +
  labs(x = "",
       y = "SES Normalised degree") +
  ggtitle("") +
  theme_classic()

# Combine
p <- (mh | pdi | nd) +
  plot_layout(ncol = 3, nrow = 1)

p + plot_annotation(tag_levels = "A", tag_suffix = ")")
