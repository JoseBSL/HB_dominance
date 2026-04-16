############################################################
# Script name: "6_1_visualization_dominance.R"
# Objective: Visualize the effect of honeybee dominance on
#            SES network metrics (ND, MH, PDI) - all pollinators pooled.
#
# Inputs:
#   - Data/fit_SES_nd_nogroup.rds
#   - Data/fit_SES_mh_nogroup.rds
#   - Data/fit_SES_pdi_nogroup.rds
#
# Outputs:
#   - Combined figure showing Dominance effects (all groups pooled)
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################

library(ggeffects)
library(ggplot2)
library(patchwork)
library(glmmTMB)

#### Load data ####
fit_SES_pdi <- readRDS("Data/fit_SES_pdi_nogroup.rds")
fit_SES_nd  <- readRDS("Data/fit_SES_nd_nogroup.rds")
fit_SES_mh  <- readRDS("Data/fit_SES_mh_nogroup.rds")

# Base theme
base_theme <- theme_bw(base_size = 7.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# Single colour for pooled effect
line_col <- "#0072B2"

# Marginal effects
gp_mh  <- ggpredict(fit_SES_mh,  terms = "Dominance[all]", type = "fixed")
gp_pdi <- ggpredict(fit_SES_pdi, terms = "Dominance[all]", type = "fixed")
gp_nd  <- ggpredict(fit_SES_nd,  terms = "Dominance[all]", type = "fixed")

# Plots
mh <- plot(gp_mh) +
  scale_color_manual(values = c(line_col)) +
  scale_fill_manual(values = c(line_col)) +
  labs(
    x = "",
    y = "SES of Resource Overlap"
  ) +
  ggtitle("")

pdi <- plot(gp_pdi) +
  scale_color_manual(values = c(line_col)) +
  scale_fill_manual(values = c(line_col)) +
  labs(
    x = "Honeybee Dominance (ln-ratio)",
    y = "SES of Specialisation"
  ) +
  ggtitle("")

nd <- plot(gp_nd) +
  scale_color_manual(values = c(line_col)) +
  scale_fill_manual(values = c(line_col)) +
  labs(
    x = "",
    y = "SES of Normalised Degree"
  ) +
  ggtitle("")

# Helper: change ribbon alpha only
set_ribbon_alpha <- function(p, alpha = 0.08) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomRibbon")) {
      p$layers[[i]]$aes_params$alpha <- alpha
    }
  }
  p
}

# Helper: change line width only
set_linewidth <- function(p, lw = 0.7) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomLine")) {
      p$layers[[i]]$aes_params$linewidth <- lw
    }
  }
  p
}

# Helper: ensure ribbons are behind lines
bring_lines_after <- function(p) {
  ribbon_idx <- which(sapply(p$layers, function(x) inherits(x$geom, "GeomRibbon")))
  line_idx   <- which(sapply(p$layers, function(x) inherits(x$geom, "GeomLine")))
  other_idx  <- setdiff(seq_along(p$layers), c(ribbon_idx, line_idx))
  
  p$layers <- c(p$layers[ribbon_idx], p$layers[line_idx], p$layers[other_idx])
  p
}

mh  <- set_linewidth(mh,  0.9)
pdi <- set_linewidth(pdi, 0.9)
nd  <- set_linewidth(nd,  0.9)

mh  <- set_ribbon_alpha(mh,  alpha = 0.07)
pdi <- set_ribbon_alpha(pdi, alpha = 0.07)
nd  <- set_ribbon_alpha(nd,  alpha = 0.07)

mh  <- bring_lines_after(mh)
pdi <- bring_lines_after(pdi)
nd  <- bring_lines_after(nd)

# Remove legends completely
mh  <- mh  + theme(legend.position = "none")
pdi <- pdi + theme(legend.position = "none")
nd  <- nd  + theme(legend.position = "none")

# Combine and apply theme
p <- (mh | pdi | nd) +
  plot_layout(ncol = 3, nrow = 1) &
  base_theme &
  theme(
    panel.spacing = grid::unit(0.1, "lines"),
    plot.margin = grid::unit(c(0, 2.5, 2, 0), "pt"),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0.35)
  )

p + plot_annotation(tag_levels = "A", tag_suffix = ")")