library(ggeffects)
library(ggplot2)
library(patchwork)
library(glmmTMB)
library(scales)

fit_SES_pdi <- readRDS("Data/fit_SES_pdi_groups.rds")
fit_SES_nd  <- readRDS("Data/fit_SES_nd_groups.rds")
fit_SES_mh  <- readRDS("Data/fit_SES_mh_groups.rds")

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

# Marginal effects
gp_mh  <- ggpredict(fit_SES_mh,  terms = c("Dominance[all]", "Group"), type = "fixed")
gp_pdi <- ggpredict(fit_SES_pdi, terms = c("Dominance[all]", "Group"), type = "fixed")
gp_nd  <- ggpredict(fit_SES_nd,  terms = c("Dominance[all]", "Group"), type = "fixed")

# Common ordered groups across all predictions
all_groups <- union(union(unique(gp_mh$group), unique(gp_pdi$group)), unique(gp_nd$group))
all_groups <- sort(all_groups)
n_groups <- length(all_groups)

# Palette for many groups
my_cols <- setNames(
  scales::hue_pal(h = c(15, 375), c = 100, l = 55)(n_groups),
  all_groups
)

# Linetypes
base_lty <- c("solid", "dashed", "dotdash", "longdash", "twodash", "dotted")
my_lty <- setNames(rep(base_lty, length.out = n_groups), all_groups)

# Helper to standardize group order
standardize_groups <- function(gp) {
  gp$group <- factor(gp$group, levels = all_groups)
  gp
}

gp_mh  <- standardize_groups(gp_mh)
gp_pdi <- standardize_groups(gp_pdi)
gp_nd  <- standardize_groups(gp_nd)

# Plots
mh <- plot(gp_mh) +
  aes(linetype = group) +
  scale_color_manual(values = my_cols, breaks = all_groups, drop = FALSE) +
  scale_fill_manual(values = my_cols, breaks = all_groups, drop = FALSE) +
  scale_linetype_manual(values = my_lty, breaks = all_groups, drop = FALSE) +
  labs(
    x = "",
    y = "SES of Resource Overlap",
    color = "Group"
  ) +
  guides(
    fill = "none",
    linetype = "none",
    colour = guide_legend(
      order = 1,
      override.aes = list(
        fill = NA,
        linetype = unname(my_lty[all_groups]),
        linewidth = 0.9
      )
    )
  ) +
  ggtitle("")

pdi <- plot(gp_pdi) +
  aes(linetype = group) +
  scale_color_manual(values = my_cols, breaks = all_groups, drop = FALSE) +
  scale_fill_manual(values = my_cols, breaks = all_groups, drop = FALSE) +
  scale_linetype_manual(values = my_lty, breaks = all_groups, drop = FALSE) +
  labs(
    x = "Honeybee Dominance (ln-ratio)",
    y = "SES of Specialisation",
    color = "Group"
  ) +
  guides(
    fill = "none",
    linetype = "none",
    colour = guide_legend(
      order = 1,
      override.aes = list(
        fill = NA,
        linetype = unname(my_lty[all_groups]),
        linewidth = 0.9
      )
    )
  ) +
  ggtitle("")

nd <- plot(gp_nd) +
  aes(linetype = group) +
  scale_color_manual(values = my_cols, breaks = all_groups, drop = FALSE) +
  scale_fill_manual(values = my_cols, breaks = all_groups, drop = FALSE) +
  scale_linetype_manual(values = my_lty, breaks = all_groups, drop = FALSE) +
  labs(
    x = "",
    y = "SES of Normalised Degree",
    color = "Group"
  ) +
  guides(
    fill = "none",
    linetype = "none",
    colour = guide_legend(
      order = 1,
      override.aes = list(
        fill = NA,
        linetype = unname(my_lty[all_groups]),
        linewidth = 0.9
      )
    )
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

set_linewidth <- function(p, lw = 0.7) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomLine")) {
      p$layers[[i]]$aes_params$linewidth <- lw
    }
  }
  p
}

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

# Combine and apply theme
p <- (mh | pdi | nd) +
  plot_layout(ncol = 3, nrow = 1, guides = "collect") &
  base_theme &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.spacing = grid::unit(0.1, "lines"),
    legend.margin = margin(t = -1, b = 0, unit = "pt"),
    plot.margin = grid::unit(c(0, 2.5, 2, 0), "pt"),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0.35),
    legend.key.width = grid::unit(1.25, "cm")
  )

p + plot_annotation(tag_levels = "A", tag_suffix = ")")