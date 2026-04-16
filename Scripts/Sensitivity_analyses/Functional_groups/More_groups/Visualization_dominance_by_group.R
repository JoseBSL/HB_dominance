library(ggeffects)
library(ggplot2)
library(patchwork)

fit_SES_pdi <- readRDS("Data/fit_SES_pdi_groups.rds")
fit_SES_nd  <- readRDS("Data/fit_SES_nd_groups.rds")
fit_SES_mh  <- readRDS("Data/fit_SES_mh_groups.rds")

gp_mh  <- ggpredict(fit_SES_mh,  terms = c("Dominance[all]", "Group"), type = "fixed")
gp_pdi <- ggpredict(fit_SES_pdi, terms = c("Dominance[all]", "Group"), type = "fixed")
gp_nd  <- ggpredict(fit_SES_nd,  terms = c("Dominance[all]", "Group"), type = "fixed")

all_groups <- union(union(unique(gp_mh$group), unique(gp_pdi$group)), unique(gp_nd$group))

pal <- setNames(scales::hue_pal()(length(all_groups)), all_groups)

mh <- plot(gp_mh) +
  labs(
    x = "",
    y = "SES Morisita–Horn",
    color = "Group"
  ) +
  scale_color_manual(values = pal, drop = FALSE) +
  theme_classic()

pdi <- plot(gp_pdi) +
  labs(
    x = "Honeybee dominance (ln-ratio)",
    y = "SES Paired difference index",
    color = "Group"
  ) +
  scale_color_manual(values = pal, drop = FALSE) +
  theme_classic()

nd <- plot(gp_nd) +
  labs(
    x = "",
    y = "SES normalised degree",
    color = "Group"
  ) +
  scale_color_manual(values = pal, drop = FALSE) +
  theme_classic()

p <- (mh | pdi | nd) +
  plot_layout(ncol = 3, nrow = 1, guides = "collect") &
  theme(legend.position = "right")

p + plot_annotation(tag_levels = "A", tag_suffix = ")")