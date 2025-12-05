############################################################
# 3 scenarios of plant–pollinator networks
# with different Apis dominance
# + network plots (Poll1 & Poll2 coloured, Apis black)
# + single summary barplot (Degree, Morisita, PDI)
############################################################

## Load libraries ----
library(dplyr)
library(ggplot2)
library(viridis)
library(data.table)
library(patchwork)
library(tidyr)
library(grid)

# Function with compute_metrics() should be here
# (must create columns: degree, pdi, morisita_horn at least)
source("Scripts/compute_metric_function.R")

############################################################
## 1. Define 3 interaction matrices (scenarios) ----------
# Rows = plants (Plant1–Plant4)
# Cols = pollinators (Apis, Poll1, Poll2)
############################################################

# LOW Apis dominance (Apis relatively weak)
M_low <- matrix(
  c(1,   0, 0,   # Plant1
    1,   2, 0,   # Plant2
    1,   2, 1,   # Plant3
    0,   0, 2),  # Plant4
  nrow = 4,
  byrow = TRUE
)

# MEDIUM Apis dominance
M_med <- matrix(
  c(2,   0, 0,
    2,   1, 0,
    1.5, 2, 1,
    0,   1, 2),
  nrow = 4,
  byrow = TRUE
)

# HIGH Apis dominance (Apis much stronger)
M_high <- matrix(
  c(5,   0,   0,
    5,   1,   0.5,
    3,   1.5, 0.5,
    1,   1.5, 1.5),
  nrow = 4,
  byrow = TRUE
)

# Give all matrices the same row/col names
set_dimnames <- function(M) {
  colnames(M) <- c("Apis", "Poll1", "Poll2")
  rownames(M) <- c("Plant1", "Plant2", "Plant3", "Plant4")
  M
}

M_low  <- set_dimnames(M_low)
M_med  <- set_dimnames(M_med)
M_high <- set_dimnames(M_high)

############################################################
## 2. Helper function to build network plot + metrics -----
############################################################

build_scenario <- function(M, title) {
  # ---- edge list for plotting ----
  edges <- as.data.frame(as.table(M)) %>%
    dplyr::rename(
      Plant      = Var1,
      Pollinator = Var2,
      Weight     = Freq
    ) %>%
    dplyr::filter(Weight > 0)
  
  # ---- nodes ----
  plants <- data.frame(
    name = rownames(M),
    x    = 0,
    y    = seq(1, nrow(M))
  )
  
  polls <- data.frame(
    name = colnames(M),
    x    = 1,
    y    = seq(1.5, 3.5, length.out = ncol(M))  # centered vertically
  )
  
  nodes <- rbind(plants, polls)
  
  # ---- edges + coordinates ----
  edges_xy <- edges %>%
    dplyr::left_join(nodes, by = c("Plant" = "name")) %>%
    dplyr::rename(x1 = x, y1 = y) %>%
    dplyr::left_join(nodes, by = c("Pollinator" = "name")) %>%
    dplyr::rename(x2 = x, y2 = y)
  
  # ---- define colour groups: Poll1 & Poll2 coloured, Apis black ----
  edges_xy <- edges_xy %>%
    dplyr::mutate(
      color_group = dplyr::case_when(
        Pollinator == "Poll1" ~ "Poll1",
        Pollinator == "Poll2" ~ "Poll2",
        TRUE                  ~ "Apis"   # Apis and others = black
      )
    )
  
  polls <- polls %>%
    dplyr::mutate(
      color_group = dplyr::case_when(
        name == "Poll1" ~ "Poll1",
        name == "Poll2" ~ "Poll2",
        TRUE            ~ "Apis"
      )
    )
  
  # ---- long format for metrics ----
  df_long <- as.data.frame(as.table(M)) %>%
    dplyr::rename(
      plant       = Var1,
      pollinator  = Var2,
      interaction = Freq
    ) %>%
    dplyr::filter(interaction > 0)
  
  # ---- compute metrics (focal = Apis) ----
  metrics <- compute_metrics(
    df              = df_long,
    focal_species   = "Apis",
    return_observed = TRUE
  )
  
  # ---- NETWORK plot ----
  p_net <- ggplot() +
    # edges: only Poll1 & Poll2 coloured; Apis edges black
    geom_segment(
      data = edges_xy,
      aes(x = x1, y = y1,
          xend = x2, yend = y2,
          size = Weight,
          color = color_group),
      alpha = 0.8
    ) +
    # manual colours: Poll1/Poll2 coloured, Apis black
    scale_color_manual(values = c(
      "Poll1" = viridis::viridis(3)[1],
      "Poll2" = viridis::viridis(3)[2],
      "Apis"  = "black"
    )) +
    
    # plant nodes (always black, no aes(color))
    geom_point(
      data = plants,
      aes(x = x, y = y),
      color = "black", size = 5
    ) +
    geom_text(
      data = plants,
      aes(x = x, y = y, label = name),
      vjust = 2,
      hjust = 0.5,
      color = "black"
    ) +
    
    # pollinator nodes (Poll1 & Poll2 coloured, Apis black)
    geom_point(
      data = polls,
      aes(x = x, y = y, color = color_group),
      size = 5
    ) +
    geom_text(
      data = polls,
      aes(x = x, y = y, label = name),
      vjust = -1.4,
      hjust = 0.5,
      color = "black"
    ) +
    
    # scales, theme, layout
    scale_size(range = c(0.5, 3)) +
    guides(colour = "none", size = "none") +  # no legends from network
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin  = margin(5, 20, 5, 20),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
    ) +
    coord_flip() +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    ggtitle(title)
  
  list(plot = p_net, metrics = metrics)
}

############################################################
## 3. Build scenarios: low / medium / high ---------------
############################################################

sc_low  <- build_scenario(M_low,  "Low")
sc_med  <- build_scenario(M_med,  "Medium")
sc_high <- build_scenario(M_high, "High")

############################################################
## 4. Combine metrics from all scenarios -----------------
############################################################

metrics_all <- rbindlist(list(
  cbind(scenario = "low",    sc_low$metrics),
  cbind(scenario = "medium", sc_med$metrics),
  cbind(scenario = "high",   sc_high$metrics)
), use.names = TRUE, fill = TRUE)

metrics_all  # optional: inspect in console

############################################################
## 5. Build single summary barplot (grouped by pollinator)
############################################################

metrics_plot <- metrics_all %>%
  dplyr::filter(pollinator %in% c("Poll1", "Poll2")) %>%
  dplyr::select(scenario, pollinator, degree, pdi, morisita_horn) %>%
  tidyr::pivot_longer(
    cols      = c(degree, morisita_horn, pdi),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    scenario = factor(scenario,
                      levels = c("low", "medium", "high"),
                      labels = c("Low", "Medium", "High")),
    pollinator = factor(pollinator, levels = c("Poll1", "Poll2")),
    metric = factor(metric,
                    levels = c("degree", "morisita_horn", "pdi"),
                    labels = c("Degree", "Morisita-Horn", "PDI")),
    scenario_label = as.character(scenario),
    x_group = interaction(pollinator, scenario_label, sep = "_", drop = TRUE)
  )

# order x-axis: Poll1 Low/Med/High, then Poll2 Low/Med/High
metrics_plot$x_group <- factor(
  metrics_plot$x_group,
  levels = c("Poll1_Low", "Poll1_Medium", "Poll1_High",
             "Poll2_Low", "Poll2_Medium", "Poll2_High")
)

# x-axis labels with pollinator indicated
x_labs2 <- c(
  "Poll1_Low"    = "Low",
  "Poll1_Medium" = "Medium",
  "Poll1_High"   = "High",
  "Poll2_Low"    = "Low",
  "Poll2_Medium" = "Medium",
  "Poll2_High"   = "High"
)

# colours for pollinators
poll_colors <- c(
  "Poll1" = viridis::viridis(3)[1],
  "Poll2" = viridis::viridis(3)[2]
)

p_summary <- ggplot(metrics_plot,
                    aes(x = x_group, y = value,
                        fill = pollinator)) +
  geom_col(width = 0.7, color = "black") +
  scale_fill_manual(values = poll_colors,
                    name   = "Pollinator") +
  scale_x_discrete(labels = x_labs2) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 0, hjust = 0.5),
    strip.background = element_rect(fill = "grey90", colour = "black")
  )

############################################################
## 6. Combine: networks on top, summary barplot below ----
############################################################

networks_row <- sc_low$plot | sc_med$plot | sc_high$plot

final_plot <- networks_row / p_summary +
  plot_layout(heights = c(2, 1))  # networks taller than barplot

final_plot


networks_row



# Create arrow row
arrow_row <- ggplot() +
  annotate("segment",
           x = 0.1, xend = 0.9,
           y = 0.5, yend = 0.5,
           arrow = arrow(type = "closed", length = unit(0.25, "cm")),
           size = 1) +
  annotate("text",
           x = 0.5, y = 0.7,
           label = "Increasing Apis dominance",
           fontface = "bold", size = 5) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 5, 0)
  )

# Combine everything
final_plot <- arrow_row / networks_row  +
  plot_layout(
    heights = c(1, 2, 1),   # adjust arrow : networks : bars
    guides = "collect"
  ) &
  theme(legend.position = "bottom")

final_plot
