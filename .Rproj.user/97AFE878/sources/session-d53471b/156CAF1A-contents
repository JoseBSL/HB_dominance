#Function to export tables

# helper to format p-values
fmt_p = function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# return a formatted data.frame for one component
tab_df = function(mod, component, label){
  df = as.data.frame(parameters::model_parameters(
    mod, component = component, effects = "fixed", ci_method = "wald"
  ))
  
  # (Your names already match, but keep these guards)
  if (!"Coefficient" %in% names(df) && "Estimate" %in% names(df)) df$Coefficient <- df$Estimate
  if (!"SE" %in% names(df) && "Std. Error" %in% names(df))       df$SE <- df[["Std. Error"]]
  if (!"z" %in% names(df) && "Statistic" %in% names(df))         df$z  <- df$Statistic
  
  df |>
    dplyr::mutate(
      Component   = label,
      Coefficient = round(Coefficient, 3),
      SE          = round(SE, 3),
      z           = round(z, 2),
      p           = fmt_p(p)
    ) |>
    dplyr::select(Component, Parameter, Coefficient, SE, CI, z, p)
}

#### Helper to build a formatted flextable ####
make_ft <- function(df_all, caption = NULL) {
  idx <- which(!duplicated(df_all$Component))
  idx <- idx[idx > 1]  # skip first row
  
  ft <- flextable(df_all) |>
    merge_v(j = "Component", part = "body") |>
    valign(j = "Component", valign = "top", part = "body") |>
    align(j = c("Component", "Parameter"), align = "left", part = "body") |>
    autofit()
  
  if (!is.null(caption)) ft <- set_caption(ft, caption)
  
  if (length(idx)) {
    ft <- hline(
      ft,
      i = idx, j = seq_len(ncol(df_all)),
      border = fp_border(color = "gray50", width = 0.75),
      part = "body"
    )
  }
  ft
}

