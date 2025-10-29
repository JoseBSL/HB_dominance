#### Load libraries ####
library(data.table)
library(dplyr)

# COMPUTE ALL METRICS ####
net_long_no012 = readRDS("Data/net_long_no0.rds")

#Create meaningful ids
net_long_no0 = net_long_no0  %>%
  dplyr::rename(null_id=type)
net_long_no0[, Type := ifelse(null_id == "observed", "observed", "null")]
net_long_no0[, final_id := ifelse(Type == "observed", id, null_id)]

###################### loop all the metrics with data.table 16 august
setDT(net_long_no0)

metrics_all_dt <- net_long_no0[, {
  res <- compute_metrics5(.SD,
                          focal_species   = "Apis mellifera",
                          return_observed = TRUE)
  if (is.null(res)) NULL else as.data.table(res)
}, by = .(final_id, id, Type)]

saveRDS(metrics_all_dt, "Data/all_metrics_500.rds")

################  Z-score ##### 

# columns to compute SES
metric_cols <- c(
  "pianka","morisita_horn","pref_pianka","pref_morisita_horn",
  "pac_from_apis", "norm_degree", "d_prime_no_apis_finite", "d_prime", "pdi"
)

# 1) Compute metrics for every network (final_id) and keep Type + id
#metrics_all <- net_long_no0 %>%
#  dplyr::group_by(id, final_id, Type) %>%
#  dplyr::group_modify(~{
#    res <- compute_metrics5(
#      df = .x,
#      focal_species   = "Apis mellifera",
#      return_observed = TRUE
#    )
#    if (is.null(res)) tibble::tibble() else tibble::as_tibble(res)
# }) %>%
#  dplyr::ungroup()

#################################################### COMPUTE SES (Z-SCORES) ####

# 1) Summary of nulls: mean, sd, IQR (+ optional: number of distinct values)
null_stats <- metrics_all_dt %>%
  dplyr::filter(Type == "null") %>%
  dplyr::group_by(id, pollinator) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(metric_cols),
      list(
        mu   = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x,   na.rm = TRUE),
        iqr  = ~IQR(.x,  na.rm = TRUE),
        uniq = ~dplyr::n_distinct(.x[is.finite(.x)])
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# 2 observed values
obs_vals <- metrics_all_dt %>%
  dplyr::filter(Type == "observed") %>%
  dplyr::select(id, final_id, pollinator, dplyr::all_of(metric_cols))

# join null summaries
ses_tbl <- dplyr::left_join(obs_vals, null_stats, by = c("id","pollinator"))

# 3) compute SES = (obs - mu_null) / sd_null for each metric (robust to dplyr versions)
for (m in metric_cols) {
  x  <- ses_tbl[[m]]
  mu <- ses_tbl[[paste0(m, "_mu")]]
  sd <- ses_tbl[[paste0(m, "_sd")]]
  z  <- (x - mu) / sd
  ses_tbl[[paste0("SES_", m)]] <- ifelse(is.finite(z), z, NA_real_)  # handles sd=0 / NA
}
dplyr::glimpse(ses_tbl)
names(ses_tbl)

saveRDS(ses_tbl, "Data/ses_tbl.rds")
