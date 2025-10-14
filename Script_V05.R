library(gamlss)
library(dplyr)
library(reshape)
library(spaa)
library(tidyr)
library(purrr)
library(tibble)
library(terra)
library(tidyverse)
library(FNN)
library(raster)
library(sf)
library(exactextractr)
library(ncdf4)
library (future)
library(DHARMa)
library(effects)
library(TMB)
library(glmmTMB)
library(gamlss)
library(brms)
library(patchwork)
library(ggplot2)
library(data.table)
library(car)
library(psych)
library(GGally)
library(performance)
library(ggeffects)
library(ggrepel)
library(modelsummary)
library(sjPlot)
library(parameters)
library(flextable)
library(officer)
library(GGally)
library(ggplot2)
library(mgcv)
library(cmdstanr)


########### CUSTOM FUNCTION TO COMPUTE ALL METRICS ##################

compute_metrics5<-function(df,
                           focal_species = "Apis mellifera",
                           return_observed = TRUE,
                           eps = 1e-12) {
  
  ## ---------- helpers (mirror bipartite::dfun for d′ finite d_min) ----------
  .d_row_kl <- function(x, q) {
    s <- sum(x); if (s <= 0) return(0)
    p <- x / s
    nz <- p > 0
    sum(p[nz] * log(p[nz] / q[nz]))
  }
  
  .dmin_bip <- function(x, q, cs = NULL) {
    n <- sum(x)
    if (n <= 0) return(0)
    expec <- floor(q * n)
    rest  <- n - sum(expec)
    xnew  <- expec
    for (j in seq_len(rest)) {
      i.vec <- if (is.null(cs)) seq_along(xnew) else which(xnew < cs)
      xsum <- sum(xnew)
      p1 <- xnew / (xsum + 1)
      dstep1 <- ifelse(xnew != 0, p1 * log(p1 / q), 0)
      dcheck <- rep(Inf, length(xnew))
      for (i in i.vec) {
        pi1 <- (xnew[i] + 1) / (xsum + 1)
        dcheck[i] <- pi1 * log(pi1 / q[i]) + sum(dstep1[-i])
      }
      i.best <- which.min(dcheck)[1]
      xnew[i.best] <- xnew[i.best] + 1L
    }
    .d_row_kl(xnew, q)
  }
  
  ## ---------- basic checks & matrix ----------
  req <- c("pollinator", "plant", "interaction")
  stopifnot(all(req %in% names(df)))
  df <- as.data.frame(df)
  df <- df[stats::complete.cases(df[req]), req, drop = FALSE]
  
  if (nrow(df) == 0) return(NULL)
  
  W <- xtabs(interaction ~ pollinator + plant, df)  # fills unobserved pairs with 0
  if (!(focal_species %in% rownames(W))) return(NULL)
  
  keep_r <- rowSums(W) > 0
  keep_c <- colSums(W) > 0
  if (!keep_r[focal_species]) return(NULL)
  W <- W[keep_r, keep_c, drop = FALSE]
  
  others <- setdiff(rownames(W), focal_species)
  if (!length(others)) return(NULL)
  
  out <- data.frame(pollinator = others, stringsAsFactors = FALSE)
  
  ## ------------ normalized degree (observed adjacency) ------------
  A <- W > 0
  deg <- rowSums(A)
  n_plants <- ncol(W)
  out$degree         <- as.numeric(deg[others])
  out$norm_degree    <- as.numeric(deg[others] / max(n_plants, 1L))
  out$n_plants_total <- n_plants
  
  ## ------------ PDI (Poisot) standardized to 0–1 ------------
  if (n_plants <= 1) {
    out$pdi <- NA_real_
  } else {
    pdi_all <- apply(W, 1, function(v) {
      if (all(v <= 0)) return(NA_real_)
      vmax <- max(v)
      if (vmax <= 0) return(NA_real_)
      vnorm <- v / vmax                      # scale by strongest link
      imax  <- which.max(v)[1]               # one instance of the max
      rest  <- vnorm[-imax]
      (sum(1 - rest)) / (length(v) - 1)      # = average (1 - v_rest)
    })
    out$pdi <- as.numeric(pdi_all[others])
  }
  
  ## ------------ observed overlaps (resource use) ------------
  if (return_observed) {
    f      <- as.numeric(W[focal_species, ])
    f_sum  <- sum(f)
    f_sq   <- sum(f^2)
    f_d    <- f_sq / (f_sum^2)
    
    O      <- W[others, , drop = FALSE]
    O_sum  <- rowSums(O)
    O_sq   <- rowSums(O^2)
    cross  <- as.numeric(O %*% f)
    
    out$pianka <- cross / sqrt(pmax(f_sq * O_sq, eps))
    out$morisita_horn <- (2 * cross) /
      ((f_d + (O_sq / pmax(O_sum^2, eps))) * f_sum * pmax(O_sum, eps))
  }
  
  ## ------------ preference overlaps (closed-form E = r*c/N) ------------
  r  <- rowSums(W); c <- colSums(W); N <- sum(W)
  E      <- (r %o% c) / max(N, eps)
  Pref   <- W / pmax(E, eps)
  Pref[E == 0] <- 0
  rs     <- rowSums(Pref)
  Pref   <- Pref / pmax(rs, eps)   # row-normalize
  
  Fp     <- as.numeric(Pref[focal_species, ])
  Op     <- Pref[others, , drop = FALSE]
  numer  <- as.numeric(Op %*% Fp)
  Sf     <- sum(Fp^2)
  So     <- rowSums(Op^2)
  
  out$pref_pianka        <- numer / sqrt(pmax(Sf * So, eps))
  out$pref_morisita_horn <- (2 * numer) / pmax(Sf + So, eps)
  
  ## ------------ Blüthgen's d' (standard, all species in q; optional keep) ------------
  P_use <- W / pmax(rowSums(W), eps)
  Q_all <- colSums(W); Q_all <- Q_all / pmax(sum(Q_all), eps); Q_all[Q_all < eps] <- eps
  log_ratio_all <- log(P_use / matrix(Q_all, nrow(P_use), ncol(P_use), byrow = TRUE))
  log_ratio_all[!is.finite(log_ratio_all)] <- 0
  D_all <- rowSums(P_use * log_ratio_all)
  dmax_all <- if (any(Q_all > 0)) log(1 / min(Q_all[Q_all > 0])) else NA_real_
  Dprime_all <- if (is.finite(dmax_all) && dmax_all > 0) pmin(pmax(D_all / dmax_all, 0), 1) else NA_real_
  out$d_prime <- as.numeric(Dprime_all[others])
  
  ## ------------ d' with Apis excluded from availability (bipartite-consistent) ------------
  nonapis_rows <- setdiff(rownames(W), focal_species)
  Wq <- if (length(nonapis_rows)) W[nonapis_rows, , drop = FALSE] else W * 0
  q  <- colSums(Wq)
  keep_cols <- q > 0
  if (!any(keep_cols)) {
    out$d_prime_no_apis_finite <- NA_real_
  } else {
    q  <- q[keep_cols]; q <- q / sum(q)
    cs <- colSums(Wq[, keep_cols, drop = FALSE])   # column caps
    m  <- sum(Wq)                                  # total (non-Apis) interactions
    X  <- W[, keep_cols, drop = FALSE]
    
    Ai    <- rowSums(X)
    draw  <- apply(X, 1, .d_row_kl, q = q)
    dmin  <- vapply(seq_len(nrow(X)), function(i) .dmin_bip(X[i, ], q, cs = cs), numeric(1))
    dmax  <- log(pmax(m / pmax(Ai, 1), 1))        # bipartite::dfun style
    
    dprime_excl <- (draw - dmin) / pmax(dmax - dmin, eps)
    dprime_excl <- pmin(pmax(dprime_excl, 0), 1)
    out$d_prime_no_apis_finite <- as.numeric(dprime_excl[others])
  }
  
  ## ------------ Müller's PAC: exposure of i to Apis (directional) ------------
  # s_apis(k) = alpha_{Apis,k} / sum_m alpha_{mk}
  s_apis <- as.numeric(W[focal_species, ]) / pmax(colSums(W), eps)
  
  # p_i(k) = alpha_{i,k} / sum_l alpha_{i,l}; PAC(i<-Apis) = sum_k p_i(k) * s_apis(k)
  O      <- W[others, , drop = FALSE]
  O_rs   <- pmax(rowSums(O), eps)
  P_i    <- O / O_rs
  out$pac_from_apis <- as.numeric(P_i %*% s_apis)
  
  rownames(out) <- NULL
  out
}

###




##### IMPORT AND FILTERING NETWORKS ####
########## load data version 1.3 10.5281/zenodo.15183272
#lanuza<-readRDS("G:\\Il mio Drive\\Articoli\\Resource Overlap EU\\Interaction_data.rds")
#sum(lanuza$Interaction)
#
######## create a study x network x date ID
#lanuza$id<-as.factor(
#  paste(lanuza$Study_id, lanuza$Network_id, lanuza$Date, sep = "//")) # id used for further analyses
#
################### remove all network ids that do not meet including criteria
#
## Convert to data.table if not already
#setDT(lanuza)
#
## Keep only ids that contain "Apis mellifera"
#df_withApis <- lanuza[lanuza[, any(Pollinator_accepted_name == "Apis mellifera"), by = id][V1 == TRUE, id], on = "id"]
#length(unique(lanuza$id)) ### all nets
#length(unique(df_withApis$id)) ### all nets with Apis
#sum(df_withApis$Interaction)
## Define the bee families
#bee_families <- c("Apidae", "Halictidae", "Andrenidae", "Colletidae", "Megachilidae", "Melittidae")
#
## Create the Pollinator_yes column based on your criteria
#df_withApis <- df_withApis %>%
#  dplyr::mutate(
#    Pollinator_yes = case_when(
#      Pollinator_order == "Coleoptera" ~ "Coleo.",
#      Pollinator_order == "Diptera" & Pollinator_family == "Syrphidae" ~ "Syrph.",
#      Pollinator_order == "Diptera" & Pollinator_family != "Syrphidae" ~ "Non-syr.",
#      Pollinator_order == "Hymenoptera" & Pollinator_family %in% bee_families ~ "Bees",
#      Pollinator_order == "Hymenoptera" & !(Pollinator_family %in% bee_families) ~ "Hymen.",
#      Pollinator_order == "Lepidoptera" ~ "Lepid.",
#      TRUE ~ "no"
#    )
#  )
#
###sum of the interactions for each plant x pollinator interaction in each site x date
#dfs <- df_withApis %>%
#  group_by(Bioregion, Country, Study_id, EuPPollNet_habitat, Network_id, Date, id,Pollinator_yes,
#           Pollinator_order, Pollinator_family, Pollinator_genus, Plant_accepted_name, Pollinator_accepted_name) %>%
#  summarise(Interactions = sum(Interaction, na.rm = TRUE), .groups = "drop")
#
###### rename variables
#dfs <- dfs %>%
#  dplyr::rename(pollinator = Pollinator_accepted_name)
#dfs <- dfs %>%
#  dplyr::rename(plant = Plant_accepted_name)
#dfs <- dfs %>%
#  dplyr::rename(interaction = Interactions)
#

# Convert to data.table if not already
#DT <- as.data.table(dfs)

## Step 1: Remove ids with interaction > 1000
#ids_high_interaction <- DT[interaction > 1000, unique(id)]
#
## Step 2: Remove ids with only one plant
#ids_one_plant <- DT[, .(n_plant = uniqueN(plant)), by = id][n_plant == 1, id]
#
## Step 3: Remove ids with only one pollinator
#ids_one_pollinator <- DT[, .(n_pollinator = uniqueN(pollinator)), by = id][n_pollinator == 1, id]
#
## Step 4: Combine all sets of problematic ids
#ids_to_remove <- unique(c(ids_high_interaction, ids_one_plant, ids_one_pollinator))
#
## Step 5: Filter them out
#DT_filtered <- DT[!id %in% ids_to_remove] ########### WORKING DATASET 
#sum(DT_filtered$interaction)
#
## Step 6: Keep only required columns
#dt_simple <- unique(DT_filtered [, .(plant, pollinator, interaction, id)])
#dt_simple[, `:=`(plant = as.character(plant),
#                 pollinator = as.character(pollinator),
#                 id = as.character(id),
#                 interaction = as.numeric(interaction))]
#

#### CREATE NULL MODELS USING PATEFIELD ####
# Create incidence matrix list (observed)
make_matrix <- function(subdt) {
  mat <- dcast(subdt, plant ~ pollinator, value.var = "interaction", fill = 0)
  rn <- mat$plant
  mat <- as.matrix(mat[, -1, with = FALSE])
  rownames(mat) <- rn
  return(mat)
}

# Build list of matrices: observed + 100 nulls per id
net_list <- list()

for (network_id in unique(dt_simple$id)) {
  sub <- dt_simple[id == network_id]
  obs_mat <- make_matrix(sub)
  
  # Store observed matrix
  net_list[[paste0(network_id, "_observed")]] <- obs_mat
  
  # Nulls using r2dtable
  rs <- rowSums(obs_mat)
  cs <- colSums(obs_mat)
  null_mats <- r2dtable(500, rs, cs) #########define number of null models
  
  for (i in seq_along(null_mats)) {
    null_name <- paste0(network_id, "_null", i)
    rownames(null_mats[[i]]) <- rownames(obs_mat)
    colnames(null_mats[[i]]) <- colnames(obs_mat)
    net_list[[null_name]] <- null_mats[[i]]
  }
}

# Convert all networks back to long format
net_long <- rbindlist(lapply(names(net_list), function(name) {
  m <- as.data.table(as.table(net_list[[name]]))
  setnames(m, c("plant", "pollinator", "interaction"))
  m[, id := sub("(_observed|_null\\d+)$", "", name)]
  m[, type := ifelse(grepl("observed$", name), "observed", name)]
  return(m)
}), use.names = TRUE)

# reorder columns
setcolorder(net_long, c("id", "type", "plant", "pollinator", "interaction"))

# Remove 0s
net_long_no0 <- net_long[interaction > 0]
saveRDS(net_long_no0, "G:\\Il mio Drive\\Articoli\\Resource Overlap EU\\Null_Pate_500.rds")

## create meaningful ids
net_long_no0 <-net_long_no0  %>%
  dplyr::rename(null_id=type)
net_long_no0[, Type := ifelse(null_id == "observed", "observed", "null")]
net_long_no0[, final_id := ifelse(Type == "observed", id, null_id)]

# COMPUTE ALL METRICS ####
net_long_no0<-readRDS("G:\\Il mio Drive\\Articoli\\Resource Overlap EU\\Null_Pate_500.rds")
net_long_no0 <-net_long_no0  %>%
  dplyr::rename(null_id=type)
net_long_no0[, Type := ifelse(null_id == "observed", "observed", "null")]
net_long_no0[, final_id := ifelse(Type == "observed", id, null_id)]

###################### loop all the metrics with data.table 16 august

library(data.table)
setDT(net_long_no0)

metrics_all_dt <- net_long_no0[, {
  res <- compute_metrics5(.SD,
                          focal_species   = "Apis mellifera",
                          return_observed = TRUE)
  if (is.null(res)) NULL else as.data.table(res)
}, by = .(final_id, id, Type)]

saveRDS(metrics_all_dt, "G:\\Il mio Drive\\Articoli\\Resource Overlap EU\\all_metrics_500.rds")

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

saveRDS(ses_tbl, "G:\\Il mio Drive\\Articoli\\Resource Overlap EU\\ses_tbl.rds")


####### IMPORT METRICS #### 

ses_tbl<-readRDS("G:\\Il mio Drive\\Articoli\\Resource Overlap EU\\ses_tbl.rds")
head(ses_tbl)
hist(ses_tbl$SES_pdi)

############## CREATE NETWORK AND SPECIES PREDICTORS ####
DT_filtered#take working dataset and rename in df to create all predictors

#### APIS METRICS: SECTION  BASED ON DF #####
# Compute total number of Apis interactions per network
apis_interactions <- DT_filtered%>%
  filter(pollinator == "Apis mellifera") %>%
  group_by(id) %>%
  summarise(
    Apis_interactions = sum(interaction, na.rm = TRUE),
    Apis_degree = n_distinct(plant)
  )

##############  Compute Apis interaction evenness per network 
apis_interaction_evenness <- DT_filtered %>%
  filter(pollinator == "Apis mellifera") %>%
  group_by(id, pollinator, plant) %>%
  summarise(w_ij = sum(interaction), .groups = "drop") %>%
  group_by(id, pollinator) %>%
  mutate(
    total_w = sum(w_ij),
    p_ij = w_ij / total_w,
    log_k = log(n_distinct(plant))  # log of partner count
  ) %>%
  summarise(
    H = -sum(p_ij * log(p_ij)),  # Shannon entropy
    log_k = first(log_k),
    Apis_interaction_evenness = ifelse(log_k > 0, H / log_k, NA),
    .groups = "drop"
  )

##############  Compute Apis Simpson dominance per network 
apis_simpson <- DT_filtered %>%
  filter(pollinator == "Apis mellifera") %>%
  group_by(id, plant) %>%
  summarise(n = sum(interaction), .groups = "drop") %>%
  group_by(id) %>%
  mutate(p = n / sum(n)) %>%
  summarise(Apis_simpson = sum(p^2))

# Compute total number of all interactions per network and total number of pollinator species
all_interactions <- DT_filtered %>%
  group_by(id) %>%
  summarise(
    All_interactions = sum(interaction, na.rm = TRUE),
    pollinator_n = n_distinct(pollinator),
    plant_n = n_distinct(plant)
  )

# Join Apis interactions, Apis interaction evenness, and all interactions
apis_summary <- left_join(apis_interactions, all_interactions, by = "id")
apis_summary1 <- left_join(apis_summary, apis_simpson, by = "id")
apis_summary2 <- left_join(apis_summary1, apis_interaction_evenness, by = "id")
apis_summary2$Dominance<-log(
  apis_summary2$Apis_interactions/(apis_summary2$All_interactions-apis_summary2$Apis_interactions) 
)
apis_summary2$Apis_ndegree<-apis_summary2$Apis_degree/apis_summary2$plant_n 
apis_summary2$Apis_1_D<-1-apis_summary2$Apis_simpson 


########### EXTRACT CLIMATE ############
###### Load temperature raster Chelsea 

r<-rast("G:\\Il mio Drive\\Articoli\\Resource Overlap EU\\CHELSA_EUR11_tas_norm_1981-2005_V1.1.nc")
Temp_summer<-(r[[5]]+r[[6]]+r[[7]]+r[[8]])/4 ## average summer temperature (May-Aug)
summary(Temp_summer)

####### Coordinates of the network
network_coords <- df_withApis %>%
  dplyr::select(Latitude, Longitude, id) %>%
  distinct()

# Convert to sf with WGS84
points_sf <- st_as_sf(network_coords, coords = c("Longitude", "Latitude"), crs = 4326)

# Create 5 km buffer (in meters)
points_buf <- st_buffer(points_sf, dist = 5000)

# Extract for all months using a stack
extracted <- as.data.frame(exact_extract(Temp_summer, points_buf, 'mean'))
names(extracted )[1] <- "Temperature"

#Add site IDs
extracted$id <- network_coords$id

temperature_summary <- extracted %>%
  group_by(id) %>%
  dplyr::summarise(
    Temperature = mean(Temperature, na.rm = TRUE)
  )

####### EXTRACT TAXONOMY PER SPECIES ############
species_data <- DT_filtered %>%
  dplyr::select(Pollinator_order, Pollinator_yes, Pollinator_family, Pollinator_genus, pollinator) %>%
  distinct()
species_data <- species_data %>%
  mutate(Pollinator_order2 = if_else(Pollinator_order %in% c("Diptera", "Lepidoptera", "Hymenoptera"), Pollinator_order, "Other"))


######### COMPUTE N PLANT PER NETWORK #############
plants_per_network <- DT_filtered %>%
  filter(!is.na(plant)) %>%
  group_by(id) %>%
  summarise(n_plants = n_distinct(plant)) %>%
  arrange(desc(n_plants))
plants_per_network_df<-as.data.frame(plants_per_network)


######### COMPUTE N INTERACTIONS PER POLLINATOR PER NETWORK #############

degree_interactions <- DT_filtered %>%
  dplyr::mutate(interaction = as.numeric(interaction)) %>%
  dplyr::group_by(id, pollinator) %>%
  dplyr::summarise(
    sp_interactions = sum(interaction, na.rm = TRUE),
    sp_degree       = dplyr::n_distinct(plant[interaction > 0]),
    .groups = "drop"
  )


############## CREATE HABITAT ###########
habitat_net<- DT_filtered %>%
  group_by(id) %>%
  summarise(habitat = first(EuPPollNet_habitat))
unique(habitat_net$habitat)

habitat_net <- habitat_net %>%
  mutate(habitat_type = case_when(
    habitat %in% c(
      "Forest/woodland",
      "Riparian vegetation",
      "Semi-natural grasslands",
      "Moors and heathland",
      "Sclerophyllous vegetation",
      "Montane to alpine grasslands",
      "Beaches, dunes, sands",
      "Riparian"
    ) ~ "Natural",
    habitat %in% c(
      "Ruderal vegetation",
      "Intensive grasslands",
      "Agricultural land",
      "Agricultural margins",
      "Green urban areas"
    ) ~ "Anthropic",
    TRUE ~ NA_character_  # in case there's any unclassified habitat
  ))


############# JOINS ###############
############### join metrics with species data 

join0 <- ses_tbl %>%
  left_join(degree_interactions, by =c("id" , "pollinator"))

join1 <- join0 %>%
  left_join(species_data, by = "pollinator")

############# join plant species 
join2 <- join1  %>%
  left_join(plants_per_network_df, by = "id")

############# join climate
join3 <- join2 %>%
  left_join(temperature_summary, by = "id")

############# join habitat
join4 <- join3 %>%
  left_join(habitat_net, by = "id")

############# join Apis 
final_join <- join4  %>%
  left_join(apis_summary2, by = "id")

str(final_join) ############ dataframe entering the modeling

##### FILTER OUT SINGLETONS, DOUBLETONS ETC. #########
names(final_join)
df_mod <- subset(final_join,
                 Pollinator_yes!="no"   &   
                 sp_interactions >4     &
                 Apis_interactions>4    &
                 Pollinator_genus!="NA" &
                 sp_degree>0            &
                 n_plants>1)##

df_mod<-df_mod%>%
  mutate(Group=as.factor(Pollinator_yes),
          Plant=plant_n,
          Habitat=habitat_type,
          id=as.factor(id),
          pollinator.x=as.factor(pollinator.x),
         .keep = "unused")

names(df_mod)

library(dplyr)
library(tidyr)

df_mod <- df_mod %>%
  separate(id, into = c("study_id", "network_id", "date"), sep = "//")

########MODELING ########### 

##### PDI ####
fit_SES_pdi<- glmmTMB(
  SES_pdi ~  Group*scale(Dominance) + Habitat + scale(Temperature) + scale(log(Plant))*scale(Dominance)+scale(I(Dominance^2))+
    (1|study_id/network_id/date),
 dispformula = ~scale(Dominance)*scale(log(Plant))+
   (1|study_id/network_id/date),
   data = subset(df_mod, pdi_iqr>0 & pdi_sd>1e-6 & pdi_uniq>3)
)

summary(fit_SES_pdi)
DHARMa::simulateResiduals(fit_SES_pdi, plot = TRUE, n=1000)
Anova(fit_SES_pdi)
plot(allEffects(fit_SES_pdi))
performance::check_collinearity(fit_SES_pdi)


######### NDEGREE #######
names(df_mod2)

fit_SES_nd <- glmmTMB(
  SES_norm_degree ~  Group*scale(Dominance) + Habitat + scale(Temperature) + scale(log(Plant))*scale(Dominance)+scale(I(Dominance^2))+
    (1|study_id/network_id/date),
  dispformula = ~scale(Dominance)+scale(log(Plant))+
    (1|study_id/network_id/date),
  data = subset(df_mod, norm_degree_iqr>0 & norm_degree_sd>1e-6 & norm_degree_uniq>3)
)

DHARMa::simulateResiduals(fit_SES_nd, plot = TRUE, n=1000)
summary(fit_SES_nd)
Anova(fit_SES_nd)
plot(allEffects(fit_SES_nd))
performance::check_collinearity(fit_SES_nd)

######### M-H #######
names(df_mod)
fit_SES_mh <- glmmTMB(
  SES_morisita_horn ~  Group*scale(Dominance)+scale(Plant) + Habitat + scale(Temperature) + scale(Plant)+scale(Dominance)+
    (1|study_id/network_id/date),
  dispformula = ~scale(Dominance)+scale(Plant)+
    (1|study_id/network_id/date),
  data = subset(df_mod, morisita_horn_iqr>0 & morisita_horn_sd>1e-6 & morisita_horn_uniq>3)
)

DHARMa::simulateResiduals(fit_SES_mh, plot = TRUE, n=1000)
summary(fit_SES_mh)
Anova(fit_SES_mh)
plot(allEffects(fit_SES_mh))
performance::check_collinearity(fit_SES_mh)


########  EXPORT TABLES  ###

# helper to format p-values
fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# return a formatted data.frame for one component
tab_df <- function(mod, component, label){
  df <- as.data.frame(parameters::model_parameters(
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

# build one combined table (same columns)
df_cond <- tab_df(fit_SES_nd, "conditional", "Conditional")
df_disp <- tab_df(fit_SES_nd, "dispersion",  "Dispersion")
df_all  <- dplyr::bind_rows(df_cond, df_disp)

# make one flextable
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("officer", quietly = TRUE))   install.packages("officer")

library(flextable); library(officer)

library(flextable)

# build the table
ft_all <- flextable(df_all) |>
  merge_v(j = "Component", part = "body") |>
  valign(j = "Component", valign = "top", part = "body") |>
  align(j = c("Component","Parameter"), align = "left", part = "body") |>
  autofit() |>
  set_caption("PAC")

# optional: add a thin rule before each new Component group (skip the first row)
idx <- which(!duplicated(df_all$Component))
idx <- idx[idx > 1]  # don’t draw before the very first row

ft_all <- hline(
  ft_all,
  i = idx, j = 1:ncol(df_all),
  border = fp_border(color = "gray50", width = 0.75),
  part = "body"
)

# write to docx
doc <- officer::read_docx() |>
  body_add_flextable(ft_all)

print(doc, target = "G:\\Il mio Drive\\Articoli\\Resource Overlap EU\\submission\\table_nd.docx")



###################### DOMINANCE X GROUP ################ 


library(ggeffects)
library(ggplot2)
# get marginal effects
    gp_mh <- ggpredict(fit_SES_mh, terms = c("Dominance[all]", "Group"), type = "fixed")
# plot with custom order
 mh<-plot(gp_mh) +
    labs(x = "", 
        y = "SES Morisita–Horn",
        color = "Group",
        main="") +
     ggtitle("")+
     theme_classic()

  
 

   gp_pdi <- ggpredict(fit_SES_pdi, terms = c("Dominance[all]", "Group"), type = "fixed")
 
   # plot with custom order
   pdi<-plot(gp_pdi) +
     labs(x = "Honeybee dominance (ln-ratio)", 
          y = "SES Paired difference index",
          color = "Group",
          main="") +
    ggtitle("")+
   theme_classic()
 
   gp_nd <- ggpredict(fit_SES_nd, terms = c("Dominance[all]", "Group"), type = "fixed")

   # plot with custom order
   nd<-plot(gp_nd) +
     labs(x = "", 
          y = "SES normalised degree",
          color = "Group",
          main="") +
     ggtitle("")+
     theme_classic()

   library(patchwork)




# One row, equal panel sizes, single legend on the right, automatic A)–D) tags
p <- (mh |pdi | nd) +
  plot_layout(ncol = 3, nrow=1, guides = "collect") &
  theme(legend.position = "right")

p + plot_annotation(tag_levels = "A", tag_suffix = ")")




################## DOMINANCE X PLANTS

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



######## FIG. 1 COVARIATION BETWEEN DOMINANCE, ND, N POLL AND N PLANT ####

# Load required packages
apis_summary2_no1<-apis_summary2%>%
  filter(Apis_interactions>1)

vars <- c("plant_n","pollinator_n","Apis_ndegree", "Dominance")
df_subset <- apis_summary2_no1[, vars, drop = FALSE]

# Custom lower panel: GAM + linear
lower_gam_lm <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "gray40", size = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed", ...) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "red", ...)
}

ggpairs(
  df_subset,
  columns = c("plant_n","pollinator_n","Apis_ndegree",  "Dominance"),
  columnLabels = c("Plant richness", "Pollinator richness", "Hb normalized degree",  "Hb dominance"),
  lower = list(continuous = wrap(lower_gam_lm)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5, fill = "lightblue")),
  upper = list(continuous = wrap("cor", size = 3))
  #diag = list(continuous = wrap(custom_density_diag))
)

#####################  FIGURE RANK-ABUNDANCE  ##########

# Frequency table
pollinator_freq <- df %>%
  distinct(id, pollinator) %>%
  count(pollinator, name = "networks_present_in") %>%
  arrange(desc(networks_present_in)) %>%
  mutate(rank = row_number())

pollinator_freq <- pollinator_freq %>%
  mutate(
    label = str_trim(pollinator),                         # remove leading/trailing spaces
    label = str_split(label, "\\s+"),                     # split into words
    label = sapply(label, function(words) {
      paste0(substr(words, 1, 3), collapse = " ")         # take first 3 letters of each word
    })
  )

# Add labels only for top 20
pollinator_freq <- pollinator_freq %>%
  mutate(label = ifelse(rank <= 50, label, NA))

# Plot
ggplot(pollinator_freq, aes(x = rank, y = log(networks_present_in))) +
  geom_point(col="grey55") +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf,         # Allow all labels to be shown
    box.padding = 0.5,          # Space around labels
    point.padding = 0.5,        # Space between label and point
    nudge_x = 50,             # Push labels upward
    segment.size = 0.2,         # Line thickness from point to label
    segment.size = 0.2,
    force_pull = 2,
    segment.color = "grey80",     # light grey color
    segment.alpha = 0.5           # Line transparency
  ) +
  scale_y_continuous(
    name = "Species count in the networks",
    breaks = log(c(1, 2, 10, 100, 500, 1500)),  # choose meaningful original values
    labels = c(1, 2, 10, 100, 500, 1500)        # show original scale on axis
  ) +
  labs(
    x = "Pollinator species rank",
    y = "Log(Frequency in networks)",
    title = ""
  ) +
  theme_minimal()





