#### Load libraries ####
library(data.table)
library(dplyr)

#### Load processed data ####
dt_simple = readRDS("Data/dt_simple.rds")

#### CREATE NULL MODELS USING PATEFIELD ####
# Create incidence matrix list (observed)
make_matrix = function(subdt) {
  mat = dcast(subdt, plant ~ pollinator, 
              value.var = "interaction", 
              fun.aggregate = function(x) sum(x, na.rm = TRUE),
              fill = 0)
  rn = mat$plant
  mat = as.matrix(mat[, -1, with = FALSE])
  rownames(mat) = rn
  return(mat)
}

# Build list of matrices: observed + 100 nulls per id
net_list = list()

for (network_id in unique(dt_simple$id)) {
  sub = dt_simple[id == network_id]
  obs_mat = make_matrix(sub)
  
  # Store observed matrix
  net_list[[paste0(network_id, "_observed")]] = obs_mat
  
  # Nulls using r2dtable
  rs = rowSums(obs_mat)
  cs = colSums(obs_mat)
  null_mats = r2dtable(500, rs, cs) #########define number of null models
  
  for (i in seq_along(null_mats)) {
    null_name = paste0(network_id, "_null", i)
    rownames(null_mats[[i]]) = rownames(obs_mat)
    colnames(null_mats[[i]]) = colnames(obs_mat)
    net_list[[null_name]] = null_mats[[i]]
  }
}


# Convert all networks back to long format
net_long = rbindlist(lapply(names(net_list), function(name) {
  m = as.data.table(as.table(net_list[[name]]))
  setnames(m, c("plant", "pollinator", "interaction"))
  m[, id := sub("(_observed|_null\\d+)$", "", name)]
  m[, type := ifelse(grepl("observed$", name), "observed", name)]
  return(m)
}), use.names = TRUE)

# reorder columns
setcolorder(net_long, c("id", "type", "plant", "pollinator", "interaction"))

# Remove 0s
net_long_no0 = net_long[interaction > 0]

#Save data
saveRDS(net_long_no0, "Data/net_long_no0.rds")





