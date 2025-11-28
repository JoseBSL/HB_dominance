############################################################
# Script name: "5_2_export_modelling_output.R"
# Objective: Create tables in word format from modelling output
#
# Inputs:
#   - Data/fit_SES_nd_over9.rds
#   - Data/fit_SES_mh_over9.rds
#   - Data/fit_SES_pdi_over9.rds
#
# Outputs:
#   - Data/Tables/table_nd_over9.docx
#   - Data/Tables/table_mh_over9.docx
#   - Data/Tables/table_pdi_over9.docx
#
# Author: L. Marini | Reviewer: J.B. Lanuza | Year: 2025
############################################################
#### Load libraries ####
library(flextable)
library(officer)
library(dplyr)

#### Load data ####
fit_SES_pdi = readRDS("Data/fit_SES_pdi_over9.rds")
fit_SES_nd  = readRDS("Data/fit_SES_nd_over9.rds")
fit_SES_mh  = readRDS("Data/fit_SES_mh_over9.rds")
#### Load function ####
source("Scripts/export_table_function.R")

#### ND ####
df_nd  = bind_rows(
  tab_df(fit_SES_nd, "conditional", "Conditional"),
  tab_df(fit_SES_nd, "dispersion",  "Dispersion"))
ft_nd  = make_ft(df_nd, caption = "ND")
doc_nd = read_docx() |> body_add_flextable(ft_nd)
print(doc_nd, target = "Data/Tables/table_nd_over9.docx")

#### MH ####
df_mh  = bind_rows(
  tab_df(fit_SES_mh, "conditional", "Conditional"),
  tab_df(fit_SES_mh, "dispersion",  "Dispersion"))
ft_mh  = make_ft(df_mh, caption = "MH")
doc_mh = read_docx() |> body_add_flextable(ft_mh)  # fixed: ft_mh
print(doc_mh, target = "Data/Tables/table_mh_over9.docx")

#### PDI ####
df_pdi  = bind_rows(
  tab_df(fit_SES_pdi, "conditional", "Conditional"),
  tab_df(fit_SES_pdi, "dispersion",  "Dispersion"))
ft_pdi  = make_ft(df_pdi, caption = "PDI")
doc_pdi = read_docx() |> body_add_flextable(ft_pdi)
print(doc_pdi, target = "Data/Tables/table_pdi_over9.docx")
