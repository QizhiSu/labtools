library(tidyverse)
df <- rio::import("/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/NQI_MSMS_database_20250523.xlsx") |> # nolint
  distinct(Name, .keep_all = TRUE) |>
  select(Name, CAS, SMILES, InChI, InChIKey, Ontology)
set.seed(123)
df <- df[sample(nrow(df), 20), ]
df_cid <- extract_cid(df, cas_col = "CAS", inchikey_col = "InChIKey", name_col = "Name", timeout = 120)
df_meta <- extract_meta(df_cid)
df_cla <- extract_classyfire(df_meta)

rio::export(df_meta, "test.xlsx")
a <- data.frame(Name =df_meta[, 1])

b <- assign_meta(a, "test.xlsx")










