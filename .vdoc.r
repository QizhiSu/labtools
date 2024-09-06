#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: read files in R
#| message: false
#| warning: false

# Reading different file formats in R
# If rio package is not installed, install it first, otherwise skip this step.
# install.packages("rio") 
  library(rio)

# Reading xlsx file. 
  df <- import("inst/chemicals.xlsx")
  df
#
#
#
#
#
#| label: write files in R
#| message: false
#| warning: false
#| echo: false

# Writing different file formats in R
  export(df, "inst/chemicals.xlsx")

# Check the documentation for more details. Using ? before the function name will display the help page.
  ?import
  ?export
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: installation
#| message: false
#| warning: false
#| eval: false

# If you haven't installed the remotes package, you can install it with:
  install.packages("remotes")

# Then you can install the labtools package with:
  remotes::install_github("QizhiSu/labtools")
#
#
#
#
#
#
#
#| label: cid extraction
#| message: false
#| warning: false
#| class: output

# Given that you have read in your data as df as shown in the Section 1, you can use the following function to extract the metadata from Pubchem:

library(labtools)
df_cid <- extract_cid(df, inchikey_col = 3, cas_col = 2, name_col = 1)
df_cid
#
#
#
#
#
#
#
#| label: manual CID assignment
#| message: false
#| warning: false

# Manually assign CID to the first row
df_cid[1, "CID"] <- 6184

# Just remember that CID is integer and not character. So don't do this:
# df_cid[1, "CID"] <- "6184"
#
#
#
#
#
#| label: meta extraction
#| message: false
#| warning: false
#| class-output: output

df_meta <- extract_meta(df_cid)
df_meta

# If you also want to retrieve CAS number, in case you haven't provided it in the input dataframe, you can use the following function, but it will be slower:
df_meta <- extract_meta(df_cid, cas = TRUE)
df_meta

# It is also possible to retrieve odour information from the Flavonet (https://www.flavornet.org/flavornet.html):
df_meta <- extract_meta(df_cid, cas = TRUE, flavornet = TRUE)
df_meta

# Once metadata are retrieved, you can export them into a new xlsx file. 
export(df_meta, "inst/chemicals_meta.xlsx")
#
#
#
#
#
#
#
#| label: class extraction
#| message: false
#| warning: false
#| class: output

df_cla <- extract_classyfire(df_meta)
df_cla
#
#
#
#
#
#
#
#| eval: false

# If rCDK and Retip are not installed, you can install them with, otherwise skip this step.
install.packages("rCDK")
install.packages("Retip")
#
#
#
#| label: descriptor calculation
#| message: false
#| warning: false

library(Retip)

df_cd <- Retip::getCD(df_meta[, c("Name", "SMILES")])
df_cd[, 1:10] # show the first 10 columns
#
#
#
#
#
#
#
#
#
#| label: structure viewer
#| message: false
#| warning: false

# subset the dataframe to show only Name, CAS, CID, and SMILES
df2view <- df_meta[, c("Name", "CAS", "CID", "SMILES")]
df2view

# navigate_chem(df2view)
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: read MS-DIAL data and extract metadata
#| message: false
#| warning: false
# read MS-DIAL data
library(labtools)

df <- read_msdial("inst/gc_msdial.txt", keep_unknown = FALSE, keep_spectrum = FALSE, keep_mean_sd = FALSE)
#
#
#
#
#
#| label: extract metadata and chemical classes
#| message: false
#| warning: false
df_meta <- df[, 1:12] # exclude peak area columns
df_meta <- extract_cid(df_meta, inchikey_col = 9, name_col = 7) 
df_meta <- extract_meta(df_meta, cas = TRUE, flavornet = TRUE)
df_meta <- extract_classyfire(df_meta)
#
#
#
#
#
#| label: calculate detection rate
df_meta$Comment <- gsub(";.*$", "", df_meta$Comment) # remove any notes in the Comment column
df_meta <- calculate_freq(df_meta, num_sample = 2, sep = ",")
df_meta[, 1:15]
#
#
#
#
#
#
#
#| label: sample code mapping
sam_code <- rio::import("inst/sample_code.xlsx")
sam_code
#
#
#
#
#
#
#
#
#| label: fcmsafety installation
#| eval: false
# remotes::install_github("QizhiSu/fcmsafety") # in case not installed
library(fcmsafety)

# export structural information for toxtree
export4toxtree(df_meta, cas_col = 19, name_col = 7, output = "inst/gc_for_toxtree.csv")
# then open Toxtree and run batch mode to calculate Cramer level

# load databases
load_databases()
df_meta <- assign_toxicity(df_meta, toxtree_result = "inst/gc_toxtree_results.csv")
df_meta
#
#
#
#
#
#
#
#| label: calculate molecular descriptors
#| message: false
#| warning: false
library(Retip)

# load reference compounds and extract metadata
ref <- rio::import("inst/reference.xlsx")
ref_meta <- extract_cid(ref, name_col = 1) 
ref_meta <- extract_meta(ref_meta)

# calculate molecular descriptors for reference compounds.
ref_cd <- Retip::getCD(ref_meta)

# Since the reference compounds may not be changed that frequently, we can calculate the descriptors once and save them as a xlsx file.
rio::export(ref_cd, "inst/ref_cd.xlsx")

# In case you do not change reference compounds, Then, we can load the xlsx file and use it for subsequent analysis. Otherwise, it is better to do it from scratch.
ref_cd <- rio::import("inst/ref_cd.xlsx")
ref_cd[, 1:15]

# calculate molecular descriptors for sample compounds
df_cd <- Retip::getCD(df_meta[, c("ID", "Name", "SMILES", "Comment")])
df_cd[, 1:15]
#
#
#
#
#
#
#
#| label: select reference standard
# select reference standard for each target compound
# do not use response for the calculation
df_con <- select_std(std_md = ref_cd, std_res_col = FALSE, std_md1_col = 14, data_md = df_cd, data_md1_col = 5)

# use the top 50% of the molecular descriptors for the calculation
# df_con <- select_std(std_md = ref_cd, std_res_col = 3, std_md1_col = 14, data_md = df_cd, data_md1_col = 4, top_npct_md = 50)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: combine data
#| message: false
#| warning: false
library(dplyr)

# combining data
df_con <- left_join(df_con, df[, c(7, 13:ncol(df))], by = "Name")
#
#
#
#
#
#
#| label: keep peak areas
# load sample code mapping file
sam_code <- rio::import("inst/sample_code.xlsx")
sam_code
colnames(df_con)

# keep peak areas of identified sample and turn others into NA
df_con <- keep_area(df_con, sam_code = sam_code, start_col = 12, keep_bk = TRUE, keep_d8 = TRUE) 
#
#
#
#
#
#| label: calculate concentration
# calculate concentration of each compound in each sample
sam_weight <- rio::import("inst/sample_weight.xlsx")
sam_weight

# concentration considering sample weight and final volumn
df_con <- calculate_con(df_con, sam_weight = sam_weight, start_col = 12)

# concentration without sample weight and final volumn adjustment
# df_con <- calculate_con(df_con, start_col = 12, sam_weight = NULL)
#
#
#
#
#
#| label: organize data
# organize data
# set NA to zero and not bind mean and sd
df_con <- organize_con(df_con, sam_code = sam_code, start_col = 12, digits = 2, na2zero = TRUE, bind_mean_sd = FALSE)

# do not set NA to zero and do not bind mean and sd
# df_con <- organize_con(df_con, sam_code = sam_code, start_col = 12, digits = 2, na2zero = FALSE, bind_mean_sd = FALSE)

# if you want to combine all data into one, then you can run
df_bind <- left_join(df_meta, select(df_con, -c("SMILES", "Comment")), by = c("ID", "Name"))
#
#
#
#
#
