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
df <- read_msdial("inst/gc_msdial.txt")
#
#
#
#
#
#| label: extract metadata and chemical classes
#| message: false
#| warning: false
df <- extract_cid(df[!grepl("Naphthalene-D8", df$Name), ], inchikey_col = 9, name_col = 7) 
df <- extract_meta(df, cas = TRUE, flavornet = TRUE)
df <- extract_classyfire(df)
#
#
#
#
#
#| label: calculate detection rate
df$Comment <- gsub(";.*$", "", df$Comment) # remove any notes in the Comment column
df <- calculate_freq(df, num_sample = 2, sep = ",")
df[, 1:30]
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
export4toxtree(df, cas_col = 19, name_col = 7, output = "inst/gc_for_toxtree.csv")
# then open Toxtree and run batch mode to calculate Cramer level

# load databases
load_databases()
df <- assign_toxicity(df, toxtree_result = "inst/gc_toxtree_results.csv")

#
#
#
