library(tidyverse)
library(rinchi)
library(rcdk)
library(labtools)
library(webchem)
library(Retip)
library(fcmsafety)

pos <- msp2df("/Users/sukis/Desktop/project/2022 NQI/data/all_1016-neg.msp")
neg <- msp2df("/Users/sukis/Desktop/project/2022 NQI/data/all_1016-pos.msp")

bind <- dplyr::bind_rows(pos, neg)
# export for manual curation
rio::export(bind, "/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_pos_and_neg_database_from_msp.xlsx") # nolint

# section -----------------------------------------------------------
# import the curated file
msp <- readxl::read_xlsx('/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_pos_and_neg_database_corrected_for_msp.xlsx', guess_max = 1000000) # nolint

# generate new InChI and InChIKey columns
msp$INCHI <- sapply(msp$SMILES, get.inchi)
msp$INCHIKEY_2 <- sapply(msp$SMILES, get.inchi.key)
msp$INCHIKEY_dif <- ifelse(msp$INCHIKEY_2 == msp$INCHIKEY, "Same", "Different")
msp <- relocate(msp, INCHIKEY_2:INCHIKEY_dif, .after = INCHIKEY)

# generate new SMILES
msp$SMILES_2 <- sapply(msp$SMILES, function(x) {
  mol <- parse.smiles(x)[[1]]
  get.smiles(mol, flavor = smiles.flavors(c("Isomeric")))
})
msp$SMILES_dif <- ifelse(msp$SMILES_2 == msp$SMILES, "Same", "Different")
msp <- relocate(msp, SMILES_2:SMILES_dif, .after = SMILES)

# generate new formula column
msp$FORMULA_2 <- sapply(msp$SMILES, function(x) {
  mol <- parse.smiles(x)[[1]]
  formula_obj <- get.mol2formula(mol)
  formula_obj@string # 直接提取 string 槽
})
msp$FORMULA_dif <- ifelse(msp$FORMULA_2 == msp$FORMULA, "Same", "Different")
msp <- relocate(msp, FORMULA_2:FORMULA_dif, .after = FORMULA)

# get exact mass
msp$ExactMass <- sapply(msp$SMILES, function(x) {
  mol <- parse.smiles(x)[[1]]
  get.exact.mass(mol)
})

rio::export(msp, "/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_pos_and_neg_database_corrected_for_msp_2.xlsx") # nolint

# section -----------------------------------------------------------
# read the curated file
msp <- readxl::read_xlsx("/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_pos_and_neg_database_corrected_for_msp_20250523_2.xlsx", guess_max = 1000000) # nolint

# msp_meta <- rio::import("/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_pos_and_neg_database_corrected_for_msp_meta.xlsx") # nolint

# # msp_cid <- msp |>
# #   select(Name, INCHIKEY) |>
# #   distinct(INCHIKEY, .keep_all = TRUE) |>
# #   extract_cid(inchikey_col = 2, timeout = 100)
# # msp_meta <- extract_meta(msp_cid, cas = TRUE)

# msp$CAS <- msp_meta$CAS_retrieved[match(msp$INCHIKEY, msp_meta$InChIKey)]
# msp$IUPACName <- msp_meta$IUPACName[match(msp$INCHIKEY, msp_meta$InChIKey)]
# msp$CID <- msp_meta$CID[match(msp$INCHIKEY, msp_meta$InChIKey)]
# msp <- relocate(msp, CID, IUPACName, .after = CAS)

# msp$AUTHORS[!is.na(msp$INSTRUMENT)] <- "Yining Xia"
# msp$AUTHORS[is.na(msp$AUTHORS)] <- "Jinxin Liang, Qizhi Su"
# msp$INSTRUMENT[is.na(msp$INSTRUMENT)] <- "Agilent 6546 QTOF"
# msp$INSTRUMENTTYPE <- "UPLC-QTOF/MS"
# msp$COMMENT <- NULL
# msp$DATE <- Sys.Date()
df2msp(msp, "/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_pos_and_neg_curated_20250523_2.msp") # nolint

# read in the curated file
msp <- msp2df("/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_pos_and_neg_curated_20250523.msp") # nolint

# # export the curated file in msp format
# pos <- msp[msp$IONMODE == "Positive", ]
# df2msp(pos, "/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_pos_curated_202505223.msp") # nolint
# neg <- msp[msp$IONMODE == "Negative", ]
# df2msp(neg, "/Users/sukis/Desktop/project/2022 NQI/data/IQTC_MSMS_neg_curated_20250523.msp") # nolint

# import guia database
guia_pos <- msp2df("/Users/sukis/Desktop/Databases/ms_spectra_libraries/LC-MS/combind/pos/GUIA_pos_435_compounds_1248_spectra_20250519_curated.msp") # nolint
# combine with
guia_neg <- msp2df("/Users/sukis/Desktop/Databases/ms_spectra_libraries/LC-MS/combind/neg/GUIA_neg_162_compounds_310_spectra_20250507_curated.msp") # nolint
guia <- dplyr::bind_rows(guia_pos, guia_neg)

# change column names of msp to allow binding
colnames(msp) <- c("Name", "RetentionTime", "PrecursorMZ", "PrecursorType", "IonMode", "CollisionEnergy", "Formula", "SMILES", "InChIKey", "InChI", "Authors", "Instrument", "InstrumentType", "CAS", "CID", "IUPACName", "ExactMass", "Ontoloty", "Date", "Num_peaks", "Spectrum") # nolint

all_msp <- dplyr::bind_rows(msp, guia) |>
  select(-c(ID, SpectrumType, MSLevel, Comment))
# all_msp <- all_msp |>
#   distinct(InChIKey, .keep_all = TRUE)
# msp_md <- Retip::getCD(msp)
# guia_md <- Retip::getCD(guia)
all_md <- dplyr::bind_rows(msp_md, guia_md) |>
  select(InChIKey, XLogP, MLogP, ALogP, TopoPSA, nB, nAtom) |>
  distinct(InChIKey, .keep_all = TRUE)
all_msp <- left_join(all_msp, all_md, by = "InChIKey") |>
  relocate(XLogP:nAtom, .before = Num_peaks)

# classifier
all_cla <- extract_classyfire(all_msp)
rio::export(all_cla, "/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/NQI_MSMS_database_classyfire.xlsx") # nolint
df2msp(msp, "/Users/sukis/Library/CloudStorage/OneDrive-unizar.es/project/2022 NQI/data/MSMS_library/IQTC_MSMS_pos_and_neg_database_combined.msp") #

all_cla <- relocate(all_cla, nAtom, .after = nB)
all_cla$RetentionTime[all_cla$RetentionTime == 0] <- NA

rio::export(all_cla, "/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/NQI_MSMS_database_classyfire_20250523.xlsx")

export4toxtree(all_msp,
  name_col = 1, cas_col = 14,
  "/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/NQI_MSMS_database_for_toxtree.csv"
)

all_qsar <- rio::import('/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/QSAR_results.xlsx')
all_qsar$InChIKey <- sapply(all_qsar$SMILES, get.inchi.key)
rio::export(all_qsar, "/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/QSAR_results_2.xlsx")

all_msp <- rio::import('/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/NQI_MSMS_database_20250523.xlsx') # nolint

df2msp(all_msp, "/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/NQI_MSMS_database_20250523.msp") # nolint

# section -----------------------------------------------------------
extract_meta <- function(data, cas = FALSE, flavornet = FALSE) {
  if (!("CID" %in% colnames(data))) {
    stop("Your data does not contain a 'CID' column. Please run extract_cid() first.")
  }

  message("Extracting metadata from Pubchem based on CID.")
  data$CID <- as.integer(data$CID)

  # Rename any conflicting old columns
  to_rename <- c("InChIKey", "SMILES", "Formula", "ExactMass")
  for (col in to_rename) {
    if (col %in% colnames(data)) {
      data <- data %>% rename(!!paste0(col, "_old") := all_of(col))
    }
  }

  # Extract properties from PubChem
  pubchem_data <- data %>%
    filter(!is.na(CID)) %>%
    distinct(CID) %>%
    pc_prop(properties = c(
      "IsomericSMILES", "InChIKey", "InChI", "MolecularWeight",
      "ExactMass", "MolecularFormula", "IUPACName"
    ))

  data <- data %>%
    left_join(pubchem_data, by = "CID") %>%
    rename(SMILES = IsomericSMILES, Formula = MolecularFormula)

  # Optional: Extract CAS
  if (cas) {
    data$CAS_retrieved <- sapply(seq_len(nrow(data)), function(i) {
      message(i)
      cas_result <- suppressWarnings(webchem::pc_sect(data$CID[i], "cas")$Result[1])
      if (length(cas_result) > 0) cas_result else NA
    })
    message("CAS extraction done.")
  }

  # Optional: Extract Flavornet perception
  if (flavornet) {
    message("Extracting flavor from Flavornet.")
    data$Flavornet <- sapply(seq_len(nrow(data)), function(i) {
      message(i)
      if (!is.na(data$CAS_retrieved[i])) {
        fn_result <- suppressWarnings(webchem::fn_percept(data$CAS_retrieved[i]))
        if (length(fn_result) > 0) fn_result else NA_character_
      } else {
        NA_character_
      }
    })
    message("Flavornet extraction done.")
  }

  # Move new columns after CID and drop *_old columns
  if ("CID" %in% colnames(data)) {
    data <- relocate(data, Formula:ncol(data), .after = CID)
  }
  data <- select(data, !contains("_old"))

  return(data)
}



# section -----------------------------------------------------------
df <- rio::import("/Users/sukis/Desktop/project/2022 NQI/data/MSMS_library/NQI_MSMS_database_20250523.xlsx") |> # nolint
  distinct(Name, .keep_all = TRUE) |>
  select(Name, CAS, InChIKey, Ontology)
a <- df[sample(nrow(df), 5), ]

a_cid <- extract_cid(a, name_col = "Name", inchikey_col = "InChIKey", timeout = 100)


library(tidyverse) # Ensure tidyverse is loaded for dplyr and purrr

extract_classyfire <- function(data) {
  # Check if classyfireR is installed
  if (!requireNamespace("classyfireR", quietly = TRUE)) {
    stop("classyfireR is not installed. Please install it first.", call. = FALSE)
  }

  # Check if InChIKey column exists
  if (!("InChIKey" %in% colnames(data))) {
    stop("Your data does not contain InChIKey column. Please ensure it's present.", call. = FALSE)
  }

  # Check if SMILES column exists, which determines if real-time classification is possible
  has_smiles <- "SMILES" %in% colnames(data)
  if (!has_smiles) {
    message("Warning: No 'SMILES' column found. New compounds (not in ClassyFire's cache) cannot be classified in real-time.")
  }

  # Get unique InChIKeys for processing
  unique_inchikeys <- unique(data$InChIKey)

  # --- Step 1: Attempt to retrieve classifications via InChIKey (faster, cached lookup) ---
  message("Attempting to retrieve classifications via InChIKey (cached lookup)...")

  classyfire_cached_results <- list()
  classified_inchikeys_cached <- character(0)

  for (key in unique_inchikeys) {
    result <- tryCatch(
      classyfireR::get_classification(key),
      error = function(e) {
        NULL
      }
    )
    if (!is.null(result)) {
      classyfire_cached_results[[key]] <- result
      classified_inchikeys_cached <- c(classified_inchikeys_cached, key)
      message(paste0("  ✔ Found classification for ", key))
    } else {
      message(paste0("  ✖ No cached classification for ", key))
    }
  }

  # Extract classified data from cached results
  classyfire_meta_cached <- classyfire_cached_results %>%
    purrr::map(classyfireR::meta) %>%
    sapply("[[", 1) %>%
    gsub("InChIKey=", "", .) %>%
    as_tibble() %>%
    rename(InChIKey = value)

  # Apply the corrected extract_cla
  classified_data_cached <- classyfire_cached_results %>%
    purrr::map(classyfireR::classification) %>%
    purrr::map_dfr(extract_cla) %>% # Use map_dfr for direct tibble binding
    cbind(classyfire_meta_cached, .) # Bind InChIKey and classification data

  message(paste0("Found ", nrow(classified_data_cached), " classifications from cache."))

  # --- Step 2: Process compounds not found in cache ---
  unclassified_inchikeys <- setdiff(unique_inchikeys, classified_inchikeys_cached)

  # Initialize all_classified_data with cached results
  all_classified_data <- classified_data_cached

  if (length(unclassified_inchikeys) > 0 && has_smiles) {
    message(paste0("Attempting to submit ", length(unclassified_inchikeys), " new compounds for real-time classification (via SMILES)... This may take some time."))

    # Filter out InChIKey and corresponding SMILES for real-time classification
    to_classify_data <- data %>%
      filter(InChIKey %in% unclassified_inchikeys) %>%
      select(InChIKey, SMILES) %>%
      distinct() # Ensure each InChIKey has only one SMILES

    smiles_for_query_map <- setNames(to_classify_data$SMILES, to_classify_data$InChIKey)

    classyfire_realtime_results <- list()
    classified_inchikeys_realtime <- character(0)

    for (key in names(smiles_for_query_map)) {
      current_smiles <- smiles_for_query_map[[key]]
      message(paste0("  Classifying InChIKey: ", key, " (SMILES: ", current_smiles, ")"))

      # FIX: Use submit_smiles_query instead of submit_query
      job_id <- tryCatch(
        classyfireR::submit_smiles_query(current_smiles),
        error = function(e) {
          warning(paste0("Submission failed for SMILES: ", current_smiles, " (Error: ", e$message, ")"), call. = FALSE)
          NULL
        }
      )
      if (is.null(job_id)) {
        message(paste0("    ✖ Failed to submit for ", key))
        next # Skip to the next InChIKey
      }

      # Polling loop for results with retries
      max_retries <- 10 # Increased retries
      retries <- 0
      result <- NULL
      while (is.null(result) && retries < max_retries) {
        Sys.sleep(3) # Wait for a few seconds
        result <- tryCatch(
          classyfireR::get_classification(job_id),
          error = function(e) {
            NULL # Continue retrying
          }
        )
        if (is.is.null(result)) { # Corrected here: is.null(result)
          message(paste0("    Still waiting for result (retry ", retries + 1, "/", max_retries, ") for job_id: ", job_id))
          retries <- retries + 1
        }
      }

      if (!is.null(result)) {
        classyfire_realtime_results[[key]] <- result
        classified_inchikeys_realtime <- c(classified_inchikeys_realtime, key)
        message(paste0("    ✔ Successfully classified ", key))
      } else {
        message(paste0("    ✖ No classification found for InChIKey: ", key, " after real-time submission and retries."))
      }
    }

    # Extract results from real-time classification
    if (length(classyfire_realtime_results) > 0) {
      classyfire_meta_realtime <- classyfire_realtime_results %>%
        purrr::map(classyfireR::meta) %>%
        sapply("[[", 1) %>%
        gsub("InChIKey=", "", .) %>%
        as_tibble() %>%
        rename(InChIKey = value)

      # Apply the corrected extract_cla
      classified_data_realtime <- classyfire_realtime_results %>%
        purrr::map(classyfireR::classification) %>%
        purrr::map_dfr(extract_cla) %>% # Use map_dfr for direct tibble binding
        cbind(classyfire_meta_realtime, .) # Bind InChIKey and classification data

      message(paste0("Classified ", nrow(classified_data_realtime), " new compounds in real-time."))

      # Combine cached and real-time classification results
      all_classified_data <- bind_rows(all_classified_data, classified_data_realtime)
    } else {
      message("No new compounds were successfully classified in real-time.")
    }

  } else if (length(unclassified_inchikeys) > 0 && !has_smiles) {
    message(paste0("Skipping real-time classification for ", length(unclassified_inchikeys), " compounds. 'SMILES' column is missing."))
  } else {
    message("No new compounds to classify.")
  }

  # --- Step 3: Merge back into the original data ---
  # Use left_join to ensure all original data rows are preserved
  final_data <- left_join(data, all_classified_data, by = "InChIKey")

  return(final_data)
}

# Corrected extract_cla Function
extract_cla <- function(classification_list) {
  # Define all expected classification levels
  expected_levels <- c("Kingdom", "Superclass", "Class", "Subclass", "DirectParent")

  # Initialize a list to store extracted names
  extracted_names <- list()

  for (level in expected_levels) {
    if (!is.null(classification_list[[level]]) && !is.null(classification_list[[level]]$name)) {
      extracted_names[[level]] <- classification_list[[level]]$name
    } else {
      extracted_names[[level]] <- NA_character_ # Assign NA if level or name is missing
    }
  }

  # Convert the list to a tibble (single row)
  # Ensure column order is consistent
  as_tibble(extracted_names)[, expected_levels]
}

my_compounds_data <- tibble(
  Name = c("CompoundA", "CompoundB", "NewPolymer1", "NewPolymer2"),
  InChIKey = c("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", # Caffeine
               "FMAONNMJHLVSHI-UHFFFAOYSA-N", # Ethanol
               "XABCCDEFGXXXXXX-UHFFFAOYSA-N", # An imaginary InChIKey
               "ZYXWVUTSRPQMLK-UHFFFAOYSA-N"), # Another imaginary InChIKey
  SMILES = c("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # Caffeine SMILES
             "CCO", # Ethanol SMILES
             "CCCCCCCCCCC(O)CCOCCOCCO", # Simple oligomer with an OH (ClassyFire might classify alcohols)
             "C(CC)CC(C)CC(C)C") # Another simple oligomer
)

classified_data <- extract_classyfire(my_compounds_data)
print(classified_data)
