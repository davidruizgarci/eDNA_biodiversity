#-------------------------------------------------------------------------------
#
# 2. Check species taxonomy
#
#-------------------------------------------------------------------------------

library(dplyr)
#install.packages("worrms")
library(purrr)
library(worrms)
library(tidyr)
library(stringr)

# 1. Load data------------------------------------------------------------------
file_path <- paste0(input_data, "/processed_df")
file <- paste0(file_path, "/merged_eDNA_spp_pa.csv")

data <- read.csv(file, sep = ";")
head(data)

data <- data %>%
  filter(!sequencing_comment %in% c("Filtration_Control", "No_target_taxa"))


# 2. Make species list----------------------------------------------------------
# Seleccionar solo las columnas de especies
species_data <- data %>%
  select(s_tursiops_truncatus:last_col())

# Especies detectadas al menos una vez
detected_species <- names(species_data)[colSums(species_data > 0, na.rm = TRUE) > 0]

# Número total de especies detectadas
length(detected_species)

# Listado de especies
detected_species

cat("Total de especies detectadas:", length(detected_species), "\n\n")
cat(paste(detected_species, collapse = "\n"))


# Clean species names
detected_species_clean <- detected_species |>
  gsub("^s_", "", x = _) |>
  gsub("_", " ", x = _)

# Capitalize genus
detected_species_clean <- sub(
  "^([a-z])",
  "\\U\\1",
  detected_species_clean,
  perl = TRUE
)

# Remove duplicates created after cleaning
detected_species_clean <- unique(detected_species_clean)

# View
cat(paste(detected_species_clean, collapse = "\n"))


# 3. Check species ID against WoRMS---------------------------------------------
# Starting from your cleaned species list:
# detected_species_clean
check_worms_one <- function(sp) {
  
  res <- tryCatch(
    wm_records_name(name = sp, marine_only = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(res) || nrow(res) == 0) {
    return(tibble(
      input_name = sp,
      worms_found = FALSE,
      aphia_id = NA,
      worms_name = NA,
      status = NA,
      valid_name = NA,
      valid_aphia_id = NA,
      deprecated_or_synonym = NA
    ))
  }
  
  # Keep first/best returned match
  res <- res[1, ]
  
  tibble(
    input_name = sp,
    worms_found = TRUE,
    aphia_id = res$AphiaID,
    worms_name = res$scientificname,
    status = res$status,
    valid_name = res$valid_name,
    valid_aphia_id = res$valid_AphiaID,
    deprecated_or_synonym = status != "accepted"
  )
}

worms_check <- map_dfr(detected_species_clean, check_worms_one)

worms_check

# Inspect problematic names:
worms_problematic <- worms_check %>%
  filter(
    worms_found == FALSE |
      deprecated_or_synonym == TRUE |
      input_name != valid_name
  )

worms_problematic

# Inspect those not found in WoRMS:
worms_not_found <- worms_check %>%
  filter(worms_found == FALSE)

worms_not_found

# Inspect those which name may have changed:
worms_deprecated <- worms_check %>%
  filter(worms_found == TRUE, status != "accepted")

worms_deprecated

# Correct those whose valid name has been detected:
accepted_species <- worms_check %>%
  mutate(final_name = if_else(!is.na(valid_name), valid_name, input_name)) %>%
  distinct(final_name) %>%
  arrange(final_name)

accepted_species

# Correct those whose name is not clear:
head(worms_not_found$input_name)

# Manual corrections/removals for species not resolved by WoRMS
species_list <- accepted_species %>%
  mutate(
    final_name = case_when(
      final_name == "Centrophorus cf uyato rc 2022" ~ "Centrophorus uyato",
      TRUE ~ final_name
    )
  ) %>%
  filter(
    !final_name %in% c(
      "Engraulis environmental sample",
      "Trachurus environmental sample",
      "Triturus marmoratus",
      "Autographa gamma",
      "Vesiculariidae gen n sp n aw 2011"
    )
  ) %>%
  distinct(final_name) %>%
  arrange(final_name)

species_list
nrow(species_list)
cat(paste(species_list$final_name, collapse = "\n"))

# Let's now filter only marine vertabrates:
# Detect AphiaID column 
get_worms_record <- function(sp) {
  
  res <- tryCatch(
    wm_records_name(name = sp, marine_only = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(res) || nrow(res) == 0) {
    return(tibble(
      final_name = sp,
      aphia_id = NA_integer_,
      worms_name = NA_character_,
      status = NA_character_,
      valid_name = sp,
      valid_aphia_id = NA_integer_
    ))
  }
  
  res <- res[1, ]
  
  tibble(
    final_name = sp,
    aphia_id = res$AphiaID,
    worms_name = res$scientificname,
    status = res$status,
    valid_name = res$valid_name,
    valid_aphia_id = res$valid_AphiaID
  )
}

species_list_worms <- species_list %>%
  mutate(final_name = as.character(final_name)) %>%
  pull(final_name) %>%
  map_dfr(get_worms_record)

species_list_worms

# Function to get WoRMS classification
get_worms_classification <- function(aphia_id) {
  
  if (is.na(aphia_id)) return(NULL)
  
  tryCatch(
    wm_classification(id = aphia_id),
    error = function(e) NULL
  )
}

extract_rank <- function(classification, rank_name) {
  
  if (is.null(classification)) return(NA_character_)
  
  out <- classification %>%
    filter(rank == rank_name) %>%
    pull(scientificname)
  
  if (length(out) == 0) NA_character_ else out[1]
}

# Apply:
species_list_taxonomy <- species_list_worms %>%
  mutate(
    final_aphia_id = if_else(!is.na(valid_aphia_id), valid_aphia_id, aphia_id),
    classification = map(final_aphia_id, get_worms_classification),
    class  = map_chr(classification, extract_rank, "Class"),
    order  = map_chr(classification, extract_rank, "Order"),
    family = map_chr(classification, extract_rank, "Family")
  ) %>%
  select(-classification)

species_list_taxonomy

#species_list_taxonomy %>%
#  filter(valid_name %in% c("Caretta caretta", "Chelonia mydas", "Dermochelys coriacea")) %>%
#  select(valid_name, class, order, family)

# Filter only marine vertebrates:
marine_vertebrates <- species_list_taxonomy %>%
  filter(
    class %in% c("Teleostei", "Elasmobranchii", "Mammalia") |
      valid_name %in% c("Caretta caretta", "Chelonia mydas", "Dermochelys coriacea") |
      class == "Aves"
  )

cat(paste(marine_vertebrates$final_name, collapse = "\n"))

terrestrial_birds_to_remove <- c(
  "Ardea cinerea",
  "Bubo bubo",
  "Columba livia",
  "Erithacus rubecula",
  "Gallus gallus",
  "Meleagris gallopavo",
  "Parus major",
  "Phoenicurus phoenicurus",
  "Scolopax rusticola",
  "Sylvia atricapilla",
  "Oryctolagus cuniculus" # conejo
)

marine_vertebrates_clean <- marine_vertebrates %>%
  filter(!valid_name %in% terrestrial_birds_to_remove) %>%
  arrange(class, order, family, valid_name)

cat("Total marine vertebrates:", nrow(marine_vertebrates_clean), "\n\n")
cat(paste(marine_vertebrates_clean$valid_name, collapse = "\n"))

## Rare species:
## Fresh water: 
#Lepomis gibbosus
#Leuciscus leuciscus
#Oncorhynchus mykiss
#Salmo salar

## From other regions:
#Calonectris borealis
#Engraulis australis
#Engraulis japonicus
#Thunnus tonggol
#Delphinapterus leucas

# Export
file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/speciesList.csv")

write.csv2(marine_vertebrates_clean, file, row.names = FALSE,  na = "")



# 4. Change species names to correct them---------------------------------------
head(data)
head(marine_vertebrates_clean)

data_renamed <- data

meta_cols <- names(data_renamed)[1:which(names(data_renamed) == "original_sample_id")]
species_cols <- setdiff(names(data_renamed), meta_cols)

new_species_names <- species_cols %>%
  str_remove("^s_") %>%
  str_replace_all("_", " ") %>%
  str_replace("^([a-z])", str_to_upper(str_sub(., 1, 1)))

new_species_names <- case_when(
  new_species_names == "Centrophorus cf uyato rc 2022" ~ "Centrophorus uyato",
  new_species_names == "Oblada melanura" ~ "Oblada melanurus",
  new_species_names == "Pomatomus saltator" ~ "Pomatomus saltatrix",
  TRUE ~ new_species_names
)

names(data_renamed)[match(species_cols, names(data_renamed))] <- new_species_names

dim(data)
dim(data_renamed)

# Export
file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/renamed_dataset.csv")

write.csv2(data_renamed, file, row.names = FALSE,  na = "")


# 5. Metadata file to make filtration later----------------------------------------
# Metadata from original species columns + WoRMS correction + taxonomy + target info
species_metadata <- tibble(
  original_colname = detected_species,
  input_name = detected_species %>%
    str_remove("^s_") %>%
    str_replace_all("_", " ") %>%
    str_replace("^([a-z])", str_to_upper(str_sub(., 1, 1)))
) %>%
  left_join(
    worms_check,
    by = "input_name"
  ) %>%
  mutate(
    corrected_name = if_else(!is.na(valid_name), valid_name, input_name),
    corrected_name = case_when(
      corrected_name == "Centrophorus cf uyato rc 2022" ~ "Centrophorus uyato",
      TRUE ~ corrected_name
    ),
    removed_manual = corrected_name %in% c(
      "Engraulis environmental sample",
      "Trachurus environmental sample",
      "Triturus marmoratus",
      "Autographa gamma",
      "Vesiculariidae gen n sp n aw 2011"
    )
  ) %>%
  left_join(
    species_list_taxonomy %>%
      select(
        final_name,
        valid_name_taxonomy = valid_name,
        final_aphia_id,
        class,
        order,
        family
      ),
    by = c("corrected_name" = "final_name")
  ) %>%
  mutate(
    TaxaTargeted = if_else(
      corrected_name %in% marine_vertebrates_clean$valid_name,
      "Yes",
      "No"
    )
  ) %>%
  select(
    original_colname,
    input_name,
    worms_found,
    aphia_id,
    worms_name,
    status,
    valid_name,
    valid_aphia_id,
    corrected_name,
    removed_manual,
    TaxaTargeted,
    final_aphia_id,
    class,
    order,
    family
  )

species_metadata

# Export
file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/metadataFiltering.csv")

write.csv2(species_metadata, file, row.names = FALSE,  na = "")
