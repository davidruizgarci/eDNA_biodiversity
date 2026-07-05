#-------------------------------------------------------------------------------
#
# 1. Prepare eDNA dataset
#
#-------------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(writexl)
library(hms)

# 1. Load data------------------------------------------------------------------
file <- paste0(input_data, "/orginal_data/final/_eDNA_data_biodiversity.xlsx")

# 1.1. Read sheets
asv <- read_excel(file, sheet = "All_ASV")
tax <- read_excel(file, sheet = "All_Taxa")
samples <- read_excel(file, sheet = "All_Samples")

# 2. Clean taxonomy table-------------------------------------------------------
tax_clean <- tax %>%
  rename(ASV_ID = ASVs) %>% #`#OTU ID`
  mutate(
    species = str_remove(S_taxonomy, "^\\s*s__"),
    species = str_squish(species) ) %>%
  filter(!is.na(species), species != "")

# Check duplicated ASV IDs in taxonomy
#tax_duplicates <- tax_clean %>%
#  count(ASV_ID, sort = TRUE) %>%
#  filter(n > 1)
#
#tax_duplicates
#
#tax_duplicate_summary <- tax_clean %>%
#  semi_join(tax_duplicates, by = "ASV_ID") %>%
#  group_by(ASV_ID) %>%
#  summarise(
#    n_rows = n(),
#    n_species = n_distinct(species),
#    species_list = paste(unique(species), collapse = "; "),
#    .groups = "drop"
#  ) %>%
#  arrange(desc(n_species), desc(n_rows))
#
#tax_duplicate_summary
#
#asv_duplicates <- asv %>%
#  count(ASVs, sort = TRUE) %>%
#  filter(n > 1)
#
#asv_duplicates

# 2.1. Convert ASV table to long format
asv_long <- asv %>%
  rename(ASV_ID = ASVs) %>%
  select(
    -Order,
    -Sequencing_Run,
    -`#OTU ID`
  ) %>%
  pivot_longer(
    cols = -ASV_ID,
    names_to = "Sample_ID",
    values_to = "reads"
  ) %>%
  mutate(
    reads = as.numeric(reads),
    present = reads > 0
  )

# 2.2. Join ASVs with species names
detections_long <- asv_long %>%
  left_join(
    tax_clean %>%
      select(ASV_ID, species),
    by = "ASV_ID") %>%
  filter(
    !is.na(species),
    present)


# 3. Create presence/absence matrix by sample and species-----------------------
species_pa <- detections_long %>%
  distinct(Sample_ID, species) %>% #just in case there is one species with more than one ASV asigned to it.
  mutate(present = 1) %>%
  pivot_wider(
    names_from = species,
    values_from = present,
    values_fill = 0)

# 3.1. Prepare sample metadata
# Sample names present in the ASV table
asv_sample_ids <- names(asv) %>%
  setdiff(c("Order", "Sequencing_Run", "ASVs", "#OTU ID"))

# Metadata IDs
meta_asv_ids <- samples$Sample_ID_ASV
meta_original_ids <- samples$Original_sample_ID

# How many ASV table samples match Sample_ID_ASV?
sum(asv_sample_ids %in% meta_asv_ids)

# How many ASV table samples match Original_sample_ID?
sum(asv_sample_ids %in% meta_original_ids)

# ASV samples with no metadata match
asv_sample_ids[!asv_sample_ids %in% meta_asv_ids]

# Metadata rows with missing Sample_ID_ASV
samples %>%
  filter(is.na(Sample_ID_ASV)) %>%
  select(Original_sample_ID, Sample_ID_ASV, Site, StudyArea, Campaign)
names(samples)

samples_clean <- samples %>%
  mutate(
    Sample_ID = case_when(
      !is.na(Sample_ID_ASV) ~ Sample_ID_ASV,
      Original_sample_ID %in% asv_sample_ids ~ Original_sample_ID,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Sample_ID)) %>% #REMOVE SAMPLES THAT DID NOT SEARCH FOR VERTABRATES OR WERE CONTROLS.
  select(
    Sample_ID,
    'Sequencing comment',
    Original_sample_ID,
    Site,
    StudyArea,
    Campaign,
    Sampling_Date,
    Year,
    Season,
    Day_Night,
    Type,
    Project,
    Latitude,
    Longitude,
    Filt_Vol
  )





# 3.2. Combine metadata + species presence/absence
edna_final <- samples_clean %>%
  left_join(
    species_pa,
    by = "Sample_ID")

species_cols <- setdiff(names(edna_final), names(samples_clean))

edna_final[species_cols] <- lapply(
  edna_final[species_cols],
  function(x) tidyr::replace_na(x, 0)
)

edna_final <- edna_final %>%
  janitor::clean_names()

# Export
file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/eDNA_spp_pa.csv")

write.csv2(edna_final, file, row.names = FALSE,  na = "")


# 5. Merge different filters into one-------------------------------------------
# Metadata columns to keep
names(edna_final)
metadata_cols <- c(
  "sample_id",
  "sequencing_comment",
  "original_sampling_id",
  "site",
  "study_area",
  "campaign",
  "sampling_date",
  "year",
  "season",
  "day_night",
  "type",
  "project",
  "latitude",
  "longitude",
  "filt_vol"
)

metadata_cols <- intersect(metadata_cols, names(edna_final))
species_cols <- setdiff(names(edna_final), metadata_cols)

# Merge columns based on naming, for example:
# 03-06-B-R1-16S & 03-06-B-R2-16Sm belong to the same sample but they are different 
# filters and this is reflected by just a change in the number after the R (R1 and R2)


edna_merged <- edna_final %>%
  mutate(sample_id = str_remove(sample_id, "_[0-9]+$")) %>%
  group_by(sample_id) %>%
  summarise(
    across(all_of(species_cols), ~ max(.x, na.rm = TRUE)),
    across(all_of(setdiff(metadata_cols, "sample_id")), ~ dplyr::first(na.omit(.x))),
    .groups = "drop"
  ) %>%
  select(
    all_of(metadata_cols),
    all_of(species_cols)
  )

# Export
file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/merged_eDNA_spp_pa.csv")

write.csv2(edna_merged, file, row.names = FALSE,  na = "")
