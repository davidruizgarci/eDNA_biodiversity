#-------------------------------------------------------------------------------
#
# 6. Prepare particle backtracking polygons
#
#-------------------------------------------------------------------------------

library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)


# 1. Set paths------------------------------------------------------------------

# eDNA presence/absence dataset
file_path <- file.path(input_data, "processed_df")
data_file <- file.path(file_path, "merged_eDNA_spp_pa.csv")

# Particle backtracking outputs
particles_dir <- "C:/Users/David/SML Dropbox/gitdata/edna-particles/output"
polygons_dir <- file.path(particles_dir, "Poligonos")
trajectories_file <- file.path(particles_dir, "Trajectories_index.csv")

# Outputs from this script
backtracking_output_dir <- file.path(output_data, "particle_backtracking")
dir.create(
  backtracking_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Species requiring a depth-distribution check-------------------------------

species_to_check <- unique(c(
  "Anguilla anguilla",
  "Chelon auratus",
  "Chelon ramada",
  "Dicentrarchus labrax",
  "Bathypterois dubius",
  "Parablennius zvonimiri",
  "Thunnus alalunga",
  "Dasyatis pastinaca",
  "Gymnura altavela",
  "Myliobatis aquila",
  "Centrophorus uyato",
  "Pteroplatytrygon violacea",
  "Chimaera monstrosa",
  "Ariosoma balearicum",
  "Balaenoptera acutorostrata",
  "Balaenoptera physalus",
  "Callanthias ruber",
  "Capros aper",
  "Centrolophus niger",
  "Coris julis",
  "Crystallogobius linearis",
  "Dentex dentex",
  "Echelus myrus",
  "Globicephala melas",
  "Gobius vittatus",
  "Lepadogaster candolii",
  "Lepidotrigla cavillone",
  "Macroramphosus scolopax",
  "Maurolicus muelleri",
  "Nerophis maculatus",
  "Pegusa lascaris",
  "Serranus scriba",
  "Sphyraena viridensis",
  "Spicara maena",
  "Stomias boa",
  "Tetragonurus cuvieri",
  "Notoscopelus kroyeri",
  "Trachipterus arcticus",
  "Trachyscorpia cristulata",
  "Pseudaphya ferreri",
  "Dalatias licha",
  "Centracanthus cirrus",
  "Brama brama"
))


# 3. Helper functions-----------------------------------------------------------

# Standardise names before matching samples or files
standardise_key <- function(x) {
  x %>%
    as.character() %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

# Find a column using a list of possible names
find_column <- function(data, candidates, label, required = TRUE) {
  
  names_clean <- names(data) %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  
  match_position <- match(candidates, names_clean, nomatch = 0)
  match_position <- match_position[match_position > 0]
  
  if (length(match_position) == 0) {
    
    if (!required) return(NA_character_)
    
    stop(
      "Could not identify the ", label, " column in Trajectories_index.csv.\n",
      "Available columns: ", paste(names(data), collapse = ", ")
    )
  }
  
  names(data)[match_position[1]]
}

# Read every layer from one GeoPackage
read_gpkg_layers <- function(gpkg_file) {
  
  layers <- st_layers(gpkg_file)$name
  
  map_dfr(
    layers,
    function(layer_name) {
      
      st_read(
        gpkg_file,
        layer = layer_name,
        quiet = TRUE
      ) %>%
        mutate(
          polygon_file = basename(gpkg_file),
          polygon_layer = layer_name,
          .before = 1
        )
    }
  )
}


# 4. Load eDNA detections--------------------------------------------------------

if (!file.exists(data_file)) {
  stop("The eDNA dataset was not found: ", data_file)
}

data <- read.csv(
  data_file,
  sep = ";",
  check.names = FALSE
) %>%
  filter(
    !sequencing_comment %in% c(
      "Filtration_Control",
      "No_target_taxa"
    )
  )

# Dataset species columns use the format s_species_name
species_lookup <- tibble(
  species = species_to_check,
  species_column = paste0(
    "s_",
    species %>%
      str_to_lower() %>%
      str_replace_all("[^a-z0-9]+", "_")
  )
)

# Check which requested species are present as columns in the dataset
species_lookup <- species_lookup %>%
  mutate(column_found = species_column %in% names(data))

species_missing_from_data <- species_lookup %>%
  filter(!column_found)

if (nrow(species_missing_from_data) > 0) {
  warning(
    "The following species columns were not found in merged_eDNA_spp_pa.csv: ",
    paste(species_missing_from_data$species, collapse = ", ")
  )
}

species_columns_found <- species_lookup %>%
  filter(column_found) %>%
  pull(species_column)

if (length(species_columns_found) == 0) {
  stop("None of the requested species were found in merged_eDNA_spp_pa.csv.")
}

# One row per species detection and sample
species_detections <- data %>%
  select(
    sample_id,
    any_of(c(
      "original_sampling_id",
      "site",
      "study_area",
      "campaign",
      "latitude",
      "longitude"
    )),
    all_of(species_columns_found)
  ) %>%
  pivot_longer(
    cols = all_of(species_columns_found),
    names_to = "species_column",
    values_to = "presence"
  ) %>%
  filter(!is.na(presence), presence > 0) %>%
  left_join(
    species_lookup %>%
      select(species, species_column),
    by = "species_column"
  ) %>%
  distinct(species, sample_id, .keep_all = TRUE) %>%
  mutate(sample_key = standardise_key(sample_id))


# 5. Load the trajectory index--------------------------------------------------

if (!file.exists(trajectories_file)) {
  stop(
    "Trajectories_index.csv was not found: ", trajectories_file, "\n",
    "Place the file in the particle output folder or update trajectories_file."
  )
}

# read_csv() automatically detects comma-separated CSV files
trajectories_index <- read_csv(
  trajectories_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

# Detect relevant columns while retaining their original names
trajectory_sample_col <- find_column(
  trajectories_index,
  candidates = c(
    "sample",
    "sample_id",
    "sampleid",
    "sampling_id",
    "original_sample_id",
    "original_sampling_id"
  ),
  label = "sample"
)

trajectory_particles_col <- find_column(
  trajectories_index,
  candidates = c(
    "particles",
    "particle",
    "particle_file",
    "particles_file",
    "trajectory_file",
    "trajectories_file",
    "netcdf",
    "nc_file",
    "file"
  ),
  label = "particle file",
  required = FALSE
)

trajectory_polygon_col <- find_column(
  trajectories_index,
  candidates = c(
    "polygon",
    "polygons",
    "polygon_file",
    "polygons_file",
    "gpkg",
    "gpkg_file"
  ),
  label = "polygon file",
  required = FALSE
)

if (
  is.na(trajectory_particles_col) &&
  is.na(trajectory_polygon_col)
) {
  stop(
    "Trajectories_index.csv must contain either a particle-file column ",
    "or a polygon-file column."
  )
}

# Standardise index fields for matching
trajectories_index_clean <- trajectories_index %>%
  mutate(
    trajectory_row = row_number(),
    trajectory_sample = .data[[trajectory_sample_col]],
    particle_file = if (
      !is.na(trajectory_particles_col)
    ) {
      as.character(.data[[trajectory_particles_col]])
    } else {
      NA_character_
    },
    polygon_file_index = if (
      !is.na(trajectory_polygon_col)
    ) {
      as.character(.data[[trajectory_polygon_col]])
    } else {
      NA_character_
    },
    sample_key = standardise_key(trajectory_sample),
    particle_key = standardise_key(particle_file),
    polygon_key_index = standardise_key(polygon_file_index)
  )


# 6. Match detections to trajectories-------------------------------------------

detection_trajectories <- species_detections %>%
  left_join(
    trajectories_index_clean %>%
      select(
        trajectory_row,
        trajectory_sample,
        particle_file,
        polygon_file_index,
        sample_key,
        particle_key,
        polygon_key_index
      ),
    by = "sample_key",
    relationship = "many-to-many"
  ) %>%
  mutate(
    trajectory_found = !is.na(trajectory_row),
    particle_file_available = !is.na(particle_file) &
      particle_file != ""
  )


# 7. Match trajectories to GeoPackage files------------------------------------

if (!dir.exists(polygons_dir)) {
  stop("The polygon folder was not found: ", polygons_dir)
}

polygon_files <- list.files(
  polygons_dir,
  pattern = "\\.gpkg$",
  full.names = TRUE,
  ignore.case = TRUE,
  recursive = TRUE
)

if (length(polygon_files) == 0) {
  stop("No GeoPackage files were found in: ", polygons_dir)
}

polygon_index <- tibble(
  polygon_path = polygon_files,
  polygon_file = basename(polygon_files),
  polygon_key = standardise_key(polygon_files)
)

# First use an explicit polygon name from Trajectories_index, when available.
# Otherwise match the polygon and particle filenames without their extensions.
detection_polygon_index <- detection_trajectories %>%
  mutate(
    expected_polygon_key = case_when(
      !is.na(polygon_file_index) &
        polygon_file_index != "" ~ polygon_key_index,
      particle_file_available ~ particle_key,
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(
    polygon_index,
    by = c("expected_polygon_key" = "polygon_key"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    polygon_found = !is.na(polygon_path),
    check_status = case_when(
      !trajectory_found ~ "No trajectory-index match",
      !particle_file_available &
        is.na(polygon_file_index) ~ "No backtracking file (e.g. ferry sample)",
      !polygon_found ~ "Trajectory found but GeoPackage missing",
      TRUE ~ "GeoPackage found"
    )
  )


# 8. Check duplicated or ambiguous matches--------------------------------------

ambiguous_sample_matches <- detection_polygon_index %>%
  filter(trajectory_found) %>%
  count(species, sample_id, name = "n_index_matches") %>%
  filter(n_index_matches > 1)

ambiguous_polygon_keys <- polygon_index %>%
  count(polygon_key, name = "n_polygon_files") %>%
  filter(n_polygon_files > 1)

if (nrow(ambiguous_sample_matches) > 0) {
  warning(
    "Some species/sample combinations matched more than one trajectory-index row. ",
    "Check ambiguous_sample_matches."
  )
}

if (nrow(ambiguous_polygon_keys) > 0) {
  warning(
    "Some standardised polygon filenames are duplicated. ",
    "Check ambiguous_polygon_keys."
  )
}


# 9. Read matched polygons------------------------------------------------------

matched_polygon_files <- detection_polygon_index %>%
  filter(polygon_found) %>%
  distinct(polygon_path) %>%
  pull(polygon_path)

if (length(matched_polygon_files) > 0) {
  
  backtracking_polygons_raw <- map_dfr(
    matched_polygon_files,
    read_gpkg_layers
  )
  
  # Add the species and sample information associated with each polygon
  backtracking_polygons <- backtracking_polygons_raw %>%
    mutate(polygon_key = standardise_key(polygon_file)) %>%
    left_join(
      detection_polygon_index %>%
        filter(polygon_found) %>%
        transmute(
          species,
          sample_id,
          trajectory_sample,
          particle_file,
          polygon_key = expected_polygon_key
        ) %>%
        distinct(),
      by = "polygon_key",
      relationship = "many-to-many"
    ) %>%
    relocate(
      species,
      sample_id,
      trajectory_sample,
      particle_file,
      polygon_file,
      polygon_layer
    )
  
  # Use a common geographic coordinate system for later IUCN comparisons
  backtracking_polygons <- st_transform(
    backtracking_polygons,
    4326
  )
  
} else {
  
  backtracking_polygons <- NULL
  
  warning(
    "No matching GeoPackage files were found. ",
    "Inspect detection_polygon_index and the filename keys."
  )
}


# 10. Summarise checks----------------------------------------------------------

backtracking_summary <- detection_polygon_index %>%
  count(species, check_status, name = "n_samples") %>%
  complete(
    species,
    check_status = c(
      "GeoPackage found",
      "No trajectory-index match",
      "No backtracking file (e.g. ferry sample)",
      "Trajectory found but GeoPackage missing"
    ),
    fill = list(n_samples = 0)
  ) %>%
  arrange(species, check_status)

backtracking_summary_wide <- backtracking_summary %>%
  pivot_wider(
    names_from = check_status,
    values_from = n_samples,
    values_fill = 0
  )

backtracking_summary_wide


# 11. Export quality-control tables and polygons--------------------------------

write.csv2(
  species_lookup,
  file.path(
    backtracking_output_dir,
    "species_column_check.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv2(
  detection_polygon_index,
  file.path(
    backtracking_output_dir,
    "species_sample_polygon_index.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv2(
  backtracking_summary_wide,
  file.path(
    backtracking_output_dir,
    "species_polygon_check_summary.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv2(
  ambiguous_sample_matches,
  file.path(
    backtracking_output_dir,
    "ambiguous_sample_matches.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv2(
  ambiguous_polygon_keys,
  file.path(
    backtracking_output_dir,
    "ambiguous_polygon_keys.csv"
  ),
  row.names = FALSE,
  na = ""
)

if (!is.null(backtracking_polygons)) {
  
  polygons_output_file <- file.path(
    backtracking_output_dir,
    "species_backtracking_polygons.gpkg"
  )
  
  if (file.exists(polygons_output_file)) {
    file.remove(polygons_output_file)
  }
  
  st_write(
    backtracking_polygons,
    polygons_output_file,
    layer = "backtracking_polygons",
    quiet = TRUE
  )
}


# 12. Final checks--------------------------------------------------------------

cat(
  "\nSpecies requested:", length(species_to_check),
  "\nSpecies columns found:", sum(species_lookup$column_found),
  "\nSpecies detections:", nrow(species_detections),
  "\nDetections with GeoPackage:", sum(detection_polygon_index$polygon_found),
  "\nDetections without GeoPackage:", sum(!detection_polygon_index$polygon_found),
  "\n\nOutputs saved in:\n", backtracking_output_dir, "\n"
)
