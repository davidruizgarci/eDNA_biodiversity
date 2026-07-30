#-------------------------------------------------------------------------------
#
# 7. Check whether particle backtracking polygons overlap IUCN ranges
#
#-------------------------------------------------------------------------------

library(sf)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)


# 1. Set paths------------------------------------------------------------------

backtracking_dir <- file.path(output_data, "particle_backtracking")
polygon_file <- file.path(
  backtracking_dir,
  "species_backtracking_polygons.gpkg"
)
polygon_index_file <- file.path(
  backtracking_dir,
  "species_sample_polygon_index.csv"
)
iucn_dir <- file.path(input_data, "IUCNmaps")

iucn_output_dir <- file.path(
  backtracking_dir,
  "IUCN_overlap"
)
dir.create(iucn_output_dir, recursive = TRUE, showWarnings = FALSE)


# 2. Load backtracking polygons and sample information--------------------------

if (!file.exists(polygon_file)) {
  stop("Run scr/6_particle_backtracking_polygons.R first.")
}

backtracking_polygons <- st_read(
  polygon_file,
  layer = "backtracking_polygons",
  quiet = TRUE
) %>%
  st_transform(4326) %>%
  filter(!is.na(species), !is.na(sample_id)) %>%
  st_make_valid()

polygon_index <- read.csv(
  polygon_index_file,
  sep = ";",
  check.names = FALSE
)

sample_information <- polygon_index %>%
  filter(polygon_found) %>%
  select(
    species,
    sample_id,
    latitude,
    longitude
  ) %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  ) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  distinct()

species_to_check <- sort(unique(backtracking_polygons$species))


# 3. Read IUCN polygons for the relevant species--------------------------------

iucn_files <- c(
  list.files(
    file.path(iucn_dir, "MAMMALS_MARINE_ONLY"),
    pattern = "\\.shp$",
    full.names = TRUE,
    ignore.case = TRUE
  ),
  list.files(
    file.path(iucn_dir, "TURTLES"),
    pattern = "\\.shp$",
    full.names = TRUE,
    ignore.case = TRUE
  ),
  list.files(
    file.path(iucn_dir, "MARINEFISH"),
    pattern = "\\.shp$",
    full.names = TRUE,
    ignore.case = TRUE
  )
)

if (length(iucn_files) == 0) {
  stop("No IUCN shapefiles were found in: ", iucn_dir)
}

read_relevant_iucn <- function(shp_file) {
  
  x <- st_read(shp_file, quiet = TRUE)
  
  species_col <- names(x)[
    str_to_lower(names(x)) %in% c(
      "sci_name",
      "binomial",
      "scientific",
      "species"
    )
  ][1]
  
  if (is.na(species_col)) return(NULL)
  
  x %>%
    rename(species = all_of(species_col)) %>%
    filter(species %in% species_to_check) %>%
    mutate(iucn_file = basename(shp_file)) %>%
    st_transform(4326)
}

iucn_polygons <- map(iucn_files, read_relevant_iucn) %>%
  compact() %>%
  bind_rows() %>%
  filter(!st_is_empty(geometry)) %>%
  st_make_valid()

iucn_species_found <- sort(unique(iucn_polygons$species))
iucn_species_missing <- setdiff(species_to_check, iucn_species_found)


# 4. Identify samples whose original point was outside the IUCN range-----------

sample_points <- sample_information %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

check_original_point <- function(sp, sample) {
  
  point_sp <- sample_points %>%
    filter(species == sp, sample_id == sample)
  
  range_sp <- iucn_polygons %>%
    filter(species == sp)
  
  tibble(
    species = sp,
    sample_id = sample,
    iucn_map_found = nrow(range_sp) > 0,
    original_point_overlaps_iucn = if (
      nrow(range_sp) == 0
    ) {
      NA
    } else {
      any(lengths(st_intersects(point_sp, range_sp)) > 0)
    }
  )
}

point_iucn_overlap <- map2_dfr(
  sample_points$species,
  sample_points$sample_id,
  check_original_point
)

samples_outside_iucn <- point_iucn_overlap %>%
  filter(
    iucn_map_found,
    original_point_overlaps_iucn == FALSE
  )


# 5. Check whether the backtracking polygon reaches the IUCN range--------------

check_backtracking_iucn <- function(sp, sample) {
  
  particle_sp <- backtracking_polygons %>%
    filter(species == sp, sample_id == sample)
  
  range_sp <- iucn_polygons %>%
    filter(species == sp)
  
  tibble(
    species = sp,
    sample_id = sample,
    n_particle_features = nrow(particle_sp),
    backtracking_overlaps_iucn = any(
      lengths(st_intersects(particle_sp, range_sp)) > 0
    )
  )
}

iucn_backtracking_overlap <- map2_dfr(
  samples_outside_iucn$species,
  samples_outside_iucn$sample_id,
  check_backtracking_iucn
) %>%
  left_join(
    sample_information,
    by = c("species", "sample_id")
  ) %>%
  arrange(species, sample_id)


# 6. Summarise and export results-----------------------------------------------

iucn_backtracking_summary <- iucn_backtracking_overlap %>%
  group_by(species) %>%
  summarise(
    n_original_points_outside = n(),
    n_backtracking_overlap = sum(backtracking_overlaps_iucn),
    any_backtracking_overlap = any(backtracking_overlaps_iucn),
    .groups = "drop"
  )

write.csv2(
  point_iucn_overlap,
  file.path(iucn_output_dir, "original_point_IUCN_overlap.csv"),
  row.names = FALSE,
  na = ""
)

write.csv2(
  iucn_backtracking_overlap,
  file.path(iucn_output_dir, "backtracking_IUCN_overlap_by_sample.csv"),
  row.names = FALSE,
  na = ""
)

write.csv2(
  iucn_backtracking_summary,
  file.path(iucn_output_dir, "backtracking_IUCN_overlap_summary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv2(
  tibble(species = iucn_species_missing),
  file.path(iucn_output_dir, "species_without_IUCN_map.csv"),
  row.names = FALSE,
  na = ""
)

cat(
  "\nOriginal detections outside IUCN:", nrow(samples_outside_iucn),
  "\nBacktracking polygons reaching IUCN:",
  sum(iucn_backtracking_overlap$backtracking_overlaps_iucn),
  "\nSpecies without IUCN map:", length(iucn_species_missing),
  "\n"
)

