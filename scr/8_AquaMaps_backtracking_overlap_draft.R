#-------------------------------------------------------------------------------
#
# 8. Draft: check backtracking polygons against AquaMaps probability > 0.4
#
#-------------------------------------------------------------------------------

library(sf)
library(sp)
library(raster)
library(dplyr)
library(purrr)
library(stringr)
library(aquamapsdata)


# 1. Set paths and parameters---------------------------------------------------

prob_threshold <- 0.4

backtracking_dir <- file.path(output_data, "particle_backtracking")
polygon_file <- file.path(
  backtracking_dir,
  "species_backtracking_polygons.gpkg"
)
polygon_index_file <- file.path(
  backtracking_dir,
  "species_sample_polygon_index.csv"
)

aquamaps_output_dir <- file.path(
  backtracking_dir,
  "AquaMaps_overlap"
)
dir.create(aquamaps_output_dir, recursive = TRUE, showWarnings = FALSE)


# 2. Load backtracking polygons and sample information--------------------------

if (!file.exists(polygon_file)) {
  stop("Run scr/6_particle_backtracking_polygons.R first.")
}

backtracking_polygons <- st_read(
  polygon_file,
  layer = "backtracking_polygons",
  quiet = TRUE
) %>%
  filter(!is.na(species), !is.na(sample_id)) %>%
  st_make_valid()

polygon_index <- read.csv(
  polygon_index_file,
  sep = ";",
  check.names = FALSE
)

sample_information <- polygon_index %>%
  filter(polygon_found) %>%
  select(species, sample_id, latitude, longitude) %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  ) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  distinct()


# 3. Find the AquaMaps raster for one species-----------------------------------

default_db("sqlite")

# Avoid searching and loading the same raster once per sample
aquamaps_raster_cache <- new.env(parent = emptyenv())

find_aquamaps_raster <- function(sp) {
  
  cache_key <- str_replace_all(sp, "[^A-Za-z0-9]+", "_")
  
  if (exists(cache_key, envir = aquamaps_raster_cache)) {
    return(get(cache_key, envir = aquamaps_raster_cache))
  }
  
  genus <- word(sp, 1)
  species_epithet <- word(sp, 2)
  
  hit <- tryCatch(
    am_search_exact(
      Genus = genus,
      Species = species_epithet
    ),
    error = function(e) NULL
  )
  
  if (is.null(hit) || nrow(hit) == 0) {
    hit <- tryCatch(
      am_search_fuzzy(sp),
      error = function(e) NULL
    )
  }
  
  if (is.null(hit) || nrow(hit) == 0) {
    assign(cache_key, NULL, envir = aquamaps_raster_cache)
    return(NULL)
  }
  
  id_candidates <- c(
    "SpeciesID", "speciesid", "SpecCode", "speccode",
    "species_id", "id", "ID", "key"
  )
  id_col <- id_candidates[id_candidates %in% names(hit)][1]
  
  if (is.na(id_col)) {
    assign(cache_key, NULL, envir = aquamaps_raster_cache)
    return(NULL)
  }
  
  candidate_ids <- unique(hit[[id_col]])
  candidate_ids <- candidate_ids[!is.na(candidate_ids)]
  
  for (id in candidate_ids) {
    
    ras <- tryCatch(
      am_raster(id),
      error = function(e) NULL
    )
    
    if (!is.null(ras)) {
      
      result <- list(
        id = as.character(id),
        raster = ras
      )
      
      assign(cache_key, result, envir = aquamaps_raster_cache)
      return(result)
    }
  }
  
  assign(cache_key, NULL, envir = aquamaps_raster_cache)
  NULL
}


# 4. Check one species and sample-----------------------------------------------

check_aquamaps_backtracking <- function(sp, sample) {
  
  message("Checking AquaMaps: ", sp, " - ", sample)
  
  aquamaps_result <- find_aquamaps_raster(sp)
  
  if (is.null(aquamaps_result)) {
    return(tibble(
      species = sp,
      sample_id = sample,
      aquamaps_found = FALSE,
      aquamaps_id = NA_character_,
      original_point_probability = NA_real_,
      original_point_above_p04 = NA,
      n_polygon_cells = NA_integer_,
      max_polygon_probability = NA_real_,
      backtracking_overlaps_p04 = NA
    ))
  }
  
  ras <- aquamaps_result$raster
  
  sample_sp <- sample_information %>%
    filter(species == sp, sample_id == sample)
  
  point_probability <- raster::extract(
    ras,
    as.matrix(sample_sp[, c("longitude", "latitude")])
  )[1]
  
  particle_sp <- backtracking_polygons %>%
    filter(species == sp, sample_id == sample) %>%
    st_transform(crs(ras))
  
  # raster::extract returns the raster-cell values touched by each polygon
  polygon_values <- raster::extract(
    ras,
    as(particle_sp, "Spatial")
  ) %>%
    unlist(use.names = FALSE)
  
  polygon_values <- polygon_values[!is.na(polygon_values)]
  
  tibble(
    species = sp,
    sample_id = sample,
    aquamaps_found = TRUE,
    aquamaps_id = aquamaps_result$id,
    original_point_probability = point_probability,
    original_point_above_p04 = ifelse(
      is.na(point_probability),
      NA,
      point_probability > prob_threshold
    ),
    n_polygon_cells = length(polygon_values),
    max_polygon_probability = ifelse(
      length(polygon_values) == 0,
      NA_real_,
      max(polygon_values)
    ),
    backtracking_overlaps_p04 = ifelse(
      length(polygon_values) == 0,
      FALSE,
      any(polygon_values > prob_threshold)
    )
  )
}


# 5. Run all species/sample combinations----------------------------------------

species_samples <- sample_information %>%
  distinct(species, sample_id) %>%
  arrange(species, sample_id)

aquamaps_backtracking_all <- map2_dfr(
  species_samples$species,
  species_samples$sample_id,
  check_aquamaps_backtracking
)

# The question of interest concerns points originally at probability <= 0.4
aquamaps_backtracking_overlap <- aquamaps_backtracking_all %>%
  filter(
    aquamaps_found,
    is.na(original_point_above_p04) |
      original_point_above_p04 == FALSE
  )


# 6. Summarise and export draft results-----------------------------------------

aquamaps_backtracking_summary <- aquamaps_backtracking_overlap %>%
  group_by(species, aquamaps_id) %>%
  summarise(
    n_original_points_outside_p04 = n(),
    n_backtracking_overlap_p04 = sum(
      backtracking_overlaps_p04,
      na.rm = TRUE
    ),
    any_backtracking_overlap_p04 = any(
      backtracking_overlaps_p04,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv2(
  aquamaps_backtracking_all,
  file.path(
    aquamaps_output_dir,
    "AquaMaps_backtracking_all_samples_draft.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv2(
  aquamaps_backtracking_overlap,
  file.path(
    aquamaps_output_dir,
    "AquaMaps_backtracking_p04_by_sample_draft.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv2(
  aquamaps_backtracking_summary,
  file.path(
    aquamaps_output_dir,
    "AquaMaps_backtracking_p04_summary_draft.csv"
  ),
  row.names = FALSE,
  na = ""
)

cat(
  "\nSamples originally outside AquaMaps p > ", prob_threshold, ": ",
  nrow(aquamaps_backtracking_overlap),
  "\nBacktracking polygons reaching p > ", prob_threshold, ": ",
  sum(
    aquamaps_backtracking_overlap$backtracking_overlaps_p04,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
