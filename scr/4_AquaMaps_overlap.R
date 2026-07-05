#-------------------------------------------------------------------------------
#
# 4. Check AquaMaps overlap for species not accepted by IUCN decision
#
#-------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(sf)
library(raster)
library(ggplot2)
# Install AquaMaps package from GitHub
#remotes::install_github("raquamaps/aquamapsdata", dependencies = TRUE)
library(aquamapsdata)

# 1. Load metadata and original data------------------------------------------------

file_path <- paste0(input_data, "/processed_df")

metadata_file <- paste0(file_path, "/metadata_iucn_decision.csv")
metadata <- read.csv(metadata_file, sep = ";")

data_file <- paste0(file_path, "/merged_eDNA_spp_pa.csv")
data <- read.csv(data_file, sep = ";") %>%
  filter(!sequencing_comment %in% c("Filtration_Control", "No_target_taxa"))

head(metadata)
head(data)


# 2. Select species needing AquaMaps validation------------------------------------------------

#species_aquamaps_check <- metadata %>%
#  filter(
#    TaxaTargeted == "Yes",
#    is.na(finalIUCN_acceptance) | finalIUCN_acceptance == FALSE
#  ) %>%
#  distinct(corrected_name, original_colname) %>%
#  arrange(corrected_name)

#species_aquamaps_check

#species_to_check <- unique(species_aquamaps_check$corrected_name)


# 3. Extract detections for these species from original data#----------------------------

#species_cols <- species_aquamaps_check$original_colname
#species_cols <- species_cols[species_cols %in% names(data)]

#aquamaps_detections <- data %>%
#  mutate(row_id = row_number()) %>%
#  dplyr::select(row_id, sample_id, latitude, longitude, all_of(species_cols)) %>%
#  pivot_longer(
#    cols = all_of(species_cols),
#    names_to = "original_colname",
#    values_to = "presence"
#  ) %>%
#  filter(presence > 0) %>%
#  left_join(
#    species_aquamaps_check,
#    by = "original_colname"
#  ) %>%
#  rename(species = corrected_name)

#aquamaps_detections

# 4. Check function working------------------------------------------------

# Test with the tiny example database included in the package
default_db("extdata")

am_search_fuzzy("trevally")

## First time only, if not done yet
#download_db()

# Use local sqlite AquaMaps database
default_db("sqlite")

# Check one download:
am_search_fuzzy("Auxis thazard")
am_search_fuzzy("Thunnus alalunga")


am_search_exact(
  Genus = "Auxis",
  Species = "thazard"
)

am_search_exact(
  Species = "thazard"
)

am_search_fuzzy("Auxis")

am_search_fuzzy("Auxis thazard")

hit <- am_search_fuzzy("Auxis thazard")
names(hit)

ras <- am_raster(hit$key[1])
plot(ras)

# 5. Function to check AquaMaps probability at detection points------------------------------------------------

# AquaMaps uses modelled probability of occurrence.
# Here we check if detections fall in cells with probability > 0.4 (following Coll et al 2010; doi: 10.1371/journal.pone.0011842).

check_aquamaps_one_species <- function(sp, detections_df, prob_threshold = 0.4) {
  
  sp <- stringr::str_squish(sp)
  message("Checking AquaMaps: ", sp)
  
  empty_result <- function(found = FALSE, id = NA_character_, n_det = NA_integer_) {
    tibble::tibble(
      species = sp,
      aquamaps_found = found,
      aquamaps_id = id,
      n_detections = n_det,
      n_overlap_aquamaps = NA_integer_,
      overlaps_aquamaps_p04 = NA,
      max_aquamaps_prob = NA_real_,
      mean_aquamaps_prob = NA_real_
    )
  }
  
  points_sp <- detections_df %>%
    dplyr::filter(stringr::str_squish(species) == sp) %>%
    dplyr::mutate(
      longitude = as.numeric(longitude),
      latitude  = as.numeric(latitude)
    ) %>%
    dplyr::filter(!is.na(longitude), !is.na(latitude))
  
  if (nrow(points_sp) == 0) {
    return(empty_result(found = NA, n_det = 0))
  }
  
  coords_sp <- points_sp %>%
    dplyr::select(longitude, latitude) %>%
    as.matrix()
  
  sp_lower <- stringr::str_to_lower(sp)
  genus_sp <- stringr::word(sp, 1)
  species_epithet <- stringr::word(sp, 2)
  
  hit_fuzzy <- tryCatch(
    am_search_fuzzy(sp),
    error = function(e) NULL
  )
  
  if (is.null(hit_fuzzy) || nrow(hit_fuzzy) == 0) {
    return(empty_result(found = FALSE, n_det = nrow(points_sp)))
  }
  
  # Collapse all character columns into one searchable text field
  char_cols <- names(hit_fuzzy)[vapply(hit_fuzzy, is.character, logical(1))]
  
  if (length(char_cols) > 0) {
    hit_fuzzy <- hit_fuzzy %>%
      dplyr::mutate(
        .search_text = apply(
          dplyr::select(., dplyr::all_of(char_cols)),
          1,
          function(x) paste(x, collapse = " ")
        ),
        .search_text = stringr::str_squish(.search_text),
        .search_text_lower = stringr::str_to_lower(.search_text)
      )
  } else {
    hit_fuzzy$.search_text <- NA_character_
    hit_fuzzy$.search_text_lower <- NA_character_
  }
  
  # Exact binomial match somewhere in the AquaMaps search result
  hit_match <- hit_fuzzy %>%
    dplyr::filter(stringr::str_detect(.search_text_lower, stringr::fixed(sp_lower)))
  
  # If there are Genus/Species columns in any capitalization, use them too
  genus_col <- names(hit_fuzzy)[stringr::str_to_lower(names(hit_fuzzy)) == "genus"]
  species_col <- names(hit_fuzzy)[stringr::str_to_lower(names(hit_fuzzy)) == "species"]
  
  if (nrow(hit_match) == 0 && length(genus_col) > 0 && length(species_col) > 0) {
    hit_match <- hit_fuzzy %>%
      dplyr::filter(
        stringr::str_to_lower(.data[[genus_col[1]]]) == stringr::str_to_lower(genus_sp),
        stringr::str_to_lower(.data[[species_col[1]]]) == stringr::str_to_lower(species_epithet)
      )
  }
  
  if (nrow(hit_match) == 0) {
    message(
      "No validated AquaMaps match for ", sp,
      ". Available columns were: ",
      paste(names(hit_fuzzy), collapse = ", ")
    )
    return(empty_result(found = FALSE, n_det = nrow(points_sp)))
  }
  
  id_candidates <- c(
    "SpeciesID", "speciesid", "SpecCode", "speccode",
    "species_id", "id", "ID", "key"
  )
  
  id_col <- id_candidates[id_candidates %in% names(hit_match)][1]
  
  if (is.na(id_col)) {
    stop(
      "No usable AquaMaps ID column found for ", sp,
      ". Available columns: ",
      paste(names(hit_match), collapse = ", ")
    )
  }
  
  candidate_ids <- unique(hit_match[[id_col]])
  candidate_ids <- candidate_ids[!is.na(candidate_ids)]
  
  if (length(candidate_ids) == 0) {
    return(empty_result(found = FALSE, n_det = nrow(points_sp)))
  }
  
  candidate_results <- purrr::map_dfr(candidate_ids, function(id) {
    
    ras_tmp <- tryCatch(
      am_raster(id),
      error = function(e) NULL
    )
    
    if (is.null(ras_tmp)) {
      return(tibble::tibble(
        aquamaps_id = as.character(id),
        n_overlap = NA_integer_,
        max_prob = NA_real_,
        mean_prob = NA_real_
      ))
    }
    
    probs_tmp <- raster::extract(ras_tmp, coords_sp)
    
    tibble::tibble(
      aquamaps_id = as.character(id),
      n_overlap = sum(probs_tmp > prob_threshold, na.rm = TRUE),
      max_prob = ifelse(all(is.na(probs_tmp)), NA_real_, max(probs_tmp, na.rm = TRUE)),
      mean_prob = ifelse(all(is.na(probs_tmp)), NA_real_, mean(probs_tmp, na.rm = TRUE))
    )
  })
  
  best_result <- candidate_results %>%
    dplyr::arrange(
      dplyr::desc(n_overlap),
      dplyr::desc(max_prob),
      dplyr::desc(mean_prob)
    ) %>%
    dplyr::slice(1)
  
  aquamaps_id <- best_result$aquamaps_id
  
  message(
    "Selected AquaMaps ID: ",
    aquamaps_id,
    " | overlap: ",
    best_result$n_overlap,
    "/",
    nrow(points_sp)
  )
  
  ras <- tryCatch(
    am_raster(aquamaps_id),
    error = function(e) NULL
  )
  
  if (is.null(ras)) {
    return(empty_result(found = TRUE, id = aquamaps_id, n_det = nrow(points_sp)))
  }
  
  probs <- raster::extract(ras, coords_sp)
  
  tibble::tibble(
    species = sp,
    aquamaps_found = TRUE,
    aquamaps_id = as.character(aquamaps_id),
    n_detections = nrow(points_sp),
    n_overlap_aquamaps = sum(probs > prob_threshold, na.rm = TRUE),
    overlaps_aquamaps_p04 = any(probs > prob_threshold, na.rm = TRUE),
    max_aquamaps_prob = ifelse(all(is.na(probs)), NA_real_, max(probs, na.rm = TRUE)),
    mean_aquamaps_prob = ifelse(all(is.na(probs)), NA_real_, mean(probs, na.rm = TRUE))
  )
}


# 6. Run AquaMaps overlap check for species which occurrence in the study area is unclear------------------------------------------------

#default_db("sqlite")

#aquamaps_overlap <- map_dfr(
#  species_to_check,
#  check_aquamaps_one_species,
#  detections_df = aquamaps_detections,
#  prob_threshold = 0.4)

#aquamaps_overlap

# 7. Check AquaMaps overlap for all species------------------------------------------------


# Species to evaluate
species_to_check <- metadata %>%
  dplyr::filter(TaxaTargeted == "Yes") %>%
  dplyr::pull(corrected_name) %>%
  unique() %>%
  sort()

# Metadata map: original dataset column -> corrected species name
aquamaps_species_map <- metadata %>%
  dplyr::filter(
    TaxaTargeted == "Yes",
    corrected_name %in% species_to_check
  ) %>%
  dplyr::select(original_colname, corrected_name) %>%
  dplyr::distinct()

# Keep only columns that actually exist in data
aquamaps_cols <- aquamaps_species_map$original_colname
aquamaps_cols <- aquamaps_cols[aquamaps_cols %in% names(data)]

length(aquamaps_cols)
aquamaps_cols


# Prepare detection dataset
aquamaps_detections <- data %>%
  mutate(row_id = row_number()) %>%
  dplyr::select(
    row_id,
    sample_id,
    latitude,
    longitude,
    all_of(aquamaps_cols)
  ) %>%
  pivot_longer(
    cols = all_of(aquamaps_cols),
    names_to = "original_colname",
    values_to = "presence"
  ) %>%
  dplyr::filter(presence > 0) %>%
  left_join(
    aquamaps_species_map,
    by = "original_colname"
  ) %>%
  rename(species = corrected_name) %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  )

aquamaps_detections


# Run AquaMaps overlap

default_db("sqlite")

aquamaps_overlap <- purrr::map_dfr(
  species_to_check,
  check_aquamaps_one_species,
  detections_df = aquamaps_detections,
  prob_threshold = 0.4
)


aquamaps_overlap

# Check rare ones:
aquamaps_overlap %>%
  dplyr::filter(species == "Caretta caretta")

aquamaps_overlap %>%
  dplyr::filter(species == "Dentex dentex")
hit <- am_search_exact(Genus = "Dentex", Species = "dentex")
hit
ras <- am_raster(hit$SpeciesID[1])
plot(ras)

sp <- "Dentex dentex"
genus_sp <- stringr::word(sp, 1)
species_epithet <- stringr::word(sp, 2)
hit_exact <- am_search_exact(
  Genus = genus_sp,
  Species = species_epithet)
candidate_ids <- unique(hit_exact$SpeciesID)
candidate_ids
#"Fis-23108"
check_aquamaps_one_species(
  sp = "Dentex dentex",
  detections_df = aquamaps_detections,
  prob_threshold = 0.4)
#Fis-23108


# 8. Add relevant AquaMaps information to metadata------------------------------------------------

metadata <- metadata %>%
  left_join(
    aquamaps_overlap %>%
      dplyr::select(
        species,
        n_overlap_aquamaps,
        overlaps_aquamaps_p04,
        max_aquamaps_prob,
        mean_aquamaps_prob
      ) %>%
      rename(corrected_name = species),
    by = "corrected_name"
  )

# Check
metadata %>%
  dplyr::select(
    corrected_name,
    n_overlap_aquamaps,
    overlaps_aquamaps_p04,
    max_aquamaps_prob,
    mean_aquamaps_prob
  ) %>%
  arrange(corrected_name)


# 9. Export updated metadata------------------------------------------------

file <- paste0(file_path, "/metadata_iucn_aquamaps_decision.csv")

write.csv2(
  metadata,
  file,
  row.names = FALSE,
  na = ""
)

# 10. Plots to check those not overlapping------------------------------------------------

default_db("sqlite")

plot_dir <- file.path(input_data, "figures", "AquaMaps_false_species")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

world_plot <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_make_valid()

plot_extent <- raster::extent(-10, 45, 25, 50)

aquamaps_false_species <- aquamaps_overlap %>%
  dplyr::filter(overlaps_aquamaps_p04 == FALSE) %>%
  dplyr::filter(aquamaps_found == TRUE) %>%
  dplyr::filter(!is.na(aquamaps_id)) %>%
  dplyr::pull(species) %>%
  unique() %>%
  sort()

for (sp in aquamaps_false_species) {
  
  message("Plotting AquaMaps: ", sp)
  
  am_id <- aquamaps_overlap %>%
    dplyr::filter(species == sp) %>%
    dplyr::pull(aquamaps_id) %>%
    unique() %>%
    .[1]
  
  ras <- tryCatch(
    am_raster(am_id),
    error = function(e) NULL
  )
  
  if (is.null(ras)) {
    message("Skipping ", sp, ": raster not available")
    next
  }
  
  extent_overlap <- raster::intersect(raster::extent(ras), plot_extent)
  
  if (is.null(extent_overlap)) {
    message("Skipping ", sp, ": raster does not overlap plot extent")
    next
  }
  
  ras_crop <- raster::crop(ras, plot_extent)
  
  ras_df <- as.data.frame(
    ras_crop,
    xy = TRUE,
    na.rm = FALSE
  )
  
  names(ras_df) <- c("longitude", "latitude", "aquamaps_prob")
  
  points_sp <- aquamaps_detections %>%
    dplyr::filter(species == sp) %>%
    dplyr::mutate(
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude)
    )
  
  if (nrow(points_sp) == 0) {
    message("Skipping ", sp, ": no detection points")
    next
  }
  
  coords_sp <- points_sp %>%
    dplyr::select(longitude, latitude) %>%
    as.matrix()
  
  points_sp$aquamaps_prob <- raster::extract(ras, coords_sp)
  
  points_sp <- points_sp %>%
    dplyr::mutate(
      aquamaps_overlap_class = dplyr::case_when(
        is.na(aquamaps_prob) ~ "No AquaMaps cell",
        aquamaps_prob > 0.4 ~ "Probability > 0.4",
        aquamaps_prob <= 0.4 ~ "Probability <= 0.4"
      )
    )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = ras_df,
      ggplot2::aes(
        x = longitude,
        y = latitude,
        fill = aquamaps_prob
      )
    ) +
    ggplot2::geom_sf(
      data = world_plot,
      fill = "grey90",
      colour = "grey40",
      linewidth = 0.2
    ) +
    ggplot2::geom_point(
      data = points_sp,
      ggplot2::aes(
        x = longitude,
        y = latitude,
        colour = aquamaps_overlap_class
      ),
      size = 1.5,
      alpha = 0.9
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Probability > 0.4" = "darkgreen",
        "Probability <= 0.4" = "red3",
        "No AquaMaps cell" = "grey30"
      ),
      name = "Point classification"
    ) +
    ggplot2::scale_fill_viridis_c(
      name = "AquaMaps\nprobability",
      limits = c(0, 1),
      na.value = NA
    ) +
    ggplot2::coord_sf(
      xlim = c(-10, 45),
      ylim = c(25, 50),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = paste0(sp, " - AquaMaps ID: ", am_id),
      subtitle = "Threshold for overlap: probability > 0.4"
    ) +
    ggplot2::theme_bw()
  
  file_name <- stringr::str_replace_all(sp, " ", "_")
  
  ggplot2::ggsave(
    filename = file.path(plot_dir, paste0(file_name, "_AquaMaps_overlap_check.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  rm(ras, ras_crop, ras_df, points_sp, p)
  gc()
}



# 11.Check rare results------------------------------------------------

aquamaps_false_species <- aquamaps_overlap %>%
  dplyr::filter(overlaps_aquamaps_p04 == FALSE) %>%
  dplyr::pull(species) %>%
  unique() %>%
  sort()

aquamaps_false_species

## Species not overlapping with p>0.4 but actually ocurring in the area:
## they are either in shallower or deeper waters:
#Ariosoma_balearicum
#Balaenoptera acutorostrata
#Balaenoptera physalus
#Bathypterois dubius
#Callanthias ruber
#Capros aper
#Centrolophus niger
#Coris julis
#Crystallogobius linearis
#Dentex dentex
#Echelus myrus
#Globicephala melas
#Gobius vittatus
#Lepadogaster candolii
#Lepidotrigla cavillone
#Macroramphosus scolopax
#Maurolicus muelleri
#Nerophis maculatus
#Notoscopelus kroyeri # rare, check if close species
#Pegusa lascaris
#Serranus scriba   
#Sphyraena viridensis          
#Spicara maena   
#Stomias boa      
#Tetragonurus cuvieri         
#Trachipterus arcticus # rare          
#Trachyscorpia cristulata  #rare 


## Does not occur in the Mediterranean:
#Cyclothone atraria
#Dagetichthys lusitanicus 
#Delphinapterus leucas
#Diaphus dumerilii    
#Engraulis australis
#Engraulis japonicus
#Foetorepus agassizii     
#Salmo salar
#Schedophilus velaini
#Scomber japonicus
#Scorpaena neglecta
#Thunnus orientalis
#Thunnus tonggol

## Freshwater
#Leuciscus leuciscus
#Oncorhynchus mykiss


# 12.Update metadata------------------------------------------------------------

head(metadata)

# Update metadata with AquaMaps visual inspection

aquamaps_true_species <- c(
  "Ariosoma balearicum",
  "Balaenoptera acutorostrata",
  "Balaenoptera physalus",
  "Bathypterois dubius",
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
  "Notoscopelus kroyeri",
  "Pegusa lascaris",
  "Serranus scriba",
  "Sphyraena viridensis",
  "Spicara maena",
  "Stomias boa",
  "Tetragonurus cuvieri",
  "Trachipterus arcticus",
  "Trachyscorpia cristulata"
)

aquamaps_false_species <- c(
  "Cyclothone atraria",
  "Dagetichthys lusitanicus",
  "Delphinapterus leucas",
  "Diaphus dumerilii",
  "Engraulis australis",
  "Engraulis japonicus",
  "Foetorepus agassizii",
  "Salmo salar",
  "Schedophilus velaini",
  "Scomber japonicus",
  "Scorpaena neglecta",
  "Thunnus orientalis",
  "Thunnus tonggol",
  "Leuciscus leuciscus",
  "Oncorhynchus mykiss"
)

metadata <- metadata %>%
  dplyr::mutate(
    corrected_name_clean = stringr::str_squish(corrected_name),
    
    visualInspection_aquaMaps = dplyr::case_when(
      corrected_name_clean %in% aquamaps_true_species ~ TRUE,
      corrected_name_clean %in% aquamaps_false_species ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  dplyr::select(-corrected_name_clean)

metadata <- metadata %>%
  dplyr::mutate(
    finalAquamaps_acceptance = dplyr::case_when(
      
      # If either automatic or manual inspection is TRUE
      overlaps_aquamaps_p04 == TRUE |
        visualInspection_aquaMaps == TRUE ~ TRUE,
      
      # Both missing
      is.na(overlaps_aquamaps_p04) &
        is.na(visualInspection_aquaMaps) ~ NA,
      
      # Everything else (FALSE/FALSE or FALSE/NA or NA/FALSE)
      TRUE ~ FALSE
    )
  )


# 13. Export updated metadata----------------------------------------------------

file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)

file <- paste0(file_path, "/metadata_iucn_aquamaps_decision.csv")

write.csv2(
  metadata,
  file,
  row.names = FALSE,
  na = ""
)
