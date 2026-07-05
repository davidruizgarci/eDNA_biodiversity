#-------------------------------------------------------------------------------
#
# 3. Download IUCN maps to check distribution
#
#-------------------------------------------------------------------------------

library(sf)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyr)
library(lwgeom)


# 1. Load data------------------------------------------------------------------
# Species list:
file_path <- paste0(input_data, "/processed_df")
file <- paste0(file_path,  "/speciesList.csv")

marine_vertebrates_clean <- read.csv(file, sep = ";")
head(marine_vertebrates_clean)

iucn_dir <- file.path(input_data, "IUCNmaps")

# Metadata dataset:
file_path <- paste0(input_data, "/processed_df")
file <- paste0(file_path, "/metadataFiltering.csv")

metadata <- read.csv(file, sep = ";")
head(metadata)

# Full dataset:
file_path <- paste0(input_data, "/processed_df")
file <- paste0(file_path, "/renamed_dataset.csv")

data <- read.csv(file, sep = ";")
head(data)



# 2. Check if marine mammal species belong to study area using IUCN polygons-----------------

# Mammals
mammals_shp <- file.path(
  iucn_dir,
  "MAMMALS_MARINE_ONLY",
  "MAMMALS_MARINE_ONLY.shp"
)

# Read mammals IUCN shapefile
iucn_mammals <- st_read(mammals_shp, quiet = TRUE) %>%
  st_transform(4326)
#names(iucn_mammals)

# filter species involved:
# Marine mammals detected in our species list
our_mammals <- marine_vertebrates_clean %>%
  filter(class == "Mammalia") %>%
  distinct(valid_name)

# Check matches:
iucn_mammals_filtered <- iucn_mammals %>%
  filter(sci_name %in% our_mammals$valid_name)
# matching
matched_mammals <- sort(unique(iucn_mammals_filtered$sci_name))
matched_mammals
# missing
mammals_missing_iucn <- setdiff(
  our_mammals$valid_name,
  matched_mammals)

mammals_missing_iucn

# summaries
cat("Marine mammals in our list:", nrow(our_mammals), "\n")
cat("Marine mammals matched in IUCN:", length(matched_mammals), "\n")
cat("Marine mammals missing from IUCN:", length(mammals_missing_iucn), "\n")
mammals_iucn_summary <- our_mammals %>%
  mutate(
    in_iucn_mammals_shp = valid_name %in% matched_mammals)
mammals_iucn_summary

# Keep only the polygons for our detected marine mammals
iucn_mammals_filtered <- iucn_mammals %>%
  filter(sci_name %in% our_mammals$valid_name)

# Remove the full shapefile from memory
rm(
  our_mammals,
  matched_mammals,
  mammals_missing_iucn)

# Run garbage collection
gc()

# World land
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot them all to check
ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey50", linewidth = 0.2) +
  geom_sf(data = iucn_mammals_filtered,
          fill = "steelblue",
          colour = NA,
          alpha = 0.6) +
  facet_wrap(~sci_name) +
  theme_bw()


# 2.1. Check if mammal species belongs to study area-----------------------------------
head(data)
head(metadata)
head(our_mammals)

summary(data$latitude)
summary(data$longitude)

# Create sample points from data
data_points <- data %>%
  mutate(row_id = row_number()) %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# Get mammal columns in data using metadata
mammal_metadata <- metadata %>%
  filter(
    TaxaTargeted == "Yes",
    class == "Mammalia",
    corrected_name %in% our_mammals$valid_name
  )

mammal_cols <- mammal_metadata$corrected_name

# Careful with the naming style
mammal_cols_data <- str_replace_all(mammal_cols, " ", ".")
mammal_cols_data <- mammal_cols_data[mammal_cols_data %in% names(data)]
mammal_cols_data

# Format
mammal_detections <- data %>%
  mutate(row_id = row_number()) %>%
  select(row_id, sample_id, latitude, longitude, all_of(mammal_cols_data)) %>%
  pivot_longer(
    cols = all_of(mammal_cols_data),
    names_to = "species_col",
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%
  mutate(
    species = str_replace_all(species_col, "\\.", " "))

mammal_detections

# Make sure points are correct for check distributions:
mammal_detection_points <- mammal_detections %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# 2.1.1. Check overlap:
# Temporarily disable S2 for validity checks in lon/lat polygons
old_s2 <- sf_use_s2()
sf_use_s2(FALSE)

# Validity report
iucn_mammals_validity <- iucn_mammals_filtered %>%
  mutate(
    valid_geom = st_is_valid(geometry),
    valid_reason = st_is_valid(geometry, reason = TRUE)
  )

# Summary
iucn_mammals_validity_summary <- iucn_mammals_validity %>%
  st_drop_geometry() %>%
  count(valid_geom, sort = TRUE)

iucn_mammals_validity_summary

# Invalid polygons only
iucn_mammals_invalid <- iucn_mammals_validity %>%
  filter(!valid_geom)

# Invalid polygon details
iucn_mammals_invalid_report <- iucn_mammals_invalid %>%
  st_drop_geometry() %>%
  select(
    sci_name,
    id_no,
    presence,
    origin,
    seasonal,
    legend,
    valid_geom,
    valid_reason
  )

iucn_mammals_invalid_report

# Invalid polygons by species/reason
iucn_mammals_invalid_report %>%
  count(sci_name, valid_reason, sort = TRUE)

# Restore previous S2 setting
sf_use_s2(old_s2)


# Asssess overlap
# Repair / standardise polygons
old_s2 <- sf_use_s2()
sf_use_s2(FALSE)

#iucn_mammals_filtered_valid <- iucn_mammals_filtered %>%
#  st_make_valid() %>%
#  st_collection_extract("POLYGON", warn = FALSE) %>%
#  st_transform(4326)

check_iucn_overlap_one_species <- function(sp) {
  
  points_sp <- mammal_detection_points %>%
    filter(species == sp)
  
  poly_sp <- iucn_mammals_filtered %>%
    filter(sci_name == sp)
  
  if (nrow(points_sp) == 0) {
    return(tibble(
      species = sp,
      n_detections = 0,
      n_overlap_iucn = 0,
      overlaps_iucn = NA
    ))
  }
  
  if (nrow(poly_sp) == 0) {
    return(tibble(
      species = sp,
      n_detections = nrow(points_sp),
      n_overlap_iucn = 0,
      overlaps_iucn = FALSE
    ))
  }
  
  overlap <- lengths(st_intersects(points_sp, poly_sp)) > 0
  
  tibble(
    species = sp,
    n_detections = nrow(points_sp),
    n_overlap_iucn = sum(overlap),
    overlaps_iucn = any(overlap)
  )
}

mammal_iucn_overlap <- purrr::map_dfr(
  unique(mammal_detection_points$species),
  check_iucn_overlap_one_species
)

sf_use_s2(old_s2)

mammal_iucn_overlap

# Check point by point:
mammal_detection_iucn_detail <- map_dfr(
  unique(mammal_detections$species),
  function(sp) {
    
    points_sp <- mammal_detection_points %>% filter(species == sp)
    poly_sp <- iucn_mammals_filtered_valid %>% filter(sci_name == sp)
    
    if (nrow(points_sp) == 0 || nrow(poly_sp) == 0) return(NULL)
    
    overlap_logical <- lengths(st_intersects(points_sp, poly_sp)) > 0
    
    points_sp %>%
      st_drop_geometry() %>%
      mutate(overlaps_iucn = overlap_logical)
  }
)

# Add IUCN overlap information to metadata
head(metadata)

metadata <- metadata %>%
  left_join(
    mammal_iucn_overlap %>%
      select(species, overlaps_iucn) %>%
      rename(
        corrected_name = species,
        overlaps_iucn = overlaps_iucn
      ),
    by = "corrected_name"
  )

# Export
file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/metadata_iucnMammals.csv")

write.csv2(metadata, file, row.names = FALSE,  na = "")


# Clean memory from mammal objects no longer needed

rm(
  mammal_metadata,
  mammal_cols,
  mammal_cols_data,
  mammal_detections,
  mammal_detection_points,
  mammal_iucn_overlap,
  mammal_detection_iucn_detail,
  iucn_mammals_validity,
  iucn_mammals_validity_summary,
  iucn_mammals_invalid,
  iucn_mammals_invalid_report
)

gc()





# 3. Check if turtle species belong to study area using IUCN polygons---------------------

# 3.1. Load IUCN turtle polygons------------------------------------------------

iucn_turtles <- st_read(
  file.path(
    iucn_dir,
    "TURTLES",
    "TURTLES.shp"
  ),
  quiet = TRUE
) %>%
  st_transform(4326)

names(iucn_turtles)
head(iucn_turtles)


# 3.2 Select turtle species from our metadata------------------------------------------------

turtle_metadata <- metadata %>%
  filter(
    TaxaTargeted == "Yes",
    corrected_name %in% c(
      "Caretta caretta",
      "Chelonia mydas",
      "Dermochelys coriacea"
    )
  )

turtle_cols <- turtle_metadata$corrected_name

# Dataset column names use dots instead of spaces
turtle_cols_data <- str_replace_all(turtle_cols, " ", ".")

# Keep only columns actually present in data
turtle_cols_data <- turtle_cols_data[turtle_cols_data %in% names(data)]

turtle_cols_data


# 3.3. Extract turtle detections from data------------------------------------------------

turtle_detections <- data %>%
  mutate(row_id = row_number()) %>%
  select(row_id, sample_id, latitude, longitude, all_of(turtle_cols_data)) %>%
  pivot_longer(
    cols = all_of(turtle_cols_data),
    names_to = "species_col",
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%
  mutate(
    species = str_replace_all(species_col, "\\.", " ")
  )

turtle_detections


# Convert turtle detections to spatial points

turtle_detection_points <- turtle_detections %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )


# 3.4. Filter IUCN turtle polygons to our detected turtle species------------------------------------------------

iucn_turtles_filtered <- iucn_turtles %>%
  filter(sci_name %in% turtle_metadata$corrected_name)

unique(iucn_turtles_filtered$sci_name)


# 3.5. Check geometry validity------------------------------------------------

old_s2 <- sf_use_s2()
sf_use_s2(FALSE)

iucn_turtles_validity <- iucn_turtles_filtered %>%
  mutate(
    valid_geom = st_is_valid(geometry),
    valid_reason = st_is_valid(geometry, reason = TRUE)
  )

iucn_turtles_validity_summary <- iucn_turtles_validity %>%
  st_drop_geometry() %>%
  count(valid_geom, sort = TRUE)

iucn_turtles_validity_summary

iucn_turtles_invalid_report <- iucn_turtles_validity %>%
  filter(!valid_geom) %>%
  st_drop_geometry() %>%
  select(
    sci_name,
    id_no,
    presence,
    origin,
    seasonal,
    legend,
    valid_geom,
    valid_reason
  )

iucn_turtles_invalid_report

sf_use_s2(old_s2)


# 3.6. Check whether turtle detections overlap their IUCN polygons------------------------------------------------

old_s2 <- sf_use_s2()
sf_use_s2(FALSE)

check_iucn_turtle_overlap_one_species <- function(sp) {
  
  points_sp <- turtle_detection_points %>%
    filter(species == sp)
  
  poly_sp <- iucn_turtles_filtered %>%
    filter(sci_name == sp)
  
  if (nrow(points_sp) == 0) {
    return(tibble(
      species = sp,
      n_detections = 0,
      n_overlap_iucn = 0,
      overlaps_iucn = NA
    ))
  }
  
  if (nrow(poly_sp) == 0) {
    return(tibble(
      species = sp,
      n_detections = nrow(points_sp),
      n_overlap_iucn = 0,
      overlaps_iucn = FALSE
    ))
  }
  
  overlap <- lengths(st_intersects(points_sp, poly_sp)) > 0
  
  tibble(
    species = sp,
    n_detections = nrow(points_sp),
    n_overlap_iucn = sum(overlap),
    overlaps_iucn = any(overlap)
  )
}

turtle_iucn_overlap <- map_dfr(
  unique(turtle_detection_points$species),
  check_iucn_turtle_overlap_one_species
)

sf_use_s2(old_s2)

turtle_iucn_overlap


# point-by-point overlap detail

old_s2 <- sf_use_s2()
sf_use_s2(FALSE)

turtle_detection_iucn_detail <- map_dfr(
  unique(turtle_detection_points$species),
  function(sp) {
    
    points_sp <- turtle_detection_points %>%
      filter(species == sp)
    
    poly_sp <- iucn_turtles_filtered %>%
      filter(sci_name == sp)
    
    if (nrow(points_sp) == 0 || nrow(poly_sp) == 0) return(NULL)
    
    overlap_logical <- lengths(st_intersects(points_sp, poly_sp)) > 0
    
    points_sp %>%
      st_drop_geometry() %>%
      mutate(overlaps_iucn = overlap_logical)
  }
)

sf_use_s2(old_s2)

turtle_detection_iucn_detail


# 3.7. Add turtle IUCN overlap information to metadata------------------------------------------------

metadata <- metadata %>%
  left_join(
    turtle_iucn_overlap %>%
      select(species, overlaps_iucn) %>%
      rename(
        corrected_name = species,
        overlaps_iucn_turtles = overlaps_iucn
      ),
    by = "corrected_name"
  ) %>%
  mutate(
    overlaps_iucn = case_when(
      !is.na(overlaps_iucn) ~ overlaps_iucn,
      !is.na(overlaps_iucn_turtles) ~ overlaps_iucn_turtles,
      TRUE ~ overlaps_iucn
    )
  ) %>%
  select(-overlaps_iucn_turtles)


# Export updated metadata

file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)

file <- paste0(file_path, "/metadata_iucnMammals_Turtles.csv")

write.csv2(
  metadata,
  file,
  row.names = FALSE,
  na = ""
)


# Clean memory from turtle objects no longer needed

rm(
  turtle_metadata,
  turtle_cols,
  turtle_cols_data,
  turtle_detections,
  turtle_detection_points,
  turtle_iucn_overlap,
  turtle_detection_iucn_detail,
  iucn_turtles_validity,
  iucn_turtles_validity_summary,
  iucn_turtles_invalid_report,
  iucn_turtles
)

gc()



# 4. Check if fish species belong to study area using IUCN polygons------------------------------------------------

# 4.1. Select fish species from metadata------------------------------------------------

fish_metadata <- metadata %>%
  filter(
    TaxaTargeted == "Yes",
    class %in% c("Teleostei", "Elasmobranchii")
  )

fish_species <- unique(fish_metadata$corrected_name)

fish_species


# 4.2. Identify fish columns in data------------------------------------------------

fish_cols_data <- str_replace_all(fish_species, " ", ".")
fish_cols_data <- fish_cols_data[fish_cols_data %in% names(data)]

fish_cols_data


# 4.3. Extract fish detections from data------------------------------------------------

fish_detections <- data %>%
  mutate(row_id = row_number()) %>%
  select(row_id, sample_id, latitude, longitude, all_of(fish_cols_data)) %>%
  pivot_longer(
    cols = all_of(fish_cols_data),
    names_to = "species_col",
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%
  mutate(
    species = str_replace_all(species_col, "\\.", " ")
  )

fish_detections


# 4.4. Convert fish detections to spatial points------------------------------------------------

fish_detection_points <- fish_detections %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )


# 4.5. Read IUCN fish shapefiles one by one and keep only our species------------------------------------------------

fish_shps <- list.files(
  file.path(iucn_dir, "MARINEFISH"),
  pattern = "\\.shp$",
  full.names = TRUE,
  ignore.case = TRUE
)

iucn_fish_filtered_list <- list()

for (i in seq_along(fish_shps)) {
  
  message("Reading fish shapefile ", i, " of ", length(fish_shps), ": ", basename(fish_shps[i]))
  
  fish_part <- st_read(fish_shps[i], quiet = TRUE) %>%
    st_transform(4326)
  
  fish_part_filtered <- fish_part %>%
    filter(sci_name %in% fish_species) %>%
    mutate(iucn_file = basename(fish_shps[i]))
  
  message("Matched polygons: ", nrow(fish_part_filtered))
  
  iucn_fish_filtered_list[[i]] <- fish_part_filtered
  
  rm(fish_part, fish_part_filtered)
  gc()
}

iucn_fish_filtered <- bind_rows(iucn_fish_filtered_list)

unique(iucn_fish_filtered$sci_name)


# 4.6. Check which fish species have IUCN polygons------------------------------------------------

fish_with_iucn <- sort(unique(iucn_fish_filtered$sci_name))

fish_missing_iucn <- setdiff(fish_species, fish_with_iucn)

fish_with_iucn
fish_missing_iucn


# 4.7. Check geometry validity------------------------------------------------

old_s2 <- sf_use_s2()
sf_use_s2(FALSE)

iucn_fish_validity <- iucn_fish_filtered %>%
  mutate(
    valid_geom = st_is_valid(geometry),
    valid_reason = st_is_valid(geometry, reason = TRUE)
  )

iucn_fish_validity_summary <- iucn_fish_validity %>%
  st_drop_geometry() %>%
  count(valid_geom, sort = TRUE)

iucn_fish_validity_summary

iucn_fish_invalid_report <- iucn_fish_validity %>%
  filter(!valid_geom) %>%
  st_drop_geometry() %>%
  select(
    sci_name,
    id_no,
    presence,
    origin,
    seasonal,
    legend,
    valid_geom,
    valid_reason
  )

iucn_fish_invalid_report

sf_use_s2(old_s2)


# 4.8. Check whether fish detections overlap their IUCN polygons------------------------------------------------

old_s2 <- sf_use_s2()
sf_use_s2(FALSE)

check_iucn_fish_overlap_one_species <- function(sp) {
  
  points_sp <- fish_detection_points %>%
    filter(species == sp)
  
  poly_sp <- iucn_fish_filtered %>%
    filter(sci_name == sp)
  
  if (nrow(points_sp) == 0) {
    return(tibble(
      species = sp,
      n_detections = 0,
      n_overlap_iucn = 0,
      overlaps_iucn = NA
    ))
  }
  
  if (nrow(poly_sp) == 0) {
    return(tibble(
      species = sp,
      n_detections = nrow(points_sp),
      n_overlap_iucn = 0,
      overlaps_iucn = FALSE
    ))
  }
  
  overlap <- lengths(st_intersects(points_sp, poly_sp)) > 0
  
  tibble(
    species = sp,
    n_detections = nrow(points_sp),
    n_overlap_iucn = sum(overlap),
    overlaps_iucn = any(overlap)
  )
}

fish_iucn_overlap <- map_dfr(
  unique(fish_detection_points$species),
  check_iucn_fish_overlap_one_species
)

sf_use_s2(old_s2)

fish_iucn_overlap

# 4.9. Point-by-point overlap detail------------------------------------------------

old_s2 <- sf_use_s2()
sf_use_s2(FALSE)

fish_detection_iucn_detail <- map_dfr(
  unique(fish_detection_points$species),
  function(sp) {
    
    points_sp <- fish_detection_points %>%
      filter(species == sp)
    
    poly_sp <- iucn_fish_filtered %>%
      filter(sci_name == sp)
    
    if (nrow(points_sp) == 0 || nrow(poly_sp) == 0) return(NULL)
    
    overlap_logical <- lengths(st_intersects(points_sp, poly_sp)) > 0
    
    points_sp %>%
      st_drop_geometry() %>%
      mutate(overlaps_iucn = overlap_logical)
  }
)

sf_use_s2(old_s2)

fish_detection_iucn_detail



# 4.10. Plot FALSE ones, to double check------------------------------------------------

# Species with at least one detection outside IUCN range
fish_false_species <- fish_detection_iucn_detail %>%
  filter(overlaps_iucn == FALSE) %>%
  distinct(species) %>%
  pull(species)

fish_false_species

# IUCN polygons for those species
iucn_fish_false <- iucn_fish_filtered %>%
  filter(sci_name %in% fish_false_species)

# Points for those species
fish_points_false_species <- fish_detection_iucn_detail %>%
  filter(species %in% fish_false_species) %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# Remove empty geometries and make valid
iucn_fish_false_plot <- iucn_fish_false %>%
  filter(!st_is_empty(geometry)) %>%
  st_make_valid()

#fish_points_false_species_plot <- fish_points_false_species %>%
#  filter(!st_is_empty(geometry))

# World land
world <- ne_countries(scale = "medium", returnclass = "sf")

#world_plot <- world %>%
#  filter(!st_is_empty(geometry)) %>%
#  st_make_valid()

# Plot

# Output folder
plot_dir <- file.path(input_data, "figures", "Overlap_IUCN")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Species with at least one point outside IUCN range
fish_false_species <- fish_detection_iucn_detail %>%
  filter(overlaps_iucn == FALSE) %>%
  distinct(species) %>%
  pull(species)

# Loop species
for (sp in fish_false_species) {
  
  message("Plotting: ", sp)
  
  # Species-specific IUCN polygon
  poly_sp <- iucn_fish_filtered %>%
    filter(sci_name == sp) %>%
    filter(!st_is_empty(geometry)) %>%
    st_make_valid()
  
  # Species-specific points
  points_sp <- fish_detection_iucn_detail %>%
    filter(species == sp) %>%
    st_as_sf(
      coords = c("longitude", "latitude"),
      crs = 4326,
      remove = FALSE
    )
  
  # Skip if no polygon or no points
  if (nrow(poly_sp) == 0 || nrow(points_sp) == 0) {
    message("Skipping ", sp, ": no polygon or no points")
    next
  }
  
  p <- ggplot() +
    geom_sf(data = world_plot, fill = "grey90", colour = "grey50", linewidth = 0.2) +
    geom_sf(
      data = poly_sp,
      fill = "steelblue",
      colour = NA,
      alpha = 0.5
    ) +
    geom_sf(
      data = points_sp,
      aes(shape = overlaps_iucn),
      size = 2
    ) +
    coord_sf(
      xlim = c(-20, 45),
      ylim = c(25, 50),
      expand = FALSE
    ) +
    labs(
      title = paste0(sp, " - IUCN range and eDNA detections"),
      shape = "Point overlaps IUCN"
    ) +
    theme_bw()
  
  # Safe filename
  file_name <- str_replace_all(sp, " ", "_")
  
  ggsave(
    filename = file.path(plot_dir, paste0(file_name, "_IUCN_overlap_check.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# 4.11. Visually inspect the plots and make a manual decision------------------------------------------------
# Some are not detected within the polygons just because the polygon is more 
# coastal, but that is fine, eDNA moves with currents

# Coastal polygon not matching species observed occurrence:
#Anguilla anguilla
#Chelon auratus
#Chelon ramada
#Dicentrarchus labrax
#Dasyatis pastinaca
#Gymnura altavela
#Myliobatis aquila


# Deeper polygon not matching species observed occurence:
#Bathypterois dubius
#Centrophorus uyato


# Species not present in the Mediterranean
#Auxis thazard # Aquamaps reports it as present here
#Cyclothone atraria #
#Diaphus dumerlii
#Engraulis australis
#Engraulis japonicus
#Foetorepus agassizii
#Maurolicus muelleri
#Oncorhynchus mykiss # North East Pacific but present in Aquaculture
#Salmo salar # Aquaculture, also I would say it occurs here (AquaMaps has it in catalonia too)
#Schedophilus velaini
#Scomber japonicus
#Thunnus orientalis
#Thunnus tonggol # Indian ocean
#Trachipterus arcticus


# Species not present in the western Mediterranean:
#Gobius bucchichi
#Millerigobius macrocephalus


# Species present in the western Med but not in the reported zone:
#Parablennius zvonimiri # TRUE, the species is very common in the study region, JUST COASTAL
#Pteroplatytrygon violacea # TRUE, the species is very common in the study region
#Thunnus alalunga # TRUE, the species is very common in the study region


# Also fish with no IUCN maps:
## fish_missing_iucn
# "Scorpaena neglecta" # Pacific ocean (No asessed by IUCN)
# "Oblada melanurus"   # Not asssessed by IUCN neither Aquamaps, but present in the Mediterranean      
# "Lepidotrigla cavillone"  # No global assessment by IUCN, but common in the area (also according AquaMaps)
# "Atherina hepsetus"       # No global assessment by IUCN, but common in the area (also according AquaMaps)
# "Oedalechilus labeo"      # No global assessment by IUCN, but common in the area (also according AquaMaps)
# "Dactylopterus volitans"  # No global map by IUCN, but common in the area (also according AquaMaps)
# "Lepomis gibbosus"        # freshwater fish from norway
# "Leuciscus leuciscus"     # freshwater fish from USA
# "Epinephelus marginatus"  # It wasnt on the list of the IUCN for some reason but it is common in the area (also according AquaMaps)
# "Dagetichthys lusitanicus" # From western central atlantic (west Africa)


# 4.12. Add fish IUCN overlap information to metadata------------------------------------------------
# Metadata dataset:
#file_path <- paste0(input_data, "/processed_df")
#file <- paste0(file_path,  "/metadata_iucnMammals_Turtles_Fish.csv")

#metadata <- read.csv(file, sep = ";")
#head(metadata)

metadata <- metadata %>%
  left_join(
    fish_iucn_overlap %>%
      select(species, overlaps_iucn) %>%
      rename(
        corrected_name = species,
        overlaps_iucn_fish = overlaps_iucn
      ),
    by = "corrected_name"
  ) %>%
  mutate(
    overlaps_iucn = case_when(
      !is.na(overlaps_iucn) ~ overlaps_iucn,
      !is.na(overlaps_iucn_fish) ~ overlaps_iucn_fish,
      TRUE ~ overlaps_iucn
    )
  ) %>%
  select(-overlaps_iucn_fish)



# 5.13. Add manual visual inspection decisions to metadata


# Species visually accepted despite IUCN polygon mismatch
visual_iucn_true <- c(
  "Anguilla anguilla",
  "Chelon auratus",
  "Chelon ramada",
  "Dicentrarchus labrax",
  "Dasyatis pastinaca",
  "Gymnura altavela",
  "Myliobatis aquila",
  "Bathypterois dubius",
  "Centrophorus uyato",
  "Parablennius zvonimiri",
  "Pteroplatytrygon violacea",  
  "Pomatomus saltatrix",
  "Thunnus alalunga"
)

# Species visually rejected after IUCN inspection
visual_iucn_false <- c(
  "Auxis thazard",
  "Cyclothone atraria",
  "Diaphus dumerilii",
  "Engraulis australis",
  "Engraulis japonicus",
  "Foetorepus agassizii",
  "Maurolicus muelleri",
  "Oncorhynchus mykiss",
  "Salmo salar",
  "Schedophilus velaini",
  "Scomber japonicus",
  "Thunnus orientalis",
  "Thunnus tonggol",
  "Trachipterus arcticus",
  "Gobius bucchichi",
  "Millerigobius macrocephalus",
  "Notoscopelus kroyeri"
)

# General visual inspection for species without IUCN maps
general_visual_true <- c(
  "Oblada melanurus",
  "Lepidotrigla cavillone",
  "Atherina hepsetus",
  "Oedalechilus labeo",
  "Dactylopterus volitans",
  "Epinephelus marginatus"
)

general_visual_false <- c(
  "Scorpaena neglecta",
  "Lepomis gibbosus",
  "Leuciscus leuciscus",
  "Dagetichthys lusitanicus"
)

metadata <- metadata %>%
  mutate(
    visualInspection_IUCN = case_when(
      corrected_name %in% visual_iucn_true ~ TRUE,
      corrected_name %in% visual_iucn_false ~ FALSE,
      TRUE ~ NA
    ),
    General_visualInspection = case_when(
      corrected_name %in% general_visual_true ~ TRUE,
      corrected_name %in% general_visual_false ~ FALSE,
      TRUE ~ NA
    )
  )

# Check
metadata %>%
  filter(
    !is.na(visualInspection_IUCN) |
      !is.na(General_visualInspection)
  ) %>%
  select(
    corrected_name,
    overlaps_iucn,
    visualInspection_IUCN,
    General_visualInspection
  ) %>%
  arrange(corrected_name)

# 4.13. Final decision combining automatic and manual inspection------------------------------------------------
metadata <- metadata %>%
  mutate(
    finalIUCN_acceptance = case_when(
      is.na(overlaps_iucn) &
        is.na(visualInspection_IUCN) &
        is.na(General_visualInspection) ~ NA,
      
      coalesce(overlaps_iucn, FALSE) |
        coalesce(visualInspection_IUCN, FALSE) |
        coalesce(General_visualInspection, FALSE) ~ TRUE,
      
      TRUE ~ FALSE
    )
  )

# Check results
metadata %>%
  select(
    corrected_name,
    overlaps_iucn,
    visualInspection_IUCN,
    General_visualInspection,
    finalIUCN_acceptance
  ) %>%
  arrange(corrected_name)

# 4.14. Export updated metadata------------------------------------------------

file_path <- paste0(input_data, "/processed_df")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)

file <- paste0(file_path, "/metadata_iucn_decision.csv")

write.csv2(
  metadata,
  file,
  row.names = FALSE,
  na = ""
)
