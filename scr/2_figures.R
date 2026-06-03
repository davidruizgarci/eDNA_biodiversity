library(tidyverse)
library(janitor)
library(readr)
library(ggtext)

# Load merged eDNA presence/absence table
file <- paste0(input_data, "/processed_df/merged_eDNA_spp_pa.csv")
edna <- read_csv2(file) %>%
  clean_names()

# Load Calpe diving observations
file <- paste0(input_data, "/original_df/calpe_diving_observations.csv")
diving <- read_csv2(file) %>%
  clean_names()

# Figure 1 — eDNA species detections across sites-------------------------------
metadata_cols <- c(
  "sample_id", "location", "site", "date_of_sampling",
  "latitude", "longitude", "time_start", "time_stop",
  "sampling_comments", "cetacean_sighting", "year", "season",
  "sequencing_highlights_16sr_rna", "comments2_16sr_rna")

species_cols <- setdiff(names(edna), metadata_cols)

# Long format
edna_long <- edna_merged %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "species",
    values_to = "presence"
  ) %>%
  mutate(
    species_clean = species %>%
      str_replace_all("_", " ") %>%
      str_replace("^mega blast cf ", "cf. ") %>%
      str_replace("^mega blast ", "") %>%
      str_squish(),
    
    species_label = str_to_sentence(species_clean),
    site_label = sample_id,
    presence = as.numeric(presence)
  )

species_order <- edna_long %>%
  group_by(species_label) %>%
  summarise(n_sites = sum(presence, na.rm = TRUE), .groups = "drop") %>%
  arrange(n_sites, species_label) %>%
  pull(species_label)

edna_long <- edna_long %>%
  mutate(
    species_label = factor(species_label, levels = species_order),
    site_label = factor(site_label, levels = c(
      "240619_CAL", "240824_CABP", "240824_CABSM", "270714_MAL"
    ))
  )

fig1 <- ggplot(edna_long, aes(x = site_label, y = species_label, fill = factor(presence))) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_manual(
    values = c("0" = "grey92", "1" = "grey20"),
    labels = c("Not detected", "Detected"),
    name = NULL
  ) +
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic('", x, "')"))
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "eDNA species detections across sampling sites"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    legend.justification = "left"
  )

fig1

file_path <- paste0(output_data, "/figures")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/spp_region.png")

ggsave(
  filename = file,
  plot = fig1,
  width = 7,
  height = 10,
  dpi = 600,
  bg = "white"
)

#Figure 2 — Calpe eDNA vs diving observations-----------------------------------
clean_species_name <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("_", " ") %>%
    str_replace("^mega blast cf ", "cf. ") %>%
    str_replace("^mega blast ", "") %>%
    str_replace_all("\\?", "") %>%
    str_squish() %>%
    str_to_sentence()
}

# eDNA detections
edna_calpe <- edna_long %>%
  filter(str_detect(site_label, "CAL")) %>%
  distinct(species) %>%
  mutate(
    species = clean_species_name(species),
    eDNA = 1
  )

# Diving detections
diving_calpe <- diving %>%
  filter(tipo == "Vertebrados (Primer's target)") %>%
  mutate(
    species = clean_species_name(species),
    Diving = 1
  ) %>%
  select(species, Diving)

# Merge detections
comparison_calpe <- full_join(
  edna_calpe,
  diving_calpe,
  by = "species"
) %>%
  mutate(
    eDNA = replace_na(eDNA, 0),
    Diving = replace_na(Diving, 0),
    
    category = case_when(
      eDNA == 1 & Diving == 1 ~ "Both",
      eDNA == 1 & Diving == 0 ~ "Only eDNA",
      eDNA == 0 & Diving == 1 ~ "Only diving"
    )
  )

# Long format for plotting
comparison_long <- comparison_calpe %>%
  pivot_longer(
    cols = c(eDNA, Diving),
    names_to = "method",
    values_to = "detected"
  ) %>%
  filter(detected == 1)

# Species order
species_order <- comparison_calpe %>%
  arrange(category, species) %>%
  pull(species)

comparison_long <- comparison_long %>%
  mutate(
    species = factor(species, levels = species_order)
  )

# Plot
fig2 <- ggplot(
  comparison_long,
  aes(x = method, y = species, fill = category)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  
  scale_fill_manual(
    values = c(
      "Only diving" = "steelblue",
      "Only eDNA"   = "#d95f02",
      "Both"        = "#1b9e77"
    )
  ) +
  
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic('", x, "')"))
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = "Calpe: species detected by eDNA and diving surveys"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 11),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

fig2


file_path <- paste0(output_data, "/figures")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/spp_method.png")

ggsave(
  filename = file,
  plot = fig2,
  width = 7,
  height = 10,
  dpi = 600,
  bg = "white"
)


#Figure 3 — Taxonomic composition of eDNA detections----------------------------
# Simple classification from species names only
# Later we can improve this by joining with eDNA_TAX taxonomy
edna_composition <- edna_long %>%
  
  # Keep only detections
  filter(presence == 1) %>%
  
  mutate(
    broad_group = case_when(
      
      # Marine turtles
      str_detect(species, "caretta") ~ "Marine turtle",
      
      # Cetaceans
      str_detect(species,
                 "stenella|tursiops|delphin|balaenoptera") ~ "Cetacean",
      
      # Chondrichthyans
      str_detect(species,
                 "dasyatis|prionace") ~ "Chondrichthyan",
      
      # Everything else
      TRUE ~ "Teleost"
    )
  ) %>%
  
  count(sample_id, broad_group, name = "n_species")

# Order sites
edna_composition <- edna_composition %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = c(
        "240619_CAL",
        "240824_CABP",
        "240824_CABSM",
        "270714_MAL"
      )
    )
  )

# Plot
fig3 <- ggplot(
  edna_composition,
  aes(x = sample_id,
      y = n_species,
      fill = broad_group)
) +
  
  geom_col(width = 0.75) +
  
  scale_fill_manual(
    values = c(
      "Teleost" = "#4C78A8",
      "Chondrichthyan" = "#F58518",
      "Cetacean" = "#54A24B",
      "Marine turtle" = "#B279A2"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Number of detected species",
    fill = NULL,
    title = "Taxonomic composition of eDNA detections"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

fig3

file_path <- paste0(output_data, "/figures")
if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
file <- paste0(file_path, "/spp_prop_site.png")

ggsave(
  filename = file,
  plot = fig3,
  width = 10,
  height = 6,
  dpi = 600,
  bg = "white"
)
