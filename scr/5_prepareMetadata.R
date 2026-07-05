#-------------------------------------------------------------------------------
#
# 5. Prepare metadata for BLAST analyses
#
#-------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# 1. Load metadata data---------------------------------------------------------

file_path <- paste0(input_data, "/processed_df")

metadata_file <- paste0(file_path, "/metadata_iucn_aquamaps_decision.csv")
metadata <- read.csv(metadata_file, sep = ";")

data_file <- paste0(file_path, "/merged_eDNA_spp_pa.csv")
data <- read.csv(data_file, sep = ";") %>%
  filter(!sequencing_comment %in% c("Filtration_Control", "No_target_taxa"))

head(metadata)
names(metadata)
head(data)




# 2. Plot differences in distribution checks------------------------------------

# Prepare output folder
plot_dir <- file.path(input_data, "figures", "distributionChecks")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare data
metadata <- metadata %>%
  filter(TaxaTargeted == "Yes")   # or == TRUE if logical

plot_df <- metadata %>%
  distinct(
    valid_name,
    overlaps_iucn,
    overlaps_aquamaps_p04,
    finalIUCN_acceptance,
    finalAquamaps_acceptance
  ) %>%
  mutate(
    disagreement_overlap = dplyr::coalesce(overlaps_iucn, FALSE) !=
      dplyr::coalesce(overlaps_aquamaps_p04, FALSE),
    disagreement_final = dplyr::coalesce(finalIUCN_acceptance, FALSE) !=
      dplyr::coalesce(finalAquamaps_acceptance, FALSE)
  ) %>%
  arrange(desc(disagreement_final), desc(disagreement_overlap), valid_name) %>%
  mutate(valid_name = factor(valid_name, levels = unique(valid_name))) %>%
  dplyr::select(
    valid_name,
    overlaps_iucn,
    overlaps_aquamaps_p04,
    finalIUCN_acceptance,
    finalAquamaps_acceptance
  ) %>%
  pivot_longer(
    cols = -valid_name,
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "overlaps_iucn",
        "overlaps_aquamaps_p04",
        "finalIUCN_acceptance",
        "finalAquamaps_acceptance"
      ),
      labels = c(
        "IUCN overlap",
        "AquaMaps overlap",
        "Final IUCN",
        "Final AquaMaps"
      )
    ),
    Value = case_when(
      Value == TRUE ~ "TRUE",
      Value == FALSE ~ "FALSE",
      is.na(Value) ~ NA_character_
    ),
    Value = factor(Value, levels = c("TRUE", "FALSE"))
  )



# Plot
p_distribution_checks <- ggplot(
  plot_df,
  aes(x = Variable, y = valid_name, fill = Value)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "TRUE" = "#2E8B57",
      "FALSE" = "#D73027"
    ),
    na.value = "grey80",
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Species",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank()
  )

p_distribution_checks

# Save figure
ggsave(
  filename = file.path(plot_dir, "IUCN_AquaMaps_distribution_checks.png"),
  plot = p_distribution_checks,
  width = 8,
  height = max(6, length(unique(metadata$valid_name)) * 0.12),
  dpi = 300,
  bg = "white"
)

# 3. List disagreements---------------------------------------------------------

# One row per species
checks_df <- metadata %>%
  distinct(
    valid_name,
    overlaps_iucn,
    overlaps_aquamaps_p04,
    finalIUCN_acceptance,
    finalAquamaps_acceptance
  )

# Species with disagreement in raw overlap comparison
disagreement_overlap_species <- checks_df %>%
  filter(
    !is.na(overlaps_iucn),
    !is.na(overlaps_aquamaps_p04),
    overlaps_iucn != overlaps_aquamaps_p04
  ) %>%
  pull(valid_name) %>%
  unique() %>%
  sort()

# Species with disagreement in final acceptance comparison
disagreement_final_species <- checks_df %>%
  filter(
    !is.na(finalIUCN_acceptance),
    !is.na(finalAquamaps_acceptance),
    finalIUCN_acceptance != finalAquamaps_acceptance
  ) %>%
  pull(valid_name) %>%
  unique() %>%
  sort()

# Differences between the two disagreement lists
only_overlap_disagreement <- setdiff(
  disagreement_overlap_species,
  disagreement_final_species
)

only_final_disagreement <- setdiff(
  disagreement_final_species,
  disagreement_overlap_species
)

in_both_disagreement_lists <- intersect(
  disagreement_overlap_species,
  disagreement_final_species
)

# Species where all four columns are NA
all_na_species <- checks_df %>%
  filter(
    is.na(overlaps_iucn),
    is.na(overlaps_aquamaps_p04),
    is.na(finalIUCN_acceptance),
    is.na(finalAquamaps_acceptance)
  ) %>%
  pull(valid_name) %>%
  unique() %>%
  sort()

# Summary table
disagreement_summary <- tibble::tibble(
  comparison = c(
    "Overlap disagreement",
    "Final acceptance disagreement",
    "Only overlap disagreement",
    "Only final disagreement",
    "In both disagreement lists",
    "All four values NA"
  ),
  n_species = c(
    length(disagreement_overlap_species),
    length(disagreement_final_species),
    length(only_overlap_disagreement),
    length(only_final_disagreement),
    length(in_both_disagreement_lists),
    length(all_na_species)
  )
)

disagreement_summary

# Print lists
disagreement_overlap_species # Species which RAW IUCN and Aquamaps result is different 
disagreement_final_species # Species which VISUALLY CHECKED IUCN and Aquamaps result is different 
only_overlap_disagreement #These species disagreed initially but no longer disagree after manual inspection
only_final_disagreement
in_both_disagreement_lists
all_na_species

# Comments on final decisions:
# Maurolicus muelleri is reported in the region but it seems to be very rare
# Notoscopelus kroyeri is reported in the region but it seems to be very rare
# Trachipterus arcticus is reported in the region but it seems to be very rare
# Taxonomic issues Auxis thazard may be Auxis rochei check carefully


# Especies no habían solapado con los mapas de la IUCN:
species_overlaps_iucn_false <- metadata %>%
  filter(
    TaxaTargeted == "Yes",
    overlaps_iucn == FALSE
  ) %>%
  distinct(valid_name) %>%
  arrange(valid_name) %>%
  pull(valid_name)

species_overlaps_iucn_false
