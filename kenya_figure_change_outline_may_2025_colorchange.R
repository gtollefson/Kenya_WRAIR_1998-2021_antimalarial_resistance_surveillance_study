# Kenya Longitudinal Antimalarial Resistance Study - Figure Generation Script
# Automated dependency management and figure generation

# Install and load pacman for package management
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
  library(pacman)
}

# Install and load all required packages
# pacman automatically installs missing packages and loads them
# No need for individual library() calls throughout the script
pacman::p_load(
  dplyr,           # Data manipulation
  ggplot2,         # Plotting
  here,            # File path management
  readr,           # Reading CSV files
  ggpubr,          # Publication ready plots
  rstatix,         # Statistical tests
  tidyr,           # Data tidying
  sf,              # Spatial data
  rnaturalearth,   # Natural earth data
  malariaAtlas,    # Malaria data
  viridis,         # Color palettes
  scales,          # Scale functions
  cowplot,         # Plot arrangements
  gridExtra,       # Grid arrangements
  RColorBrewer,    # Color palettes
  forcats,         # Factor handling
  miplicorn,       # MIP analysis tools
  scico,           # Scientific color palettes
  stringr,         # String manipulation
  terra            # Spatial data analysis
)

# All packages are now loaded - no need for library() calls below

# ==============================================================================
# FIGURE GENERATION FOR KENYA LONGITUDINAL ANTIMALARIAL RESISTANCE STUDY
# ==============================================================================
#
# This script generates all figures and supplementary tables for the manuscript:
# "Temporal trends in antimalarial resistance markers in Kenya, 1998-2021"
#
# Outputs:
# - Figure 1: Sample collection overview and mutation prevalence comparison
# - Figure 2: MDR1 mutations and haplotype analysis  
# - Figure 3: CRT mutations and chloroquine resistance reversion
# - Figure 4: DHFR mutations and pyrimethamine resistance
# - Figure 5: DHPS mutations and sulfadoxine resistance
# - Supplementary Tables: Sample summaries and mutation prevalences
#
# ==============================================================================

# #TODO: 

# Supplemental materials outline
# S1. Sample collection summary table (number by year and site).
# S2. Details about probe target regions and reagents used for MIP capture and sequencing, target mutations
# S3. K13 propeller domain missense mutation incidence table
# DONE S4. Prevalence table of 77 individual key antimalarial resistance markers
# S5. Table of prevalence of artemisinin resistance mutations in UBP1, ATP6, AP2MU
# S6. Individual level within-sample allele frequency data for all DR2 panel target mutations

# SF1: K13 propeller domain missense mutation incidence plot (across protein)

# S1:Sample collection summary table:

# Load the table
variant_table <- read.csv(here("data", "kenya_all_DR2_variant_table.csv"), header = TRUE)
# Convert Year_start to numeric if needed
#variant_table$Year_start <- as.numeric(variant_table$Year_start)

# Get unique samples only
unique_samples <- variant_table %>%
  select(Sampleid, Year_start, site_uid) %>%
  distinct()

# Summarize samples by year and site
samples_by_year_site <- unique_samples %>%
  group_by(Year_start, site_uid) %>%
  summarise(Sample_Count = n(), .groups = "drop")

# Summarize total samples per year
samples_by_year <- unique_samples %>%
  group_by(Year_start) %>%
  summarise(Sample_Count = n(), .groups = "drop")

# Save the summary table
write.csv(samples_by_year, "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_supplement/Table_S1_samples_by_year.csv", row.names = FALSE)

# Plot
plot <- ggplot(samples_by_year, aes(x = factor(Year_start), y = Sample_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Total Samples Sequenced Per Year",
       x = "Year",
       y = "Number of Unique Samples") +
  theme_minimal()

# Save the plot
base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_supplement/Supp_figures/"

library(ggplot2)
# Save plots in multiple formats
formats <- c("svg", "png", "pdf")

for (ext in formats) {
  ggsave(paste0(base_path,"samples_by_year_barplot.", ext), 
         plot = plot, width = 15, height = 10, units = "in", dpi = 600)
}

# S2: Probe target regions tables
# Done in py script 

# S3: K13 propeller domain missense mutation incidence table

library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)

# Load necessary libraries
library(dplyr)
library(readr)

# Define the file path
file_path <- here("data", "kenya_all_DR2_prevalence_results_table.csv")

# Read in the prevalence table
all_prev_table <- read_csv(file_path)

# Print the first few rows to verify
print(head(all_prev_table))

# Remove Eastern kenya sites with too few samples spanning entire collection period
# Define sites to exclude
excluded_sites <- c("Kilifi", "Malindi")

# Filter out the unwanted rows
filtered_table <- all_prev_table[!(all_prev_table$site_uid %in% excluded_sites), ]


# Filter filtered_table for k13 mutations with positions >= 460, valid prevalence, no underscores or fs
k13_mutations <- c("Cys473Phe", "Ala569Ser", "Ala578Ser")
# 
# k13_filtered <- filtered_table %>%
#   filter(gene_name == "k13", !is.na(prevalence_percentage)) %>%  # Filter k13 with valid prevalence
#   mutate(
#     aa_position = as.numeric(gsub("[^0-9]", "", mutation))  # Extract numeric position
#   ) %>%
#   filter(aa_position >= 460) %>%  # Filter for positions >= 460
#   filter(!grepl("_", mutation) & !grepl("fs", mutation))  # Exclude invalid mutations

k13_filtered <- filtered_table %>%
  filter(gene_name == "k13") %>%
  filter(mutation %in% k13_mutations)

library(dplyr)
library(ggplot2)

# Add time period column
k13_filtered <- k13_filtered %>%
  mutate(
    time_period = case_when(
      Year_start >= 1999 & Year_start <= 2005 ~ "SP",
      Year_start >= 2006 & Year_start <= 2021 ~ "ACT",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(time_period))  # Remove any rows that don’t fit the time periods

# Calculate average prevalence per mutation per time period
k13_summary <- k13_filtered %>%
  group_by(mutation, time_period) %>%
  summarise(avg_prevalence = mean(prevalence_percentage, na.rm = TRUE), .groups = 'drop')

# # Define the color palette
# palette <- c("SP" = "#00AFBB", "ACT" = "#E7B800")
# 
# # Plot the barplot
# ggplot(k13_summary, aes(x = mutation, y = avg_prevalence, fill = time_period)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
#   scale_fill_manual(values = palette) +
#   labs(
#     title = "Average Prevalence of K13 Propeller Domain Mutations Across Eras",
#     x = "Mutation",
#     y = "Average Prevalence (%)",
#     fill = "Treatment Era"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#     axis.text.y = element_text(size = 10),
#     axis.title = element_text(size = 12, face = "bold"),
#     plot.title = element_text(size = 14, face = "bold")
#   )

# get novel k13 counts
samples_with_mutation_counts <- k13_filtered %>%
  mutate(Year_start = as.integer(Year_start)) %>%
  group_by(site_uid, Year_start, mutation) %>%
  summarise(samples_with_mutation = sum(samples_with_mutation, na.rm = TRUE), .groups = "drop") %>%
  filter(samples_with_mutation > 0) %>%
  arrange(site_uid, Year_start, mutation)

# View result
print(samples_with_mutation_counts)
# Assuming your tibble is called samples_with_mutation_counts
# And "figureoutput" is a directory in your working directory

base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_supplement/"
write.csv(samples_with_mutation_counts, paste0(base_path, "Table_S3_novel_k13_mutation_counts.csv"), row.names = FALSE)

# S4: 77 prevalences
file_path <- here("data", "kenya_all_DR2_prevalence_results_table_with_reference_resistant.csv")
supp_base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_supplement/"
all_prev_table <- read_csv(file_path)
excluded_sites <- c("Kilifi", "Malindi")
filtered_table <- all_prev_table[!(all_prev_table$site_uid %in% excluded_sites), ]
targets_file_path <- here("data", "targets_k13updated2023.tsv")
targets_table <- read_tsv(targets_file_path)
targets_table <- targets_table %>%
  mutate(mutation_name_match = paste0(gene_name, "-", aminoacid_change))
filtered_matched <- filtered_table %>%
  filter(mutation_name %in% targets_table$mutation_name_match)

write.table(filtered_matched,paste0(supp_base_path,"Table_S4_of_77_key_marker_prevalences.csv"),sep=",",quote=F)

# S5: Table of art mutations in pfap2mu, atp6, pfubp1
library(dplyr)
library(readr)

# Load the prevalence table (already includes reference_resistant column)
file_path <- here("data", "kenya_all_DR2_prevalence_results_table_with_reference_resistant.csv")
all_prev_table <- read_csv(file_path)

# Exclude unwanted sites
excluded_sites <- c("Kilifi", "Malindi")
filtered_table <- all_prev_table %>%
  filter(!(site_uid %in% excluded_sites))

# Ensure gene names are lowercase
filtered_table <- filtered_table %>%
  mutate(gene_name = tolower(gene_name))

# Filter for target genes
genes_of_interest <- c("pfap2mu", "pfubp1", "atp6")
filtered_gene_subset <- filtered_table %>%
  filter(gene_name %in% genes_of_interest)

# Save filtered table
supp_base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_supplement/"
write_csv(filtered_gene_subset, paste0(supp_base_path, "Table_S5_ubp1_atp6_ap2mu_marker_prevalences.csv"))

# S6: Individual level AF and coverage table for all missense SNPs detected

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(rstatix)
library(miplicorn)
library(tidyr)
library(viridis)

# Define file paths for metadata and miptools tables
metadata_table_path <- here("data", "metadata_DR_cleaned_no_dups_renamed_cols_from_hack.csv")
#old NO CRT 73, 74 #miptools_tables_dir <- "/Users/george/Bailey_lab/kenya_longitudinal/2024_analyses/miptools_tables"
# Variant run with miptools 4.10 on march 3,2025 
miptools_tables_dir <- here("data")

# Load metadata
metadata_table <- read.csv(metadata_table_path, header = TRUE)
met <- metadata_table %>% rename(Sampleid = sample_uid)

# Define paths for miptools reference, alternate, and coverage tables
ref_file_path <- paste0(miptools_tables_dir, "/reference_AA_table.csv")
alt_file_path <- paste0(miptools_tables_dir, "/alternate_AA_table.csv")
cov_file_path <- paste0(miptools_tables_dir, "/coverage_AA_table.csv")
genotype_file_path <- paste0(miptools_tables_dir, "/genotypes_AA_table.csv")

# Read miptools tables using custom functions
ref_file <- read_tbl_reference(ref_file_path)
alt_file <- read_tbl_alternate(alt_file_path)
cov_file <- read_tbl_coverage(cov_file_path)
genotype_file <-read_tbl_genotype(genotype_file_path) # for filtering out any samples with more than 1 heterozygous call (accomodate some low frequency seq/calling error)

# Function to filter files
filter_mutations <- function(df) {
  df[!grepl("ins|del|dup", df$mutation_name, ignore.case = TRUE) &
       df$exonic_func == "missense_variant", ]
}

# Apply to all files
ref_file      <- filter_mutations(ref_file)
alt_file      <- filter_mutations(alt_file)
cov_file      <- filter_mutations(cov_file)
genotype_file <- filter_mutations(genotype_file)

# Assemble the variant table with ref, alt, and coverage data
# This records presence/absence of all MIPtools reported variants for each sample with ref and alt allele coverage and sample metadata
variant_table <- data.frame(
  Sampleid = ref_file$sample,
  gene_id = ref_file$gene_id,
  gene_name = ref_file$gene,
  mutation_name = ref_file$mutation_name,
  mutation = ref_file$aa_change,
  umi_ref = ref_file$ref_umi_count,
  umi_alt = alt_file$alt_umi_count,
  umi_cov = cov_file$coverage
)

# GT analysis specific - rename samples for compatibility between metadata sampleid nomenclature and miptools table sampleids
variant_table$Sampleid <- gsub("S-","S", variant_table$Sampleid)

# Add total coverage (UMI ref + UMI alt)
variant_table <- variant_table %>% mutate(total_cov = umi_ref + umi_alt)

# Filter by coverage threshold of at least 2 UMIs
variant_table_cov_filtered <- variant_table %>% filter(total_cov >= 2)

# Merge metadata with the variant table
variant_table_cov_filtered_with_met <- merge(met, variant_table_cov_filtered, by = "Sampleid", all.x = TRUE)

# filter to remove east sites
filtered_variant_table <- variant_table_cov_filtered_with_met %>%
  filter(!(site_uid %in% c("Kilifi", "Malindi")))
# get total sample size after removing eastern sites
length(unique(filtered_variant_table$Sampleid))

# get sample size by era
filtered_variant_table %>%
  distinct(Sampleid, Year_start) %>%  # Count each sample once per year
  mutate(
    era = case_when(
      Year_start %in% c("1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "PRE-ACT") ~ "1998–2005 or PRE-ACT",
      Year_start %in% as.character(2006:2021) ~ "2006–2021",
      TRUE ~ "Other"
    )
  ) %>%
  count(era)

# Save final variant table to CSV so you dont need to rerun miptools tables loading and variant table formation
write.csv(filtered_variant_table, here("data", "kenya_all_DR2_variant_table_missense_SNPs_only.csv"), row.names = FALSE)

filtered_variant_table <- read.csv(here("data", "kenya_all_DR2_variant_table_missense_SNPs_only.csv"))

# Add allele_frequency column
filtered_variant_table <- filtered_variant_table %>%
  mutate(allele_frequency = umi_alt / (umi_alt + umi_ref))

# Optional: write the new table to file
write.csv(
  filtered_variant_table,
  "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_supplement/Table_S6_Kenya_20_year_surveillance_full_cohort_all_variant_calls.csv",
  row.names = FALSE
)
# SF1: K13 propeller domain missense mutation plot





################################################################################
#
# Main figure start
#
################################################################################

# make 1A collection markres collored by SP and ACT period and make same colros as used in figures yellow nad purple
# timeline increase font size
# single lettter dual coding for 1C
#Color sceme for 1C based on 2A
# use new font size setting for 1C

library(scico)
scico_palette_show()

# examine palettes (choose 8 and 2 reserved)
# library(scico)
# library(scales)

# # Corrected list of palettes
# palette_names <- c("navia", "batlow", "lajolla", "bilbao", "berlin", "managua", "lipari")
# reserved_indices <- c(5, 15)
# n_colors_needed <- 8
# 
# par(mfrow = c(4, 2), mar = c(2, 1, 3, 1))  # adjust layout
# 
# for (pname in palette_names) {
#   full_palette <- scico(n = 20, palette = pname)
#   palette_for_other_figures_pre <- full_palette[-reserved_indices]
#   
#   step <- floor(length(palette_for_other_figures_pre) / n_colors_needed)
#   selected_indices <- seq(1, length(palette_for_other_figures_pre), by = step)[1:n_colors_needed]
#   palette_for_other_figures <- palette_for_other_figures_pre[selected_indices]
#   
#   barplot(rep(1, n_colors_needed), col = palette_for_other_figures, border = NA,
#           space = 0, axes = FALSE, main = pname)
# }


# make color pallete for boxplots and other figures (reserve boxplot colors for boxplot only)

library(scico)
# make palettes for figures with scrambled order
library(scico)
library(scales)

# Step 1: Create master palette
full_palette <- scico(n = 20, palette = "managua")

# Step 2: Reserve colors
reserved_indices <- c(5, 15)
reserved_colors <- full_palette[reserved_indices]

treatment_periods <- c("SP Period (1998-2005)","ACT Period (2006-2021)")

# Step 3: Assign reserved colors to time levels
time_levels <- treatment_periods
boxplot_colors <- setNames(reserved_colors[1:length(time_levels)], time_levels)

# Step 4: Remove reserved colors
palette_for_other_figures_pre <- full_palette[-reserved_indices]

# Step 5: Select evenly spaced colors
n_colors_needed <- 8
step <- floor(length(palette_for_other_figures_pre) / n_colors_needed)
selected_indices <- seq(1, length(palette_for_other_figures_pre), by = step)[1:n_colors_needed]
palette_for_other_figures <- palette_for_other_figures_pre[selected_indices]

# Step 6: Scramble order to avoid gradient appearance
scrambled_order <- c(1, 8, 3, 6, 2, 7, 4, 5)  # You can tweak this
palette_for_other_figures <- palette_for_other_figures[scrambled_order]

# Step 7: View result
scales::show_col(palette_for_other_figures)
# # Step 1: Create master palette
# full_palette <- scico(n = 20, palette = "managua")
# 
# # Step 2: Reserve colors for boxplot
# reserved_indices <- c(5, 15)
# reserved_colors <- full_palette[reserved_indices]
# 
# # Step 3: Assign reserved colors to your time_period levels
# time_levels <- levels(gene_data$time_period)
# boxplot_colors <- setNames(reserved_colors[1:length(time_levels)], time_levels)
# 
# # Step 4: Remove reserved colors for all other figures
# # Assuming palette_for_other_figures has 18 colors
# palette_for_other_figures_pre <- full_palette[-reserved_indices]
# scales::show_col(palette_for_other_figures_pre)
# 
# # Select 7 reasonably spaced colors
# n_colors_needed <- 8
# step <- floor(length(palette_for_other_figures_pre) / n_colors_needed)
# 
# # This picks every `step`th color starting from the first
# selected_indices <- seq(1, length(palette_for_other_figures_pre), by = step)[1:n_colors_needed]
# 
# palette_for_other_figures <- palette_for_other_figures_pre[selected_indices]
# 
# # View result
# scales::show_col(palette_for_other_figures)

library(scico)
library(viridis)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

axis.title.setting_fig1 <- 30
axis.text.setting_fig1 <- 24
plot.title.setting_fig1 <- 32
legend.title.setting_fig1 <- 26
legend.text.setting_fig1 <- 26

axis.title.setting <- 40
axis.text.setting <- 40
plot.title.setting <- 42
legend.title.setting <- 40
legend.text.setting <- 40

base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data"

convert_aaa_to_a_with_gene_name <- function(x) {
  aa3 <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His",
           "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
           "Tyr", "Val")
  aa1 <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K",
           "M", "F", "P", "S", "T", "W", "Y", "V")
  ref <- setNames(aa1, aa3)
  
  sapply(x, function(label) {
    parts <- strsplit(label, "-", fixed = TRUE)[[1]]
    
    if (length(parts) != 2) return(NA)
    
    gene <- parts[1]
    mutation <- parts[2]
    
    pattern <- "^([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2})$"
    match <- regexec(pattern, mutation)
    m <- regmatches(mutation, match)[[1]]
    
    if (length(m) != 4) return(NA)
    
    from <- m[2]
    pos  <- m[3]
    to   <- m[4]
    
    from1 <- ref[[from]]
    to1   <- ref[[to]]
    
    if (is.null(from1) || is.null(to1)) return(NA)
    
    paste0(gene, "-", from1, pos, to1)
  }, USE.NAMES = FALSE)
}


convert_aaa_to_a <- function(x) {
  aa3 <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His",
           "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
           "Tyr", "Val")
  aa1 <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K",
           "M", "F", "P", "S", "T", "W", "Y", "V")
  ref <- setNames(aa1, aa3)
  
  sapply(x, function(label) {
    pattern <- "^([A-Z][a-z]{2})?([0-9]+)([A-Z][a-z]{2})?$"
    match <- regexec(pattern, label)
    parts <- regmatches(label, match)[[1]]
    
    if (length(parts) == 0) return(NA)
    
    from <- parts[2]
    pos <- parts[3]
    to <- parts[4]
    
    from1 <- if (!is.na(from) && from %in% names(ref)) ref[[from]] else ""
    to1 <- if (!is.na(to) && to %in% names(ref)) ref[[to]] else ""
    
    paste0(from1, pos, to1)
  }, USE.NAMES = FALSE)
}

# Function to save plots and data
save_results <- function(gene_name, plot, data, stat_results,fig_num) {
  base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
  
  # Convert list columns in `stat_results` to character (if any)
  if (any(sapply(stat_results, is.list))) {
    stat_results <- stat_results %>%
      mutate(across(where(is.list), ~ sapply(., toString)))
  }
  
  # Save results as CSV
  write.csv(stat_results, paste0(base_path, gene_name, "_stat_results.csv"), row.names = FALSE)
  write.csv(data, paste0(base_path, gene_name, "_data.csv"), row.names = FALSE)
  
  # Save plots in multiple formats
  formats <- c("svg", "png", "pdf")
  for (ext in formats) {
    ggsave(paste0(base_path,fig_num,"boxplot_", gene_name, ".", ext), 
           plot = plot, width = 15, height = 10, units = "in", dpi = 600)
  }
}


# Define the function
transform_haplotype_labels <- function(haplotype_labels) {
  
  # Convert three-letter amino acid codes to one-letter codes
  aa_conversion <- c(
    "Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D",
    "Cys" = "C", "Gln" = "Q", "Glu" = "E", "Gly" = "G",
    "His" = "H", "Ile" = "I", "Leu" = "L", "Lys" = "K",
    "Met" = "M", "Phe" = "F", "Pro" = "P", "Ser" = "S",
    "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V"
  )
  
  transformed_labels <- haplotype_labels  # Initialize with original values
  
  # Process each haplotype label
  for (i in seq_along(haplotype_labels)) {
    
    label <- haplotype_labels[i]
    
    if (label == "WT") next  # Skip WT
    
    # Find the first valid mutation (Three-letter AA -> Number -> Three-letter AA)
    first_match <- str_extract(label, "\\w*-(\\w{3}\\d+\\w{3})")
    
    if (is.na(first_match)) next  # Skip if no valid mutation is found
    
    # Remove everything before the first valid mutation
    label <- sub("^.*?-", "", label)
    
    # Extract all mutations in the format "Three-letter AA -> Number -> Three-letter AA"
    mutations <- str_extract_all(label, "(\\w{3})(\\d+)(\\w{3})")[[1]]
    
    if (length(mutations) == 0) next  # Skip if no matches
    
    formatted_mutations <- c()
    
    for (m in mutations) {
      parts <- unlist(str_match(m, "(\\w{3})(\\d+)(\\w{3})"))  # Extract AA, Number, AA
      if (length(parts) < 4) next  # Ensure valid match
      
      res_num <- parts[3]  # Residue number
      to_aa <- parts[4]    # Mutant amino acid
      
      # Convert to one-letter amino acid code
      one_letter_to <- aa_conversion[to_aa]
      if (is.na(one_letter_to)) next  # Skip if conversion fails
      
      formatted_mutations <- c(formatted_mutations, paste0(res_num, one_letter_to))
    }
    
    # Sort mutations by residue number
    formatted_mutations <- formatted_mutations[order(as.numeric(str_extract(formatted_mutations, "\\d+")))]
    
    # Create final formatted string
    if (length(formatted_mutations) > 0) {
      mutation_string <- paste(formatted_mutations, collapse = "-")
      aa_only_string <- paste0("(", paste(str_extract(formatted_mutations, "[A-Z]$"), collapse = ""), ")")
      transformed_labels[i] <- paste(mutation_string, aa_only_string)
    }
  }
  
  return(transformed_labels)
}

################
#
#Data prep for all figures:
#
#############


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(rstatix)
library(miplicorn)
library(tidyr)
library(viridis)

# Define file paths for metadata and miptools tables
metadata_table_path <- here("data", "metadata_DR_cleaned_no_dups_renamed_cols_from_hack.csv")
#old NO CRT 73, 74 #miptools_tables_dir <- "/Users/george/Bailey_lab/kenya_longitudinal/2024_analyses/miptools_tables"
# Variant run with miptools 4.10 on march 3,2025 
miptools_tables_dir <- here("data")

# Load metadata
metadata_table <- read.csv(metadata_table_path, header = TRUE)
met <- metadata_table %>% rename(Sampleid = sample_uid)

# Define paths for miptools reference, alternate, and coverage tables
ref_file_path <- paste0(miptools_tables_dir, "/reference_AA_table.csv")
alt_file_path <- paste0(miptools_tables_dir, "/alternate_AA_table.csv")
cov_file_path <- paste0(miptools_tables_dir, "/coverage_AA_table.csv")
genotype_file_path <- paste0(miptools_tables_dir, "/genotypes_AA_table.csv")

# Read miptools tables using custom functions
ref_file <- read_tbl_reference(ref_file_path)
alt_file <- read_tbl_alternate(alt_file_path)
cov_file <- read_tbl_coverage(cov_file_path)
genotype_file <-read_tbl_genotype(genotype_file_path) # for filtering out any samples with more than 1 heterozygous call (accomodate some low frequency seq/calling error)

# Assemble the variant table with ref, alt, and coverage data
# This records presence/absence of all MIPtools reported variants for each sample with ref and alt allele coverage and sample metadata
variant_table <- data.frame(
  Sampleid = ref_file$sample,
  gene_id = ref_file$gene_id,
  gene_name = ref_file$gene,
  mutation_name = ref_file$mutation_name,
  mutation = ref_file$aa_change,
  umi_ref = ref_file$ref_umi_count,
  umi_alt = alt_file$alt_umi_count,
  umi_cov = cov_file$coverage
)

# GT analysis specific - rename samples for compatibility between metadata sampleid nomenclature and miptools table sampleids
variant_table$Sampleid <- gsub("S-","S", variant_table$Sampleid)

# Add total coverage (UMI ref + UMI alt)
variant_table <- variant_table %>% mutate(total_cov = umi_ref + umi_alt)

# Filter by coverage threshold of at least 2 UMIs
variant_table_cov_filtered <- variant_table %>% filter(total_cov >= 2)

# Merge metadata with the variant table
variant_table_cov_filtered_with_met <- merge(met, variant_table_cov_filtered, by = "Sampleid", all.x = TRUE)

# filter to remove east sites
filtered_variant_table <- variant_table_cov_filtered_with_met %>%
  filter(!(site_uid %in% c("Kilifi", "Malindi")))
# get total sample size after removing eastern sites
length(unique(filtered_variant_table$Sampleid))

# get sample size by era
filtered_variant_table %>%
  distinct(Sampleid, Year_start) %>%  # Count each sample once per year
  mutate(
    era = case_when(
      Year_start %in% c("1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "PRE-ACT") ~ "1998–2005 or PRE-ACT",
      Year_start %in% as.character(2006:2021) ~ "2006–2021",
      TRUE ~ "Other"
    )
  ) %>%
  count(era)

# Save final variant table to CSV so you dont need to rerun miptools tables loading and variant table formation
write.csv(filtered_variant_table, "/Users/george/Bailey_lab/kenya_longitudinal/2024_analyses/kenya_all_DR2_variant_table.csv", row.names = FALSE)

variant_table_cov_filtered_with_met <- filtered_variant_table

# Define reference-resistant alleles of interest
reference_resistant_labels <- c("mdr1-Asn86Tyr", "mdr1-Asp1246Tyr", "crt-Cys72Ser", "crt-Val73Leu")

# Calculate prevalence by site and year while preserving mutation column
all_prev_table_with_ref_resistant <- variant_table_cov_filtered_with_met %>%
  group_by(site_uid, Year_start, gene_name, mutation_name, mutation) %>%  # Ensure mutation column is kept
  summarise(
    total_samples = n_distinct(Sampleid),
    samples_with_adequate_coverage = n_distinct(Sampleid[total_cov >= 2]),
    samples_with_mutation = n_distinct(Sampleid[total_cov >= 2 & umi_alt > 0]),
    samples_with_reference_resistant = n_distinct(Sampleid[total_cov >= 2 & umi_ref > 0 & mutation_name %in% reference_resistant_labels])
  ) %>%
  ungroup() %>%
  mutate(
    prevalence_percentage = (samples_with_mutation / samples_with_adequate_coverage) * 100,
    reference_resistant_prevalence = (samples_with_reference_resistant / samples_with_adequate_coverage) * 100
  )

# Print first few rows to verify
print(head(all_prev_table_with_ref_resistant))

# Save final prevalence table to CSV
write.csv(all_prev_table_with_ref_resistant, here("data", "kenya_all_DR2_prevalence_results_table_with_reference_resistant.csv"), row.names = FALSE)

################################################################################
# 
# Figure 1) Collection summary and barplot
# 
################################################################################
# 1A

# get sample size per size summary
site_sample_counts <- metadata_table %>%
  group_by(site_uid) %>%
  summarise(sample_count = n()) %>%
  arrange(desc(sample_count))

print(site_sample_counts)

site_sample_counts_west <- site_sample_counts[-c(4,10),]
sum(site_sample_counts_west$sample_count)

# setup map plotting
# Load required libraries
library(ggplot2)
library(sf)
library(readr)
library(dplyr)
library(malariaAtlas)
library(rnaturalearth)

# Load site data
sites_data_pre <- read_csv(here("data", "region_mdr_sites_samples.csv"))

# exckude Kilifi and Malindi as Eastern Kenya is not sampled sufficiently over time to include
sites_data_pre_west <- sites_data_pre[-(which(sites_data_pre$name=="Kilifi (4)" | sites_data_pre$name == "Malindi (58)")),]

# take lat and lon points
sites_data <- sites_data_pre_west[,4:5]

# Ensure latitude and longitude are numeric
sites_data$lat <- as.numeric(sites_data$lat)
sites_data$lon <- as.numeric(sites_data$lon)
colnames(sites_data)[1:2] <- c("latitude", "longitude")

# Download and load Kenya boundary and regional shapefiles from malariaAtlas
# Load Admin 0 and Admin 1 shapefiles for Kenya
sf_adm0 <- malariaAtlas::getShp(ISO = "KEN", admin_level = "admin0") %>% st_as_sf()
sf_adm1 <- malariaAtlas::getShp(ISO = "KEN", admin_level = "admin1") %>% st_as_sf()

# Load surrounding countries for context (Admin 0 shapefiles)
context_countries <- malariaAtlas::getShp(ISO = c("TZA", "UGA", "SOM", "ETH", "SSD"), admin_level = "admin0") %>% st_as_sf()

# Load major lakes and rivers using rnaturalearth
lakes <- ne_download(scale = "medium", type = 'lakes', category = 'physical', returnclass = "sf")
rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = "physical", returnclass = "sf")

# Filter for Lake Victoria
lake_victoria <- lakes[lakes$name == "Lake Victoria", ]

# Filter for Lake Victoria
lake_turkana <- lakes[lakes$name == "Lake Turkana", ]

# Intersect rivers with Kenya's boundary for major waterways within Kenya
rivers_in_kenya <- st_intersection(rivers, sf_adm0)

# Zoom in on Kenya by setting coordinate limits to the bounding box of Kenya
kenya_bbox <- st_bbox(sf_adm0)

# Define coordinates for the country labels
country_labels <- data.frame(
  country = c("Somalia", "Uganda", "Ethiopia", "Tanzania", "South\nSudan"),
  lon = c(41.95, 33.9, 39, 36, 34),
  lat = c(1.5, 2, 5.5, -4, 5.35)
)

country_labels_for_zoom <- data.frame(
  country = c("Somalia", "Uganda", "Ethiopia", "Tanzania", "South\nSudan"),
  lon = c(41.9, 33.9, 39, 34.25, 34),
  lat = c(1.5, 1.5, 5.5, -1.75, 5.5)
)

# Figure 1: Collection Sites

# malaria incidence collection site plot:
# to do: add vegetation or nataional park map overlap to cover missing incidence regions in natuional parks
# make dots proportiaonate to sample size at that site
# make dots different color from malaria incidence

# Load required libraries
library(terra)
library(ggplot2)
library(sf)
library(dplyr)

# to do
# make dot a different color than red (so its not confused with red map)
# dot size proportional to sample size (over all years)

# malaria atlas project incidence rate raster
# Define the path to TIFF files
tiff_dir <- "/Users/george/Downloads/clippedlayers/"

# Load all
# Define the path to your TIFF files
tiff_dir <- "/Users/george/Downloads/clippedlayers/"

# List all .tiff files in the directory
tiff_files <- list.files(tiff_dir, pattern = "\\.tiff$", full.names = TRUE)

# Load and stack all TIFF files, then calculate the average
rasters <- rast(tiff_files)
average_raster <- mean(rasters, na.rm = TRUE)  # Calculate the mean raster layer

# Convert the raster to a data frame for ggplot
average_raster_df <- as.data.frame(average_raster, xy = TRUE)
colnames(average_raster_df) <- c("x", "y", "malaria_rate")

# plot the map
collection_map <- ggplot() +
  # Ocean background
  geom_rect(aes(xmin = kenya_bbox["xmin"] - 5, xmax = kenya_bbox["xmax"] + 5,
                ymin = kenya_bbox["ymin"] - 5, ymax = kenya_bbox["ymax"] + 5),
            fill = "lightblue", color = NA) +
  
  # Light grey background for Kenya (added below context countries)
  geom_sf(data = context_countries, fill = "lightgrey", color = "black", size = 0.25) +
  geom_sf(data = sf_adm0, fill = "#A9B889", color = "black", size = 0.75) +  # <--- Added here
  
  # Malaria incidence raster layer (plotted on top of grey Kenya layer)
  geom_raster(data = average_raster_df, aes(x = x, y = y, fill = malaria_rate)) +
  scale_fill_gradient(low = "lightyellow", high = "red", name = "Malaria\nIncidence Rate") +
  
  # Other map layers (Admin boundaries, lakes, rivers, site points)
  geom_sf(data = sf_adm1, fill = NA, color = "white", size = 0.5, alpha = 0.6) +
  geom_sf(data = lake_victoria, fill = "lightblue", color = "lightblue", size = 0.5) +
  geom_sf(data = lake_turkana, fill = "lightblue", color = "lightblue", size = 0.5) +
  geom_sf(data = rivers_in_kenya, color = "lightblue", size = 0.5) +
  geom_point(data = sites_data, aes(x = longitude, y = latitude),  
             color = "white", size = 3.5, shape = 21, fill = palette_for_other_figures[4], stroke = 0.5, alpha = 0.8,
             position = position_jitter(width = 0.07, height = 0.07)) +
  # Coordinate limits for Kenya's bounding box
  coord_sf(xlim = c(kenya_bbox["xmin"] - 1, kenya_bbox["xmax"] + 1),
           ylim = c(kenya_bbox["ymin"] - 1, kenya_bbox["ymax"] + 1),
           expand = FALSE) +
  
  ## If want to ZOOM to western kenya (NOT RECOMMENDED DUE TO LOW RES INCIDENCE LAYER) Set coordinate limits to zoom in on Western Kenya
  #coord_sf(xlim = c(33.5, 37.5), ylim = c(-2, 2), expand = FALSE) +
  
  # Country labels
  geom_text(data = country_labels, aes(x = lon, y = lat, label = country), 
            size = 6, color = "black") +
  
  # Theme adjustments
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"), 
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    plot.title = element_blank(),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  labs(x = "Longitude", y = "Latitude", fill = "Malaria Incidence Rate")

# Save the plot
plot_title <- "1B_collection_map_kenya_malaria.svg"
ggsave(paste0("/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/", plot_title), 
       plot = collection_map, width = 14, height = 7, units = "in", dpi = 600)

# Save standardized data to /data folder
write.csv(sites_data, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig1B_collection_map_data.csv", row.names = FALSE)

# 1C)
################################################################################
#
# barplot of prevalence of key drug resistance haplotypes
# mdr1 NFD (Lmefantrine resistance) before and after ACT 
# dhfr-dhps IFN and triple.quintuple hap? dhps-437 and dhps-540 - SP resistance core allele before and after ACT
# crt resistance core haplotype (76T) before after ACT

################################################################################

# Convert era to a factor
# Load required packages
library(readr)
library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)

# Define mutation-gene mapping (modify if necessary)
mutation_gene_map <- c("Ala437Gly" = "DHPS",
                       "Lys540Glu" = "DHPS",
                       "Ser436His" = "DHPS",
                       "Lys76Thr"  = "CRT",
                       "Tyr184Phe" = "MDR1",
                       "Cys59Arg" = "DHFR")

ref_resistance_map <- c("Asn86" = "MDR1",
                        "Asp1246" = "MDR1")

# Define eras
SP_era <- 1998:2005
ACT_era <- 2006:2021

# Define mutations of interest
mutations_of_interest <- c("Ala437Gly", "Lys540Glu", "Lys76Thr", "Tyr184Phe", "Ser436His", "Cys59Arg")
reference_alleles_of_interest <- c("Asn86Tyr","Asp1246Tyr")
reference_allele_labels <- c("Asn86","Asp1246")

# Define colors
palette <- c("SP Period (1998-2005)" = "#00AFBB", "ACT Period (2006-2021)" = "#E7B800")

# Define the file path
file_path <- here("data", "kenya_all_DR2_prevalence_results_table_with_reference_resistant.csv")

# define results path
base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"

# Read in the prevalence table
all_prev_table <- read_csv(file_path)

# Print the first few rows to verify
print(head(all_prev_table))

# Remove Eastern kenya sites with too few samples spanning entire collection period
# Define sites to exclude
excluded_sites <- c("Kilifi", "Malindi")

# Filter out the unwanted rows
filtered_table <- all_prev_table[!(all_prev_table$site_uid %in% excluded_sites), ]

# --- Process mutations ---
filtered_mutations <- filtered_table %>%
  filter(mutation %in% mutations_of_interest) %>%
  mutate(era = case_when(
    as.numeric(Year_start) %in% SP_era ~ "SP Period (1998-2005)",
    as.numeric(Year_start) %in% ACT_era ~ "ACT Period (2006-2021)",
    Year_start == "PRE-ACT" ~ "SP Period (1998-2005)",  # treat PRE-ACT as SP
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(era)) %>%
  mutate(mutation_label = paste0(mutation_gene_map[mutation], "-", mutation))  # Append gene name

# --- Process reference-resistant alleles ---
filtered_reference <- filtered_table %>%
  filter(mutation %in% reference_alleles_of_interest) %>%
  mutate(era = case_when(
    as.numeric(Year_start) %in% SP_era ~ "SP Period (1998-2005)",
    as.numeric(Year_start) %in% ACT_era ~ "ACT Period (2006-2021)",
    Year_start == "PRE-ACT" ~ "SP Period (1998-2005)",  # treat PRE-ACT as SP
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(era)) %>%
  mutate(
    mutation_label = paste0(ref_resistance_map[reference_allele_labels[match(mutation, reference_alleles_of_interest)]], "-", reference_allele_labels[match(mutation, reference_alleles_of_interest)])
  )

# --- Compute mean prevalence for plotting ---
summary_mutations <- filtered_mutations %>%
  group_by(mutation, era, mutation_label) %>%
  summarise(mean_prevalence = mean(prevalence_percentage, na.rm = TRUE), .groups = "drop")

summary_reference <- filtered_reference %>%
  group_by(mutation, era, mutation_label) %>%
  summarise(mean_prevalence = mean(reference_resistant_prevalence, na.rm = TRUE), .groups = "drop")

summary_data <- bind_rows(summary_mutations, summary_reference)

# --- Wilcoxon Test on Sample-Level Data ---
wilcox_mutations <- filtered_mutations %>%
  group_by(mutation_label) %>%
  summarise(p_value = wilcox.test(prevalence_percentage ~ era, data = cur_data())$p.value) %>%
  ungroup()

wilcox_reference <- filtered_reference %>%
  group_by(mutation_label) %>%
  summarise(p_value = wilcox.test(reference_resistant_prevalence ~ era, data = cur_data())$p.value) %>%
  ungroup()

# Merge Wilcoxon results from both datasets
wilcox_results <- bind_rows(wilcox_mutations, wilcox_reference)

# Apply Holm correction for multiple comparisons
wilcox_results <- wilcox_results %>%
  mutate(p_value = p.adjust(p_value, method = "holm"))  # Adjusted p-values

# Define significance labels based on adjusted p-values
wilcox_results <- wilcox_results %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

# Print results
print(wilcox_results)

# Define y-position for significance brackets (higher than max mean prevalence)
bracket_positions <- summary_data %>%
  group_by(mutation_label) %>%
  summarise(max_y = max(mean_prevalence) + 5)  # Adds space above bars

# Merge significance results with y-positions
wilcox_results <- wilcox_results %>%
  left_join(bracket_positions, by = "mutation_label")

# Merge significance results with summary_data to ensure consistent grouping
summary_with_significance <- summary_data %>%
  left_join(wilcox_results, by = "mutation_label") %>%
  filter(significance != "ns")  # Hide non-significant brackets

# Define y-position for significance brackets (above highest prevalence)
summary_with_significance <- summary_with_significance %>%
  group_by(mutation_label) %>%
  mutate(max_y = max(mean_prevalence) + 5)  # Adjusted for better visibility

# Ensure SP bars appear first
summary_data <- summary_data %>%
  mutate(era = factor(era, levels = c("SP Period (1998-2005)", "ACT Period (2006-2021)")))  # Orders SP before ACT

# Define the desired order of mutation labels
mutation_order <- c("CRT-Lys76Thr", "DHFR-Cys59Arg", "DHPS-Ser436His", 
                    "DHPS-Ala437Gly", "DHPS-Lys540Glu", 
                    "MDR1-Asn86", "MDR1-Tyr184Phe", "MDR1-Asp1246")

# Convert mutation_label into a factor with the specified order
summary_data <- summary_data %>%
  mutate(mutation_label = factor(mutation_label, levels = mutation_order))

# define act and sp era colors from viridis palette so it matches other panels:
#viridis_colors <- viridis(2, option = "D", direction = 1)

# fix labels for x axis
summary_data <- summary_data %>%
  mutate(mutation_label = if_else(
    str_detect(mutation_label, "^[A-Za-z0-9]+-[A-Z][a-z]{2}[0-9]+$"),  # e.g., "MDR1-Asn86"
    paste0(mutation_label, str_extract(mutation_label, "[A-Z][a-z]{2}")),  # Append "Asn"
    mutation_label
  ))

summary_data <- summary_data %>%
  mutate(mutation_label = convert_aaa_to_a_with_gene_name(mutation_label))

# fix labels for bracket matching

summary_with_significance <- summary_with_significance %>%
  mutate(mutation_label = if_else(
    str_detect(mutation_label, "^[A-Za-z0-9]+-[A-Z][a-z]{2}[0-9]+$"),  # e.g., "MDR1-Asn86"
    paste0(mutation_label, str_extract(mutation_label, "[A-Z][a-z]{2}")),  # Append "Asn"
    mutation_label
  ))

summary_with_significance <- summary_with_significance %>%
  mutate(mutation_label = convert_aaa_to_a_with_gene_name(mutation_label))

# Generate the bar plot
p <- ggplot(summary_data, aes(x = mutation_label, y = mean_prevalence, fill = era)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.6) +
  theme_minimal() +
  scale_fill_manual(values = setNames(boxplot_colors, levels(summary_data$era))) +
  labs(title = "Mutation Prevalence in SP vs. ACT Eras",
       x = "Mutation",
       y = "Mean Prevalence (%)",
       fill = "Treatment Era") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting_fig1), 
    axis.text.y = element_text(size = axis.text.setting_fig1),
    axis.title.x = element_text(size = axis.title.setting_fig1),
    axis.title.y = element_text(size = axis.title.setting_fig1),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = legend.text.setting_fig1),
    legend.position = "bottom"
  )

# Manually add significance brackets & labels for significant results only
p <- p + geom_segment(data = summary_with_significance,
                      aes(x = mutation_label, xend = mutation_label, 
                          y = max_y, yend = max_y - 3),  # Bracket height
                      color = "black", size = 0.5) +
  geom_text(data = summary_with_significance,
            aes(x = mutation_label, y = max_y - 5, label = significance),  # Lowered asterisks
            vjust = -0.2, size = 8)  # Slightly adjusted vertical positioning

# write results to table for results section
write.table(summary_with_significance,sep=",",file="/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data/fig_1C_data_mutation_prevalence_SP_ACT.csv",row.names = F,col.names = T)

# Save standardized data and stats to /data folder
write.csv(summary_data, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig1C_mutation_prevalence_SP_ACT_data.csv", row.names = FALSE)
write.csv(wilcox_results, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig1C_mutation_prevalence_SP_ACT_stats.csv", row.names = FALSE)

# Save the plot
output_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
#ggsave(paste0(output_path,"mutation_prevalence_SP_ACT.svg"), plot = p, width = 8, height = 6, dpi = 300)

plot_title = "1C_mutation_prevalence_SP_ACT.svg"
ggsave(paste0(output_path, plot_title), 
       plot = p, width = 14, height = 8, units = "in", dpi = 600)

################################################################################
# 
# # 2A) mdr1 boxplot 
# 
################################################################################

library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)

# Load necessary libraries
library(dplyr)
library(readr)

gene_name="mdr1"
# Define the file path
file_path <- here("data", "kenya_all_DR2_prevalence_results_table_with_reference_resistant.csv")

# Read in the prevalence table
all_prev_table <- read_csv(file_path)

# Print the first few rows to verify
print(head(all_prev_table))

# Remove Eastern kenya sites with too few samples spanning entire collection period
# Define sites to exclude
excluded_sites <- c("Kilifi", "Malindi")

# Filter out the unwanted rows
filtered_table <- all_prev_table[!(all_prev_table$site_uid %in% excluded_sites), ]

# Define standardized time periods
define_treatment_periods <- function(data) {
  data %>%
    mutate(
      time_period = case_when(
        Year_start >= 1998 & Year_start <= 2005 ~ "SP Period (1998-2005)",
        Year_start >= 2006 & Year_start <= 2021 ~ "ACT Period (2006-2021)",
        Year_start == "PRE-ACT" ~ "SP Period (1998-2005)",  # treat PRE-ACT as SP
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(time_period)) %>%
    mutate(time_period = factor(time_period, levels = c("SP Period (1998-2005)", "ACT Period (2006-2021)")))
}

mdr1_alleles <- c("Asp1246Tyr", "Asn86Tyr", "Tyr184Phe")
mdr1_ref_resistant_alleles <- c("Asn86Tyr", "Asp1246Tyr")

gene_name_to_plot="mdr1"
alleles=mdr1_alleles
data_table = filtered_table
plot_title = "2A_MDR1"

# Filter relevant data
gene_data <- data_table %>%
  filter(gene_name == gene_name_to_plot & mutation %in% alleles) %>%
  define_treatment_periods()

# Create prevalence_plot based on whether mutation is reference-resistant or variant-resistant
# gene_data <- gene_data %>%
#   mutate(
#     prevalence_plot = ifelse(
#       mutation %in% mdr1_ref_resistant_alleles,
#       100 - prevalence_percentage,  # Flip prevalence for reference-resistant mutations
#       prevalence_percentage         # Keep original for variant-resistant mutations
#     )
#   )

gene_data <- gene_data %>%
  mutate(
    prevalence_plot = ifelse(
      mutation %in% mdr1_ref_resistant_alleles,
      reference_resistant_prevalence,  # Use reference_resistant_prevalence for reference-resistant mutations
      prevalence_percentage            # Keep original for variant-resistant mutations
    )
  )

# Order mutations by amino acid position and rename reference-resistant mutations
unique_mutations <- gene_data %>%
  mutate(
    # Extract numeric position from mutation name
    aa_position = as.numeric(gsub("[^0-9]", "", mutation)),
    # Remove the trailing three-letter amino acid code only for reference-resistant mutations
    mutation_label = ifelse(mutation %in% mdr1_ref_resistant_alleles, 
                            gsub("(?<=[0-9])[A-Za-z]{3}$", "", mutation, perl = TRUE), 
                            mutation)
  ) %>%
  distinct(mutation, aa_position, mutation_label) %>%
  arrange(aa_position) %>%
  pull(mutation_label)

# Apply the ordered factor levels to mutation_label
gene_data <- gene_data %>%
  mutate(
    # Extract amino acid position for ordering
    aa_position = as.numeric(gsub("[^0-9]", "", mutation)),
    # Remove the second amino acid entirely if mutation is reference-resistant
    mutation_label = ifelse(mutation %in% mdr1_ref_resistant_alleles, 
                            gsub("(?<=[0-9])[A-Za-z]{3}$", "", mutation, perl = TRUE), 
                            mutation),
    # Convert mutation_label to a factor with levels ordered by AA position
    mutation_label = factor(mutation_label, levels = unique_mutations)
  )

# Perform pairwise Wilcoxon tests for each mutation_label
gene_stat_results <- gene_data %>%
  group_by(mutation_label) %>%
  wilcox_test(prevalence_plot ~ time_period) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  filter(p.adj < 0.05) %>%  # Keep only significant results
  add_xy_position(x = "mutation_label", dodge = 0.8)

# Adjust y-position for significant p-values based on number of comparisons
gene_stat_results <- gene_stat_results %>%
  group_by(mutation_label) %>%
  mutate(y.position = 105 + row_number() * 5) %>%  # Start above 100% and space by 5
  ungroup()

print(gene_data)

gene_data <- gene_data %>%
  mutate(mutation_label = if_else(
    str_detect(mutation_label, "^[A-Z][a-z]{2}[0-9]+$"),  # Matches 'Asn86' pattern
    paste0(mutation_label, str_sub(mutation_label, 1, 3)), # Append first 3 letters again
    mutation_label                                        # Else keep as is
  ))

gene_data <- gene_data %>%
  mutate(mutation_label = convert_aaa_to_a(mutation_label))

# Get unique mutation labels ordered by position
ordered_levels <- gene_data %>%
  distinct(mutation_label, aa_position) %>%
  arrange(aa_position) %>%
  pull(mutation_label)

# Apply the levels
gene_data <- gene_data %>%
  mutate(mutation_label = factor(mutation_label, levels = ordered_levels))

# define act and sp era colors from viridis palette so it matches other panels:
# viridis_colors <- viridis(2, option = "D", direction = 1)

# Create the box plot using `mutation_label` as the x-axis
boxplot_gene <- ggboxplot(
  gene_data, x = "mutation_label", y = "prevalence_plot",
  fill = "time_period") +
  scale_fill_manual(values = boxplot_colors) +  # Use reserved scico colors
  labs(
    title = "MDR1",
    x = "Mutation",
    y = "Prevalence (%)",
    fill = "Treatment Era"
  ) +
  theme_minimal() +
  # scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
    axis.text.y = element_text(size = axis.text.setting),
    axis.title = element_text(size = axis.title.setting),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = legend.text.setting),
    legend.position = "bottom",
    legend.key.size = unit(1.5, "cm") 
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),  # Set y-axis breaks every 20%
    expand = expansion(mult = c(0, 0.1)),  # Add headspace for brackets
    labels = function(x) ifelse(x > 100, "", x)  # Blank labels above 100%
  )

gene_stat_results <- gene_stat_results %>%
  mutate(
    significance_label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Add significant p-values
boxplot_with_significant_pvalues <- boxplot_gene + 
  stat_pvalue_manual(gene_stat_results, label = "significance_label", tip.length = 0.01, bracket.size = 1, size = 11)

print(gene_data)

# save text results for manuscript results text
gene_summary <- gene_data %>%
  group_by(mutation_label, time_period) %>%
  summarise(
    prevalence_plot = mean(prevalence_plot, na.rm = TRUE),
    total_samples = sum(total_samples, na.rm = TRUE)
  ) %>%
  ungroup()
print(gene_summary)

final_data <- gene_summary %>%
  left_join(gene_stat_results %>% select(mutation_label, p.adj), by = "mutation_label")

# Save plots and results
save_results(plot_title, boxplot_with_significant_pvalues, gene_data, gene_stat_results,"2A_")

# save numbers for manuscript results section
base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
write.csv(final_data, paste0(base_path, gene_name, "_results_for_manuscript.csv"), row.names = FALSE)

# Save standardized data and stats to /data folder
write.csv(gene_data, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig2A_MDR1_boxplot_data.csv", row.names = FALSE)
write.csv(gene_stat_results, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig2A_MDR1_boxplot_stats.csv", row.names = FALSE)

################################################################################
#
# # 2B) mdr1 stacked barplot
#
################################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
#library(viridis)


mdr1_alleles <- c("mdr1-Asn86Tyr","mdr1-Tyr184Phe","mdr1-Ser1034Cys","mdr1-Asn1042Asp","mdr1-Asp1246Tyr")
reference_resistant <- c("mdr1-Asn86Tyr","mdr1-Tyr184Phe", "mdr1-Asp1246Tyr")
reference_resistant_labels <- c("Asn86Asn","Tyr184Tyr", "Asp1246Asp")

gene_to_plot <- "mdr1"  # Example gene (can contain numbers or dashes)
treatment_era_to_plot <- "All_Eras"
output_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
mutations_of_interest <- mdr1_alleles

# Convert gene name to uppercase for plot title
gene_title <- toupper(gene_to_plot)

# Define time bins (1998-2021)
time_bins <- c(1998, 2002, 2006, 2010, 2014, 2018, 2022)
time_bin_labels <- c("1998-2001", "2002-2005", "2006-2009", "2010-2013", "2014-2017", "2018-2021")

################################################################################
# Step 1: Filter & Format the Genotype Table
################################################################################

filtered_genotype <- genotype_file %>%
  filter(mutation_name %in% mutations_of_interest) %>%
  select(sample, mutation_name, genotype)

wide_genotype_table <- filtered_genotype %>%
  pivot_wider(names_from = mutation_name, values_from = genotype)

# valid_samples <- wide_genotype_table %>%
#   rowwise() %>%
#   filter(all(c_across(where(is.numeric)) != 1)) %>%
#   ungroup() %>%
#   pull(sample)
valid_samples <- wide_genotype_table %>%
  rowwise() %>%
  filter(all(!c_across(where(is.numeric)) %in% c(1, -1))) %>%
  ungroup() %>%
  pull(sample)

valid_samples_clean <- gsub("S-", "S", valid_samples)

variant_table_homozygous <- variant_table_cov_filtered_with_met %>%
  filter(Sampleid %in% valid_samples_clean)


# check out 03MDH37-HAKpf-1
#variant_table_homozygous[which(variant_table_homozygous$Sampleid=="03MDH37-HAKpf-1" & variant_table_homozygous$mutation == "Asp1246Tyr"),]
# has alt_ref 6 but homozygous alt genotype for 1246 (2) # question - is homo/hetero dtermined by percent allele frequency or string presence of at least 1 alt position (probably within sample freq)
# use a minimum within sample AF to determine haplotype prevalence
# wide_genotype_table[wide_genotype_table$sample=="03MDH37-HAKpf-1",]
# just filter out anything with alt count over 1 for the reference af using count of 3 filter)
################################################################################
# Step 2: Filter & Format the Variant Table, Include Reference-Resistant Alleles
################################################################################

# Include only mutations where umi_alt >= 1 and umi_ref < 1
variant_mutations <- variant_table_homozygous %>%
  filter(mutation_name %in% mutations_of_interest & umi_alt >= 1 & umi_ref < 1) %>%
  mutate(present = 1) %>%
  select(Sampleid, mutation_name, present) %>%
  pivot_wider(names_from = mutation_name, values_from = present, values_fill = list(present = 0))

# Initialize the upset table with variant mutations
upset_table <- variant_mutations

# Add reference-resistant alleles where umi_ref >= 1 and umi_alt < 1
for (i in seq_along(reference_resistant)) {
  allele <- reference_resistant[i]
  label <- reference_resistant_labels[i]
  
  reference_present <- variant_table_homozygous %>%
    filter(mutation_name == allele & umi_ref >= 1 & umi_alt < 1) %>%
    select(Sampleid) %>%
    mutate(!!label := 1)  # Create new column with reference-resistant label
  
  upset_table <- upset_table %>%
    left_join(reference_present, by = "Sampleid") %>%
    mutate(!!label := ifelse(is.na(.data[[label]]), 0, .data[[label]]))  # Fill missing values with 0
}

library(dplyr)
library(stringr)

upset_table <- upset_table %>%
  rename_with(~ str_remove(., "^mdr1-"), .cols = starts_with("mdr1-"))

# Keep only mutations with at least one occurrence
columns_to_keep <- c(1, which(colSums(upset_table[, -1]) > 1) + 1)
upset_table_df_filtered <- upset_table[, columns_to_keep]

sample_year_data <- variant_table_homozygous %>%
  select(Sampleid, Year_start) %>%
  distinct(Sampleid, .keep_all = TRUE) %>%
  mutate(Year_start = as.numeric(Year_start))

upset_table_df_filtered <- upset_table_df_filtered %>%
  left_join(sample_year_data, by = "Sampleid") %>%
  mutate(year_bin = cut(Year_start, breaks = time_bins, include.lowest = TRUE, right = FALSE, labels = time_bin_labels))

# convert to amino acid single letter code
# First, apply convert_aaa_to_a to colnames
new_colnames <- colnames(upset_table_df_filtered)[-c(1,length(colnames(upset_table_df_filtered))-1,length(colnames(upset_table_df_filtered)))]
new_colnames <- sapply(new_colnames, convert_aaa_to_a, USE.NAMES = FALSE)

# Assign new column names back to upset_table
colnames(upset_table_df_filtered)[2:(length(colnames(upset_table_df_filtered)) - 2)] <- new_colnames

################################################################################
# Step 3: Compute Haplotype Counts & Prevalence
################################################################################

# remove any samples with low frequency alt reads (resulting in missing call for the ref resistant allele using ref_cov =0)
# Recalculate row sums for allele presence across columns 2 to 7
row_sums <- rowSums(upset_table_df_filtered[, 2:7], na.rm = TRUE)
# Find rows where the sum is not 3
incomplete_haplotype_rows_to_remove <- which(row_sums != 3)

upset_table_df_filtered <- upset_table_df_filtered[-incomplete_haplotype_rows_to_remove, ]

haplotype_counts <- upset_table_df_filtered %>%
  select(-Sampleid, -Year_start) %>%
  group_by(year_bin, across(everything())) %>%
  summarise(count = n(), .groups = "drop")

total_samples_per_bin <- upset_table_df_filtered %>%
  group_by(year_bin) %>%
  summarise(total_samples = n_distinct(Sampleid), .groups = "drop")

haplotype_counts <- haplotype_counts %>%
  left_join(total_samples_per_bin, by = "year_bin") %>%
  mutate(prevalence = ifelse(is.na(count), 0, (count / total_samples) * 100))

haplotype_counts <- full_join(tibble(year_bin = time_bin_labels), haplotype_counts, by = "year_bin") %>%
  replace_na(list(count = 0, prevalence = 0, total_samples = 0))

################################################################################
# Step 4: Process Haplotype Labels
################################################################################

# haplotype_counts <- haplotype_counts %>%
#   mutate(haplotype_label = apply(select(., -count, -prevalence, -year_bin, -total_samples), 1,
#                                  function(x) paste(names(x)[which(x == 1)], collapse = "-")))

haplotype_counts <- haplotype_counts %>%
  mutate(haplotype_label = apply(select(., -count, -prevalence, -year_bin, -total_samples), 1, function(x) {
    present_mutations <- names(x)[which(x == 1)]
    
    # Extract numeric positions from names like N86Y, Y184F, D1246Y
    residue_positions <- as.numeric(str_extract(present_mutations, "\\d+"))
    
    # Sort by position
    sorted_mutations <- present_mutations[order(residue_positions)]
    
    paste(sorted_mutations, collapse = "-")
  }))

plot_data <- haplotype_counts %>%
  select(year_bin, haplotype_label, prevalence, count, total_samples)

plot_data$haplotype_label[plot_data$haplotype_label == ""] <- "WT"

transform_haplotype_labels_simple <- function(haplotype_labels) {
  sapply(haplotype_labels, function(label) {
    if (label == "WT") return("WT")
    
    # Split the haplotype string into individual mutations
    mutations <- unlist(strsplit(label, "-"))
    
    # Extract the final amino acid letter from each mutation
    aa_letters <- substr(mutations, nchar(mutations), nchar(mutations))
    
    # Combine original label with trailing letters in parentheses
    paste0(label, " (", paste(aa_letters, collapse = ""), ")")
  }, USE.NAMES = FALSE)
}


plot_data <- plot_data %>%
  mutate(haplotype_label = transform_haplotype_labels_simple(haplotype_label)) %>%
  filter(!is.na(year_bin))

################################################################################
# Step 5: Plot Stacked Bar Chart for Haplotype Prevalences with Increased Font Sizes
################################################################################

hap_prev_plot <- ggplot(plot_data, aes(x = year_bin, y = prevalence, fill = haplotype_label)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.95) +
  scale_fill_manual(values = palette_for_other_figures) +  # Use reserved scico colors
  labs(x = "Time Period", y = "Haplotype Prevalence (%)", fill = "Haplotype") +
  ggtitle(gene_title) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
    axis.text.y = element_text(size = axis.text.setting),
    axis.title = element_text(size = axis.title.setting),
    plot.title = element_blank(),
    legend.title = element_text(size = legend.title.setting),
    legend.text = element_text(size = legend.text.setting)
  )

################################################################################
# Step 6: Save the Plot
################################################################################

plot_title <- paste0("2B_",gene_to_plot, "_", treatment_era_to_plot, "_hap_prev_stacked_barplot_time_bins.svg")
ggsave(file.path(output_path, plot_title),
       plot = hap_prev_plot, width = 16, height = 10, units = "in", dpi = 600)

# save the data for manuscript results
base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data/"
write.csv(plot_data, paste0(base_path, gene_to_plot, "_stacked_barplot_results_for_manuscript.csv"), row.names = FALSE)

# Save standardized data to /data folder  
write.csv(plot_data, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig2B_MDR1_haplotype_stacked_data.csv", row.names = FALSE)

print(paste("Plot saved to:", file.path(output_path, plot_title)))

################################################################################
#
# 2C) mdr1 lineplot
#
################################################################################

# Load necessary libraries
library(readr)
library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)

# Define the file path
file_path <- here("data", "kenya_all_DR2_prevalence_results_table_with_reference_resistant.csv")

# define results path
base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"

# Read in the prevalence table
all_prev_table <- read_csv(file_path)

# Print the first few rows to verify
print(head(all_prev_table))

# Remove Eastern kenya sites with too few samples spanning entire collection period
# Define sites to exclude
excluded_sites <- c("Kilifi", "Malindi")

# Filter out the unwanted rows
filtered_table <- all_prev_table[!(all_prev_table$site_uid %in% excluded_sites), ]

key_mutations_aaa <- c("Asp1246Tyr", "Asn86Tyr", "Tyr184Phe")
mdr1_reference_resistant_alleles <- c("Asn86Tyr", "Asp1246Tyr")
mdr1_variant_resistant_alleles <- c("Tyr184Phe")

mdr1_data <- filtered_table %>%
  filter(gene_name == "mdr1" & mutation %in% c(mdr1_reference_resistant_alleles, mdr1_variant_resistant_alleles)) %>%
  mutate(
    prevalence_custom = case_when(
      mutation %in% mdr1_reference_resistant_alleles ~
        ((samples_with_adequate_coverage - samples_with_mutation) / samples_with_adequate_coverage) * 100,
      mutation %in% mdr1_variant_resistant_alleles ~ prevalence_percentage
    ),
    # Create display labels, removing "Tyr" from reference alleles
    display_label = case_when(
      mutation %in% mdr1_reference_resistant_alleles ~ sub("Tyr$", "", mutation),  # Remove "Tyr" from end
      TRUE ~ mutation  # Keep full name for variant alleles
    )
  ) %>%
  filter(!is.na(prevalence_custom))  # Remove rows with missing prevalence

mdr1_data <- mdr1_data %>%
  filter(Year_start != "PRE-ACT" & Year_start != "POST-ACT") %>%  # Exclude non-numeric years
  mutate(Year_start = as.numeric(Year_start),
         Year_group = cut(Year_start, breaks = seq(1998, 2022, by = 2),
                          labels = paste(seq(1998, 2020, by = 2), seq(1999, 2021, by = 2), sep = "-"),
                          right = FALSE)) %>%  # Right = FALSE for left-inclusive intervals
  group_by(display_label, Year_group) %>%
  summarize(avg_prevalence = mean(prevalence_custom, na.rm = TRUE)) %>%
  ungroup()

library(dplyr)
library(stringr)

mdr1_data <- mdr1_data %>%
  mutate(display_label = if_else(
    str_detect(display_label, "^[A-Z][a-z]{2}[0-9]+$"),  # Matches 'Asn86' pattern or anything with only the ref codon listed
    paste0(display_label, str_sub(display_label, 1, 3)), # Append first 3 letters again
    display_label                                        # Else keep as is
  ))

mdr1_data_single_a <- mdr1_data %>%
  mutate(display_label = convert_aaa_to_a(display_label))  # Extract numbers

mdr1_data_single_a <- mdr1_data_single_a %>%
  mutate(residue_number = as.numeric(gsub("[A-Za-z]", "", display_label)))  # Extract numbers

# # Plot for mdr1 with modified labels in the legend
mdr1_plot <- ggplot(mdr1_data_single_a, aes(x = Year_group, y = avg_prevalence, color = reorder(display_label, residue_number), group = display_label)) +
  geom_line(size = 3) +
  geom_point(size = 5, alpha = 0.7) +
  labs(title = "MDR1", x = "Year (2-Year Intervals)", y = "Average Prevalence (%)", color = "Allele") +
  theme_minimal() +
  #scale_color_brewer(palette = "Paired") +  # Use "Paired" palette for more colors
  scale_color_manual(values = palette_for_other_figures) +  # Use reserved scico colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
    axis.text.y = element_text(size = axis.text.setting),
    axis.title = element_text(size = axis.title.setting),
    plot.title = element_blank(),
    legend.title = element_text(size = legend.title.setting),
    legend.text = element_text(size = legend.text.setting)
  )

# write results to table for results section
write.table(mdr1_data,sep=",",file="/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data/fig_2A_mdr1_lineplot_data.csv",row.names = F,col.names = T)

# Save standardized data to /data folder (corrected figure number)
write.csv(mdr1_data_single_a, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig2C_MDR1_lineplot_data.csv", row.names = FALSE)

# Save the plot
plot_title = "2C_lineplot_mdr1.svg"
ggsave(paste0(base_path, plot_title), 
       plot = mdr1_plot, width = 33, height = 10, units = "in", dpi = 600)

################################################################################
#
# Figures 3-5) CRT, DHFR, DHPS Figures
#
################################################################################
################################################################################
#
# A) boxplot
#
################################################################################

library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)

# Load necessary libraries
library(dplyr)
library(readr)

gene_name="mdr1"
# Define the file path
file_path <- here("data", "kenya_all_DR2_prevalence_results_table_with_reference_resistant.csv")

# Read in the prevalence table
all_prev_table <- read_csv(file_path)

# Print the first few rows to verify
print(head(all_prev_table))

# Remove Eastern kenya sites with too few samples spanning entire collection period
# Define sites to exclude
excluded_sites <- c("Kilifi", "Malindi")

# Filter out the unwanted rows
filtered_table <- all_prev_table[!(all_prev_table$site_uid %in% excluded_sites), ]

# Define standardized time periods
define_treatment_periods <- function(data) {
  data %>%
    mutate(
      time_period = case_when(
        Year_start >= 1998 & Year_start <= 2005 ~ "SP Period (1998-2005)",
        Year_start >= 2006 & Year_start <= 2021 ~ "ACT Period (2006-2021)",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(time_period)) %>%
    mutate(time_period = factor(time_period, levels = c("SP Period (1998-2005)", "ACT Period (2006-2021)")))
}



# Generalized function for gene mutation analysis
analyze_gene_mutations <- function(gene_name, alleles, data_table, plot_title,fig_num) {
  
  # Filter relevant data
  gene_data <- data_table %>%
    filter(gene_name == gene_name & mutation %in% alleles) %>%
    define_treatment_periods()
  
  if (gene_name == "mdr1"){
    
    mdr1_ref_resistant_alleles <- c("Asn86Tyr", "Asp1246Tyr")
    
    # Create prevalence_plot based on whether mutation is reference-resistant or variant-resistant
    gene_data <- gene_data %>%
      mutate(
        prevalence_plot = ifelse(
          mutation %in% mdr1_ref_resistant_alleles,
          100 - prevalence_percentage,  # Flip prevalence for reference-resistant mutations
          prevalence_percentage         # Keep original for variant-resistant mutations
        )
      )
    
    # Order mutations by amino acid position and rename reference-resistant mutations
    unique_mutations <- gene_data %>%
      mutate(
        # Extract numeric position from mutation name
        aa_position = as.numeric(gsub("[^0-9]", "", mutation)),
        # Remove the trailing three-letter amino acid code only for reference-resistant mutations
        mutation_label = ifelse(mutation %in% mdr1_ref_resistant_alleles, 
                                gsub("(?<=[0-9])[A-Za-z]{3}$", "", mutation, perl = TRUE), 
                                mutation)
      ) %>%
      distinct(mutation, aa_position, mutation_label) %>%
      arrange(aa_position) %>%
      pull(mutation_label)
    
    # Apply the ordered factor levels to mutation_label
    gene_data <- gene_data %>%
      mutate(
        # Extract amino acid position for ordering
        aa_position = as.numeric(gsub("[^0-9]", "", mutation)),
        # Remove the second amino acid entirely if mutation is reference-resistant
        mutation_label = ifelse(mutation %in% mdr1_ref_resistant_alleles, 
                                gsub("(?<=[0-9])[A-Za-z]{3}$", "", mutation, perl = TRUE), 
                                mutation),
        # Convert mutation_label to a factor with levels ordered by AA position
        mutation_label = factor(mutation_label, levels = unique_mutations)
      )
    
    # Perform pairwise Wilcoxon tests for each mutation_label
    gene_stat_results <- gene_data %>%
      group_by(mutation_label) %>%
      wilcox_test(prevalence_plot ~ time_period) %>%
      adjust_pvalue(method = "holm") %>%
      add_significance("p.adj") %>%
      filter(p.adj < 0.05) %>%  # Keep only significant results
      add_xy_position(x = "mutation_label", dodge = 0.8)
    
    # Adjust y-position for significant p-values based on number of comparisons
    gene_stat_results <- gene_stat_results %>%
      group_by(mutation_label) %>%
      mutate(y.position = 105 + row_number() * 5) %>%  # Start above 100% and space by 5
      ungroup()
    
    gene_data <- gene_data %>%
      mutate(mutation_label = if_else(
        str_detect(mutation_label, "^[A-Z][a-z]{2}[0-9]+$"),  # Matches 'Asn86' pattern
        paste0(mutation_label, str_sub(mutation_label, 1, 3)), # Append first 3 letters again
        mutation_label                                        # Else keep as is
      ))
    
    gene_data <- gene_data %>%
      mutate(mutation_label = convert_aaa_to_a(mutation_label))
    
    # Get unique mutation labels ordered by position
    ordered_levels <- gene_data %>%
      distinct(mutation_label, aa_position) %>%
      arrange(aa_position) %>%
      pull(mutation_label)
    
    # Apply the levels
    gene_data <- gene_data %>%
      mutate(mutation_label = factor(mutation_label, levels = ordered_levels))
    
    # Create the box plot using `mutation_label` as the x-axis
    boxplot_gene <- ggboxplot(
      gene_data, x = "mutation_label", y = "prevalence_plot",
      fill = "time_period") +
      labs(
        title = "MDR1",
        x = "Mutation",
        y = "Prevalence (%)",
        fill = "Treatment Era"
      ) +
      scale_fill_manual(values = boxplot_colors) +  # Use reserved scico colors
      theme_minimal() +
      # scale_fill_brewer(palette = "Set2") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
        axis.text.y = element_text(size = axis.text.setting),
        axis.title = element_text(size = axis.title.setting),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = legend.text.setting),
        legend.position = "bottom",
        legend.key.size = unit(1.5, "cm") 
      ) +
      scale_y_continuous(
        breaks = seq(0, 100, by = 20),  # Set y-axis breaks every 20%
        expand = expansion(mult = c(0, 0.1)),  # Add headspace for brackets
        labels = function(x) ifelse(x > 100, "", x)  # Blank labels above 100%
      )
    
    gene_stat_results <- gene_stat_results %>%
      mutate(
        significance_label = case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01 ~ "**",
          p.adj < 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    
    # Add significant p-values
    boxplot_with_significant_pvalues <- boxplot_gene + 
      stat_pvalue_manual(gene_stat_results, label = "significance_label", tip.length = 0.01, bracket.size = 1, size = 11)
    
    print(gene_data)
    
    # save text results for manuscript results text
    gene_summary <- gene_data %>%
      group_by(mutation_label, time_period) %>%
      summarise(
        prevalence_plot = mean(prevalence_plot, na.rm = TRUE),
        total_samples = sum(total_samples, na.rm = TRUE)
      ) %>%
      ungroup()
    print(gene_summary)
    
    
    final_data <- gene_summary %>%
      left_join(gene_stat_results %>% select(mutation_label, p.adj), by = "mutation_label")
    
  } else {
    
    # Order mutations by AA position
    unique_mutations <- gene_data %>%
      mutate(aa_position = as.numeric(gsub("[^0-9]", "", mutation))) %>%
      distinct(mutation, aa_position) %>%
      arrange(aa_position) %>%
      pull(mutation)
    
    gene_data <- gene_data %>%
      mutate(
        aa_position = as.numeric(gsub("[^0-9]", "", mutation)),
        mutation = factor(mutation, levels = unique_mutations)
      )
    
    # Perform Wilcoxon test for each mutation
    gene_stat_results <- gene_data %>%
      group_by(mutation) %>%
      wilcox_test(prevalence_percentage ~ time_period) %>%
      adjust_pvalue(method = "holm") %>%
      add_significance("p.adj") %>%
      filter(p.adj < 0.05) %>%
      add_xy_position(x = "mutation", dodge = 0.8) %>%
      group_by(mutation) %>%
      mutate(y.position = 105 + row_number() * 5) %>%
      ungroup()
    
    gene_data <- gene_data %>%
      mutate(mutation = if_else(
        str_detect(mutation, "^[A-Z][a-z]{2}[0-9]+$"),  # Matches 'Asn86' pattern
        paste0(mutation, str_sub(mutation, 1, 3)), # Append first 3 letters again
        mutation                                        # Else keep as is
      ))
    
    gene_data <- gene_data %>%
      mutate(mutation = convert_aaa_to_a(mutation))
    
    # Get unique mutation labels ordered by position
    ordered_levels <- gene_data %>%
      distinct(mutation, aa_position) %>%
      arrange(aa_position) %>%
      pull(mutation)
    
    # Apply the levels
    gene_data <- gene_data %>%
      mutate(mutation = factor(mutation, levels = ordered_levels))
    
    # define act and sp era colors from viridis palette so it matches other panels:
   # viridis_colors <- viridis(2, option = "D", direction = 1)
    
    
    # Create boxplot
    boxplot_gene <- ggboxplot(
      gene_data, x = "mutation", y = "prevalence_percentage",
      fill = "time_period"
    ) +
      labs(title = plot_title, x = "Mutation", y = "Prevalence (%)", fill = "Treatment Era") +
      theme_minimal() +
      scale_fill_manual(values = boxplot_colors) +  # Use reserved scico colors
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
        axis.text.y = element_text(size = axis.text.setting),
        axis.title = element_text(size = axis.title.setting),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = legend.text.setting),
        legend.position = "bottom",
        legend.key.size = unit(1.5, "cm") 
      ) +
      scale_y_continuous(breaks = seq(0, 100, by = 20), expand = expansion(mult = c(0, 0.1)))
    
    gene_stat_results <- gene_stat_results %>%
      mutate(
        significance_label = case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01 ~ "**",
          p.adj < 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    
    # Add significant p-values
    boxplot_with_significant_pvalues <- boxplot_gene + 
      stat_pvalue_manual(gene_stat_results, label = "significance_label", tip.length = 0.01, bracket.size = 1, size = 11)
    
    # save text results for manuscript results text
    gene_summary <- gene_data %>%
      group_by(mutation, time_period) %>%
      summarise(
        prevalence_plot = mean(prevalence_percentage, na.rm = TRUE),
        total_samples = sum(total_samples, na.rm = TRUE)
      ) %>%
      ungroup()
    
    final_data <- gene_summary %>%
      left_join(gene_stat_results %>% select(mutation, p.adj), by = "mutation")
    
  }
  
  # Save plots and results
  save_results(plot_title, boxplot_with_significant_pvalues, gene_data, gene_stat_results,fig_num)
  
  # save numbers for manuscript results section
  base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
  write.csv(final_data, paste0(base_path,gene_name, "_results_for_manuscript.csv"), row.names = FALSE)
  
  # Save standardized data and stats to /data folder
  standardized_data_filename <- paste0("/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/", 
                                      fig_num, toupper(gene_name), "_boxplot_data.csv")
  standardized_stats_filename <- paste0("/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/", 
                                       fig_num, toupper(gene_name), "_boxplot_stats.csv")
  write.csv(gene_data, standardized_data_filename, row.names = FALSE)
  
  # Convert list columns in gene_stat_results to character before saving
  gene_stat_results_clean <- gene_stat_results
  if (any(sapply(gene_stat_results_clean, is.list))) {
    gene_stat_results_clean <- gene_stat_results_clean %>%
      mutate(across(where(is.list), ~ sapply(., toString)))
  }
  write.csv(gene_stat_results_clean, standardized_stats_filename, row.names = FALSE)
  
  return(list(gene_data = gene_data, gene_stat_results = gene_stat_results, plot = boxplot_with_significant_pvalues))
}

dhps_alleles <- c("Ser436His", "Lys540Glu", "Ala437Gly", "Ala581Gly", "Ser436Ala", "Ala613Ser")
dhfr_alleles <- c("Asn51Ile", "Cys59Arg", "Ser108Asn", "Ile164Leu")
crt_alleles <- c("Cys72Ser","Val73Leu", "Met74Ile","crt-Asn75Glu", "Lys76Thr", "Ala220Ser", "Arg371Ile", "Asn326Ser", "Gln271Glu", "Asp24Tyr")

dhps_results <- analyze_gene_mutations("dhps", dhps_alleles, filtered_table, "DHPS","5A_")
dhfr_results <- analyze_gene_mutations("dhfr-ts", dhfr_alleles, filtered_table, "DHFR","4A_")
crt_results <- analyze_gene_mutations("crt", crt_alleles, filtered_table, "CRT","3A_")


################################################################################
#
# B) Stacked barplot
#
################################################################################

dhps_alleles <- c("dhps-Ile431Val", "dhps-Ser436Ala", 
                  "dhps-Ser436His", "dhps-Ala437Gly", "dhps-Lys540Glu", "dhps-Ala581Gly", 
                  "dhps-Ala613Thr", "dhps-Ala613Ser")

dhfr_alleles <- c("dhfr-ts-Ala16Val", "dhfr-ts-Asn51Ile", "dhfr-ts-Cys59Arg", "dhfr-ts-Ser108Asn", 
                  "dhfr-ts-Ser108Thr", "dhfr-ts-Ile164Leu")

crt_alleles <- c("Val73Leu", "Cys72Ser", "Met74Ile","crt-Asn75Glu", "crt-Lys76Thr", 
                 "crt-Thr93Ser", "crt-His97Tyr", "crt-His97Leu", "crt-Cys101Phe", "crt-Phe145Ile", 
                 "crt-Ile218Phe", "crt-Ala220Ser", "crt-Gln271Glu", "crt-Asn326Ser", "crt-Met343Leu", 
                 "crt-Gly353Val", "crt-Ile356Thr", "crt-Arg371Ile")


# Load necessary libraries
library(dplyr)
library(stringr)

generate_haplotype_prevalence_plot <- function(gene_to_plot, treatment_era_to_plot, output_path, 
                                               genotype_file, variant_table_cov_filtered_with_met, 
                                               mutations_of_interest,fig_num) {
  
  # Convert gene name to uppercase for plot title
  gene_title <- toupper(gene_to_plot)
  
  # Define time bins (1998-2021)
  time_bins <- c(1998, 2002, 2006, 2010, 2014, 2018, 2022)
  time_bin_labels <- c("1998-2001", "2002-2005", "2006-2009", "2010-2013", "2014-2017", "2018-2021")
  
  ################################################################################
  # Step 1: Filter & Format the Genotype Table
  ################################################################################
  
  filtered_genotype <- genotype_file %>%
    filter(mutation_name %in% mutations_of_interest) %>%
    select(sample, mutation_name, genotype)
  
  wide_genotype_table <- filtered_genotype %>%
    pivot_wider(names_from = mutation_name, values_from = genotype)
  
  valid_samples <- wide_genotype_table %>%
    rowwise() %>%
    filter(all(!c_across(where(is.numeric)) %in% c(1, -1))) %>%
    ungroup() %>%
    pull(sample)
  
  valid_samples_clean <- gsub("S-", "S", valid_samples)
  
  variant_table_homozygous <- variant_table_cov_filtered_with_met %>%
    filter(Sampleid %in% valid_samples_clean)
  
  ################################################################################
  # Step 2: Filter & Format the Variant Table
  ################################################################################
  
  upset_table <- variant_table_homozygous %>%
    filter(mutation_name %in% mutations_of_interest) %>%
    mutate(present = ifelse(umi_alt > 1, 1, 0)) %>%
    select(Sampleid, mutation_name, present) %>%
    pivot_wider(names_from = mutation_name, values_from = present, values_fill = list(present = 0))
  
  upset_table <- upset_table %>%
    rename_with(~ str_replace(., ".*-", ""), .cols = everything())
  
  columns_to_keep <- c(1, which(colSums(upset_table[, -1]) > 1) + 1)
  upset_table_df_filtered <- upset_table[, columns_to_keep]
  
  sample_year_data <- variant_table_homozygous %>%
    select(Sampleid, Year_start) %>%
    distinct(Sampleid, .keep_all = TRUE) %>%
    mutate(Year_start = as.numeric(Year_start))
  
  upset_table_df_filtered <- upset_table_df_filtered %>%
    left_join(sample_year_data, by = "Sampleid") %>%
    mutate(year_bin = cut(Year_start, breaks = time_bins, include.lowest = TRUE, right = FALSE, labels = time_bin_labels))
  
  # convert to amino acid single letter code
  # First, apply convert_aaa_to_a to colnames
  new_colnames <- colnames(upset_table_df_filtered)[-c(1,length(colnames(upset_table_df_filtered))-1,length(colnames(upset_table_df_filtered)))]
  new_colnames <- sapply(new_colnames, convert_aaa_to_a, USE.NAMES = FALSE)
  
  # Assign new column names back to upset_table
  colnames(upset_table_df_filtered)[2:(length(colnames(upset_table_df_filtered)) - 2)] <- new_colnames
  
  ################################################################################
  # Step 3: Compute Haplotype Counts & Prevalence
  ################################################################################
  
  # dont remove rows with rowsum under 3 for this (thats only necessary for wildtype reference resistant addition I think)
  
  haplotype_counts <- upset_table_df_filtered %>%
    select(-Sampleid, -Year_start) %>%
    group_by(year_bin, across(everything())) %>%
    summarise(count = n(), .groups = "drop")
  
  total_samples_per_bin <- upset_table_df_filtered %>%
    group_by(year_bin) %>%
    summarise(total_samples = n_distinct(Sampleid), .groups = "drop")
  
  haplotype_counts <- haplotype_counts %>%
    left_join(total_samples_per_bin, by = "year_bin") %>%
    mutate(prevalence = ifelse(is.na(count), 0, (count / total_samples) * 100))
  
  haplotype_counts <- full_join(tibble(year_bin = time_bin_labels), haplotype_counts, by = "year_bin") %>%
    replace_na(list(count = 0, prevalence = 0, total_samples = 0))
  
  ################################################################################
  # Step 4: Process Haplotype Labels
  ################################################################################
  
  # Define columns containing mutation markers
  mutation_cols <- setdiff(colnames(haplotype_counts), c("year_bin", "count", "total_samples", "prevalence", "haplotype_label"))
  
  # Rebuild haplotype_label using the actual column names
  haplotype_counts <- haplotype_counts %>%
    mutate(haplotype_label = apply(select(., all_of(mutation_cols)), 1, function(x) {
      present_mutations <- mutation_cols[which(x == 1)]
      
      if (length(present_mutations) == 0) return("WT")
      
      residue_positions <- as.numeric(str_extract(present_mutations, "\\d+"))
      sorted_mutations <- present_mutations[order(residue_positions)]
      
      paste(sorted_mutations, collapse = "-")
    }))
  
  # haplotype_counts <- haplotype_counts %>%
  #   mutate(haplotype_label = apply(select(., -count, -prevalence, -year_bin, -total_samples), 1, 
  #                                  function(x) paste(names(x)[which(x == 1)], collapse = "-")))
  
  plot_data <- haplotype_counts %>%
    select(year_bin, haplotype_label, prevalence)
  
  plot_data$haplotype_label[plot_data$haplotype_label == ""] <- "WT"
  
  transform_haplotype_labels_simple <- function(haplotype_labels) {
    sapply(haplotype_labels, function(label) {
      if (label == "WT") return("WT")
      
      # Split the haplotype string into individual mutations
      mutations <- unlist(strsplit(label, "-"))
      
      # Extract the final amino acid letter from each mutation
      aa_letters <- substr(mutations, nchar(mutations), nchar(mutations))
      
      # Combine original label with trailing letters in parentheses
      paste0(label, " (", paste(aa_letters, collapse = ""), ")")
    }, USE.NAMES = FALSE)
  }
  
  plot_data <- plot_data %>%
    mutate(haplotype_label = transform_haplotype_labels_simple(haplotype_label)) %>%
    filter(!is.na(year_bin))
  
  # plot_data <- plot_data %>%
  #   mutate(haplotype_label = transform_haplotype_labels(haplotype_label)) %>%
  #   filter(!is.na(year_bin))
  
  ################################################################################
  # Step 5: Plot Stacked Bar Chart for Haplotype Prevalences with Increased Font Sizes
  ################################################################################
  
  hap_prev_plot <- ggplot(plot_data, aes(x = year_bin, y = prevalence, fill = haplotype_label)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = palette_for_other_figures) +  # Use reserved scico colors
    labs(x = "Time Period", y = "Haplotype Prevalence (%)", fill = "Haplotype") +
    ggtitle(gene_title) +  # Uppercased gene name only
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
      axis.text.y = element_text(size = axis.text.setting),
      axis.title = element_text(size = axis.title.setting),
      plot.title = element_blank(),
      legend.title = element_text(size = legend.title.setting),
      legend.text = element_text(size = legend.text.setting)
    )
  
  ################################################################################
  # Step 6: Save the Plot
  ################################################################################
  
  plot_title <- paste0(fig_num,gene_to_plot, "_", treatment_era_to_plot, "_hap_prev_stacked_barplot_time_bins.svg")
  ggsave(file.path(output_path, plot_title), 
         plot = hap_prev_plot, width = 16, height = 10, units = "in", dpi = 600)
  
  # save the data for manuscript results
  base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data"
  write.csv(plot_data, paste0(base_path, gene_to_plot, "_stacked_barplot_results_for_manuscript.csv"), row.names = FALSE)
  
  # Save standardized data to /data folder
  standardized_filename <- paste0("/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/", 
                                 fig_num, toupper(gene_to_plot), "_haplotype_stacked_data.csv")
  write.csv(plot_data, standardized_filename, row.names = FALSE)
  
  print(paste("Plot saved to:", file.path(output_path, plot_title)))
}

################################################################################
#
# DHFR
#
################################################################################

gene_to_plot <- "dhfr"  # Example gene (can contain numbers or dashes)
treatment_era_to_plot <- "All_Eras"
output_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
mutations_of_interest <- dhfr_alleles

generate_haplotype_prevalence_plot(gene_to_plot, treatment_era_to_plot, output_path, 
                                   genotype_file, variant_table_cov_filtered_with_met, mutations_of_interest,"4B_")

################################################################################
#
# DHPS
#
################################################################################

gene_to_plot <- "dhps"  # Example gene (can contain numbers or dashes)
treatment_era_to_plot <- "All_Eras"
output_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
mutations_of_interest <- dhps_alleles

generate_haplotype_prevalence_plot(gene_to_plot, treatment_era_to_plot, output_path, 
                                   genotype_file, variant_table_cov_filtered_with_met, mutations_of_interest,"5B_")

################################################################################
#
# CRT
#
################################################################################
# ################################################################################
#
# CRT with reference alleles added to show CVIET (NOT USING SINCE 73V and 74I are not sequenced or called successfully with mips/miptools. Show up as 73V missing from genotype table and 74I all heterozygous calls (Abebe confirmed happens in all his data sets.))
#
################################################################################

# mdr1_alleles <- c("mdr1-Asn86Tyr","mdr1-Tyr184Phe","mdr1-Ser1034Cys","mdr1-Asn1042Asp","mdr1-Asp1246Tyr")
# reference_resistant <- c("mdr1-Asn86Tyr", "mdr1-Asp1246Tyr")
# reference_resistant_labels <- c("Asn86Asn", "Asp1246Asp")
# 
# #mdr1_alleles <- c("mdr1-Asn86Tyr","mdr1-Tyr184Phe","mdr1-Ser1034Cys","mdr1-Asn1042Asp","mdr1-Asp1246Tyr")
# crt_alleles <- c("crt-Cys72Ser", "crt-Val73Leu", "crt-Met74Ile","crt-Asn75Glu", "crt-Lys76Thr", 
#                  "crt-Ala220Ser")
# 
# reference_resistant <- c("crt-Cys72Ser", "crt-Val73Leu")
# reference_resistant_labels <- c("Cys72Cys", "crt-Val73Val")
# # reference_resistant <- c("crt-Cys72Ser", "crt-Val73Leu")
# # reference_resistant_labels <- c("Cys72Cys", "Val73Val")
# 
# gene_to_plot <- "crt"  # Example gene (can contain numbers or dashes)
# treatment_era_to_plot <- "All_Eras"
# output_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
# mutations_of_interest <- crt_alleles

# Run the function
# generate_haplotype_prevalence_plot_ref_resistant(gene_to_plot, treatment_era_to_plot, output_path,
#                                                  genotype_file, variant_table_cov_filtered_with_met,
#                                                  mutations_of_interest, reference_resistant, reference_resistant_labels)


# do not filter out heterozygous calls as 73/74 are always heterozygous in miptools. Use major allele for building haplotypes instead

generate_haplotype_prevalence_plot_crt_ref <- function(gene_to_plot, treatment_era_to_plot, output_path,
                                                       genotype_file, variant_table_cov_filtered_with_met,
                                                       mutations_of_interest, crt_reference_resistant,
                                                       crt_reference_resistant_labels,fig_num) {
  
  # Convert gene name to uppercase for plot title
  gene_title <- toupper(gene_to_plot)
  
  # Define time bins (1998-2021)
  time_bins <- c(1998, 2002, 2006, 2010, 2014, 2018, 2022)
  time_bin_labels <- c("1998-2001", "2002-2005", "2006-2009", "2010-2013", "2014-2017", "2018-2021")
  
  ################################################################################
  # Step 1: Filter & Format the Genotype Table
  ################################################################################
  
  filtered_genotype <- genotype_file %>%
    filter(mutation_name %in% mutations_of_interest) %>%
    select(sample, mutation_name, genotype)
  
  wide_genotype_table <- filtered_genotype %>%
    pivot_wider(names_from = mutation_name, values_from = genotype)
  
  valid_samples <- wide_genotype_table %>%
    rowwise() %>%
    filter(all(!c_across(where(is.numeric)) %in% c(1, -1))) %>%
    ungroup() %>%
    pull(sample)

  valid_samples_clean <- gsub("S-", "S", valid_samples)
  
  variant_table_homozygous <- variant_table_cov_filtered_with_met# %>%
  #filter(Sampleid %in% valid_samples_clean)
  
  ################################################################################
  # Step 2: Filter & Format the Variant Table, Include Reference-Resistant Alleles
  ################################################################################
  
  #00KS13-HAKpf-1
  variant_table_cov_filtered_with_met %>%
    filter(mutation_name %in% mutations_of_interest)
  
  upset_table <- variant_table_homozygous %>%
    filter(mutation_name %in% mutations_of_interest) %>%
    #mutate(present = ifelse(umi_alt > 2, 1, 0)) %>%
    mutate(present = ifelse(umi_alt > umi_ref, 1, 0)) %>%
    select(Sampleid, mutation_name, present) %>%
    pivot_wider(names_from = mutation_name, values_from = present, values_fill = list(present = 0))
  
  # Add reference-resistant alleles based on umi_ref > 1 condition
  for (i in seq_along(crt_reference_resistant)) {
    allele <- crt_reference_resistant[i]
    label <- crt_reference_resistant_labels[i]
    
    reference_present <- variant_table_homozygous %>%
      filter(mutation_name == allele & umi_ref > umi_alt) %>%
      select(Sampleid) %>%
      mutate(!!label := 1)  # Create new column with reference-resistant label
    
    upset_table <- upset_table %>%
      left_join(reference_present, by = "Sampleid") %>%
      mutate(!!label := ifelse(is.na(.data[[label]]), 0, .data[[label]]))  # Fill missing values with 0
  }
  
  # Keep only mutations with at least one occurrence
  columns_to_keep <- c(1, which(colSums(upset_table[, -1]) > 1) + 1)
  upset_table_df_filtered <- upset_table[, columns_to_keep]
  
  sample_year_data <- variant_table_homozygous %>%
    select(Sampleid, Year_start) %>%
    distinct(Sampleid, .keep_all = TRUE) %>%
    mutate(Year_start = as.numeric(Year_start))
  
  upset_table_df_filtered <- upset_table_df_filtered %>%
    left_join(sample_year_data, by = "Sampleid") %>%
    mutate(year_bin = cut(Year_start, breaks = time_bins, include.lowest = TRUE, right = FALSE, labels = time_bin_labels))
  
  ################################################################################
  # Step 3: Compute Haplotype Counts & Prevalence
  ################################################################################
  
  haplotype_counts <- upset_table_df_filtered %>%
    select(-Sampleid, -Year_start) %>%
    group_by(year_bin, across(everything())) %>%
    summarise(count = n(), .groups = "drop")
  
  total_samples_per_bin <- upset_table_df_filtered %>%
    group_by(year_bin) %>%
    summarise(total_samples = n_distinct(Sampleid), .groups = "drop")
  
  haplotype_counts <- haplotype_counts %>%
    left_join(total_samples_per_bin, by = "year_bin") %>%
    mutate(prevalence = ifelse(is.na(count), 0, (count / total_samples) * 100))
  
  haplotype_counts <- full_join(tibble(year_bin = time_bin_labels), haplotype_counts, by = "year_bin") %>%
    replace_na(list(count = 0, prevalence = 0, total_samples = 0))
  
  ################################################################################
  # Step 4: Process Haplotype Labels
  ################################################################################
  
  haplotype_counts <- haplotype_counts %>%
    mutate(haplotype_label = apply(select(., -count, -prevalence, -year_bin, -total_samples), 1,
                                   function(x) paste(names(x)[which(x == 1)], collapse = "-")))
  
  plot_data <- haplotype_counts %>%
    select(year_bin, haplotype_label, prevalence)
  
  
  # Ensure empty labels are set to "WT"
  plot_data$haplotype_label[plot_data$haplotype_label == ""] <- "WT"
  
  # Define exact haplotypes to merge into WT
  haplotypes_to_merge <- c("WT", "crt-Cys72Cys-crt-Val73Val")
  
  # Sum prevalence for selected haplotypes per year_bin
  summed_prevalence_data <- plot_data %>%
    filter(haplotype_label %in% haplotypes_to_merge) %>%
    group_by(year_bin) %>%
    summarise(prevalence = sum(prevalence, na.rm = TRUE), .groups = "drop") %>%
    mutate(haplotype_label = "WT")  # Rename merged rows as WT
  
  # Remove only WT and "crt-Cys72Cys-crt-Val73Val" from plot_data
  filtered_plot_data <- plot_data %>%
    filter(!(haplotype_label %in% haplotypes_to_merge))
  
  # Append the new summed-prevalence WT rows
  final_plot_data <- bind_rows(filtered_plot_data, summed_prevalence_data) %>%
    arrange(year_bin)  # Ensure correct chronological order
  
  # Display final dataset
  print(final_plot_data)
  
  plot_data <- final_plot_data
  plot_data <- plot_data %>%
    mutate(haplotype_label = transform_haplotype_labels(haplotype_label)) %>%
    filter(!is.na(year_bin))
  
  ################################################################################
  # Step 5: Plot Stacked Bar Chart for Haplotype Prevalences with Increased Font Sizes
  ################################################################################
  
  hap_prev_plot <- ggplot(plot_data, aes(x = year_bin, y = prevalence, fill = haplotype_label)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = palette_for_other_figures) +  # Use reserved scico colors
    labs(x = "Time Period", y = "Haplotype Prevalence (%)", fill = "Haplotype") +
    ggtitle(gene_title) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
      axis.text.y = element_text(size = axis.text.setting),
      axis.title = element_text(size = axis.title.setting),
      plot.title = element_blank(),
      legend.title = element_text(size = legend.title.setting),
      legend.text = element_text(size = legend.text.setting)
    )
  
  ################################################################################
  # Step 6: Save the Plot
  ################################################################################
  
  plot_title <- paste0(fig_num,gene_to_plot, "_", treatment_era_to_plot, "_hap_prev_stacked_barplot_time_bins.svg")
  ggsave(file.path(output_path, plot_title),
         plot = hap_prev_plot, width = 22, height = 10, units = "in", dpi = 600)
  
  # save the data for manuscript results
  base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data"
  write.csv(plot_data, paste0(base_path, gene_to_plot, "_stacked_barplot_results_for_manuscript.csv"), row.names = FALSE)
  
  # Save standardized data to /data folder
  standardized_filename <- paste0("/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/", 
                                 fig_num, toupper(gene_to_plot), "_haplotype_stacked_data.csv")
  write.csv(plot_data, standardized_filename, row.names = FALSE)
  
  print(paste("Plot saved to:", file.path(output_path, plot_title)))
}


crt_alleles <- c("crt-Cys72Ser","crt-Val73Leu","crt-Met74Ile","crt-Asn75Glu", "crt-Lys76Thr",
                 "crt-Ala220Ser","crt-Asn326Ser")

crt_reference_resistant <- c("crt-Cys72Ser","crt-Val73Leu")
crt_reference_resistant_labels <- c("crt-Cys72Cys", "crt-Val73Val")

gene_to_plot <- "crt"  # Example gene (can contain numbers or dashes)
treatment_era_to_plot <- "All_Eras"
output_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
mutations_of_interest <- crt_alleles

# Run the function
generate_haplotype_prevalence_plot_crt_ref(gene_to_plot, treatment_era_to_plot, output_path,
                                           genotype_file, variant_table_cov_filtered_with_met,
                                           mutations_of_interest, crt_reference_resistant, crt_reference_resistant_labels,"3B_")

# Save standardized data to /data folder (assumes plot_data is available from function)
# Note: Would need to modify function to return plot_data or save it within function with standardized name

################################################################################
#
# Figure 3-5) Lineplots for CRT, DHPS, DHFR
#
################################################################################

################################################################################
#
# 2C) mdr1 lineplot
#
################################################################################
# Define the alleles for each gene
dhps_alleles <- c("Ser436His", "Lys540Glu", "Ala437Gly", "Ala581Gly", "Ser436Ala","Ala613Ser")
dhfr_alleles <- c("Asn51Ile", "Cys59Arg", "Ser108Asn", "Ile164Leu")
crt_alleles <- c("Cys72Ser","Val73Leu","Met74Ile", "Asn75Glu", "Lys76Thr", "Ala220Ser")
aat1_alleles <- c("Ser258Leu")


# Filter data for relevant genes and alleles, excluding "PRE-ACT" and "POST-ACT"
filtered_data <- filtered_table %>%
  filter(
    (gene_name == "dhps" & mutation %in% dhps_alleles) |
      (gene_name == "dhfr-ts" & mutation %in% dhfr_alleles) |
      (gene_name == "crt" & mutation %in% crt_alleles)
  ) %>%
  filter(Year_start != "PRE-ACT" & Year_start != "POST-ACT") %>%  # Exclude non-numeric years
  mutate(Year_start = as.numeric(Year_start),
         Year_group = cut(Year_start, breaks = seq(1998, 2022, by = 2), 
                          labels = paste(seq(1998, 2020, by = 2), seq(1999, 2021, by = 2), sep = "-"), 
                          right = FALSE)) %>%  # `right = FALSE` ensures each interval includes the lower year
  group_by(gene_name, mutation, Year_group) %>%
  summarize(avg_prevalence = mean(prevalence_percentage, na.rm = TRUE)) %>%
  ungroup()

# Run dhfr
# Filter data for the specified gene
dhfr_data <- filtered_data %>% filter(gene_name == 'dhfr-ts')

dhfr_data <- dhfr_data %>%
  mutate(mutation = if_else(
    str_detect(mutation, "^[A-Z][a-z]{2}[0-9]+$"),  # Matches 'Asn86' pattern or anything with only the ref codon listed
    paste0(mutation, str_sub(mutation, 1, 3)), # Append first 3 letters again
    mutation                                        # Else keep as is
  ))

dhfr_data <- dhfr_data %>%
  mutate(residue_number = as.numeric(gsub("[A-Za-z]", "", mutation)))  # Extract numbers

dhfr_data <- dhfr_data %>%
  mutate(mutation = convert_aaa_to_a(mutation))  # Extract numbers

# Generate the plot
dhfr_lineplot <- ggplot(dhfr_data, aes(x = Year_group, y = avg_prevalence, color = reorder(mutation, residue_number), group = mutation)) +
  geom_line(size = 3) +
  geom_point(size = 5, alpha = 0.7) +
  labs(title = "DHFR", x = "Year (2-Year Intervals)", y = "Average Prevalence (%)", color = "Allele") +
  theme_minimal() +
  scale_color_manual(values = palette_for_other_figures) +  # Use reserved scico colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
    axis.text.y = element_text(size = axis.text.setting),
    axis.title = element_text(size = axis.title.setting),
    plot.title = element_blank(),
    legend.title = element_text(size = legend.title.setting),
    legend.text = element_text(size = legend.text.setting)
  )


# write results to table for results section
write.table(dhfr_data,sep=",",file="/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data/fig_4A_dhfr_lineplot_data.csv",row.names = F,col.names = T)

# Save standardized data to /data folder (corrected figure number)
write.csv(dhfr_data, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig4C_DHFR_lineplot_data.csv", row.names = FALSE)

# Save the plot in SVG format with double-wide dimensions and 600 dpi
plot_title = "4C_lineplot_dhfr.svg"
ggsave(paste0(base_path, plot_title), 
  plot = dhfr_lineplot, width = 33, height = 10, units = "in", dpi = 600)

# Run dhps
# Filter data for the specified gene
dhps_data <- filtered_data %>% filter(gene_name == 'dhps')

dhps_data <- dhps_data %>%
  mutate(mutation = if_else(
    str_detect(mutation, "^[A-Z][a-z]{2}[0-9]+$"),  # Matches 'Asn86' pattern or anything with only the ref codon listed
    paste0(mutation, str_sub(mutation, 1, 3)), # Append first 3 letters again
    mutation                                        # Else keep as is
  ))

dhps_data <- dhps_data %>%
  mutate(residue_number = as.numeric(gsub("[A-Za-z]", "", mutation)))  # Extract numbers

dhps_data <- dhps_data %>%
  mutate(mutation = convert_aaa_to_a(mutation))  # Extract numbers

# Generate the plot
dhps_lineplot <- ggplot(dhps_data, aes(x = Year_group, y = avg_prevalence, color = reorder(mutation, residue_number), group = mutation)) +
  geom_line(size = 3) +
  geom_point(size = 5, alpha = 0.7) +
  labs(title = "DHPS", x = "Year (2-Year Intervals)", y = "Average Prevalence (%)", color = "Allele") +
  theme_minimal() +
  scale_color_manual(values = palette_for_other_figures) +  # Use reserved scico colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
    axis.text.y = element_text(size = axis.text.setting),
    axis.title = element_text(size = axis.title.setting),
    plot.title = element_blank(),
    legend.title = element_text(size = legend.title.setting),
    legend.text = element_text(size = legend.text.setting)
  )

# write results to table for results section
write.table(dhps_data,sep=",",file="/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data/fig_5A_dhps_lineplot_data.csv",row.names = F,col.names = T)

# Save standardized data to /data folder (corrected figure number)
write.csv(dhps_data, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig5C_DHPS_lineplot_data.csv", row.names = FALSE)

# Save the plot in SVG format with double-wide dimensions and 600 dpi
plot_title = "5C_lineplot_dhps.svg"
ggsave(paste0(base_path, plot_title), 
       plot = dhps_lineplot, width = 33, height = 10, units = "in", dpi = 600)

# Run crt

# Filter data for the specified gene
crt_data <- filtered_data %>% filter(gene_name == 'crt')


crt_data <- crt_data %>%
  mutate(mutation = if_else(
    str_detect(mutation, "^[A-Z][a-z]{2}[0-9]+$"),  # Matches 'Asn86' pattern or anything with only the ref codon listed
    paste0(mutation, str_sub(mutation, 1, 3)), # Append first 3 letters again
    mutation                                        # Else keep as is
  ))

crt_data <- crt_data %>%
  mutate(residue_number = as.numeric(gsub("[A-Za-z]", "", mutation)))  # Extract numbers

crt_data <- crt_data %>%
  mutate(mutation = convert_aaa_to_a(mutation))  # Extract numbers

# Generate the plot
crt_lineplot <- ggplot(crt_data, aes(x = Year_group, y = avg_prevalence, color = reorder(mutation, residue_number), group = mutation)) +
  geom_line(size = 3, position = position_jitter(width = 0, height = 0.5)) +  # Add vertical jitter to lines
  geom_point(size = 5, alpha = 0.7, position = position_jitter(width = 0, height = 0.5)) +  # Apply jitter to points
  labs(title = "CRT", x = "Year (2-Year Intervals)", y = "Average Prevalence (%)", color = "Allele") +
  theme_minimal() +
  scale_color_manual(values = palette_for_other_figures) +  # Use reserved scico colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
    axis.text.y = element_text(size = axis.text.setting),
    axis.title = element_text(size = axis.title.setting),
    plot.title = element_blank(),
    legend.title = element_text(size = legend.title.setting),
    legend.text = element_text(size = legend.text.setting)
  )


# write results to table for results section
write.table(crt_data,sep=",",file="/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data/fig_3A_crt_lineplot_data.csv",row.names = F,col.names = T)

# Save standardized data to /data folder (corrected figure number)
write.csv(crt_data, "/Users/george/Bailey_lab/kenya_longitudinal/manuscript_submission_2025/data/Fig3C_CRT_lineplot_data.csv", row.names = FALSE)

# Save the plot in SVG format with double-wide dimensions and 600 dpi
plot_title = "3C_lineplot_crt.svg"
ggsave(paste0(base_path, plot_title), 
       plot = crt_lineplot, width = 33, height = 10, units = "in", dpi = 600)




# 
# 
# # Load necessary libraries
# library(readr)
# library(ggpubr)
# library(dplyr)
# library(rstatix)
# library(ggplot2)
# 
# # Define the file path
# file_path <- here("data", "kenya_all_DR2_prevalence_results_table_with_reference_resistant.csv")
# 
# # define results path
# base_path <- "/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/"
# 
# # Read in the prevalence table
# all_prev_table <- read_csv(file_path)
# 
# # Print the first few rows to verify
# print(head(all_prev_table))
# 
# # Remove Eastern kenya sites with too few samples spanning entire collection period
# # Define sites to exclude
# excluded_sites <- c("Kilifi", "Malindi")
# 
# # Filter out the unwanted rows
# filtered_table <- all_prev_table[!(all_prev_table$site_uid %in% excluded_sites), ]
# 
# key_mutations_aaa <- c("Asp1246Tyr", "Asn86Tyr", "Tyr184Phe")
# mdr1_reference_resistant_alleles <- c("Asn86Tyr", "Asp1246Tyr")
# mdr1_variant_resistant_alleles <- c("Tyr184Phe")
# 
# mdr1_data <- filtered_table %>%
#   filter(gene_name == "mdr1" & mutation %in% c(mdr1_reference_resistant_alleles, mdr1_variant_resistant_alleles)) %>%
#   mutate(
#     prevalence_custom = case_when(
#       mutation %in% mdr1_reference_resistant_alleles ~
#         ((samples_with_adequate_coverage - samples_with_mutation) / samples_with_adequate_coverage) * 100,
#       mutation %in% mdr1_variant_resistant_alleles ~ prevalence_percentage
#     ),
#     # Create display labels, removing "Tyr" from reference alleles
#     display_label = case_when(
#       mutation %in% mdr1_reference_resistant_alleles ~ sub("Tyr$", "", mutation),  # Remove "Tyr" from end
#       TRUE ~ mutation  # Keep full name for variant alleles
#     )
#   ) %>%
#   filter(!is.na(prevalence_custom))  # Remove rows with missing prevalence
# 
# mdr1_data <- mdr1_data %>%
#   filter(Year_start != "PRE-ACT" & Year_start != "POST-ACT") %>%  # Exclude non-numeric years
#   mutate(Year_start = as.numeric(Year_start),
#          Year_group = cut(Year_start, breaks = seq(1998, 2022, by = 2),
#                           labels = paste(seq(1998, 2020, by = 2), seq(1999, 2021, by = 2), sep = "-"),
#                           right = FALSE)) %>%  # Right = FALSE for left-inclusive intervals
#   group_by(display_label, Year_group) %>%
#   summarize(avg_prevalence = mean(prevalence_custom, na.rm = TRUE)) %>%
#   ungroup()
# 
# library(dplyr)
# library(stringr)
# 
# mdr1_data <- mdr1_data %>%
#   mutate(display_label = if_else(
#     str_detect(display_label, "^[A-Z][a-z]{2}[0-9]+$"),  # Matches 'Asn86' pattern or anything with only the ref codon listed
#     paste0(display_label, str_sub(display_label, 1, 3)), # Append first 3 letters again
#     display_label                                        # Else keep as is
#   ))
# 
# mdr1_data_single_a <- mdr1_data %>%
#   mutate(display_label = convert_aaa_to_a(display_label))  # Extract numbers
# 
# mdr1_data_single_a <- mdr1_data_single_a %>%
#   mutate(residue_number = as.numeric(gsub("[A-Za-z]", "", display_label)))  # Extract numbers
# 
# # # Plot for mdr1 with modified labels in the legend
# mdr1_plot <- ggplot(mdr1_data_single_a, aes(x = Year_group, y = avg_prevalence, color = reorder(display_label, residue_number), group = display_label)) +
#   geom_line(size = 3) +
#   geom_point(size = 5, alpha = 0.7) +
#   labs(title = "MDR1", x = "Year (2-Year Intervals)", y = "Average Prevalence (%)", color = "Allele") +
#   theme_minimal() +
#   #scale_color_brewer(palette = "Paired") +  # Use "Paired" palette for more colors
#   scale_color_manual(values = palette_for_other_figures) +  # Use reserved scico colors
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.setting),
#     axis.text.y = element_text(size = axis.text.setting),
#     axis.title = element_text(size = axis.title.setting),
#     plot.title = element_blank(),
#     legend.title = element_text(size = legend.title.setting),
#     legend.text = element_text(size = legend.text.setting)
#   )
# 
# # write results to table for results section
# write.table(mdr1_data,sep=",",file="/Users/george/Bailey_lab/kenya_longitudinal/2025_organization/figure_prep/final_figures/data/fig_2A_mdr1_lineplot_data.csv",row.names = F,col.names = T)
# 
# # Save standardized data to /data folder (this appears to be duplicate code - probably for different analysis)
# # Note: This section appears to be duplicate code and should probably be removed or clarified
# 
# # Save the plot
# plot_title = "lineplot_mdr1.svg"
# ggsave(paste0(base_path, plot_title), 
#        plot = mdr1_plot, width = 33, height = 10, units = "in", dpi = 600)
# 
