
# HADEG analysis of Thalassolituus proteomes
# Using the HADEG R script developed by Rojas-Vargas et al. (2023):
# https://github.com/jarojasva/HADEG/tree/main/Scripts

setwd("/Users/somabarawi/OneDrive - Dalhousie University/Bioinformatics/Chapter-1-Analyses/HDE-analysis/HADEG-results")

library(stringr)
library(dplyr)
library(readr)
library(scales)
library(svglite)
library(ggplot2)
library(reshape)

# Read the original file with the first row as header
df <- read.table("Results_HADEG.proteinortho.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "", fill = TRUE)

# Rename the first three columns as ColumnA, ColumnB, and ColumnC
colnames(df)[1:3] <- c("ColumnA", "ColumnB", "ColumnC")

# Adjust the header to have the 'HADEG_protein_database_231119.faa' column at the beginning
df <- df[, c(grep("HADEG_protein_database_231119.faa", colnames(df)), setdiff(seq_along(df), grep("HADEG_protein_database_231119.faa", colnames(df))))]

# Remove specific columns
columns_to_remove <- c("ColumnA", "ColumnB", "ColumnC")
df_selected <- df[, !(colnames(df) %in% columns_to_remove)]

# Split the 'HADEG_protein_database_231119.faa' column into two new columns
split_columns <- str_split_fixed(df_selected$HADEG_protein_database_231119.faa, "\\|", 2)

# Assign names to the new columns
colnames(split_columns) <- c("Protein_ID", "Gene_name")

# Combine the new columns with the existing dataframe
df_selected <- cbind(df_selected, split_columns)

# Remove the original column
df_selected <- df_selected[, !grepl("HADEG_protein_database_231119.faa", colnames(df_selected))]

# Move the columns "Protein_ID" and "Gene_name" as the first two columns
df_selected <- df_selected[, c("Protein_ID", "Gene_name", setdiff(names(df_selected), c("Protein_ID", "Gene_name")))]

# Remove lines containing "*" in the "Protein_ID" column
df_selected <- df_selected[!grepl("\\*", df_selected$Protein_ID), ]

# Split the "Gene_name" column in cases where there are commas
split_genes <- str_split(df_selected$Gene_name, ",")

# Take only the first code before the commas (or the complete code if there are no commas)
df_selected$Gene_name <- sapply(split_genes, function(x) x[1])

# Save the result in a new file
write.table(df_selected, file = "1_table_HADEG_codes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Read the modified final table
df_codes <- read.table("1_table_HADEG_codes.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "", fill = TRUE)

# Create a new table with the first two columns and the corresponding headers
df_counts <- data.frame(Protein_ID = df_codes$Protein_ID, Gene_name = df_codes$Gene_name)

# Function to count the number of codes in each cell
count_codes <- function(x) {
  ifelse(grepl("\\*", x), 0, str_count(x, ",") + 1)
}

# Apply the function to all remaining columns and keep the headers
df_counts <- cbind(df_counts, sapply(df_codes[, 3:ncol(df_codes)], count_codes))

# Name the new columns
colnames(df_counts)[3:ncol(df_counts)] <- colnames(df_codes)[3:ncol(df_codes)]

# Order by the second column "Gene_name"
df_counts <- df_counts[order(df_counts$Gene_name), ]

# Save the result in a new file
write.table(df_counts, file = "2_table_HADEG_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Read the All_pathways.csv file
all_pathways <- read.csv("7_All_pathways.csv", stringsAsFactors = FALSE)

# Ensure that there are no leading or trailing whitespaces in the Protein_ID column in both data frames
df_counts$Protein_ID <- trimws(df_counts$Protein_ID)
all_pathways$Protein_ID <- trimws(all_pathways$Protein_ID)

# Merge df_counts with all_pathways based on Protein_ID
df_counts_merged <- merge(df_counts, all_pathways, by = "Protein_ID", all.x = TRUE)

# Keep the 74 focal genomes used in this study
keep_genomes <- c(
  "GCA_002294235","GCA_002298335","GCA_002314145","GCA_002331085","GCA_002450785",
  "GCA_040118715","GCA_040118895","GCA_040119005","GCA_041247055","GCA_905478475",
  "GCA_913043175","GCA_913047515","GCA_913052225","GCA_937895325","GCA_937902495",
  "GCA_963998775","GCA_963999285","GCA_964007845","GCA_964012525",
  "HBD_gamma01","HBD_gamma04","HBD_gamma08","HBD_gamma11","HBD_gamma6","HBD_gamma9",
  "New_ARC_Gamma03",
  "GCF_000355675","GCF_000582555","GCF_001932535","GCF_002297765","GCF_002324005",
  "GCF_002324035","GCF_002331105","GCF_002335285","GCF_002336415","GCF_002339775",
  "GCF_002377545","GCF_002433655","GCF_007785795","GCF_015665225","GCF_015665495",
  "GCF_018066875","GCF_020320395","GCF_020320405","GCF_020340755","GCF_020447205",
  "GCF_020447245","GCF_022832715","GCF_025318105","GCF_025345565","GCF_026626885",
  "GCF_027018755","GCF_029210095","GCF_030867545","GCF_034111725","GCF_040119275",
  "GCF_040372045","GCF_041222825","GCF_041230555","GCF_041234595","GCF_041250145",
  "GCF_041258785","GCF_041259185","GCF_041259835","GCF_046582215","GCF_046582255",
  "GCF_900156675","GCF_913019975","GCF_913058015","GCF_937896655","GCF_950066645",
  "GCF_950069745","GCF_963922525","GCF_964005015"
)

keep_cols <- paste0(keep_genomes, ".faa")
missing_genomes <- setdiff(keep_cols, colnames(df_counts_merged))

if (length(missing_genomes) > 0) {
  warning("These genomes were not found in the Proteinortho results: ",
          paste(missing_genomes, collapse = ", "))
}

keep_cols <- intersect(keep_cols, colnames(df_counts_merged))

df_counts_merged <- df_counts_merged[, c(
  "Protein_ID", "Gene_name", "Mechanism", "Compound",
  "Pathway", "Subpathway", "Gene", keep_cols
)]

# Drop columns "Protein_ID", "Gene_name", "code_mechanism", "code_compound", and "code_subpathway"
df_counts_merged <- df_counts_merged[, !(colnames(df_counts_merged) %in% c("Protein_ID", "Gene_name", "code_mechanism", "code_compound", "code_subpathway"))]

# Reorder the columns
df_counts_merged <- df_counts_merged[, c("Mechanism", "Compound", "Pathway", "Subpathway", "Gene", setdiff(names(df_counts_merged), c("Mechanism", "Compound", "Pathway", "Subpathway", "Gene")))]

# Organize df_counts_merged by Gene alphabetically
df_counts_merged <- df_counts_merged %>%
  arrange(Gene)

# Group by Gene and sum the counts across all columns ending with ".faa"
df_counts_aggregated <- df_counts_merged %>%
  group_by(Gene) %>%
  summarise(across(ends_with(".faa"), sum, na.rm = TRUE)) %>%
  distinct()

# Create a summary data frame for all_pathways
all_pathways_summary <- all_pathways %>%
  select(-Protein_ID) %>%
  distinct()

# Merge df_counts_aggregated with all_pathways_summary based on Gene
df_counts_aggregated_2 <- merge(df_counts_aggregated, all_pathways_summary, by = "Gene", all.x = TRUE)

# Remove duplicates based on Gene
df_counts_aggregated_2 <- df_counts_aggregated_2 %>% distinct(Gene, .keep_all = TRUE)

# Reorder the columns
df_counts_aggregated_2 <- df_counts_aggregated_2[, c("Mechanism", "Compound", "Pathway", "Subpathway", "Gene", setdiff(names(df_counts_aggregated_2), c("Mechanism", "Compound", "Pathway", "Subpathway", "Gene")))]

# Drop columns "code_mechanism", "code_compound", and "code_subpathway"
df_counts_aggregated_2 <- df_counts_aggregated_2[, !(colnames(df_counts_aggregated_2) %in% c("code_mechanism", "code_compound", "code_subpathway"))]

# Save the result in a new file
write.table(df_counts_aggregated_2, file = "3_table_HADEG_final.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Get all unique names in the "Compound" column
unique_compounds <- unique(df_counts_aggregated_2$Compound)

# Iterate over the list of unique compounds
for (compound in unique_compounds) {
  # Filter the DataFrame to get only rows with the current compound
  compound_df <- df_counts_aggregated_2[df_counts_aggregated_2$Compound == compound, ]
  
  # Create the file name with the desired format and save the DataFrame to a separate file
  file_name <- paste0("4_", gsub("/", "_", compound), ".tsv")
  write.table(compound_df, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}



## BASED ON GENES

# Group by Gene and sum the counts across all columns ending with ".faa"
df_counts_aggregated_by_gene <- df_counts_merged %>%
  group_by(Gene) %>%
  summarise(across(ends_with(".faa"), sum, na.rm = TRUE)) %>%
  distinct()

# Create a summary data frame for all_pathways
all_pathways_summary <- all_pathways %>%
  select(-Protein_ID) %>%
  distinct()

# Merge df_counts_aggregated with all_pathways_summary based on Gene
df_counts_aggregated_by_gene_2 <- merge(df_counts_aggregated_by_gene, all_pathways_summary, by = "Gene", all.x = TRUE)

# Remove duplicates based on Gene
df_counts_aggregated_by_gene_2 <- df_counts_aggregated_by_gene_2 %>% distinct(Gene, .keep_all = TRUE)

# Reorder the columns
df_counts_aggregated_by_gene_2 <- df_counts_aggregated_by_gene_2[, c("Mechanism", "Compound", "Pathway", "Subpathway", "Gene", setdiff(names(df_counts_aggregated_by_gene_2), c("Mechanism", "Compound", "Pathway", "Subpathway", "Gene")))]

# Drop columns "code_mechanism", "code_compound", and "code_subpathway"
df_counts_aggregated_by_gene_2 <- df_counts_aggregated_by_gene_2[, !(colnames(df_counts_aggregated_by_gene_2) %in% c("code_mechanism", "code_compound", "code_subpathway"))]

# Save the result in a new file
write.table(df_counts_aggregated_by_gene_2, file = "3_table_HADEG_final.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Get all unique names in the "Compound" column
unique_compounds <- unique(df_counts_aggregated_by_gene_2$Compound)

# Iterate over the list of unique compounds
for (compound in unique_compounds) {
  # Filter the DataFrame to get only rows with the current compound
  compound_df <- df_counts_aggregated_by_gene_2[df_counts_aggregated_by_gene_2$Compound == compound, ]
  
  # Create the file name with the desired format and save the DataFrame to a separate file
  file_name <- paste0("4_", gsub("/", "_", compound), ".tsv")
  write.table(compound_df, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}


# Get all unique names in the "Compound" column
unique_compounds <- unique(df_counts_aggregated_by_gene_2$Compound)

# List to store the read DataFrames
df_list <- list()

# Iterate over the list of unique compounds
for (compound in unique_compounds) {
  # Create the file name with the desired format
  file_name <- paste0("4_", gsub("/", "_", compound), ".tsv")
  
  # Read the specific TSV file for the current compound
  compound_df <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "", fill = TRUE)
  
  # Add the DataFrame to the list
  df_list[[file_name]] <- compound_df
}

# Apply melt to each DataFrame in the list
melted_df_list <- lapply(df_list, function(df) melt(df, id.vars = c("Gene", "Mechanism", "Compound", "Pathway", "Subpathway")))

# Access the DataFrames in the list using the file names as keys
# Example: df_list$"4_A_Alkanes.tsv"

# You can access the melted DataFrames through melted_df_list
# melted_df_list[[1]] contains the first melted DataFrame, and so on

# Define the new column names
new_column_names <- c("Gene", "Mechanism", "Compound", "Pathway", "Subpathway", "Genome", "Number")

# Iterate over the list of melted DataFrames
for (i in seq_along(melted_df_list)) {
  # Change the column names
  colnames(melted_df_list[[i]]) <- new_column_names
}

# Combine all melted DataFrames into one
combined_df <- bind_rows(melted_df_list, .id = "group")

# Define the color names and legend titles
colors <- c("#a9d6e5", "#012a4a")  # You can add more colors as needed

# Generate the heatmap plot
svg("Heatmap_Genes_HADEG-test.svg", width=15, height=6)
ggplot(combined_df, aes(x = Gene, y = Genome, fill = Number, group = group)) +
  geom_tile(color = "white", lwd = 1) +
  scale_fill_gradientn(colors = colors, name = "", labels = scales::comma) +
  facet_grid(~ group, scales = "free_x", space = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(hjust = 1, size = 9),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.text.x = element_text(size = 10, face = "bold")
  )
dev.off()




## BASED ON SUBPATHWAYS

# Group by Gene and sum the counts across all columns ending with ".faa"
df_counts_aggregated_by_subpathway <- df_counts_merged %>%
  group_by(Subpathway) %>%
  summarise(across(ends_with(".faa"), sum, na.rm = TRUE)) %>%
  distinct()

# Create a summary data frame for all_pathways
all_pathways_summary <- all_pathways %>%
  select(-Protein_ID) %>%
  distinct()

# Merge df_counts_aggregated with all_pathways_summary based on Subpathway
df_counts_aggregated_by_subpathway_2 <- merge(df_counts_aggregated_by_subpathway, all_pathways_summary, by = "Subpathway", all.x = TRUE)

# Remove specific columns
columns_to_remove <- c("Gene")
df_counts_aggregated_by_subpathway_2 <- df_counts_aggregated_by_subpathway_2[, !(colnames(df_counts_aggregated_by_subpathway_2) %in% columns_to_remove)]

# Remove duplicates based on Subpathway
df_counts_aggregated_by_subpathway_2 <- df_counts_aggregated_by_subpathway_2 %>% distinct(Subpathway, .keep_all = TRUE)

# Reorder the columns
df_counts_aggregated_by_subpathway_2 <- df_counts_aggregated_by_subpathway_2[, c("Mechanism", "Compound", "Pathway", "Subpathway", setdiff(names(df_counts_aggregated_by_subpathway_2), c("Mechanism", "Compound", "Pathway", "Subpathway")))]

# Drop columns "code_mechanism", "code_compound", and "code_subpathway"
df_counts_aggregated_by_subpathway_2 <- df_counts_aggregated_by_subpathway_2[, !(colnames(df_counts_aggregated_by_subpathway_2) %in% c("code_mechanism", "code_compound", "code_subpathway"))]


# Generate the bubble plot

# Apply melt to obtain a general DataFrame with subtpathways
df_counts_aggregated_by_subpathway_2_melted <- melt(df_counts_aggregated_by_subpathway_2, id.vars = c("Subpathway", "Mechanism", "Compound", "Pathway"))

# Rename the first three columns as ColumnA, ColumnB, and ColumnC
colnames(df_counts_aggregated_by_subpathway_2_melted) <- c("Subpathway", "Mechanism", "Compound", "Pathway", "Genome", "Hits_number")

clean_genome <- function(x) {
  x <- sub("\\.faa$", "", as.character(x))
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub("^ARC_Gamma033$", "New_ARC_Gamma03", x)
  x <- sub("^ARC_Gamma03$", "New_ARC_Gamma03", x)
  x <- sub("^HBD_Gamma_?", "HBD_gamma", x)
  x <- sub("^HBD_gamma06$", "HBD_gamma6", x)
  x <- sub("^HBD_gamma09$", "HBD_gamma9", x)
  x
}

tree <- read.tree("/Users/somabarawi/Desktop/Trees/species.tree")
tree$tip.label <- clean_genome(tree$tip.label)
names(genome_labels) <- clean_genome(names(genome_labels))

CLADE_COLORS <- c("Clade 1" = "#f4a582", "Clade 2" = "#92c5de", "Clade 3" = "#0571b0", "Clade 4" = "#ca0020")

CLADE_TIPS <- list(
  "Clade 1" = c("ARC_Gamma_03", "HBD_Gamma_11", "HBD_Gamma_08", "Thalassolituus haligoni BB40", "Parathalassolituus penaei G-43", "uncultured Thalassolituus sp. 10"),
  "Clade 2" = c("HBD_Gamma_09", "HBD_Gamma_06", "HBD_Gamma_01", "Thalassolituus sp. UBA6592", "Thalassolituus sp. UBA1505", "uncultured Thalassolituus sp. 7", "uncultured Thalassolituus sp. 5", "uncultured Thalassolituus sp. 3", "Thalassolituus sp. UBA6928", "Thalassolituus sp. UBA1087", "Thalassolituus sp. UBA885"),
  "Clade 3" = c("HBD_Gamma_04", "uncultured Thalassolituus sp. 9", "uncultured Thalassolituus sp. 6", "uncultured Thalassolituus sp. 4", "Thalassolituus sp. UBA2109", "Thalassolituus sp. UBA674"),
  "Clade 4" = c("Thalassolituus sp. MAG485", "Thalassolituus oleivorans SRR3933287")
)

tip_clade <- setNames(rep(names(CLADE_TIPS), lengths(CLADE_TIPS)), unlist(CLADE_TIPS, use.names = FALSE))

plot_data <- df_counts_aggregated_by_subpathway_2_melted %>%
  mutate(Genome = clean_genome(Genome), Hits_number = as.numeric(Hits_number)) %>%
  filter(!is.na(Subpathway), !is.na(Compound), !is.na(Hits_number), is.finite(Hits_number), Hits_number > 0)

# Order genomes according to the species tree
plot.phylo(tree, type = "phylogram", direction = "rightwards", show.tip.label = FALSE, plot = FALSE)
tree_layout <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- as.numeric(tree_layout$yy[seq_len(Ntip(tree))])
tree_order <- rev(tree$tip.label[order(tip_y)])
tree_order <- tree_order[tree_order %in% unique(plot_data$Genome)]

missing_from_tree <- setdiff(unique(plot_data$Genome), tree_order)
genome_order <- c(tree_order, missing_from_tree)
plot_data$Genome <- factor(plot_data$Genome, levels = genome_order)

# Replace accessions with simplified genome names
display_names <- unname(genome_labels[genome_order])
display_names[is.na(display_names)] <- genome_order[is.na(display_names)]

# Color putative diazotroph labels by clade
display_clades <- unname(tip_clade[display_names])
display_colors <- unname(CLADE_COLORS[display_clades])
display_labels <- ifelse(is.na(display_clades), display_names, paste0("<span style='color:", display_colors, "'>", display_names, "</span>"))
names(display_labels) <- genome_order

# HADEG compound colors
colors <- c("#003566", "#fcbf49", "#f77f00", "#d62828", "#caf0f8", "#03045e", "#a9d6e5")
size_breaks <- seq(min(plot_data$Hits_number, na.rm = TRUE), max(plot_data$Hits_number, na.rm = TRUE), by = 1)

p <- ggplot(plot_data, aes(x = Genome, y = Subpathway, size = Hits_number, fill = Compound)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_manual(values = colors) +
  scale_size_continuous(breaks = size_breaks, range = c(2, 9)) +
  scale_x_discrete(labels = display_labels, drop = FALSE) +
  labs(x = NULL, y = NULL, size = "Number of hits", fill = "Compound") +
  theme_bw(base_size = 16, base_family = "Arial") +
  theme(text = element_text(family = "Arial", size = 16), axis.text.x = element_markdown(angle = 45, vjust = 1, hjust = 1, size = 16), axis.text.y = element_text(size = 16, hjust = 1), panel.grid = element_blank(), legend.position = "right", legend.text = element_text(size = 16), legend.title = element_text(size = 16))

p
ggsave("Bubbles_Subpathways_HADEG_Thalassolituus-final.svg", plot = p, width = 25, height = 10, device = svglite::svglite)
