source("art/dna_art.R")
library(rentrez)
library(svglite)
library(ggplot2)
library(stringr)

text_path <- "art/cftr_deltaf508.txt"
letter_plot_path <- "art/dna_letters_cftr_deltaf508.svg"
dot_plot_path <- "art/dna_dots_cftr_deltaf508.svg"

# Inputs ----
# Deletion of phenylalanine at position 508 in the protein
# Define the positions for the deletion (1-based indexing)
# Using current NCBI annotation: CDS starts at position 71
# F508 codon is the 508th codon in the CDS
# CDS position: (508-1)*3 + 1 = 1522 (start of codon)
# mRNA position: 71 + 1522 - 1 = 1592

cftr_id <- "NM_000492.4"
deletion_start <- 1522
deletion_end <- 1524
n_per_row <- 80
plot_title <- "CFTR ΔF508 Coding DNA Sequence (Location 7q31.2)"
colors <- c(
  A = "forestgreen",
  T = "firebrick",
  C = "royalblue",
  G = "goldenrod"
)

# Get sequence ----

cftr_dna <- scrape_dna(cftr_id)
cds_positions <- get_cds_positions(cftr_id)
cftr_cds <- str_sub(
  cftr_dna,
  cds_positions[1],
  cds_positions[2]
)
cftr_cds

# Apply mutation ----
# Extract the nucleotides to be deleted for verification
# should be TTT or TTC for phenylalanine
delta_f508 <- delete_dna(cftr_cds, text_path, deletion_start, deletion_end)

# Plot 1: Letters ----

letter_plot(
  dna_seq = delta_f508,
  svg_file_name = letter_plot_path,
  n_per_row = n_per_row,
  colors = colors,
  title = TRUE,
  title_text = plot_title
)

# Plot 2: DNA as colored dots ----

dot_plot(
  dna_seq = delta_f508,
  svg_file_name = dot_plot_path,
  n_per_row = n_per_row,
  colors = colors
)
