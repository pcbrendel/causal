library(rentrez)
library(svglite)
library(ggplot2)
library(stringr)

source("art/dna_art.R")
text_path <- "art/cftr_deltaf508.txt"
letter_plot_path <- "art/cftr_dna_letters_deltaf508.svg"
dot_plot_path <- "art/cftr_dna_dots_deltaf508.svg"

# Inputs ----

cftr_id <- "NM_000492.4"
# Deletion of phenylalanine at position 508 in the protein
# Define the positions for the deletion (1-based indexing)
# Using current NCBI annotation: CDS starts at position 71
# F508 codon is the 508th codon in the CDS
# CDS position: (508-1)*3 + 1 = 1522 (start of codon)
# mRNA position: 71 + 1522 - 1 = 1592
deletion_start <- 1522
deletion_end <- 1524
n_per_row <- 80
title <- "CFTR ΔF508 Coding DNA Sequence (Location 7q31.2)"
colors <- c(
  A = "forestgreen",
  T = "firebrick",
  C = "royalblue",
  G = "goldenrod"
)

# Get sequence ----

cftr_dna <- scrape_dna(cftr_id)
cds <- get_cds(cftr_id)
cftr_cds <- str_sub(
  cftr_dna,
  cds[1],
  cds[2]
)
cftr_cds

# Apply mutation ----

# Extract the nucleotides to be deleted for verification
# should be TTT or TTC for phenylalanine
nucleotides_to_delete <- str_sub(
  cftr_cds, deletion_start, deletion_end
)
message(
  paste(
    "\nNucleotides to be deleted (positions",
    deletion_start,
    "-",
    deletion_end,
    "):",
    nucleotides_to_delete
  )
)

# Create the deltaF508 sequence by removing the three nucleotides
delta_f508 <- paste0(
  str_sub(cftr_cds, 1, deletion_start - 1),
  ".", ".", ".",
  str_sub(cftr_cds, deletion_end + 1, -1)
)

# Show the region around the deletion for comparison
context_size <- 20
healthy_context <- str_sub(
  cftr_cds,
  deletion_start - context_size,
  deletion_end + context_size
)
mutant_context <- str_sub(
  delta_f508,
  deletion_start - context_size,
  deletion_end + context_size
)

message(
  paste(
    "\nHealthy sequence around deletion (positions",
    deletion_start - context_size,
    "-",
    deletion_end + context_size,
    "):"
  )
)
message(healthy_context)
message(paste("\nDeltaF508 sequence around deletion site:"))
message(mutant_context)

# Save and format ----

# Save the deltaF508 sequence
file_conn <- file(text_path)
writeLines(delta_f508, file_conn)
close(file_conn)

# format
delta_f508_spaces <- str_replace_all(
  delta_f508,
  paste0("(.{", letters_per_row, "})(?=.)"), "\\1 "
)
delta_f508_wrapped <- str_wrap(
  delta_f508_spaces,
  width = letters_per_row
)
delta_f508_wrapped <- str_replace_all(
  delta_f508_wrapped,
  "\\.",
  " "
)

cat(delta_f508_wrapped)

# Convert sequence to vector
delta_f508_vec <- strsplit(delta_f508, "")[[1]]

# Plot 1: Letters ----

svglite(file = letter_plot_path, width = 12, height = 12)

par(
  mar = c(0, 0, 0, 0),
  omi = c(0, 0, 0, 0),
  ann = FALSE,
  bty = "n"
)

# initialize
plot(
  NULL,
  xlim = c(0, 1),
  ylim = c(0, 1.1), # increased upper limit for title
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "",
  sub = "",
  frame.plot = FALSE
)

# optional title
text(
  x = 0.005,
  y = 1.02, # just above the top, reduced space
  labels = title,
  cex = 1.2, # larger font for title
  font = 2, # bold
  adj = c(0, 0) # left alignment
)

lines_to_plot <- unlist(strsplit(delta_f508_wrapped, "\n"))

text_cex <- 0.8
vfont_type <- c("serif", "plain")
line_height <- strheight("Mg", cex = text_cex, font = 1) * 1.5
char_width <- 0.01

current_y <- 0.995
start_x <- 0.005

for (line_idx in seq_along(lines_to_plot)) {
  line_text <- lines_to_plot[line_idx]
  current_x <- start_x

  for (char_idx in seq_along(strsplit(line_text, "")[[1]])) {
    char_to_plot <- str_sub(line_text, char_idx, char_idx)

    text(
      x = current_x,
      y = current_y,
      labels = char_to_plot,
      vfont = vfont_type,
      cex = text_cex,
      adj = c(0, 0.5)
    )

    current_x <- current_x + char_width
  }

  current_y <- current_y - line_height
}

dev.off()

# Plot 2: DNA as colored dots ----

dot_plot(
  dna_seq = delta_f508_vec,
  svg_file_name = dot_plot_path,
  n_per_row = n_per_row,
  colors = colors
)