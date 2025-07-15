library(rentrez)
library(svglite)
library(ggplot2)
library(stringr)

# get sequence ----
# https://www.ncbi.nlm.nih.gov/nuccore/NM_000492.4
# entrez = a search engine and data retrieval system
# that allows users to search across dozens of NCBI databases
# nuccore = "Nucleotide Core" database at NCBI
# (National Center for Biotechnology Information)
# accession = unique, stable identifier for a sequence

# Define the specific RefSeq mRNA accession for human CFTR
# This is NM_000492.4 for human CFTR mRNA
cftr_mrna <- "NM_000492.4"

# Fetch the FASTA sequence directly using the accession ID
message(
  paste("Attempting to fetch mRNA sequence for accession:", cftr_mrna)
)
cftr_fasta_record <- tryCatch(
  {
    entrez_fetch(
      db = "nuccore",
      id = cftr_mrna,
      rettype = "fasta",
      retmode = "text"
    )
  },
  error = function(e) {
    message(paste("Error fetching sequence:", e$message))
    return(NULL)
  }
)

if (!is.null(cftr_fasta_record)) {
  lines <- strsplit(cftr_fasta_record, "\n")[[1]]
  cftr_dna_sequence <- paste(lines[-1], collapse = "")

  message("\nCFTR DNA (mRNA/cDNA) Sequence (first 100 characters):\n")
  print(substring(cftr_dna_sequence, 1, 100))
  message(paste("\nTotal length:", nchar(cftr_dna_sequence), "base pairs"))
} else {
  message("Failed to retrieve CFTR mRNA sequence using direct accession.")
}
# this is the coding DNA sequence
# it includes:
# - 5' UTR (non-coding before the CDS)
# - CDS (coding sequence)
# - 3' UTR (non-coding after the CDS)

# Create deltaF508 mutation ----
# Deletion of phenylalanine at position 508 in the protein

# Define the positions for the deletion (1-based indexing)
# Using current NCBI annotation: CDS starts at position 71
# F508 codon is the 508th codon in the CDS
# CDS position: (508-1)*3 + 1 = 1522 (start of codon)
# mRNA position: 71 + 1522 - 1 = 1592
deletion_start <- 1592
deletion_end <- 1594

# Extract the nucleotides to be deleted for verification
# should be TTT or TTC for phenylalanine
nucleotides_to_delete <- str_sub(
  cftr_dna_sequence, deletion_start, deletion_end
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
cftr_deltaF508_sequence <- paste0(
  str_sub(cftr_dna_sequence, 1, deletion_start - 1),
  ".", ".", ".",
  str_sub(cftr_dna_sequence, deletion_end + 1, -1)
)

# Show the region around the deletion for comparison
context_size <- 20
healthy_context <- str_sub(
  cftr_dna_sequence,
  deletion_start - context_size,
  deletion_end + context_size
)
mutant_context <- str_sub(
  cftr_deltaF508_sequence,
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

# save, modify sequence ----

# Save the deltaF508 sequence
file_conn <- file("art/cftr_deltaF508.txt")
writeLines(cftr_deltaF508_sequence, file_conn)
close(file_conn)

# Format for display
cftr_deltaF508_with_spaces <- str_replace_all(
  cftr_deltaF508_sequence,
  paste0("(.{", 100, "})(?=.)"), "\\1 "
)
cftr_deltaF508_sequence_wrapped <- str_wrap(
  cftr_deltaF508_with_spaces,
  width = 100
)
cftr_deltaF508_sequence_wrapped <- str_replace_all(
  cftr_deltaF508_sequence_wrapped,
  "\\.",
  " "
)

cat(cftr_deltaF508_sequence_wrapped)

# CDS only
cds_start <- 71
cds_end <- 4513

cftr_deltaF508_sequence_cds <- str_sub(
  cftr_deltaF508_sequence,
  cds_start,
  cds_end
)

cftr_deltaF508_cds_with_spaces <- str_replace_all(
  cftr_deltaF508_sequence_cds,
  paste0("(.{", 80, "})(?=.)"), "\\1 "
)
cftr_deltaF508_sequence_cds_wrapped <- str_wrap(
  cftr_deltaF508_cds_with_spaces,
  width = 80
)
cftr_deltaF508_sequence_cds_wrapped <- str_replace_all(
  cftr_deltaF508_sequence_cds_wrapped,
  "\\.",
  " "
)

cat(cftr_deltaF508_sequence_cds_wrapped)

# plot ----

svg_file_name <- "art/cftr_dna_deltaF508.svg"
svglite(file = svg_file_name, width = 12, height = 12)

par(mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0), ann = FALSE, bty = "n")

# initialize
plot(
  NULL,
  xlim = c(0, 1),
  ylim = c(0, 1),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "",
  sub = "",
  frame.plot = FALSE
)

lines_to_plot <- unlist(strsplit(cftr_deltaF508_sequence_cds_wrapped, "\n"))

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
