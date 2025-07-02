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
  # Parse the FASTA record
  lines <- strsplit(cftr_fasta_record, "\n")[[1]]
  # Remove the header line and concatenate the sequence lines
  cftr_dna_sequence <- paste(lines[-1], collapse = "")

  message("\nCFTR DNA (mRNA/cDNA) Sequence (first 100 characters):\n")
  print(substring(cftr_dna_sequence, 1, 100))
  message(paste("\nTotal length:", nchar(cftr_dna_sequence), "base pairs"))

} else {
  message("Failed to retrieve CFTR mRNA sequence using direct accession.")
}

# save, modify sequence ----

cftr_dna_sequence

file_conn <- file("art/cftr_healthy.txt")
writeLines(cftr_dna_sequence, file_conn)
close(file_conn)

cftr_with_spaces <- str_replace_all(
  cftr_dna_sequence,
  paste0("(.{", 100, "})(?=.)"), "\\1 "
)
cftr_dna_sequence_wrapped <- str_wrap(
  cftr_with_spaces,
  width = 100
)
cat(cftr_dna_sequence_wrapped)

# plot ----

svg_file_name <- "art/cftr_dna_test.svg"
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

lines_to_plot <- unlist(strsplit(cftr_dna_sequence_wrapped, "\n"))

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

    char_to_plot <- substr(line_text, char_idx, char_idx)

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

# lines_to_plot <- unlist(strsplit(cftr_dna_sequence_wrapped, "\n"))

# text_cex <- 0.8
# vfont_type <- c("serif", "plain")
# line_height <- strheight("Mg", cex = text_cex, font = 1) * 1.2
# start_y <- 1

# for (i in seq_along(lines_to_plot)) {
#   current_y <- start_y - (i - 1) * line_height
#   text(
#     x = 0,
#     y = current_y,
#     labels = lines_to_plot[i],
#     vfont = vfont_type,
#     cex = text_cex,
#     adj = c(0, 1)
#   )
# }

dev.off()


# OLD ----

my_string <- substr(cftr_dna_sequence, 1, 200)
my_string <- cftr_dna_sequence

chars_per_row <- 100

# 2. Split the string into individual characters
all_chars <- unlist(strsplit(my_string, split = ""))
n_chars <- length(all_chars)

# 3. Create indices for columns and rows
# Column index repeats 1 to chars_per_row
col_indices <- rep(1:chars_per_row, ceiling(n_chars / chars_per_row))[1:n_chars]

# Row index increments every chars_per_row characters
row_indices <- ceiling((1:n_chars) / chars_per_row)

# 4. Create a data frame for plotting
plot_data <- data.frame(
  x = col_indices,
  y = -row_indices, # Use negative rows so the plot starts from the top (row 1)
  char = all_chars
)

# 5. Plot the characters
ggplot(plot_data, aes(x = x, y = y, label = char)) +
  geom_text(size = 4, family = "mono") + # Use a monospace font for consistent character width
  scale_x_continuous(breaks = 1:chars_per_row) + # Show a tick for each character position
  labs(title = "CFTR cDNA", x = "", y = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  scale_y_continuous(
    breaks = unique(plot_data$y), # Show breaks for each row
    labels = abs(unique(plot_data$y)) # Make labels positive row numbers
  )
