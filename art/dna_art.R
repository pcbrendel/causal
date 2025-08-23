library(rentrez)
library(svglite)
library(ggplot2)
library(stringr)

# Get sequence ----
# https://www.ncbi.nlm.nih.gov/nuccore/NM_000492.4
# entrez = a search engine and data retrieval system
# that allows users to search across dozens of NCBI databases
# nuccore = "Nucleotide Core" database at NCBI
# (National Center for Biotechnology Information)
# accession = unique, stable identifier for a sequence

# this is the coding DNA sequence
# it includes:
# - 5' UTR (non-coding before the CDS)
# - CDS (coding sequence)
# - 3' UTR (non-coding after the CDS)

scrape_dna <- function(id) {
  # Fetch the FASTA sequence directly using the accession ID
  message(
    paste("Attempting to fetch mRNA sequence for accession:", id)
  )
  fasta_record <- tryCatch(
    {
      entrez_fetch(
        db = "nuccore",
        id = id,
        rettype = "fasta",
        retmode = "text"
      )
    },
    error = function(e) {
      message(paste("Error fetching sequence:", e$message))
      return(NULL)
    }
  )

  if (!is.null(fasta_record)) {
    # Parse the FASTA record
    lines <- strsplit(fasta_record, "\n")[[1]]
    # Remove the header line and concatenate the sequence lines
    dna_sequence <- paste(lines[-1], collapse = "")

    message("\nDNA (mRNA/cDNA) Sequence (first 100 characters):\n")
    print(substring(dna_sequence, 1, 100))
    message(paste("\nTotal length:", nchar(dna_sequence), "base pairs"))
    return(dna_sequence)
  } else {
    message("Failed to retrieve cDNA sequence using direct accession.")
  }
}

# example
# cftr_id <- "NM_000492.4"
# test <- scrape_dna(cftr_id)

# Get CDS ----

get_cds <- function(id) {
  # Extract CDS coordinates from GenBank record
  message("Fetching CDS coordinates from GenBank record...")
  gb_record <- entrez_fetch(
    db = "nuccore",
    id = id,
    rettype = "gb",
    retmode = "text"
  )

  # Initialize default values
  cds_start <- NA
  cds_end <- NA

  # Parse CDS coordinates from GenBank format
  cds_lines <- grep(
    "^\\s+CDS\\s+",
    strsplit(gb_record, "\n")[[1]],
    value = TRUE
  )

  if (length(cds_lines) > 0) {
    # Extract coordinates from the first CDS line
    cds_line <- cds_lines[1]
    message(paste("Found CDS line:", cds_line))

    # Extract numbers from format like "join(71..4513)" or "71..4513"
    # Use a more robust regex pattern
    coords_pattern <- "(\\d+)\\.\\.(\\d+)"
    coords_match <- regexpr(coords_pattern, cds_line)

    if (coords_match != -1) {
      # Extract the matched substring
      matched_text <- substr(
        cds_line,
        coords_match,
        coords_match + attr(coords_match, "match.length") - 1
      )

      # Parse the coordinates using str_match
      coords <- str_match(matched_text, coords_pattern)
      if (!is.na(coords[1, 2]) && !is.na(coords[1, 3])) {
        cds_start <- as.numeric(coords[1, 2])
        cds_end <- as.numeric(coords[1, 3])
        message(
          paste(
            "CDS coordinates from GenBank:", cds_start, "to", cds_end
          )
        )
      } else {
        message("Could not parse coordinates from matched text")
      }
    } else {
      message("Could not find coordinate pattern in CDS line")
    }
  } else {
    message("No CDS lines found in GenBank record")
  }

  return(c(cds_start, cds_end))
}

# example
# cftr_id <- "NM_000492.4"
# test <- get_cds(cftr_id)

# Apply deletion ----

# Letter plot ----

# Dot plot ----

dot_plot <- function(
  dna_seq, svg_file_name, n_per_row, colors
) {
  n <- length(dna_seq)
  x <- rep(1:n_per_row, length.out = n)
  y <- rep(seq(1, ceiling(n / n_per_row)), each = n_per_row)[1:n]

  plot_data <- data.frame(
    x = x,
    y = -y,
    base = dna_seq,
    color = colors[dna_seq]
  )

  # Remove rows with "." (deletion sites) - no circles for these
  plot_data <- plot_data[!is.na(plot_data$color), ]

  # fallback for any unexpected chars
  plot_data$color[is.na(plot_data$color)] <- "black"

  # Split by color
  color_groups <- split(plot_data, plot_data$color)

  # Create SVG file
  svglite(file = svg_file_name, width = 12, height = 12)

  par(mar = c(1, 1, 2, 1))

  # Initialize empty plot
  plot(
    NULL,
    xlim = range(plot_data$x),
    ylim = range(plot_data$y),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )

  # Plot each color group separately to group them in SVG
  for (color_name in names(color_groups)) {
    group_data <- color_groups[[color_name]]
    points(
      group_data$x,
      group_data$y,
      col = color_name,
      pch = 1, # Open circles with strokes (pen plotter compatible)
      cex = 1.2,
      lwd = 1.5 # Stroke width for better visibility
    )
  }

  dev.off()
}

# example
# dot_plot(
#   dna_seq = delta_f508_vec,
#   svg_file_name = "art/cftr_dna_dots_deltaf508.svg",
#   n_per_row = 80,
#   colors = c(
#     A = "forestgreen",
#     T = "firebrick",
#     C = "royalblue",
#     G = "goldenrod"
#   )
# )

# old ----

# svg_file_name <- "art/cftr_dna_test.svg"
# svglite(file = svg_file_name, width = 12, height = 12)

# par(mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0), ann = FALSE, bty = "n")

# initialize
# plot(
#   NULL,
#   xlim = c(0, 1),
#   ylim = c(0, 1),
#   axes = FALSE,
#   xlab = "",
#   ylab = "",
#   main = "",
#   sub = "",
#   frame.plot = FALSE
# )

# lines_to_plot <- unlist(strsplit(cftr_dna_sequence_wrapped, "\n"))

# text_cex <- 0.8
# vfont_type <- c("serif", "plain")
# line_height <- strheight("Mg", cex = text_cex, font = 1) * 1.5
# char_width <- 0.01

# current_y <- 0.995
# start_x <- 0.005

# for (line_idx in seq_along(lines_to_plot)) {
#   line_text <- lines_to_plot[line_idx]
#   current_x <- start_x

#   for (char_idx in seq_along(strsplit(line_text, "")[[1]])) {
#     char_to_plot <- substr(line_text, char_idx, char_idx)

#     text(
#       x = current_x,
#       y = current_y,
#       labels = char_to_plot,
#       vfont = vfont_type,
#       cex = text_cex,
#       adj = c(0, 0.5)
#     )

#     current_x <- current_x + char_width
#   }

#   current_y <- current_y - line_height
# }

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

# dev.off()


# my_string <- substr(cftr_dna_sequence, 1, 200)
# my_string <- cftr_dna_sequence

# chars_per_row <- 100

# # 2. Split the string into individual characters
# all_chars <- unlist(strsplit(my_string, split = ""))
# n_chars <- length(all_chars)

# # 3. Create indices for columns and rows
# # Column index repeats 1 to chars_per_row
# col_indices <- rep(1:chars_per_row, ceiling(n_chars / chars_per_row))[1:n_chars]

# # Row index increments every chars_per_row characters
# row_indices <- ceiling((1:n_chars) / chars_per_row)

# # 4. Create a data frame for plotting
# plot_data <- data.frame(
#   x = col_indices,
#   y = -row_indices, # Use negative rows so the plot starts from the top (row 1)
#   char = all_chars
# )

# # 5. Plot the characters
# ggplot(plot_data, aes(x = x, y = y, label = char)) +
#   geom_text(size = 4, family = "mono") + # Use a monospace font for consistent character width
#   scale_x_continuous(breaks = 1:chars_per_row) + # Show a tick for each character position
#   labs(title = "CFTR cDNA", x = "", y = "") +
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(), # Remove major grid lines
#     panel.grid.minor = element_blank(), # Remove minor grid lines
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
#   ) +
#   scale_y_continuous(
#     breaks = unique(plot_data$y), # Show breaks for each row
#     labels = abs(unique(plot_data$y)) # Make labels positive row numbers
#   )
