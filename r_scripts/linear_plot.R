suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggiraph)
  library(htmlwidgets)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5) stop("Not enough arguments. Need: links_file lengths_file annots_file output_file show_intra_links")

links_file <- args[1]
lengths_file <- args[2]
annots_file <- args[3]
output_file <- args[4]
show_intra_links <- as.logical(args[5])
pandoc_path <- if (length(args) >= 6) args[6] else NULL

# Set pandoc path if provided
if (!is.null(pandoc_path) && pandoc_path != "" && file.exists(pandoc_path)) {
  # Add pandoc directory to PATH
  pandoc_dir <- dirname(pandoc_path)
  current_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste0(pandoc_dir, ":", current_path))
  
  # Set RSTUDIO_PANDOC (used by rmarkdown/htmlwidgets)
  Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
  
  print(paste("Using pandoc at:", pandoc_path))
} else {
  # Try to find pandoc in common locations
  common_paths <- c("/usr/bin/pandoc", "/usr/local/bin/pandoc", "/opt/conda/bin/pandoc", "/home/adminuser/.conda/bin/pandoc")
  found_pandoc <- NULL
  for (path in common_paths) {
    if (file.exists(path)) {
      found_pandoc <- path
      pandoc_dir <- dirname(path)
      Sys.setenv(PATH = paste0(pandoc_dir, ":", Sys.getenv("PATH")))
      Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
      print(paste("Found pandoc at:", path))
      break
    }
  }
  if (is.null(found_pandoc)) {
    print("Warning: pandoc path not found. HTML widget generation may fail.")
  }
}

# Load Data
links_df <- read_csv(links_file, show_col_types = FALSE)
sector_df <- read_csv(lengths_file, show_col_types = FALSE)
annotations_bed <- read_csv(annots_file, show_col_types = FALSE)

# Sanity checks
if (!("P1_clean" %in% colnames(links_df))) {
  stop("Error: Input CSV is missing 'P1_clean' column.")
}
if (!("P2_clean" %in% colnames(links_df))) {
  stop("Error: Input CSV is missing 'P2_clean' column.")
}
if (!("name" %in% colnames(sector_df)) || !("end" %in% colnames(sector_df))) {
  stop("Error: lengths_file must have 'name' and 'end' columns.")
}

# Determine score column for alpha mapping
if ("NumPSMs" %in% colnames(links_df)) {
  links_df <- links_df %>% mutate(AlphaValue = NumPSMs)
} else if ("Score" %in% colnames(links_df)) {
  links_df <- links_df %>% mutate(AlphaValue = Score)
} else {
  links_df <- links_df %>% mutate(AlphaValue = 1.0)
}

# Ensure start column exists (default to 0 if missing)
if (!("start" %in% colnames(sector_df))) {
  sector_df$start <- 0
}

# Order proteins by length (descending)
protein_order <- sector_df %>% arrange(desc(end)) %>% pull(name)
plot_df <- sector_df %>% 
  mutate(name = factor(name, levels = protein_order)) %>%
  select(name, start, end)  # Ensure we have the required columns

# Process annotations if available
annot_df <- NULL
if (nrow(annotations_bed) > 0 && "chr" %in% colnames(annotations_bed)) {
  annot_df <- annotations_bed %>% mutate(chr = factor(chr, levels = protein_order))
  annot_df <- annot_df %>% filter(!is.na(chr))  # Remove annotations for proteins not in plot
}

# Process links
links_df_processed <- links_df %>% 
  mutate(
    P1_clean = factor(P1_clean, levels = protein_order),
    P2_clean = factor(P2_clean, levels = protein_order)
  ) %>%
  filter(!is.na(P1_clean), !is.na(P2_clean))  # Remove links for proteins not in plot

# Separate inter-links and intra-links
inter_links <- links_df_processed %>% 
  filter(as.character(P1_clean) != as.character(P2_clean))

intra_links <- links_df_processed %>% 
  filter(as.character(P1_clean) == as.character(P2_clean), LinkPos1 != LinkPos2)

# Calculate max length for x-axis
max_len <- max(plot_df$end, na.rm = TRUE)

# Define color palette (using fixed palette from original R code)
fixed_palette <- c(
  "SMARCA4"="#558B2F", "SMCA4"="#558B2F", "BRG1"="#558B2F",
  "SMARCA2"="#6B8E23", "SMCA2"="#6B8E23", "BRM"="#6B8E23",
  "ARID1A"="#A44A4A", "ARI1A"="#A44A4A",
  "ARID1B"="#CD5C5C", "ARI1B"="#CD5C5C",
  "ARID2"="#F08080",
  "BCL7A"="#483D8B", "BCL7B"="#6A5ACD", "BCL7C"="#7B68EE",
  "SMARCB1"="#4682B4", "SNF5"="#4682B4",
  "SMARCC1"="#5F9EA0", "SMRC1"="#5F9EA0", "BAF155"="#5F9EA0",
  "SMARCC2"="#87CEEB", "SMRC2"="#87CEEB", "BAF170"="#87CEEB",
  "SMARCD1"="#B8860B", "SMRD1"="#B8860B",
  "SMARCD2"="#DAA520", "SMRD2"="#DAA520",
  "SMARCD3"="#FFD700", "SMRD3"="#FFD700",
  "SMARCE1"="#8B4513", "SMCE1"="#8B4513", "BAF57"="#8B4513",
  "SS18"="#A0522D", "SS18L1"="#D2691E",
  "ACTB"="#2E8B57",
  "ACTL6A"="#3CB371", "ACL6A"="#3CB371",
  "ACTL6B"="#90EE90",
  "DPF1"="#C71585", "DPF2"="#DB7093", "DPF3"="#FF69B4",
  "PBRM1"="#663399",
  "PHF10"="#8A2BE2",
  "BRD7"="#9400D3",
  "BICRA"="#FF4500", "BICRAL"="#FF6347",
  "BRD9"="#FF7F50"
)

# Extend palette for any missing proteins
all_proteins <- unique(c(as.character(plot_df$name), as.character(links_df_processed$P1_clean), as.character(links_df_processed$P2_clean)))
missing_proteins <- setdiff(all_proteins, names(fixed_palette))
if (length(missing_proteins) > 0) {
  # Use a default color scheme for missing proteins
  default_colors <- rainbow(length(missing_proteins))
  names(default_colors) <- missing_proteins
  color_palette <- c(fixed_palette, default_colors)
} else {
  color_palette <- fixed_palette
}

# Determine if output is HTML (interactive) or static (PNG/PDF)
is_html_output <- grepl("\\.html$", output_file, ignore.case = TRUE)

# Build the plot - use interactive geoms for HTML, static for PNG/PDF
if (is_html_output) {
  # Interactive version with tooltips
  p <- ggplot() +
    # Protein tracks (black rectangles) - interactive
    geom_rect_interactive(data = plot_df,
             aes(xmin = start, xmax = end, 
                 ymin = as.numeric(name) - 0.2, ymax = as.numeric(name) + 0.2,
                 tooltip = as.character(name), data_id = as.character(name)),
             fill = "black") +
    
    # Annotations (colored blocks on top of tracks) - interactive
    {
      if (!is.null(annot_df) && nrow(annot_df) > 0) {
        geom_rect_interactive(data = annot_df,
                  aes(xmin = start, xmax = end,
                      ymin = as.numeric(chr) - 0.2, ymax = as.numeric(chr) + 0.2,
                      fill = chr,
                      tooltip = paste0(as.character(chr), ": ", start, " - ", end, " (", name, ")"),
                      data_id = paste(as.character(chr), start, end)),
                  alpha = 0.6)
      } else {
        NULL
      }
    } +
    
    # Inter-links (straight lines between different proteins) - interactive
    {
      if (nrow(inter_links) > 0) {
        geom_segment_interactive(data = inter_links,
                     aes(x = LinkPos1, xend = LinkPos2,
                         y = as.numeric(P1_clean) + 0.2, yend = as.numeric(P2_clean) + 0.2,
                         color = P1_clean,
                         alpha = AlphaValue,
                         tooltip = paste0(as.character(P1_clean), "(", LinkPos1, ") - ", as.character(P2_clean), "(", LinkPos2, ")"),
                         data_id = paste0(as.character(P1_clean), LinkPos1, as.character(P2_clean), LinkPos2)),
                     linewidth = 0.8)
      } else {
        NULL
      }
    } +
    
    # Intra-links (curved arcs for same protein, if enabled) - interactive
    {
      if (show_intra_links && nrow(intra_links) > 0) {
        geom_curve_interactive(data = intra_links,
                   aes(x = LinkPos1, xend = LinkPos2,
                       y = as.numeric(P1_clean) + 0.2, yend = as.numeric(P1_clean) + 0.2,
                       color = P1_clean,
                       alpha = AlphaValue,
                       tooltip = paste0(as.character(P1_clean), " (", LinkPos1, " - ", LinkPos2, ")"),
                       data_id = paste0(as.character(P1_clean), LinkPos1, as.character(P1_clean), LinkPos2)),
                   curvature = -0.4,
                   linewidth = 0.8)
      } else {
        NULL
      }
    }
} else {
  # Static version (PNG/PDF)
  p <- ggplot() +
    # Protein tracks (black rectangles)
    geom_rect(data = plot_df,
             aes(xmin = start, xmax = end, 
                 ymin = as.numeric(name) - 0.2, ymax = as.numeric(name) + 0.2),
             fill = "black") +
    
    # Annotations (colored blocks on top of tracks)
    {
      if (!is.null(annot_df) && nrow(annot_df) > 0) {
        geom_rect(data = annot_df,
                aes(xmin = start, xmax = end,
                    ymin = as.numeric(chr) - 0.2, ymax = as.numeric(chr) + 0.2,
                    fill = chr),
                alpha = 0.6)
      } else {
        NULL
      }
    } +
    
    # Inter-links (straight lines between different proteins)
    {
      if (nrow(inter_links) > 0) {
        geom_segment(data = inter_links,
                   aes(x = LinkPos1, xend = LinkPos2,
                       y = as.numeric(P1_clean) + 0.2, yend = as.numeric(P2_clean) + 0.2,
                       color = P1_clean,
                       alpha = AlphaValue),
                   linewidth = 0.8)
      } else {
        NULL
      }
    } +
    
    # Intra-links (curved arcs for same protein, if enabled)
    {
      if (show_intra_links && nrow(intra_links) > 0) {
        geom_curve(data = intra_links,
                 aes(x = LinkPos1, xend = LinkPos2,
                     y = as.numeric(P1_clean) + 0.2, yend = as.numeric(P1_clean) + 0.2,
                     color = P1_clean,
                     alpha = AlphaValue),
                 curvature = -0.4,
                 linewidth = 0.8)
      } else {
        NULL
      }
    }
}

# Add scales and theme (common for both versions)
p <- p +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  scale_alpha_continuous(name = "Score/PSMs", range = c(0.25, 1.0)) +
  scale_y_continuous(breaks = 1:length(protein_order), labels = protein_order, expand = c(0.1, 0.1)) +
  scale_x_continuous(name = "Amino Acid Position", limits = c(0, max_len), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.y = element_text(size = 10, hjust = 1),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

# Save the plot
if (is_html_output) {
  # Create interactive widget
  widget <- girafe(ggobj = p, width_svg = 12, height_svg = 8,
                   options = list(
                     opts_tooltip(css = "background-color:wheat;color:black;padding:5px;border-radius:5px;"),
                     opts_hover(css = "stroke:black;stroke-width:3px;")
                   ))
  # Save as HTML widget
  saveWidget(widget, file = output_file, selfcontained = TRUE)
  print(paste("Interactive linear plot saved to:", output_file))
} else {
  # Save as static image
  ggsave(output_file, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  print(paste("Linear plot saved to:", output_file))
}

