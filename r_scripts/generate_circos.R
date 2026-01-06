suppressPackageStartupMessages({
  library(circlize)
  library(dplyr)
  library(RColorBrewer)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) stop("Not enough arguments. Need: links_file sectors_file annots_file output_file")

links_file <- args[1]
sectors_file <- args[2]
annots_file <- args[3]
output_file <- args[4]

# Load Data
links_df <- read_csv(links_file, show_col_types = FALSE)
sector_df <- read_csv(sectors_file, show_col_types = FALSE)
annotations_bed <- read_csv(annots_file, show_col_types = FALSE)

# Sanity check: Verify required columns exist
if (!("P1_clean" %in% colnames(links_df))) {
  stop("Error: Input CSV is missing 'P1_clean' column.")
}
if (!("P2_clean" %in% colnames(links_df))) {
  stop("Error: Input CSV is missing 'P2_clean' column.")
}
print(paste("Plotting", nrow(links_df), "links..."))
print(paste("Sectors:", nrow(sector_df), "proteins"))
print(paste("Unique P1_clean values:", length(unique(links_df$P1_clean))))
print(paste("Unique P2_clean values:", length(unique(links_df$P2_clean))))
print(paste("Sector names:", paste(sector_df$name, collapse=", ")))

# --- 1. Define Color Palette (Exact Match) ---
fixed_palette <- c("SMARCA4"="#558B2F", "SMCA4"="#558B2F", "BRG1"="#558B2F",
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
                   "ACTB"="#2E8B57", "ACTL6A"="#3CB371", "ACL6A"="#3CB371", "ACTL6B"="#90EE90",
                   "DPF1"="#C71585", "DPF2"="#DB7093", "DPF3"="#FF69B4",
                   "PBRM1"="#663399", "PHF10"="#8A2BE2", "BRD7"="#9400D3",
                   "BICRA"="#FF4500", "BICRAL"="#FF6347", "BRD9"="#FF7F50")

# Handle dynamic colors for unknown proteins
all_proteins <- unique(sector_df$name)
other_proteins <- all_proteins[!all_proteins %in% names(fixed_palette)]

if (length(other_proteins) > 0) {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  color_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if (length(color_vector) < length(other_proteins)) { 
    color_vector <- rep(color_vector, length.out = length(other_proteins)) 
  }
  new_colors <- setNames(color_vector[1:length(other_proteins)], other_proteins)
  color_palette <- c(fixed_palette, new_colors)
} else {
  color_palette <- fixed_palette
}

# --- 2. Generate Plot ---
png(filename = output_file, width = 2000, height = 2000, res = 300)

tryCatch({
  circos.clear()
  circos.par(start.degree = 90, gap.after = setNames(rep(0.5, nrow(sector_df)), sector_df$name), track.margin = c(1e-05, 1e-05))
  
  # Initialize
  circos.initialize(factors = sector_df$name, xlim = sector_df[, c("start", "end")])
  
  # Track 1: Axis
  circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(x, y) { 
    circos.axis(h = "top", major.at = seq(0, 5000, by = 100), labels.cex = 0.6) 
  })
  
  # Track 2: Annotations (Conditional)
  if (nrow(annotations_bed) > 0) {
    circos.genomicTrack(annotations_bed, ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(region, value, ...) {
      sector.name <- get.cell.meta.data("sector.index")
      col_val <- color_palette[sector.name]
      if(is.na(col_val)) col_val <- "grey"
      bright_color <- adjustcolor(col_val, alpha.f = 0.4)
      circos.genomicRect(region, value, col = bright_color, border = "grey60") 
    })
  }
  
  # Track 3: Labels & Colored Bands
  circos.track(ylim = c(0, 1), track.height = 0.1, bg.border = NA, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    col_val <- color_palette[sector.name]
    if(is.na(col_val)) col_val <- "grey"
    circos.rect(get.cell.meta.data("xlim")[1], 0, get.cell.meta.data("xlim")[2], 1, col = col_val, border = NA)
    circos.text(mean(get.cell.meta.data("xlim")), 0.5, sector.name, col = "white", cex = 0.8, facing = "inside", niceFacing = TRUE) 
  })
  
  # Draw Links
  if (nrow(links_df) > 0) {
    # Calculate transparency logic from R code
    shading_values <- if ("NumPSMs" %in% names(links_df)) links_df$NumPSMs else links_df$Score
    min_val <- min(shading_values, na.rm = TRUE)
    max_val <- max(shading_values, na.rm = TRUE)
    
    # Calculate alpha
    if (is.infinite(min_val) || is.infinite(max_val)) { 
      alpha_vals <- rep(0.1, nrow(links_df))
    } else if (max_val > min_val) { 
      alpha_vals <- 0.1 + 0.9 * ((shading_values - min_val)/(max_val - min_val)) 
    } else { 
      alpha_vals <- rep(0.8, nrow(links_df))
    }
    
    track_inner_radius <- get.cell.meta.data("cell.bottom.radius", track.index = ifelse(nrow(annotations_bed) > 0, 3, 2))
    
    for (i in 1:nrow(links_df)) {
      p1 <- links_df$P1_clean[i]
      p2 <- links_df$P2_clean[i]
      # Color based on Protein 1
      col_base <- color_palette[p1]
      if(is.na(col_base)) col_base <- "#808080"
      
      link_col <- adjustcolor(col_base, alpha.f = alpha_vals[i])
      
      circos.link(p1, links_df$LinkPos1[i], 
                  p2, links_df$LinkPos2[i], 
                  rou = track_inner_radius, col = link_col, lwd = 2, border = NA)
    }
  }
  
}, error = function(e) {
  message("Error in plotting:")
  message(e)
})

dev.off()

