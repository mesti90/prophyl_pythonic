# Load required libraries
library(devtools)
library(treeio)
library(ggplot2)
library(ggtree)
library(dplyr)
library(tidyr)
library(optparse)
library(rlang)

# Load external scripts
script_dir <- "."
source(file.path(script_dir, "custom_plot_tree_long.R"))
source(file.path(script_dir, "get_colors.R"))
source(file.path(script_dir, "read_df.R"))

#' Parse collection_date column into full Date objects
#'
#' Adds a new column `collection_day` to `tree_tbl` with proper Date class entries.
#' Handles three formats: yyyy-mm-dd, yyyy-mm, and yyyy.
#'
#' @param tree_tbl A tibble with a character column `collection_date`.
#' @param seed Integer; optional random seed for reproducibility.
#' @return The same tibble, with a new column `collection_day` of class Date.
parse_collection_dates <- function(tree_tbl, seed = 598234) {
	set.seed(seed)

	tree_tbl$collection_day <- sapply(tree_tbl$Date, function(date_str) {
		if (is.na(date_str)) return(NA)
		if (grepl("^\\d{4}-\\d{2}-\\d{2}$", date_str)) return(as.Date(date_str))
		if (grepl("^\\d{4}-\\d{2}$", date_str)) {
			first_day <- as.Date(paste0(date_str, "-01"))
			last_day <- seq(first_day, length = 2, by = "1 month")[2] - 1
			return(sample(seq(first_day, last_day, by = "day"), 1))
		}
		if (grepl("^\\d{4}$", date_str)) {
			year <- as.numeric(date_str)
			first_day <- as.Date(paste0(year, "-01-01"))
			last_day <- as.Date(paste0(year, "-12-31"))
			return(sample(seq(first_day, last_day, by = "day"), 1))
		}
	return(NA)
	})
	
	tree_tbl$collection_day <- as.Date(tree_tbl$collection_day, origin = "1970-01-01")
	
	return(tree_tbl)
}

#' Generate the heatmap plot for a phylogenetic tree
#'
#' Reads a tree and metadata, processes inputs, and saves a PDF plot.
#' @param tree_file Path to the Newick tree file.
#' @param heatmap_file Path to the metadata TSV file.
#' @param output_file Path where the PDF output will be saved.
#' @return NULL (plot is saved to file)
generate_tree_plot <- function(tree_file, heatmap_file, output_file, heatmap_vars) {
	# Read tree and convert to tibble
	tree_tbl <- as_tibble(read.tree(tree_file))
	
	# Read metadata and join
	gentypes <- read.delim(heatmap_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	gentypes$highlight <- FALSE
	tree_tbl <- left_join(tree_tbl, gentypes, by = c("label" = "strain"))
	
	# Replace NA values in heatmap vars
	for (col in heatmap_vars) tree_tbl[[col]][is.na(tree_tbl[[col]])] <- "NA" 

	# Prepare country color map
	country_colors <- data.frame(country = sort(unique(tree_tbl$country)))
	country_colors <- get_colors(country_colors, var = "country")

	# Extract and clean collection dates
	tree_tbl <- parse_collection_dates(tree_tbl)
	dates <- as.Date(tree_tbl$collection_day)

	# Generate tree plot
	g <- custom_plot_tree_long(
		tree_tbl,
		heatmap_opts = list(var = heatmap_vars, colors = list(country = country_colors), offset = 40),
		mrsd = max(dates, na.rm = TRUE),
		linewidth = 0.1,
		verbose = TRUE
	)

	# --- Get tip order (handle patchwork output) ---
	tree_part <- if ("patchwork" %in% class(g)) g[[1]] else g
	df <- tree_part$data
	tip_order <- df$label[df$isTip][order(df$y[df$isTip], decreasing = TRUE)]

	# --- Reorder the input data by tip order ---
	gentypes_ordered <- gentypes[match(tip_order, gentypes$strain), ]
	
	# Print (or save) the reordered table
	message("Reordered metadata (first few rows):")
	#print(head(gentypes_ordered))
	
	# Optionally save the reordered file alongside the PDF
	tip_table_file <- sub("\\.pdf$", "_tips.tsv", output_file)
	write.table(gentypes_ordered, tip_table_file, sep = "\t", quote = FALSE, row.names = FALSE)
	message(paste0("Tip-ordered metadata written to: ", tip_table_file))

	# --- Save the tree plot ---
	ggsave(
		filename = output_file,
		plot = g,
		height = 0.05 * nrow(tree_tbl),
		width = 50,
		limitsize = FALSE
	)
	
	message(paste0("Ready with ", output_file))
}



if (interactive() || length(commandArgs(trailingOnly = TRUE)) > 0) {
	option_list <- list(
		make_option(c("-t", "--tree"), type = "character", help = "Path to Newick tree file", metavar = "FILE"),
		make_option(c("-m", "--meta"), type = "character", help = "Path to metadata TSV file", metavar = "FILE"),
		make_option(c("-o", "--out"), type = "character", default = "tree_with_heatmap.pdf", help = "Output PDF file [default: %default]", metavar = "FILE"),
		make_option(c("-c", "--columns"), type = "character", help = "Comma-separated list of metadata column names to include", metavar = "col1,col2,...", default="H.type,O.type,country,Source,clb,fimH.N70S,fimH.S78N")
	) # OR: "H.type,O.type,country,Source,Titer.E2COLE7,EOP.E2COLE7"


	args <- parse_args(OptionParser(option_list = option_list))
	
	tree_file <- args$tree %||% "trees/prophyl_ST1193.final_tree.nwk"
	heatmap_file <- args$meta %||% "input/ST1193.tsv"
	output_file <- args$out
	heatmap_vars <- unlist(strsplit(args$columns, ","))
	
	generate_tree_plot(tree_file, heatmap_file, output_file, heatmap_vars)
}



##tree_file <- "work/ST105.final_dated_tree.nwk"
##heatmap_file <- "tree_input/ST105.tsv"
##output_file <- "work/ST105.tree.pdf"
##heatmap_vars <- "spatyper,Capsule.type,country"
##heatmap_vars <- unlist(strsplit(heatmap_vars, ","))
##
##generate_tree_plot(tree_file, heatmap_file, output_file, heatmap_vars)
##
