#' Validate user inputs for custom tree plotting
#'
#' Ensures that input types and lengths are appropriate for the tree and heatmap plotting.
#'
#' @param tree_tbl A tibble with class "tbl_tree", representing a phylogenetic tree.
#' @param mrsd A single Date object indicating the most recent sampling date.
#' @param heatmap_width A numeric scalar or vector specifying heatmap widths.
#' @param heatmap_var A character vector of heatmap variable names.
#'
#' @return Throws an error if validation fails; otherwise, returns invisibly.
validate_inputs <- function(tree_tbl, mrsd, heatmap_width, heatmap_var) {
  # Validate tree_tbl class
  if (!inherits(tree_tbl, "tbl_tree")) {
    stop("'tree_tbl' must be a tibble with class 'tbl_tree'.")
  }

  # Validate mrsd
  if (!is.null(mrsd)) {
    if (!inherits(mrsd, "Date")) {
      stop("Argument 'mrsd' must be of class 'Date'.")
    }
    if (length(mrsd) != 1 || is.na(mrsd)) {
      stop("Argument 'mrsd' must be a single, non-NA Date.")
    }
  }

  # Validate heatmap_width if set
  if (!is.null(heatmap_width)) {
    if (length(heatmap_width) > 1 && length(heatmap_width) != length(heatmap_var)) {
      stop("'heatmap_width' must have length 1 or match length of 'heatmap_var'.")
    }
  }
}

#' Estimate default heatmap widths based on label length
#'
#' Calculates a suitable heatmap width for each variable by examining the maximum string length
#' of values in each corresponding column. Width is scaled so that approximately 4 characters
#' fit per unit, with a minimum width of 2.
#'
#' @param tree_tbl A tibble containing the heatmap variable columns.
#' @param heatmap_var A character vector of column names to use for heatmaps.
#'
#' @return A numeric vector of estimated widths, one per variable.
estimate_heatmap_width <- function(tree_tbl, heatmap_var, tiplab_size) {
  sapply(heatmap_var, function(var) {
    max_len <- max(nchar(as.character(tree_tbl[[var]])), na.rm = TRUE)
    scale_factor = tiplab_size / 30
    max(0.15, max_len * scale_factor)
  })
}


#' Determine if internal nodes should be highlighted
#'
#' Checks whether all internal nodes have non-NA values for the highlight variable.
#'
#' @param tree_tbl A "tbl_tree" tibble.
#' @param highlight_var Character string for the variable used in highlighting.
#'
#' @return Logical value; TRUE if all internal nodes can be highlighted, FALSE otherwise.
should_highlight_ancestors <- function(tree_tbl, highlight_var) {
  root_index <- which(is.na(tree_tbl$branch.length))
  if (length(root_index) != 1) {
    stop("Could not find tree root.")
  }

  all(!is.na(tree_tbl[[highlight_var]][root_index:nrow(tree_tbl)]))
}


#' Build base ggtree plot for a phylogenetic tree
#'
#' Constructs a ggtree object from a tree in treedata format, optionally colored by a highlight
#' variable, and optionally labeled with tip text. Supports internal node coloring and temporal layout.
#'
#' @param tree_db A tree object in `treeio::treedata` format.
#' @param tree_tbl A tibble with tree metadata (must match `tree_db`).
#' @param tree A phylo object derived from the same data.
#' @param show_tiplab Logical; whether to show tip labels.
#' @param tiplab_size Numeric; size of tip labels.
#' @param tiplab_offset Numeric; horizontal offset of tip labels.
#' @param highlight_var Character; column name for node coloring (optional).
#' @param highlight_ancestors Logical; whether to color internal nodes.
#' @param highlight_colors Data frame of custom colors (optional).
#' @param mrsd A Date object for time scaling (optional).
#' @param legend_show Character vector of variables to show in the legend.
#' @param linewidth Numeric; width of tree lines.
#'
#' @return A `ggtree` plot object.
build_base_tree_plot <- function(tree_db, tree_tbl, tree, 
                                 show_tiplab, tiplab_size, tiplab_offset,
                                 highlight_var, highlight_ancestors,
                                 highlight_colors, mrsd,
                                 legend_show, linewidth) {
  options(ignore.negative.edge = TRUE)

  # Determine aesthetic mapping
  mapping <- if (!is.null(highlight_var) && highlight_ancestors) {
    aes(color = .data[[highlight_var]])
  } else {
    NULL
  }

  # Build ggtree with/without mrsd
  p <- if (is.null(mrsd)) {
    ggtree(tree_db, mapping = mapping, size = linewidth)
  } else {
    ggtree(tree_db, mapping = mapping, size = linewidth, mrsd = mrsd) + 
    theme_tree2()
  }
  
  
  p <- p + theme(axis.ticks.x = element_line(linewidth = 1))

  # Hide legend if not explicitly requested
  #if (!is.null(highlight_var) && !highlight_var %in% legend_show) {
    p <- p + guides(color = "none")
  #}

  # Add tip labels if requested
  if (show_tiplab) {
    p <- add_tip_labels(p, highlight_var, highlight_colors, tiplab_size, tiplab_offset)
  }
  return(p + coord_cartesian(clip = "off"))
}



#' Add tip labels to a ggtree plot
#'
#' Adds aligned tip labels to a ggtree object, optionally colored by a highlight variable.
#' If a highlight variable is provided and `highlight_colors` is specified, it uses
#' manual colors; otherwise, default colors are applied.
#'
#' @param p A ggtree plot object.
#' @param highlight_var Character; name of the variable used for coloring tip labels. Can be NULL.
#' @param highlight_colors Data frame with `group` and `color` columns (optional).
#' @param tiplab_size Numeric; text size for tip labels.
#' @param tiplab_offset Numeric; horizontal offset for tip labels.
#'
#' @return A ggtree plot with tip labels added.
add_tip_labels <- function(p, highlight_var, highlight_colors, tiplab_size, tiplab_offset) {
  label_args <- list(
    align = TRUE,
    linesize = 0.1,
    geom = "text",
    hjust = 0
  )

  if (is.null(highlight_var)) {
    label_args$size <- tiplab_size
    p <- p + do.call(geom_tiplab, label_args)
  } else {
    label_args$aes_params <- aes(col = .data[[highlight_var]])
    label_args$size <- tiplab_size
    label_args$offset <- tiplab_offset
    label_args$width <- 10

    p <- p + do.call(geom_tiplab, label_args)

    if (!is.null(highlight_colors)) {
      hlt_colors <- highlight_colors$color
      names(hlt_colors) <- highlight_colors$group
      p <- p + scale_color_manual(values = hlt_colors)
    }
  }

  return(p)
}


#' Annotate ambiguous and state-change nodes in a ggtree plot
#'
#' Adds visual annotations for ambiguous internal nodes and nodes where the state
#' changes from parent to child, based on a highlight variable. Ambiguous nodes
#' are detected using the presence of "|" in the highlight value.
#'
#' @param p A ggtree plot object.
#' @param tree_tbl A "tbl_tree" tibble containing the tree structure and metadata.
#' @param highlight_var Character; name of the column used for state annotations.
#' @param tiplab_size Numeric; font size for annotation labels.
#'
#' @return A ggtree plot with point and label annotations added.
annotate_nodes_with_states <- function(p, tree_tbl, highlight_var, tiplab_size) {
  # Ambiguous nodes: values containing "|"
  index_ambiguous_nodes <- grep("\\|", tree_tbl[[highlight_var]])

  # Map nodes to parents
  node_to_parent <- match(tree_tbl$parent, tree_tbl$node)
  highlight_vals <- tree_tbl[[highlight_var]]

  # State-change: child ≠ parent
  state_changed <- highlight_vals != highlight_vals[node_to_parent]
  state_changed[is.na(state_changed)] <- FALSE
  index_state_change_nodes <- which(state_changed)

  # Build state-change labels
  collection_dates <- as.Date(tree_tbl$collection_day)
  sclabel <- rep(NA_character_, nrow(tree_tbl))
  for (i in index_state_change_nodes) {
    parent_idx <- node_to_parent[i]
    if (!is.na(parent_idx)) {
      sclabel[i] <- paste0(
        "CHANGE: \n",
        collection_dates[parent_idx], "/", collection_dates[i]
      )
    }
  }

  # Add annotations
  p <- p +
    geom_point2(
      aes(subset = node %in% index_ambiguous_nodes),
      shape = 21, size = 10, fill = "orange", alpha = 0.5
    ) +
    geom_label2(
      aes(
        x = branch,
        subset = node %in% index_state_change_nodes,
        label = sclabel
      ),
      size = tiplab_size,
      col = "black",
      fill = "red",
      alpha = 0.5
    )

  return(p)
}

#' Get tip labels in top-to-bottom order from a ggtree plot
#'
#' @param tree_plot A ggtree plot object.
#'
#' @return A character vector of tip labels in plotting order (top to bottom).
get_tip_order <- function(tree_plot) {
  df <- tree_plot$data
  with(df, {
    i <- order(y, decreasing = TRUE)
    label[i][isTip[i]]
  })
}

#' Prepare heatmap data frame ordered by tree tip positions
#'
#' Subsets the tip rows from the full tree tibble and reorders them to match
#' the vertical layout of a ggtree plot (top to bottom).
#'
#' @param tree_tbl A "tbl_tree" tibble with full node and metadata info.
#' @param tree A phylo object used to extract tip labels.
#' @param tips_from_top A character vector of tip labels in plotting order.
#'
#' @return A tibble of tip-level data ordered for heatmap plotting, with the `label`
#'         column converted to a factor to preserve row order.
prepare_heatmap_data <- function(tree_tbl, tree, tips_from_top) {
  # Subset to rows corresponding to tip labels
  tip_rows <- tree_tbl[tree_tbl$label %in% tree$tip.label, ]

  # Match order of tip labels to ggtree's plotted order
  match_idx <- match(tips_from_top, tip_rows$label)
  if (any(is.na(match_idx))) {
    stop("Some tips from `tips_from_top` are not found in `tree_tbl$label`.")
  }

  tipdf <- tip_rows[match_idx, ]
  tipdf$label <- factor(tipdf$label, levels = rev(tipdf$label))

  return(tipdf)
}



#' Resolve or generate color mapping for a heatmap variable
#'
#' Returns a named vector of colors for a categorical heatmap variable. If a color
#' table is provided in `heatmap_colors`, it uses those. Otherwise, it generates a
#' distinct palette automatically.
#'
#' @param tree_tbl A tibble containing the heatmap variable.
#' @param heatmap_colors A named list of data frames specifying manual color mappings (optional).
#'                        Each data frame must have one column named after the variable and one
#'                        named `"color"`.
#' @param heatmap_var_i Character; the name of the variable to resolve colors for.
#' @param verbose Logical; whether to print status messages.
#'
#' @return A named character vector mapping each category to a color.
resolve_heatmap_colors <- function(tree_tbl, heatmap_colors, heatmap_var_i, verbose = FALSE) {
  values <- as.character(tree_tbl[[heatmap_var_i]])
  levels <- unique(values)

  # Use custom color table if available
  if (!is.null(heatmap_colors) && heatmap_var_i %in% names(heatmap_colors)) {
    if (verbose) message("Using custom colors for variable: ", heatmap_var_i)

    df <- heatmap_colors[[heatmap_var_i]]

    # Basic validation
    if (!all(c("color", heatmap_var_i) %in% names(df))) {
      stop("Color table for '", heatmap_var_i, "' must contain columns: '", heatmap_var_i, "' and 'color'.")
    }

    cols <- df$color
    names(cols) <- df[[heatmap_var_i]]

  } else {
    if (verbose) message("Auto-generating colors for variable: ", heatmap_var_i)

    if (length(levels) == 1) {
      cols <- "grey"
      names(cols) <- levels
    } else {
      cols <- qualpalr::qualpal(length(levels), colorspace = "pretty")$hex
      names(cols) <- levels
    }
  }

  return(cols)
}


#' Build a heatmap ggplot object for a single variable
#'
#' Constructs a heatmap from tip-level data aligned with a phylogenetic tree,
#' with customizable fill colors, legends, and optional value labels.
#'
#' @param tipdf A tibble containing reordered tip-level data (`label` must be a factor).
#' @param var Character; the column name to visualize in the heatmap.
#' @param hmdf_colors Named character vector of fill colors (names = levels).
#' @param legend_breaks Optional list of breaks (levels) to display in each variable's legend.
#' @param legend_guides Optional list of additional guide_legend() settings per variable.
#' @param show_text Logical; whether to display heatmap values as text.
#' @param tiplab_size Numeric; text size for heatmap values.
#'
#' @return A ggplot object representing the heatmap.
build_heatmap_plot <- function(tipdf, var, hmdf_colors, legend_breaks, legend_guides,
                               show_text = TRUE, tiplab_size = 1) {

  p <- ggplot(tipdf, aes(x = "", y = label)) +
    geom_tile(aes(fill = .data[[var]])) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) + 
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    xlab("") + ylab("")

  if (show_text) {
    p <- p + geom_text(aes(label = .data[[var]]), size = tiplab_size)
  }

  fill_breaks <- if (var %in% names(legend_breaks)) legend_breaks[[var]] else names(hmdf_colors)

  scale <- scale_fill_manual(
    name = var,
    values = hmdf_colors,
    breaks = fill_breaks,
    na.translate = FALSE
  )

  if (var %in% names(legend_guides)) {
    scale <- scale + guides(fill = rlang::exec(guide_legend, !!!legend_guides[[var]]))
  }

  return(p + scale)
}



#' Build and combine heatmaps aligned to a phylogenetic tree
#'
#' Constructs individual heatmaps for the given variables and combines them
#' with a ggtree plot using patchwork layout. Automatically resolves colors,
#' adjusts tip ordering, and configures legend behavior.
#'
#' @param tree_tbl A "tbl_tree" tibble containing tree metadata.
#' @param tree A phylo object corresponding to the tree.
#' @param tree_plot A ggtree plot object to which heatmaps will be added.
#' @param heatmap_var Character vector of variable names to plot as heatmaps.
#' @param heatmap_colors Optional list of data frames specifying color maps per variable.
#' @param heatmap_offset Not currently used (reserved for spacing).
#' @param heatmap_width Numeric; width(s) of heatmap bands. Either a single value or a vector.
#' @param heatmap_text Logical; whether to show heatmap values as text.
#' @param tiplab_size Numeric; size of heatmap text.
#' @param legend_breaks List of break vectors (levels) for each variable legend.
#' @param legend_guides List of lists passed to guide_legend() per variable.
#' @param verbose Logical; print color resolution info if TRUE.
#'
#' @return A patchwork plot combining the tree and all heatmaps.
build_all_heatmaps <- function(tree_tbl, tree, tree_plot,
                               heatmap_var, heatmap_colors = NULL,
                               heatmap_offset = 0, heatmap_width = 5,
                               heatmap_text = TRUE, tiplab_size = 1,
                               legend_breaks = list(), legend_guides = list(),
                               verbose = FALSE) {
  tips_from_top <- get_tip_order(tree_plot)
  tipdf <- prepare_heatmap_data(tree_tbl, tree, tips_from_top)

  # Generate individual heatmaps
  hm_list <- lapply(heatmap_var, function(var) {
    colors <- resolve_heatmap_colors(tree_tbl, heatmap_colors, var, verbose)
    build_heatmap_plot(
      tipdf = tipdf,
      var = var,
      hmdf_colors = colors,
      legend_breaks = legend_breaks,
      legend_guides = legend_guides,
      show_text = heatmap_text,
      tiplab_size = tiplab_size
    ) + ggtitle(var) + theme(plot.title = element_text(hjust = 0.5, size = 3, margin = margin(b = 2)))

  })

  # Validate/expand heatmap_width
  
	tree_width <- 8
	if (length(heatmap_width) == 1) {
	  heatmap_width <- rep(heatmap_width, length(hm_list))
	}
	if (length(heatmap_width) != length(hm_list)) {
	  stop("Length of 'heatmap_width' must match number of heatmaps.")
	}
	widths <- c(tree_width, heatmap_width)


  combined <- tree_plot + hm_list +
    patchwork::plot_layout(nrow = 1, widths = widths, guides = "collect")

  return(combined)
}


#' Plot a phylogenetic tree with optional heatmaps (refactored interface)
#'
#' Same functionality as `custom_plot_tree_long()`, but groups related parameters into lists.
#'
#' @param tree_tbl A tibble with class "tbl_tree".
#' @param tiplab_opts List with `show`, `size`, and `offset` (e.g., list(show = TRUE, size = 1, offset = 2)).
#' @param highlight_opts List with `var`, `colors`, `ancestors`.
#' @param heatmap_opts List with `var`, `colors`, `text`, `width`, `offset`.
#' @param legend_opts List with `show`, `breaks`, `guides`, `position`, `box`.
#' @param mrsd Date of the most recent sample (optional).
#' @param linewidth Numeric; line width for the tree.
#' @param verbose Logical; verbose output?
#'
#' @return A patchwork plot combining the tree and heatmaps.
#' @export
custom_plot_tree_long <- function(
  tree_tbl,
  tiplab_opts = list(show = TRUE, size = 1, offset = 2),
  highlight_opts = list(var = NULL, colors = NULL, ancestors = NULL),
  heatmap_opts = list(var = NULL, colors = NULL, text = TRUE, width = NULL, offset = 10),
  legend_opts = list(show = NA, breaks = list(), guides = list(), position = "right", box = "horizontal"),
  mrsd = NULL,
  linewidth = 0.5,
  verbose = getOption("verbose")
) {
  # Validate legend box
  legend_opts$box <- match.arg(legend_opts$box, c("horizontal", "vertical"))

  validate_inputs(tree_tbl, mrsd, heatmap_opts$width, heatmap_opts$var)

  if (is.null(heatmap_opts$width)) {
    heatmap_opts$width <- estimate_heatmap_width(tree_tbl, heatmap_opts$var, tiplab_opts$size)
  }
  

  tree_db <- treeio::as.treedata(tree_tbl)
  tree <- ape::as.phylo(tree_db)

  if (is.null(highlight_opts$ancestors)) {
    highlight_opts$ancestors <- if (!is.null(highlight_opts$var)) {
    should_highlight_ancestors(tree_tbl, highlight_opts$var)
  } else {
    FALSE
  }
  }

  p <- build_base_tree_plot(
    tree_db = tree_db,
    tree_tbl = tree_tbl,
    tree = tree,
    show_tiplab = tiplab_opts$show,
    tiplab_size = tiplab_opts$size,
    tiplab_offset = tiplab_opts$offset,
    highlight_var = highlight_opts$var,
    highlight_ancestors = highlight_opts$ancestors,
    highlight_colors = highlight_opts$colors,
    mrsd = mrsd,
    legend_show = legend_opts$show,
    linewidth = linewidth
  )

  if (highlight_opts$ancestors) {
    p <- annotate_nodes_with_states(p, tree_tbl, highlight_opts$var, tiplab_opts$size)
  }

  if (!is.null(heatmap_opts$var)) {
    p <- build_all_heatmaps(
      tree_tbl = tree_tbl,
      tree = tree,
      tree_plot = p,
      heatmap_var = heatmap_opts$var,
      heatmap_colors = heatmap_opts$colors,
      heatmap_offset = heatmap_opts$offset,
      heatmap_width = heatmap_opts$width,
      heatmap_text = TRUE,
      tiplab_size = tiplab_opts$size,
      #legend_breaks = legend_opts$breaks,
      #legend_guides = legend_opts$guides,
      verbose = verbose
    )
  }

  #p <- p &
  #  theme(
  #    legend.position = legend_opts$position,
  #    legend.box = legend_opts$box
  #  )
  return(p)
}
