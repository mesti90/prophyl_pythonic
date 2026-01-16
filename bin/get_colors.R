get_colors <- function(
    df,
    var,
    colors = NULL,
    colorspace = "pretty",
    na.color = "gray",
    other.color = "#808080",
    verbose = getOption("verbose")
    ) {
  testthat::expect_true("data.frame" %in% class(df))
  testthat::expect_true(class(var) == "character")
  testthat::expect_true(class(colors) == "NULL" | "data.frame" %in% class(colors))
  testthat::expect_true(class(na.color) == "character")
  testthat::expect_true(class(other.color) == "character")
  testthat::expect_true(class(verbose) == "logical")
  var <- match.arg(var, choices = names(df))
  
  vars <- sort(unique(df[[var]]), na.last = TRUE)
  if (is.null(colors)) {
    colors <- data.frame(var = vars)
    names(colors) <- var
    colors$color <- qualpalr::qualpal(nrow(colors), colorspace = colorspace)$hex
  } else {
    if (!is.data.frame(colors)) {
      stop(paste("Colors must be provided as a data frame."))
    }
    if (any(table(colors[[var]])) > 1) {
      stop("Colors data frame must contain a single row for each ", var, ".")
    }
    if (any(!vars %in% colors[[var]])) {
      if (verbose) message("Adding new colors.")
      new_vars <- vars[which(!vars %in% colors[[var]])]
      new_colors <- data.frame(var = new_vars)
      names(new_colors) <- var
      if (nrow(new_colors) == 1) {
        new_colors$color <- qualpalr::qualpal(
          nrow(new_colors)+1, colorspace = colorspace)$hex[1]
      } else {
        new_colors$color <- qualpalr::qualpal(
          nrow(new_colors), colorspace = colorspace)$hex
      } 
      colors <- dplyr::bind_rows(
        colors,
        new_colors
      )
    }
    if (any(colors[[var]] %in% vars == FALSE)) {
      if (verbose) message("Removing unused colors.")
      colors <- colors[which(colors[[var]] %in% vars),]
    }
  }
  colors <- colors[order(colors[[var]]), ]
  if ("Other" %in% vars) {
    colors$color[which(colors[[var]] == "Other")] <- other.color
    index_other <- which(colors[[var]] == "Other")
    if (nrow(colors) > 1) {
      colors <- dplyr::bind_rows(colors[-index_other,], colors[index_other,])
    }
  }
  is_na_like <- is.na(colors[[var]]) | colors[[var]] == "" | colors[[var]] == "NA"
  colors$color[which(is_na_like)] <- na.color
  
  if (any(is_na_like)) {
	index_na <- which(is_na_like)
	if (nrow(colors) > 1) {
		colors <- dplyr::bind_rows(colors[-index_na, ], colors[index_na, ])
	}
}
  return(colors)
}
