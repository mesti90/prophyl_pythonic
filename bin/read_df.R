#' Read table
#' 
#' This is a convenience function for importing table formats that are commonly
#' used in this package, i.e. csv, tsv, rds using a single function.
#' @param df character; a file path
#' @return a data frame
read_df <- function(df) {
  na_strings <- c(
    naniar::common_na_strings, "Missing"
  )
  if (grepl("csv$", df)) {
    df <- read.csv(df, na.strings = na_strings)
  } else if (grepl(".tsv$", df)) {
    df <- read.csv(df, sep = "\t", na.strings = na_strings)
  } else if (grepl("rds$", df)) {
    df <- readRDS(df)
  } else {
    stop("File format not supported.")
  }
  return(df)
}
