rm(list=ls())

if (!interactive()) {
  library(optparse)
  args_list <- list(
    make_option(
      "--project_dir",
      type = "character",
      help = "Path to the project directory",
      default = "prophyl-priv"
    ),
    make_option(
      "--store_dir",
      type = "character",
      help = "Path to the storage directory for the process",
      default = "results"
    ),
    make_option(
      "--assemblies",
      type = "character",
      help = "Path to a tbl of assemblies",
      default = "assemblies.tsv"
    )
  )
  args_parser  <- OptionParser(option_list = args_list)
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "prophyl-priv",
    store_dir = "results",
    assemblies = "assemblies.tsv"
  )
}

library(devtools)
library(dplyr)
library(R.utils)
load_all(args$project_dir)

df <- read_df(args$assemblies)

df$mode = NA
df$error = NA
for (i in 1:nrow(df)) {
  if (!is.na(df$R1_path[i])) {
    if (!is.na(df$R2_path[i])) {
      # paired_reads
      R1_exists <- file.exists(df$R1_path[i])
      R2_exists <- file.exists(df$R2_path[i])
      if (all(R1_exists, R2_exists)) {
        df$mode[i] <- "paired"
      } else {
        df$error[i] <- "R1 or R2 file not found."
      }
    } else {
      # append single_reads
      R1_exists <- file.exists(df$R1_path[i])
      if (R1_exists) {
        df$mode[i] <- "single"
      } else {
        df$error[i] <- "R1 file not found."
      }
    }
  } else if (!is.na(df$R2_path[i])) {
    df$error[i] <- "Syntax error. For single end reads use R1."
  } else if (!is.na(df$assembly_path[i])) {
    # append contigs
    assembly_exists <- file.exists(df$assembly_path[i])
    if (assembly_exists) {
      df$mode[i] <- "contigs"
    } else {
      df$error[i] <- "Assembly file not found."
    }
  } else {
    df$error[i] <- "No read or assembly files found."
  }
}

df <- df %>%
  dplyr::select(
    assembly,
    R1_path,
    R2_path,
    assembly_path,
    mode,
    error
  )

unzipped_assemblies_dir <- "unzipped_assemblies"

if (!dir.exists("unzipped_assemblies")) {
  dir.create(unzipped_assemblies_dir)
}

index <- which(df$mode == "contigs")

if (length(index) > 0) {
  for (i in index) {
    assembly_path_original <- df$assembly_path[i]
    if (grepl(".gz$", assembly_path_original)) {
      # copy file to a new location
      assembly_path_copy <- paste0(
        getwd(), "/",
        unzipped_assemblies_dir, "/",
        basename(assembly_path_original)
      )
      file.copy(from = assembly_path_original, to = assembly_path_copy)
      # decompress
      R.utils::gunzip(filename = assembly_path_copy, remove = TRUE)
      # update file path in dataframe
      assembly_path_store <- paste0(
        args$store_dir, "/",
        unzipped_assemblies_dir, "/",
        basename(assembly_path_copy)
      )
      df$assembly_path[i] <- gsub(".gz$", "", assembly_path_store)
    }
  }
}

write_tsv(df, "snippy_input.tsv")

if (any(!is.na(df$error))) {
  stop("Some files were not found. See error.tsv for details.")
}
