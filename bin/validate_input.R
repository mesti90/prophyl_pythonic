library(optparse)
rm(list=ls())

args_list <- list(
  make_option(
    "--project_dir",
    type = "character",
    help = "The project directory.",
    default = "~/Methods/prophyl"
  ),
  make_option(
    "--assemblies",
    type = "character",
    help = "A list of assemblies in tsv format.",
    default = "assemblies.tsv"
  )
)

args_parser  <- OptionParser(option_list = args_list)
args  <- parse_args(args_parser)

library(devtools)
load_all(args$project_dir)

df <- read_df(args$assemblies)

# Data frame must contain a number of columns.

varnames <- c(
  "assembly", 
  "R1_path",
  "R2_path",
  "assembly_path",
  "collection_date",
  "genome_size",
  "longest_contig"
)

index <- which(!varnames %in% colnames(df))

if (length(index) > 0) {
  varnames_missing_collapsed <- paste(varnames[index], collapse = ", ")
  msg <- paste0(
    "The following required variables are missing from the input table: ",
    varnames_missing_collapsed,
    "."
  )
  stop(msg)
}

# Validate values in the assembly column.

if (length(unique(df$assembly)) != length(df$assembly)) {
  stop("Assembly names must be unique. Check.")
}

if (sum(grepl("root", df$assembly, ignore.case = TRUE)) > 0) {
  stop("Assembly name cannot be 'root'. Specify another name.")
}

if (sum(grepl("^node_", df$assembly, ignore.case = TRUE)) > 0) {
  stop("Assembly name cannot start with 'Node_'. Specify another name.")
}

if (sum(!is.na(suppressWarnings(as.numeric(df$assembly)))) > 0) {
  stop("Assembly name cannot be a number. Specify another name.")
}

# All genomes must have collection dates.

if (any(is.na(df$collection_date))) {
  stop("All assemblies must contain collection dates. Check.")
}

# non-na paths must exist

R1_paths <- df$R1_path[!is.na(df$R1_path)]
R2_paths <- df$R2_path[!is.na(df$R2_path)]
assembly_paths <- df$assembly_path[!is.na(df$assembly_path)]

all_paths <- c(R1_paths, R2_paths, assembly_paths)
path_exists <- file.exists(all_paths)

if (any(!path_exists)) {
  stop("Some of the specified file paths do not exist. Check.")
}

write.table(
  df,
  file = "assemblies.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
