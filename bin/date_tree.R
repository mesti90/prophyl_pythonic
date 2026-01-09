# This script takes a list of rooted phylogentic trees as a single rds file and
# other input files and calculates dated phylogenetic trees for each tree in
# the list.
#
# If the collection date of a sample is known by day, it will be used as is. If
# the collection date of the sample is known by month or year, a matching range
# will be used as input and the middle value as the starting value. If the
# collection date of the sample is in another format or not known, the sample
# will be dropped from the analysis. This is different from the approach used
# by treedater authors in the article https://doi.org/10.1093/ve/vex025 where
# they seem to have interpolated unknown samples to fall between the boundaries
# of the rest of the genomes.

library(devtools)
library(dplyr)
library(magrittr)
library(optparse)
library(treedater)
rm(list = ls())

args_list <- list(
  make_option(
    "--project_dir",
    type = "character",
    help = "Path to project directory.",
    default = "~/Methods/prophyl"
  ),
  make_option(
   "--trees",
   type = "character",
   help = "A list of rooted trees in rds format.",
   default = "rooted_trees.rds"
  ),
  make_option(
    "--snps",
    type = "character",
    help = "Path to a fasta file containing snps.",
    default = "chromosomes.nodup.filtered_polymorphic_sites.fasta"
  ),
  make_option(
    "--assemblies",
    type = "character",
    help = "Path to assemblies file.",
    default = "assemblies.tsv"
  ),
  make_option(
    "--threads",
    type = "integer",
    help = "Number of threads to use.",
    default = 10
  ),
  make_option(
    "--branch_dimension",
    type = "character",
    help = "Dimension of branch lengths in tree. Either 'snp_per_genome' or 'snp_per_site'.",
    default = "snp_per_genome"
  ),
  make_option(
    "--clock",
    type = "character",
    help = "The choice of molecular clock model. Choices are 'uncorrelated', 'additive', or 'strict'.",
    default = "strict"
  ),
  make_option(
    "--reroot",
    type = "logical",
    help = "Whether to reroot the tree using treedater's standard rerooting functionality.",
    default = FALSE
  )
)

args_parser  <- OptionParser(option_list = args_list)
args  <- parse_args(args_parser)

load_all(args$project_dir)

if (!interactive()) {
  # create log file and start logging
  con <- file("log.txt")
  sink(con, split = TRUE)
}

# Input validation
branch_dimension <- match.arg(
  args$branch_dimension,
  choices = c("snp_per_genome", "snp_per_site")
)

clock <- match.arg(
  args$clock,
  choices = c("uncorrelated", "additive", "strict")
)

# import trees
trees <- readRDS(args$trees)

# consistency checks

# check that all trees have the same number of tips
ntips <- sapply(trees, function(x) length(x$tip.label))
testthat::expect_equal(length(unique(ntips)), 1)

# import snps
f <- seqinr::read.fasta(args$snps)

# rescale branch lengths if necessary
alignment_length <- length(f[[1]])
if (branch_dimension == "snp_per_genome") {
  # rescale branch lengths from per genome to per site (required for treedater)
  trees <- lapply(trees, function(tree) {
    tree$edge.length <- tree$edge.length / alignment_length
    return(tree)
  })
}

# import assemblies
assemblies <- read.csv(args$assemblies, sep = "\t", header = TRUE)

# drop tips which cannot be found in assembly table and give a warning
index <- which(trees[[1]]$tip.label %in% assemblies$assembly == FALSE)
if (length(index) > 0) {
  tips_to_drop <- trees[[1]]$tip.label[index]
  tips_to_drop_collapsed <- paste(tips_to_drop, collapse = ", ")
  trees <- lapply(trees, function(tree) {
    tree <- ape::drop.tip(tree, tips_to_drop)
    return(tree)
  })
  msg <- paste0(
    "One or more tips could not be found in assembly table and were dropped: ",
    tips_to_drop_collapsed,
    "."
  )
  warning(msg)
}

# filter to assemblies that are included in the tree
index <- which(assemblies$assembly %in% trees[[1]]$tip.label == FALSE)
if (length(index) > 0) {
  assemblies <- assemblies[-index, ]
}

# filter to relevant columns
assemblies <- assemblies[, which(names(assemblies) %in% c(
  "assembly", "collection_date"
))]

# add new temporal variables
assemblies$date <- date_middle(assemblies$collection_date)
assemblies$lower <- date_lower(assemblies$collection_date)
assemblies$upper <- date_upper(assemblies$collection_date)

# define range for dates which are not known exactly
# note, the which function automatically removes NA values
# so this step will only include dates which are known to some extent
# but not exactly
index_uncertain <- which(assemblies$upper > assemblies$lower)
if (length(index_uncertain) > 0) {
  uncertain_dates <- assemblies[index_uncertain, c("lower", "upper")]
  rownames(uncertain_dates) <- assemblies$assembly[index_uncertain]
} else {
  uncertain_dates <- NULL
}

# drop tips with unknown sampling dates
index <- which(is.na(assemblies$date))
tips_to_drop <- assemblies$assembly[index]
if (length(index) > 0) {
  # drop from tree
  trees <- lapply(trees, function(tree) {
    tree <- ape::drop.tip(tree, tips_to_drop)
    return(tree)
  })
  # drop from assemblies
  assemblies <- assemblies[-which(assemblies$assembly %in% tips_to_drop), ]
  # print warning message
  tips_to_drop_collapsed <- paste(tips_to_drop, collapse = ", ")
  msg <- paste0(
    "One or more tips were dropped because no data on sampling data was found: ",
    tips_to_drop_collapsed,
    "."
  )
  warning(msg)
}

# rename tips to include dates
for (i in 1:length(trees)) {
  tree <- trees[[i]]
  tree$tip.label <- sapply(tree$tip.label, function(x) {
  index <- which(assemblies$assembly == x)
  paste(x, assemblies$date[index], sep = "|")
  }, USE.NAMES = FALSE)
  trees[[i]] <- tree
}

if (!is.null(uncertain_dates)) {
  # rename rownames in uncertain dates to include dates
  row.names(uncertain_dates) <- sapply(row.names(uncertain_dates), function(x) {
    index <- which(assemblies$assembly == x)
    # this will include the middle value of the range
    paste(x, assemblies$date[index], sep = "|")
  }, USE.NAMES = FALSE)
}

# extract dates from tip labels in appropriate format
sts <- sampleYearsFromLabels(trees[[1]]$tip.label, delimiter = "|")

# if the reroot argument is set to TRUE, and the tree is rooted, unroot the tree
if (as.logical(args$reroot)) {
  trees <- lapply(trees, function(tree) {
    if (ape::is.rooted(tree)) {
      tree <- ape::unroot(tree)
    }
    return(tree)
  })
}

# date trees
# if the tree is rooted, the function will not reroot the tree
# if the tree is unrooted, the function will root the tree
dtr <- list()
for (i in 1:length(trees)) {
  dtr[[i]] <- treedater::dater(
    trees[[i]],
    sts,
    s = alignment_length,
    estimateSampleTimes = uncertain_dates,
    clock = clock,
    ncpu =  args$threads)
}

# link dated trees with rooting methods
names(dtr) <- names(trees)

# rescale non-dated branch lengths from per site to per genome (more intuitive)
dtr <- lapply(dtr, function(tree) {
  tree$intree$edge.length <- tree$intree$edge.length * alignment_length
  return(tree)
})

# rename internal nodes
dtr <- lapply(dtr, function(tree) {
  tree %<>% makeNodeLabel(., method = "number", prefix = "Node_")
  return(tree)
})

# add root method to dated tree objects
dtr_class <- class(dtr[[1]])

for (i in 1:length(dtr)) {
  dtr[[i]]$root_method <- names(dtr)[i]
  class(dtr[[i]]) <- dtr_class
}

# export RTT plots
if (!dir.exists("rtt_plots")) dir.create("rtt_plots")

for (i in 1:length(dtr)) {
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  pdf(file = paste0("rtt_plots/",names(dtr)[i],".pdf"))
  rootToTipRegressionPlot(dtr[[i]])
  dev.off()
}

for (i in 1:length(dtr)) {
  dtr[[i]]$tip.label <- unname(sapply(dtr[[i]]$tip.label, function(x) {
    strsplit(x, "\\|")[[1]][1]
  }))
}

# export dated trees
if (!dir.exists("dated_trees")) dir.create("dated_trees")

for (i in 1:length(dtr)) {
  write.tree(dtr[[i]], file = paste0("dated_trees/",names(dtr)[i],".tre"))
}

# export dated_trees object as rds
saveRDS(dtr, file = "dated_trees.rds")

# end logging
if (!interactive()) {
  # end logging
  sink(con)
}
