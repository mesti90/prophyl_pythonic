# The script identifies and removes small outgroup lineages that branch directly 
# from the root of a phylogenetic tree. These outgroups, if left in the tree, 
# can disproportionately affect the estimation of the Most Recent Common 
# Ancestor (MRCA). The script operates by:
#   1. Identifying the immediate descendants (splits) at the MRCA/root node.
#   2. Calculating the fraction of tips in each descendant lineage.
#   3. Flagging and removing lineages where the fraction of tips falls below a 
#      user-defined threshold (1 percent of the tree, by default).
#   4. Iteratively repeating the process until no further lineages meet the 
#      pruning criteria. Note, by default we also do not allow the overall 
#      fraction of tips removed to exceed 1 percent of the original tree. Note,
#      if cgLIN codes are available, we do not prune if the outgroup lineage
#      shares cgLIN codes with the ingroup. NA cgLIN codes are ignored.
#   5. Outputting a pruned tree with the outgroup(s) removed for downstream 
#      analyses.
# This process helps ensure that rare, deeply branching lineages do not skew 
# phylogenetic inferences.
# Note, the tree returned by Gubbins is a rooted tree, probably midpoint rooted.
# Earlier version of this script used the rooted tree as is but I was concerned
# that multiple rooting strategies within a single pipeline may be odd. This new
# version of the script uses rtt-rms rooting instead of midpoint rooting.

library(devtools)
library(optparse)
library(dplyr)
rm(list = ls())

args_list <- list(
  make_option(
    "--project_dir",
    type = "character",
    help = "The project directory.",
    default = "prophyl"
  ),
  make_option(
    "--tree",
    type = "character",
    help = "A phylogenetic tree",
    default = "results/tree_input_ST147/build_tree/chromosomes.nodup.node_labelled.final_tree.tre"
  ),
  make_option(
    "--gentypes",
    type = "character",
    help = "gentype summary table",
    default = "data/tree_input/tree_input_ST147.tsv"
  ),
  make_option(
    "--step_threshold",
    type = "numeric",
    help = "Maximum ratio of tips to prune in a single step.",
    default = 0.01
  ),
  make_option(
    "--overall_threshold",
    type = "numeric",
    help = "Maximum ratio of tips to prune overall.",
    default = 0.01
  ),
  make_option(
    "--threads",
    type = "numeric",
    help = "Number of CPU threads to use for rooting",
    default = 10
  ),
  make_option(
    "--outtree",
    type = "character",
    help = "Pruned tree",
    default = "mrca_pruned_tree.tre"
  ),
  make_option(
    "--dropped_tips",
    type = "character",
    help = "Tsv file for dropped tips",
    default = "pruning.dropped_tips.tsv"
  )

  
)

args_parser  <- OptionParser(option_list = args_list)
args  <- parse_args(args_parser)

load_all(args$project_dir)
tree <- ape::read.tree(args$tree) 

# if tree is rooted, unroot
if (ape::is.rooted(tree)) {
  tree <- ape::unroot(tree)
}

# reroot tree using rtt-rms rooting
gentypes <- read.csv(args$gentypes, sep = "\t")


date_mid <- function(x) {
	x <- as.character(x)
	as.Date(dplyr::case_when(
		grepl("^\\d{4}-\\d{2}-\\d{2}$", x) ~ x,
		grepl("^\\d{4}-\\d{2}$", x)        ~ paste0(x, "-15"),
		grepl("^\\d{4}$", x)               ~ paste0(x, "-07-01"),
		TRUE                               ~ NA_character_
	))
}


tip_dates_raw <- as.numeric(date_mid(gentypes$Date))
names(tip_dates_raw) <- gentypes$strain
tip_dates <- tip_dates_raw[tree$tip.label]

## collect tip dates in the same order as tree$tip.label
#tip_dates <- unname(sapply(tree$tip.label, function(x) {
#	index <- which(gentypes$assembly == x)
#	out <- try(date_middle(gentypes$collection_date[index]), silent = TRUE)
#	if (inherits(out, "try-error")) stop(x)
#	out
#}))


# root tree
rtree <- root_rtt(
  t = tree,
  tip.dates = tip_dates,
  topx = 1, 
  ncpu = args$threads,
  objective = "rms",
  objective_fn = NULL
)[[1]]

count_tips <- function(tree, node) {
  children <- tree$edge[tree$edge[,1] == node, 2]
  if (length(children) == 0) {
    return(1)  # this node is a tip
  } else {
    return(sum(sapply(children, function(child) count_tips(tree, child))))
  }
}

# Function to get tip labels under a node
get_tips <- function(tree, node) {
  children <- tree$edge[tree$edge[,1] == node, 2]
  if (length(children) == 0) {
    return(tree$tip.label[node])
  } else {
    return(unlist(lapply(children, function(child) get_tips(tree, child))))
  }
}

get_lin <- function(lin, level) {
  foo <- function(x, level) {
    out <- strsplit(x, split = ",")[[1]]
    out[which(out == "*")] <- NA
    out <- as.integer(out)
    if (all(is.na(out[1:level]))) {
      out <- NA
    } else {
      out <- paste(out[1:level], collapse = "_")
    }
    return(out)
  }
  unname(sapply(lin, function(x) foo(x, level)))
}

if (!"LIN_code" %in% names(gentypes)) {
  warning("cgLIN codes are not available. Pruning without cgLIN codes.")
}

tree <- rtree
dropped_tips <- vector()
drop = TRUE
while (drop == TRUE) {
  # Root node index in a rooted tree
  root_node <- ape::Ntip(tree) + 1
  # Get the two daughter nodes of the root
  daughters <- tree$edge[tree$edge[,1] == root_node, 2]
  # Count tips on each side of the root
  tip_counts <- sapply(daughters, function(node) count_tips(tree, node))
  # Find daughter with the smaller number of tips
  min_side <- daughters[which.min(tip_counts)]
  # Extract tip labels on the branch with fewer tips
  tip_labels_min_branch <- get_tips(tree, min_side)
  # Extract tip labels on the branch with more tips
  tip_labels_max_branch <- tree$tip.label[which(!tree$tip.label %in% tip_labels_min_branch)]
  # If LIN codes are available and non-NA level 3 LIN codes are isolated, drop
  # Also overwrite initial tree to reset thresholds
  if ("LIN_code" %in% names(gentypes)) {
    min_lins <- gentypes$LIN_code[which(gentypes$assembly %in% tip_labels_min_branch)]
    min_lins <- get_lin(min_lins, level = 3) |> unique() |> na.omit()
    max_lins <- gentypes$LIN_code[which(gentypes$assembly %in% tip_labels_max_branch)]
    max_lins <- get_lin(max_lins, level = 3) |> unique() |> na.omit()
    if (any(min_lins %in% max_lins) == FALSE) {
      tree <- ape::drop.tip(tree, tip_labels_min_branch)
      rtree <- tree
      next()
    }
  }
  # If split is not distorted enough, break
  drop_length <- length(tip_labels_min_branch) 
  drop_threshold <- args$step_threshold * length(tree$tip.label)
  if (drop_length > drop_threshold) break()
  # If LIN codes are available and non-NA LIN codes are not isolated, break
  if ("LIN_code" %in% names(gentypes)) {
    min_lins <- gentypes$LIN_code[which(gentypes$assembly %in% tip_labels_min_branch)]
    min_lins <- get_lin(min_lins, level = 4) |> unique() |> na.omit()
    max_lins <- gentypes$LIN_code[which(gentypes$assembly %in% tip_labels_max_branch)]
    max_lins <- get_lin(max_lins, level = 4) |> unique() |> na.omit()
    if (any(min_lins %in% max_lins)) break()
  }
  # If tree is too small after drop, break, keep the previous iteration
  tree_staged <- ape::drop.tip(tree, tip_labels_min_branch)
  if (length(tree_staged$tip.label) < (1-args$overall_threshold) * length(rtree$tip.label)) break()
  # Update tree and vector of dropped tips
  tree <- tree_staged
  dropped_tips <-  c(dropped_tips, tip_labels_min_branch)
}

ape::write.tree(tree, file = args$outtree)

write_tsv(as.data.frame(dropped_tips), file = args$dropped_tips)
