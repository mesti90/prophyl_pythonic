# This script takes a phylogenetic tree and attempts to root it using three
# aproaches:
#
# 1. Midpoint rooting
# 2. Minimum Ancestor Deviation (MAD)
# 3. Root-to-tip regression
#
# Root-to-tip regression uses a custom function that was built on the
# non-exported .multi.rtt() function from the treedater package which in turn
# was built on the rtt() function from the ape package. For more details, check
# the function documentation ?root_rtt().

library(devtools)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    "--project_dir",
    type = "character",
    help = "Path to project directory.",
    default = "~/Methods/prophyl"
  ),
  make_option(
    "--tree",
    type = "character",
    help = "Path to tree file.",
    default = "treeshrink.tre"
  ),
  make_option(
    "--assemblies",
    type = "character",
    help = "Path to assemblies file.",
    default = "assemblies.tsv"
  ),
  make_option(
    "--root_method",
    type = "character",
    help = "Rooting method to use. Options are 'midpoint', 'rtt_correlation', 
    'rtt_rms', 'rtt_rsquared', 'all'.",
    default = "rtt_rms"
  ),
  make_option(
    "--root_topn",
    type = "integer",
    help = "Number of top trees to return for the selected root method. Only
    used for root-to-tip regression methods.",
    default = 1
  ),
  make_option(
    "--threads",
    type = "integer",
    help = "Number of threads to use.",
    default = 10
  ),
  make_option(
    "--outprefix",
    type="character",
    help="Prefix of output files",
    default="output"
  )
)

args_parser  <- OptionParser(option_list = args_list)

args  <- parse_args(args_parser)

load_all(args$project_dir)

args$rooted_trees_rds = paste0(args$outprefix, ".trees.rds")
args$rtt_metrics_rds=paste0(args$outprefix,".rtt_metrics.rds")
args$rtt_plots_pdf=paste0(args$outprefix,".rtt_plots.pdf")

date_mid <- function(x) {
	x <- as.character(x)
	as.Date(dplyr::case_when(
		grepl("^\\d{4}-\\d{2}-\\d{2}$", x) ~ x,
		grepl("^\\d{4}-\\d{2}$", x)        ~ paste0(x, "-15"),
		grepl("^\\d{4}$", x)               ~ paste0(x, "-07-01"),
		TRUE                               ~ NA_character_
	))
}


# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

# read tree
tree <- ape::read.tree(args$tree)

# if tree is rooted, unroot
if (ape::is.rooted(tree)) {
  tree <- ape::unroot(tree)
}

rooted_trees <- list()

# OPTION 1: MIDPOINT ROOTING

if (args$root_method %in% c("midpoint", "all")) {
  # root tree
  rooted_trees[["midpoint"]] <- phytools::midpoint.root(tree)
}

# Note: OPTION 2  is removed until Issue #73 is fixed.
# # OPTION 2: MINIMUM ANCESTOR DEVIATION (MAD)

# # root tree
# mad <- root_mad(
#   tree,
#   output_mode = "full",
#   cache = TRUE,
#   threads = args$threads,
#   verbose = TRUE
# )

# # add tips that were collapsed

# mad_tree <- mad[[3]]
# collapsed_tips <- mad[[7]]

# for (i in 1:nrow(collapsed_tips)) {
#   mad_tree <- TreeTools::AddTip(
#     mad_tree,
#     where = collapsed_tips$keep[i],
#     label = collapsed_tips$drop[i],
#     edgeLength = 0
#   )
# }

# rooted_trees[["mad"]] <- mad_tree

# OPTION 3: ROOT-TO-TIP REGRESSION

if (args$root_method %in% c(
  "rtt_correlation", "rtt_rms", "rtt_rsquared", "all")) {
  
  # read assemblies
  assemblies <- read.csv(args$assemblies, sep = "\t")
  
  # collect tip dates in the same order as tree$tip.label
  tip_dates_raw <- as.numeric(date_mid(assemblies$Date))
  names(tip_dates_raw) <- assemblies$strain
  tip_dates <- tip_dates_raw[tree$tip.label]
  
  # TODO: look for better objectives
  objective_rlm_slope <- function(x,y) MASS::rlm(y ~ x)$coef[2]
  objective_rlm_rms <- function(x,y) -summary(MASS::rlm(y ~ x))$sigma^2
  
  # Remove custom objectives for now
  # objective <- list(
  #   "correlation" = NULL,
  #   "rsquared" = NULL,
  #   "rms" = NULL,
  #   "rlm_slope" = objective_rlm_slope,
  #   "rlm_rms" = objective_rlm_rms
  # )
  
  if (args$root_method == "rtt_correlation") {
    objective <- list("correlation" = NULL)
  } else if (args$root_method == "rtt_rsquared") {
    objective <- list("rsquared" = NULL)
  } else if (args$root_method == "rtt_rms") {
    objective <- list("rms" = NULL)
  } else {
    objective <- list(
      "correlation" = NULL,
      "rsquared" = NULL,
      "rms" = NULL
    )
  }
  
  # return the top_n trees for each objective
  top_n <- args$root_topn
  cat(top_n)
  cat("\n")
  for (i in seq_along(objective)) {
    rtree <- root_rtt(
      t = tree,
      tip.dates = tip_dates,
      topx = top_n, 
      ncpu = args$threads,
      objective = names(objective)[[i]],
      objective_fn = objective[[i]]
    )
    names(rtree) <- paste0("rtt_", names(objective)[i], "_", 1:top_n)
    index_from = length(rooted_trees) + 1
    index_to = length(rooted_trees) + top_n
    rooted_trees[index_from:index_to] <- rtree
    names(rooted_trees)[index_from:index_to] <- names(rtree)
  }
}

# CALCULATE ROOT TO TIP METRICS FOR EACH ROOTED TREE

# calculate snps for each rooted tree
snp <- lapply(rooted_trees, function(x) {
  ape::node.depth.edgelength(x)[1:ape::Ntip(tree)]
})

# rescale tip_dates to calendar dates
tip_dates <- as.Date(tip_dates, origin = "1970-01-01")

# recalculate root-to-tip regression using calendar dates
fit <- lapply(snp, function(x) lm(x~tip_dates))

# calculate metrics for each fit
results <- data.frame(
  r.squared = sapply(fit, function(x) summary(x)$r.squared),
  adj.r.squared = sapply(fit, function(x) summary(x)$adj.r.squared),
  rse = sapply(fit, function(x) summary(x)$sigma),
  ssr = sapply(fit, function(x) sum((summary(x)$residuals)^2)),
  mrca = sapply(fit, function(x) -x$coef[1]/x$coef[2]),
  first = min(tip_dates, na.rm = TRUE)[1]
)
results$first <- as.Date(results$first, origin = "1970-01-01")
results$mrca <- as.Date(results$mrca, origin = "1970-01-01")

df <- data.frame()
for (i in seq_along(fit)) {
  new_df <- data.frame(
    name = names(rooted_trees)[i],
    snp = fit[[i]]$model$x,
    date = fit[[i]]$model$tip_dates
  )
  df <- dplyr::bind_rows(df, new_df)
}

g <- ggplot(df, aes(date, snp)) + 
  geom_point() + 
  facet_grid(name~.) + 
  geom_smooth(method = "lm")

# export rooted tree object
saveRDS(rooted_trees, file = args$rooted_trees_rds)
# export root to tip metrics
saveRDS(results, file = args$rtt_metrics_rds)
# export root to tip regression plots
ggsave(
  filename = args$rtt_plots_pdf,
  plot = g,
  width = 10,
  height = 5 * length(rooted_trees),
  limitsize = FALSE
)

# if (!dir.exists("rooted_trees")) dir.create("rooted_trees")

# for (i in seq_along(rooted_trees)) {
#   ape::write.tree(
#     rooted_trees[i],
#     file = paste0("rooted_trees/rooted_tree_", names(rooted_trees)[i], ".tre")
#   )
# }

# consistency checks

# check that all trees have the same number of tips
ntips <- sapply(rooted_trees, function(x) length(x$tip.label))
testthat::expect_equal(length(unique(ntips)), 1)

# end logging
if (!interactive()) {
  sink(con)
}
