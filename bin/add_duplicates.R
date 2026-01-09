library(dplyr)
library(ggplot2)
library(optparse)
library(TreeTools)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("-l", "--launch_dir"),
    type = "character",
    help = "Directory from which the pipeline was launched."
  ),
  make_option(
    c("-t", "--tree"),
    type = "character",
    help = "A dated tree in rds format."
  ),
  make_option(
    c("-d", "--duplicates"),
    type = "character",
    help = "A text file which contains tip labels for identical tips."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/prophyl",
    launch_dir = getwd(),
    tree = "final_dated_tree.rds",
    duplicates = "duplicates.txt"
  )
}

library(devtools)
load_all(args$project_dir)

# set cache directory for R.cache
R.cache::setCacheRootPath(paste0(
  args$launch_dir,
  "/.cache/R/R.cache"
))

# import tree
tree <- readRDS(args$tree)

# import duplicates
if (grepl("\\.txt$", args$duplicates)) {
  duplicates <- parse_duplicates(args$duplicates)
} else if (grepl("\\.rds", args$duplicates)) {
  duplicates <- readRDS(args$duplicates)
}

if (length(duplicates) > 0) {
  
  # add duplicates to tree
  for (i in seq_along(duplicates)) {
    for (j in 1:length(duplicates[[i]])) {
      if (duplicates[[i]][j] != names(duplicates)[i]) {
        tree <- TreeTools::AddTip(
          tree,
          where = names(duplicates)[i],
          label = duplicates[[i]][j],
          edgeLength = 0
        )
      }
    }
  }
  
  # remove tips which were used as reference for duplicates but were not in tree
  # this is relevant when tree is a subset tree.
  for (i in seq_along(duplicates)) {
    # if the name of the list entry (reference) is not identical with its first
    # element then the reference was not part of the subset and should be
    # removed. 
    if (names(duplicates)[i] != duplicates[[i]][1]) {
      tree <- ape::drop.tip(tree, names(duplicates)[i])
    }
  }
  
  # consistency check
  duplicate_tips <- duplicates %>% unlist() %>% unique()
  testthat::expect_true(all(duplicate_tips %in% tree$tip.label))
  
}

# Rename internal nodes. This is needed because TreeTools::AddTips() does not
# alter node.label when adding tips. Raised an issue here:
# https://github.com/ms609/TreeTools/issues/149
tree$node.label <- paste0("Node_", 1:tree$Nnode)

# Note that tips are added only to the dated tree. Also the internal
# node labels are regenerated only for the dated tree. The
# rests of the list elements within the dated tree object remain unchanged.

tree_tbl <- treeio::as_tibble(tree)

# CALCULATE TREE METRICS

dT <- ape::node.depth.edgelength(tree)
sts <- (tree$timeOfMRCA + dT[1:ape::Ntip(tree)])

dG <- ape::node.depth.edgelength(tree$intree)
snps <- dG[1:ape::Ntip(tree)]

rtt_df <- data.frame(snp = snps, date = sts)

g <- ggplot(rtt_df, aes(date, snp)) + 
  geom_point() +
  geom_smooth(
    method = "lm", 
    formula = y ~ I(x - tree$timeOfMRCA) - 1
  ) +
  xlab("Date") +
  ylab("Number of substitutions")

fit <- lm(snp ~ date, data = rtt_df)

# calculate metrics
results <- data.frame(
  ntips = ape::Ntip(tree),
  root_method = tree$root_method,
  r.squared = round(summary(fit)$r.squared, 3),
  adj.r.squared = round(summary(fit)$adj.r.squared, 3),
  rse = round(summary(fit)$sigma, 0),
  ssr = round(sum((summary(fit)$residuals)^2), 0),
  first = lubridate::date_decimal(min(sts, na.rm = TRUE)[1]),
  mrca = lubridate::date_decimal(tree$timeOfMRCA),
  rate = round(unname(fit$coefficients[2]), 1)
)
results$first <- as.Date(results$first, origin = "1970-01-01")
results$mrca <- as.Date(results$mrca, origin = "1970-01-01")

# export tree
saveRDS(tree, file = "final_tree.rds")
ape::write.tree(tree, file = "final_tree.nwk")
write_tsv(tree_tbl, "final_tree.tsv")

# export rtt plot
ggsave(g, file = "rtt_plot.png")

# export tree metrics
write_tsv(results, file = "metrics.tsv")
