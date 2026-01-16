library(ggtree)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    "--tree",
    type = "character",
    help = "Path to newick tree before pruning",
    default = "chromosomes.nodup.node_labelled.final_tree.tre"
  ),
  make_option(
    "--pruned_tree",
    type = "character",
    help = "Path to a newick tree after pruning",
    default = "treeshrink.tre"
  )
)

args_parser  <- OptionParser(option_list = args_list)
args  <- parse_args(args_parser)

tree <- ape::read.tree(args$tree)
stree <- ape::read.tree(args$pruned_tree)

pruned_tips <- setdiff(tree$tip.label, stree$tip.label)

write.table(
  pruned_tips,
  file = "pruned_tips.txt"
)

# note, even though the "pruned" variable does not exist at this point, it will
# be defined later and because of this, the plot will work.

if (ape::is.rooted(tree) == FALSE) tree <- ape::unroot(tree)

p <- ggtree(tree, aes(color = pruned), layout = "equal_angle")

p <- p %<+% 
  # this is where a data.frame containing the "pruned" variable is added.
  data.frame(
    label = tree$tip.label,
    pruned = tree$tip.label %in% pruned_tips
  ) +
  geom_tiplab(aes(color = pruned), size = 1) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))

ggsave(
  p,
  file = "marked_tree.pdf",
  width = sqrt(ape::Ntip(tree)),
  height = sqrt(ape::Ntip(tree)),
  limitsize = FALSE
)
