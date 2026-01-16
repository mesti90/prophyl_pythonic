#' Root a tree using root to tip distances
#' 
#' This functions roots a phylogenetic tree using root to tip distances.
#' @param t phylo; a phylogentic tree.
#' @param tip.dates	character; a vector of sampling times associated to the tips
#' of t, in the same order as t$tip.label.
#' @param ncpu integer; number of cores to use.
#' @param objective character; the name of the objective to use for finding the
#' best root. See Details for more information.
#' @param objective_fn function; either \code{NULL}, or a function that takes
#' two numeric vectors x then y, and returns a numeric value. This function will
#' be used to calculate the objective value for a given root. See Details for
#' more information.
#' @param opt.tol tolerance for optimization precision.
#' @return The function returns a list of rooted phylogenetic trees.
#' @details There are two ways to define an objective. The simplest approach is
#' to set \code{objective} to \code{"correlation"}, \code{"rsquared"}, or
#' \code{"rms"} and leave \code{objective_fn} as \code{NULL}. This way the
#' function will use built in functions for finding the best root position:
#' \itemize{
#'   \item \code{"correlation"}: look for the tree with the largest
#'     \code{cor.test(y, x)$estimate}
#'   \item \code{"rsquared"}: look for the tree with the largest
#'     \code{summary(lm(y ~ x))$r.squared}
#'   \item \code{"rms"}: look for the tree with the smallest
#'     \code{summary(lm(y ~ x))$sigma^2}
#' }
#' @details If none of these builtin objective functions are appropriate, you
#' can define a custom objective function in \code{objective_fn}. This function
#' must take two numeric vectors x then y and return a numeric value.
#' @note At the time of writing this documentation \code{ape::rtt()} cannot
#' handle missing dates, but \code{treedater:::.multi.rtt()} can. This
#' flexibility is preserved in \code{root_rtt()} as well.
#' @references Based on the non-exported \code{.multi.rtt()} function in the
#' \code{treedater} package.
#' https://github.com/emvolz/treedater/blob/master/R/multiRtT.R.
#' @references Based on the \code{rtt()} function in the \code{ape} package.
#' https://cran.r-project.org/web/packages/ape/index.html.

root_rtt <-function (
    t,
    tip.dates,
    topx = 1,
    ncpu = 1,
    objective = NULL,
    objective_fn = NULL,
    opt.tol = .Machine$double.eps^0.25) {
    topx <- max(1, topx)
    if (!is.null(objective)) {
      if (is.null(objective_fn)) {
        if (objective %in% c("correlation", "rsquared", "rms")) {
          if (objective == "correlation") 
            objective <- function(x, y) cor.test(y, x)$estimate
          else if (objective == "rsquared") 
            objective <- function(x, y) summary(lm(y ~ x))$r.squared
          else if (objective == "rms") 
            objective <- function(x, y) -summary(lm(y ~ x))$sigma^2
        } else {
          msg <- paste0(
            "To use built in objectives, 'objective' must be 'correlation', ",
            "'rsquared' or 'rms'."
          )
          stop(msg)
        }
      } else {
        if (objective %in% c("correlation", "rsquared", "rms")) {
          msg <- paste0(
            "The following 'objective' names are reserved: ",
            "'correlation', 'rsquared', 'rms'. Either use a different ",
            "'objective' name for your custom function or set 'objective_fn' ",
            "to NULL."
          )
          stop(msg)
        } else {
          objective <- objective_fn
        }
      }
    } else {
      if (is.null(objective_fn)) {
        stop("Please define an objective. See ?root.rtt() for more information.")
      } else {
        objective <- objective_fn
      }
    }
    t <- ape::multi2di(t)
    ut <- ape::unroot(t)
    dist <- dist.nodes(ut)[, 1:(ut$Nnode + 2)]
    f <- function(x, parent, child) {
      edge.dist <- x * dist[parent, ] + (1 - x) * dist[child, 
      ]
      objective(tip.dates, edge.dist)
    }
    obj.edge <- if (ncpu > 1) 
      unlist(parallel::mclapply(1:nrow(ut$edge), function(e) {
        opt.fun <- function(x) f(x, ut$edge[e, 1], ut$edge[e, 
                                                           2])
        optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
      }, mc.cores = ncpu))
    else apply(ut$edge, 1, function(e) {
      opt.fun <- function(x) f(x, e[1], e[2])
      optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
    })
    obj.edge[ ut$edge.length<=0 ]  <- -Inf # excludes nodes with zero branch parents 
    
    best.edges <- order( obj.edge, decreasing=TRUE)[1:topx]
    #for (best.edge in best.edges )
    lapply(best.edges , function(best.edge){
      best.edge.parent <- ut$edge[best.edge, 1]
      best.edge.child <- ut$edge[best.edge, 2]
      best.edge.length <- ut$edge.length[best.edge]
      opt.fun <- function(x) f(x, best.edge.parent, best.edge.child)
      best.pos <- optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$maximum
      new.root <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new.root", 
                       edge.length = 1, Nnode = 1L, root.edge = 1)
      class(new.root) <- "phylo"
      ut2 <- ape::bind.tree(ut, new.root, where = best.edge.child, position = best.pos * 
                         best.edge.length)
      ut2 <- ape::collapse.singles(ut2)
      ut2 <- ape::root(ut2, "new.root")
      x <- drop.tip(ut2, "new.root")
      if (!is.rooted(x)) return(NULL)
      x
    })-> tres
    tres[ !sapply( tres, is.null) ]
    }
