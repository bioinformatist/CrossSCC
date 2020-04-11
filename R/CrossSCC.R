#' CrossSCC: An Single-Cell Clustering analysis framework Crossing batches of data
#'
#' @docType package
#' @name CrossSCC
#' @import data.tree R.utils visNetwork data.table
#' @importFrom graphics plot
#' @importFrom methods is new
#' @importFrom stats dnorm pnorm setNames
#' @importFrom utils capture.output
#' @importFrom R.oo isVisible
#' @importFrom crayon green bold
#' @importFrom org.Hs.eg.db org.Hs.eg.db
NULL

#' CrossSCC
#'
#' @param m a matrix of single cell expression values.
#' @param ncores number of CPU cores used.
#' @param verbose verbose level. By default, CrossSCC will output all logs as well as progress bars.
#' @param show.progress.bar Set to FALSE if you don't want to see progress bar.
#' @param log.file a test/debug use option. Messages will be directed to this new text file instead of stderr().
#' Note that verbose threshold will be set to -1 with this parameter.
#' @param markers a dataframe contains marker information. Try examples for details.
#'
#' @return a data.tree object.
#' @export
#'
#' @examples
#' library(data.table)
#' markers <- fread(system.file("extdata", "markers.csv", package = "CrossSCC"))
#' data(immu)
#' handsome.zuo <- CrossSCC(immu, markers, ncores = 16)
CrossSCC <- function(m, markers, ncores = 16, 
                     verbose = R.utils::Verbose(threshold = -1, timestamp = TRUE), show.progress.bar = TRUE,
                     log.file = NULL) {
  message('Note: if you met error message as: 
          "Error in serialize(data, node$con) : error writing to connection"
          Please restart R session and try again.
          It may caused by interrupt during previous running since CrossSCC processes data concurrently by default.')
  # Progress bar should be turned off if verbose is set to FALSE
  if (!verbose) {
    show.progress.bar <- FALSE
  }
  
  # If logs is directed to a file, prgress bar should be closed
  if (!is.null(log.file)) {
    verbose <- R.utils::Verbose(con = file(log.file, open = 'a+'), threshold = -1, timestamp = TRUE, removeFile = FALSE)
    show.progress.bar <- FALSE
  }

  verbose && newline(verbose)
  verbose && ruler(verbose, char = emojifont::emoji("dash"), length = 40)
  verbose && newline(verbose)

  verbose && header(verbose, 'Starting tree splitting',
                    char = emojifont::emoji("sunflower"))

  assigned <- rate_sample(m, ncores = ncores)
  
  is.root <- TRUE
  
  for (type in colnames(markers)) {
    assigned.T <- assigned[markers$type]
    assigned.T <- Filter(function(x) !is.null(x), assigned.T)
    snowfall::sfRemoveAll(except = 'm')
    snowfall::sfExport('assigned.T')
    assigned.T <- pbapply::pblapply(colnames(m),
                                    function(x) vapply(assigned.T,
                                                       function(z) vapply(z,
                                                                          function(y) x %in% y,
                                                                          logical(1)),
                                                       logical(2)), cl = cl)
    names(assigned.T) <- colnames(m)
    ownership <- vapply(assigned.T, function(x) names(which.max(apply(x, 1, sum))), character(1))
    if (is.root) {
      result <- Node$new('All', sampleNames = colnames(m))
      is.root <- FALSE
      result$AddChild(type, sampleNames = colnames(m)[which(ownership == "comp.expressed")])
      no.type <- paste0('non', type)
      result$AddChild(no.type, sampleNames = colnames(m)[which(ownership == "comp.no")])
      m <- m[, colnames(m)[which(ownership == "comp.no")]]
    } else {
      pointer <- Traverse(result, filterFun = function(x) x$name == no.type)[[1]]
      pointer$AddChild(type, sampleNames = colnames(m)[which(ownership == "comp.expressed")])
      no.type <- paste0('non', type)
      result$AddChild(no.type, sampleNames = colnames(m)[which(ownership == "comp.no")])
      m <- m[, colnames(m)[which(ownership == "comp.no")]]
    }
  }
  
  snowfall::sfStop()
  verbose && newline(verbose)
  verbose && header(verbose, 'CrossSCC finished!
                    Try with another member in "Cross" family:
                    CrossICC: An Interactive Consensus Clustering framework for Cross-platform data analysis
                    See https://bioconductor.org/packages/CrossICC',
                    char = emojifont::emoji("sunflower"))

  # Prune only returns removed tree number
  result
}

# https://stackoverflow.com/a/39926819
# To avoid using <<-, I use function parameter instead of "global variables"
render_tree <- function(m, ncores, 
                         nnode = 1, decision.node = NULL, pending.node = 0,
                         result, verbose = FALSE, show.progress.bar = TRUE) {
  # To return with sample names instead of indices
  names(rated.feature) <- rownames(m)
  rated.feature <- Filter(function(x) !is.null(x), rated.feature)
  if (!length(rated.feature)) {
    verbose && enter(verbose, 'Trying split decision node: ', emphasize(decision.node),
            ', working on its child (decision/terminal) node: ', emphasize('#'), emphasize(nnode), indent = 0)
    verbose && enter(verbose, '\tCannot find proper representative feature', indent = 0)
    return(result)
  }
  good.feature <- Filter(function(x) x[['type']] == 'good', rated.feature)

  if (length(good.feature)) {
    scores <- vapply(good.feature, function(x) x[['score']], 2333.2333)
    best.name <- names(which.max(scores))
    best.feature <- good.feature[[best.name]]
    # First node is root
    if (is.null(decision.node)) {
      if (is.null(min.group.size)) {
        min.group.size <- min.group.ratio * NCOL(m)
      }

      verbose && enter(verbose, emphasize(best.name), ' was chosen as root node', indent = 0)
      # Initialize whole tree here
      result <- Node$new(best.name, sampleNames = list(best.feature[['comp.1']], best.feature[['comp.2']]))

      if (length(best.feature[['comp.1']]) > min.group.size & length(best.feature[['comp.2']]) > min.group.size) {
        rank_feature(m[, best.feature[['comp.1']]], ncores = ncores,
                     decision.node = best.name, nnode = 1, result = result$root,
                     mean.posterior.cutoff = mean.posterior.cutoff, ovl.cutoff = ovl.cutoff,
                     mean.posterior.weight= mean.posterior.weight, min.group.ratio = min.group.ratio,
                     ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff, min.group.size = min.group.size,
                     verbose = verbose, show.progress.bar = show.progress.bar)
        
        rank_feature(m[, best.feature[['comp.2']]], ncores = ncores,
                     decision.node = best.name, nnode = 2, result = result$root,
                     mean.posterior.cutoff = mean.posterior.cutoff, ovl.cutoff = ovl.cutoff,
                     mean.posterior.weight= mean.posterior.weight, min.group.ratio = min.group.ratio,
                     ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff, min.group.size = min.group.size,
                     verbose = verbose, show.progress.bar = show.progress.bar)
      }
    } else {
      verbose && enter(verbose, 'Trying split decision node: ', emphasize(decision.node),
              ', working on its child (decision/terminal) node: ', emphasize('#'), emphasize(nnode), indent = 0)
      verbose && enter(verbose, '\tDetermined as ', emphasize('decision node'), '. ',
                       emphasize(best.name), ' was chosen as representative feature', indent = 0)
      # Should use Traverse() here for FindNode() can only return the 1st node who matches
      # CAN NOT use Do() for all nodes kept after filtering will be performed same operation
      pointer <- Traverse(result$root, traversal = "post-order",
                          filterFun = function(x) x$name == decision.node & length(x$children) < 2)[[1]]
      
      pointer$AddChild(best.name,
                       sampleNames = list(best.feature[['comp.1']],
                                          best.feature[['comp.2']]))
      

      if (length(best.feature[['comp.1']]) > min.group.size & length(best.feature[['comp.2']]) > min.group.size) {
        pending.node <- pending.node + 1
        # Must use Traverse() here, for result is a reference
        rank_feature(m[, best.feature[['comp.1']]], ncores = ncores,
                     decision.node = best.name, nnode = 1,
                     # Traverse() will return all matched nodes, we use the 1st one in according to "post-order"
                     # post-order means from bottom to up, but we still need "right-to-left" by rev()
                     result = result,
                     mean.posterior.cutoff = mean.posterior.cutoff, ovl.cutoff = ovl.cutoff,
                     mean.posterior.weight= mean.posterior.weight, min.group.ratio = min.group.ratio,
                     ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff, min.group.size = min.group.size,
                     verbose = verbose, show.progress.bar = show.progress.bar)
        pending.node <- pending.node - 1

        pending.node <- pending.node + 1
        rank_feature(m[, best.feature[['comp.2']]], ncores = ncores,
                     decision.node = best.name, nnode = 2,
                     result = result,
                     mean.posterior.cutoff = mean.posterior.cutoff, ovl.cutoff = ovl.cutoff,
                     mean.posterior.weight= mean.posterior.weight, min.group.ratio = min.group.ratio,
                     ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff, verbose = verbose,
                     min.group.size = min.group.size,
                     show.progress.bar = show.progress.bar)
        pending.node <- pending.node - 1
      } else {
        terminal.pointer <- Traverse(pointer, filterFun = function(x) x$name == best.name)[[1]]
        lapply(seq_along(terminal.pointer$sampleNames),
               function(i) terminal.pointer$AddChild(paste0(terminal.pointer$name, '\nsub-component ', i),
                                                     sampleNames = terminal.pointer$sampleNames[[i]]))
      }

      if (pending.node == 0) {
        # Should NOT call Prune() here, for it may not the last decision node
        return(result)
      }
    }
  } else {
    # Those have no good features at beginning
    if (is.null(decision.node)) {
      stop('Got too few feature for clustering. Please provide a bigger matrix with more features!')
    }
    only.feature <- Filter(function(x) x[['type']] == 'only', rated.feature)
    scores <- vapply(only.feature, function(x) x[['score']], 2333.2333)
    representative.name <- names(which.max(scores))
    verbose && enter(verbose, 'Trying split decision node: ', emphasize(decision.node),
            ', working on its child (decision/terminal) node: ', emphasize('#'), emphasize(nnode), indent = 0)
    verbose && enter(verbose, '\tDetermined as ', emphasize('terminal node'), '. ',
            emphasize(representative.name), ' was chosen as representative feature', indent = 0)
    representative.feature <- only.feature[[representative.name]]
    Traverse(result$root, traversal = "post-order",
             filterFun = function(x) x$name == decision.node & length(x$children) < 2)[[1]]$AddChild(representative.name, sampleNames = colnames(m))
    result
  }
}
