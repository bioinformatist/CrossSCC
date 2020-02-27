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
#' @param var.cutoff cutoff value when filtering features.
#' @param mapping database package name according to organism. e.g. "org.Hs.eg.db".
#' @param mean.posterior.cutoff cutoff value for mean posterior of each component when rating Gaussion Mixture Model.
#' @param ovl.cutoff cutoff value for OVL (overlapping coefficient) when rating Gaussion Mixture Model.
#' @param mean.posterior.weight the weight for mean.posterior.cutoff when assessing the performance of the model.
#' @param ovl.weight the weight for ovl.cutoff when assessing the performance of the model.
#' @param lambda.cutoff Gaussian components with lambda over this cutoff value
#' will be classified as a model representing a subtype of samples.
#' @param ontos ontology terms used when converting features to meta-gene features.
#' Default is "BP", could be arbitrary combination of c('BP', 'MF', 'CC').
#' @param min.group.ratio minimal final group ratio (proportion of all sample number).
#' At each decision node, only groups has more samples than this threshold will be further splitted.
#' @param verbose verbose level. By default, CrossSCC will output all logs as well as progress bars.
#' @param show.progress.bar Set to FALSE if you don't want to see progress bar.
#' @param min.group.size minimal final group size. Note: this parameter will overwrite parameter min.group.ratio.
#' @param play.leaves a test/debug use option. Default is TRUE.
#' @param log.file a test/debug use option. Messages will be directed to this new text file instead of stderr().
#' Note that verbose threshold will be set to -1 with this parameter.
#'
#' @return a data.tree object.
#' @export
#'
#' @examples
#' data('cl.b1')
#' handsome.zuo <- CrossSCC(cl.b1[1:500,], ncores = 10, mean.posterior.cutoff = 0.18)
CrossSCC <- function(m, ncores = 4, var.cutoff = 0.9, mapping = "org.Hs.eg.db",
                     mean.posterior.cutoff = 0.3,
                     ovl.cutoff = 0.05, mean.posterior.weight = 0.5, min.group.ratio = 0.1,
                     ovl.weight = 0.5, lambda.cutoff = 0.9, ontos = 'BP', min.group.size = NULL,
                     verbose = R.utils::Verbose(threshold = -1), show.progress.bar = TRUE,
                     play.leaves = TRUE, log.file = NULL) {
  # Progress bar should be turned off if verbose is set to FALSE
  if (!verbose) {
    show.progress.bar <- FALSE
  }
  
  # If logs is directed to a file, prgress bar should be closed
  if (log.file) {
    verbose <- R.utils::Verbose(con = file(log.file), threshold = -1)
    show.progress.bar <- FALSE
  }

  m <- as_go(m, ncores = ncores, var.cutoff = var.cutoff, ontos = ontos,
             mapping = mapping, verbose = verbose, show.progress.bar = show.progress.bar)

  verbose && newline(verbose)
  verbose && ruler(verbose, char = emojifont::emoji("dash"), length = 40)
  verbose && newline(verbose)

  verbose && header(verbose, 'Starting tree splitting',
                    char = emojifont::emoji("sunflower"))

  result <- rank_feature(m, ncores = ncores, mean.posterior.cutoff = mean.posterior.cutoff,
                         ovl.cutoff = ovl.cutoff, mean.posterior.weight= mean.posterior.weight,
                         min.group.ratio = min.group.ratio, min.group.size = min.group.size,
                         ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff,
                         result = NULL, verbose = verbose, show.progress.bar = show.progress.bar)

  if (play.leaves) {
    # Cut redundant nodes
    Prune(result, function(x) length(x$siblings))
    # Expand leaves as two components
    result$Do(function(x) lapply(seq_along(x$sampleNames),
                                 function(i) x$AddChild(paste0(x$name, '\nsub-component ', i),
                                                        sampleNames = x$sampleNames[[i]])),
              filterFun = function(x) isLeaf(x) & is.list(x$sampleNames))
  }
  
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
rank_feature <- function(m, ncores, mean.posterior.cutoff, var.cutoff, ovl.cutoff,
                         mean.posterior.weight, ovl.weight, lambda.cutoff, min.group.ratio,
                         nnode = 1, decision.node = NULL, pending.node = 0, min.group.size,
                         result, verbose = FALSE, show.progress.bar = TRUE) {

  # Matrix with too few features need NOT paralleled computing
  # Must use NROW() instead of nrow() for latter may produce NULL
  # https://stackoverflow.com/questions/27674937/why-do-ncol-and-nrow-only-yield-null-when-i-do-have-data
  if (NROW(m) < 400 | ncores == 1) {
    rated.feature <- lapply(split(m, seq(NROW(m))), rate_feature, sample.names = colnames(m),
                            mean.posterior.cutoff = mean.posterior.cutoff,
                            ovl.cutoff = ovl.cutoff, mean.posterior.weight= mean.posterior.weight,
                            ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff)
  } else {
    if (show.progress.bar) {
      pbo <- pbapply::pboptions(type = "timer", use_lb = TRUE,
                         char = emojifont::emoji('rocket'),
                         txt.width = 30)
    } else {
      pbo <- pbapply::pboptions(type = "none", use_lb = TRUE)
    }

    on.exit(pbapply::pboptions(pbo))

    expected.cores <- floor(NROW(m) / 400) + 1
    if (expected.cores > ncores) {
      expected.cores <- ncores
    }
    verbose && enter(verbose, 'Rating features on ', emphasize(expected.cores), ' cores', indent = 0)
    capture.output(suppressMessages(snowfall::sfInit(parallel = TRUE, cpu = expected.cores)))
    cl <- snowfall::sfGetCluster()

    rate_feature <- rate_feature
    snowfall::sfExport('m', 'rate_feature')
    rated.feature <- pbapply::pblapply(split(m, seq(NROW(m))), rate_feature, sample.names = colnames(m),
                                       mean.posterior.cutoff = mean.posterior.cutoff,
                                       ovl.cutoff = ovl.cutoff, mean.posterior.weight= mean.posterior.weight,
                                       ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff, cl = cl)
    suppressMessages(snowfall::sfStop())
  }

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
      Traverse(result$root,traversal = "post-order",
               filterFun = function(x) x$name == decision.node)[[1]]$AddChild(best.name,
                                                                              sampleNames = list(best.feature[['comp.1']],
                                                                                                 best.feature[['comp.2']]))

      if (length(best.feature[['comp.1']]) > min.group.size & length(best.feature[['comp.2']]) > min.group.size) {
        pending.node <- pending.node + 1
        # Must use Traverse() here, for result is a reference
        rank_feature(m[, best.feature[['comp.1']]], ncores = ncores,
                     decision.node = best.name, nnode = 1,
                     # Traverse() will return all matched nodes, we use the 1st one in according to "post-order"
                     result = Traverse(result$root,traversal = "post-order",
                                       filterFun = function(x) x$name == best.name)[[1]],
                     mean.posterior.cutoff = mean.posterior.cutoff, ovl.cutoff = ovl.cutoff,
                     mean.posterior.weight= mean.posterior.weight, min.group.ratio = min.group.ratio,
                     ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff, min.group.size = min.group.size,
                     verbose = verbose, show.progress.bar = show.progress.bar)
        pending.node <- pending.node - 1

        pending.node <- pending.node + 1
        rank_feature(m[, best.feature[['comp.2']]], ncores = ncores,
                     decision.node = best.name, nnode = 2,
                     result = Traverse(result$root,traversal = "post-order",
                                       filterFun = function(x) x$name == best.name)[[1]],
                     mean.posterior.cutoff = mean.posterior.cutoff, ovl.cutoff = ovl.cutoff,
                     mean.posterior.weight= mean.posterior.weight, min.group.ratio = min.group.ratio,
                     ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff, verbose = verbose,
                     min.group.size = min.group.size,
                     show.progress.bar = show.progress.bar)
        pending.node <- pending.node - 1
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
    result$AddChild(representative.name, sampleNames = colnames(m))
    result
  }
}
