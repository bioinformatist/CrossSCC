assign_sample <- function(x, sample.names) {
  min.index <- which.min(x$sigma)
  max.comp <- apply(x$posterior, 1, which.max)
  list(
    comp.no = sample.names[which(max.comp == min.index)],
    comp.expressed = sample.names[which(max.comp != min.index)]
  )
}

rate_sample <- function(m, ncores) {
  # capture.output(suppressMessages(snowfall::sfInit(parallel = TRUE, cpu = ncores)))
  snowfall::sfInit(parallel = TRUE, cpu = ncores)
  cl <- snowfall::sfGetCluster()
  snowfall::sfExport('m')
  snowfall::sfExport('normalmixEM', namespace = 'mixtools')
  rated.features <- pbapply::pbapply(m, 1, function(x) tryCatch(normalmixEM(x, k = 2), error=function(e) NULL), cl = cl)
  snowfall::sfRemove('normalmixEM')
  filter.rules <- function(x) !is.null(x) & sum(between(x$sigma, 0, 1)) == 1 & sum(between(x$mu, 0, 1)) == 1
  filtered.features <- Filter(filter.rules, rated.features)
  snowfall::sfExport('filtered.features', 'assign_sample')
  (assigned <- pbapply::pblapply(filtered.features, assign_sample, sample.names = colnames(m), cl = cl))
}