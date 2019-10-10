rate_feature <- function(feature, sample.names, mean.posterior.cutoff = 0.1, ovl.cutoff = 0.05, mean.posterior.weight = 0.5, ovl.weight = 0.5, lambda.cutoff = 0.9) {
  # We only use GMM with 2 components as more components may cause overfitting
  # 2 components may not be the best fitting model, thus one will cause variance closing to 0 and Inf log-likelihood
  # Must set a small-enough epsilon for robustness
  capture.output(a.optim <- tryCatch(mixtools::normalmixEM(feature, k = 2, epsilon = 1e-16), error=function(e) NULL))
  # No good boy
  if (is.null(a.optim)) {
    return()
  }

  mean.posterior <- apply(a.optim$posterior, 2, function(x) .Internal(mean(x)))
  mean.posterior.check <- all(mean.posterior > mean.posterior.cutoff)

  ovl <- get_overlap_coef(a.optim$mu[1], a.optim$mu[2], a.optim$sigma[1], a.optim$sigma[2])
  ovl.check <- ovl < ovl.cutoff

  # Check if all components are good boys
  if (all(mean.posterior.check, ovl.check)) {
    # Determine sample assignment
    max.comp <- apply(a.optim$posterior, 1, which.max)
    comp.1 <- sample.names[which(max.comp == 1)]
    comp.2 <- sample.names[which(max.comp == 2)]

    mean.posterior.score <- sum(mean.posterior - mean.posterior.cutoff)
    ovl.score <- ovl.cutoff - ovl
    score <- mean.posterior.score * mean.posterior.weight + ovl.score * ovl.weight

    list(type = 'good',
         score = score,
         comp.1 = comp.1,
         comp.2 = comp.2)
    # For those have only one good boy
  } else if (a.optim$lambda[1] > lambda.cutoff) {
    max.comp <- apply(a.optim$posterior, 1, which.max)
    comp.1 <- sample.names[which(max.comp == 1)]
    score <- a.optim$lambda[1] - lambda.cutoff
    list(type = 'only', score = score, comp = comp.1)
  } else if (a.optim$lambda[2] > lambda.cutoff) {
    max.comp <- apply(a.optim$posterior, 1, which.max)
    comp.2 <- sample.names[which(max.comp == 2)]
    score <- a.optim$lambda[2] - lambda.cutoff
    list(type = 'only', score = score, comp = comp.2)
    # No good boy
  } else {
    return()
  }
}
