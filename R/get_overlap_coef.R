get_overlap_coef <- function(mu1, mu2, sd1, sd2){
  xs  <- seq(min(mu1 - 4 * sd1, mu2 - 4 * sd2),
             max(mu1 + 4 * sd1, mu2 + 4 * sd2),
             length.out = 500)
  f1  <- dnorm(xs, mean = mu1, sd = sd1)
  f2  <- dnorm(xs, mean = mu2, sd = sd2)
  int <- xs[which.max(pmin(f1, f2))]
  l   <- pnorm(int, mu1, sd1, lower.tail = mu1 > mu2)
  r   <- pnorm(int, mu2, sd2, lower.tail = mu1 < mu2)
  l + r
}
