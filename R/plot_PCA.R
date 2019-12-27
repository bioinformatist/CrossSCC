#' Perform PCA analysis and plot result
#'
#' @param x CrossSCC result.
#'
#' @return a PCA plot
#' @export
#'
#' @examples
#' data('cl.b1')
#' handsome.zuo <- CrossSCC(cl.b1[1:500,], ncores = 10, mean.posterior.cutoff = 0.18)
#' plot_PCA(handsome.zuo)
plot_PCA <- function(x) {
  tmp.list <- x$Get('sampleNames', filterFun = isNotLeaf, simplify = FALSE)
  sample.names <- unique(unlist(tmp.list, use.names = FALSE))
  tmp.mat <- matrix(nrow = length(sample.names), ncol = length(tmp.list),
                    dimnames = list(unique(sample.names), names(tmp.list)))
  for (i in names(tmp.list)) {
    for (j in seq_len(length(tmp.list[[i]]))) {
      for (k in tmp.list[[i]][[j]]) {
        tmp.mat[k, i] <- ifelse(j == 2, 0, 1)
      }
    }
  }

  tmp.mat[is.na(tmp.mat)] <- 9
  # https://www.jianshu.com/p/f15625700b3b
  pca <- prcomp(tmp.mat)

  percentage <- round(pca$sdev / sum(pca$sdev) *100,2)
  percentage <- paste0(colnames(pca),"(", paste0(as.character(percentage),"%",")", sep=""))

  ggplot(data.frame(pca$x),aes(x = PC1,y = PC2)) + geom_point() + geom_jitter() + xlab(paste0('PC1 ', percentage[1])) +ylab(paste0('PC2 ', percentage[2]))
}
