# @todo Clean up existed GO terms and build a database: remove terms if it has too many genes (>= 20% of all genes); Use Jaccard similarity coefficient to remove duplicated terms (> 0.9, keep the smaller one).
# @todo Perform PCA on real data first: only genes in PCs with significant variance (PC1 + PC2 + ... + PCn >= 90%) will be kept.
# @todo Support user-defined GO terms and geneset relationship.
# @todo Geneset score: Mean, median, weighted PCA score and GSVA enrichment score.
as_go <- function(m, var.cutoff = 0.85, ontos = 'BP', mapping = "org.Hs.eg.db", verbose = FALSE, show.progress.bar = TRUE) {
  if (is(m, "ExpressionSet")) {
    m <- m
  } else if (is(m, "matrix")) {
    # m <- new("ExpressionSet", m = as.matrix(m), annotation = 'org.Hs.eg.db')
    m <- Biobase::ExpressionSet(m)
  } else {
    stop('Input data of CrossSCC must be ExpressionSet or matrix!')
  }
  
  m <- Biobase::exprs(genefilter::varFilter(m, var.cutoff = var.cutoff))

  # Perform filtering based on PCA result
  # m.pca <- prcomp(t(m))
  # m.summary <- summary(m.pca)
  # loadings <- m.pca$rotation[, seq_len(which(m.summary$importance[3,] > 0.9)[1])]
  # loadings.logic <- apply(abs(loadings) > 0.05, 1, any)
  # m <- m[names(loadings.logic[loadings.logic == TRUE]),]

  verbose && header(verbose, 'Performing Gene Ontology based meta-gene annotation using ontology level(s): ',
                    emphasize(ontos), char = emojifont::emoji("sunflower"))

  GO.uniq <- unique(hGO, by = 'GO')[, GO]
  setkey(hGO, GO, Onto)
  geneset <- unlist(lapply(ontos,
                           function(o) lapply(GO.uniq,
                                              function(x) tryCatch(m[hGO[.(x, o), GID], ],
                                                                   error = function(e) NULL))),
                    recursive = FALSE)
  names(geneset) <- GO.uniq
  geneset <- Filter(function(x) !is.null(x), geneset)
  geneset <- lapply(geneset, function(x) apply(x, 2, function(y) .Internal(mean(y))))
  symbols <- suppressMessages(AnnotationDbi::mapIds(get(mapping), keys = rownames(m), keytype = "ENTREZID", column="SYMBOL", multiVals = 'first'))
  rownames(m) <- unname(symbols)
  rbind(do.call(rbind, geneset), m)
}
