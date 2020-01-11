# @todo Clean up existed GO terms and build a database: remove terms if it has too many genes (>= 20% of all genes); Use Jaccard similarity coefficient to remove duplicated terms (> 0.9, keep the smaller one).
# @todo Perform PCA on real data first: only genes in PCs with significant variance (PC1 + PC2 + ... + PCn >= 90%) will be kept.
# @todo Support user-defined GO terms and geneset relationship.
# @todo Geneset score: Mean, median, weighted PCA score and GSVA enrichment score.
as_go <- function(m, ncores = 4, var.cutoff = 0.9, ontos = 'BP', simple.mapping = 'org.HsSimple.eg.db', mapping = "org.Hs.eg.db", verbose = FALSE, show.progress.bar = TRUE) {
  if (is(m, "ExpressionSet")) {
    m <- m
  } else if (is(m, "matrix")) {
    m <- new("ExpressionSet", m = as.matrix(m), annotation = 'org.Hs.eg.db')
  } else {
    stop('Input data of CrossSCC must be ExpressionSet or matrix!')
  }

  ontos <- ontos

  if (show.progress.bar) {
    pbo <- pbapply::pboptions(type = "timer", use_lb = TRUE,
                            char = emojifont::emoji('rocket'),
                            txt.width = 30)
  } else {
    pbo <- pbapply::pboptions(type = "none", use_lb = TRUE)
  }

  on.exit(pbapply::pboptions(pbo))

  if (ncores == 'all') {
    ncores <- parallel::detectCores(all.tests = TRUE, logical = TRUE) - 1
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
  verbose && enter(verbose, 'Using ', emphasize(ncores), ' core(s)...\nInitializing... DO NOT STOP AT THIS STEP', indent = 0, suffix = '!')
  capture.output(suppressMessages(snowfall::sfInit(parallel = TRUE, cpu = ncores)))
  # For cl parameter of pbapply (see source code for details)
  cl <- snowfall::sfGetCluster()

  annot_onto_Wrapper <- function(x, onto) {
    names(topGO::annFUN.org(onto, feasibleGenes = x, mapping = simple.mapping, ID = "entrez"))
  }

  snowfall::sfExport('ontos', 'm', 'annot_onto_Wrapper', 'mapping')

  verbose && enter(verbose, 'Initializing finished', indent = 0, suffix = '!')
  gene2go <- lapply(ontos, function(x) pbapply::pblapply(rownames(m), function(y) annot_onto_Wrapper(y, x), cl = cl))

  suppressMessages(snowfall::sfStop())

  for (i in seq_len(length(ontos))) {
    names(gene2go[[i]]) <- rownames(m)
  }

  gene2go <- unlist(gene2go, recursive = FALSE)

  # https://stackoverflow.com/a/35163845
  gene2go <- setNames(unlist(gene2go, use.names = FALSE), rep(names(gene2go), lengths(gene2go)))

  go2gene <- list()

  # lapply + <<- usually implies you should be using a for loop. – hadley
  # Ref: https://stackoverflow.com/q/14321517
  for (i in seq_len(length(gene2go))) {
    go2gene[[gene2go[i]]] <- c(go2gene[[gene2go[i]]], names(gene2go)[i])
  }

  m.go <- t(vapply(go2gene, function(x) apply(subset(m, rownames(m) %in% x), 2, function(x) .Internal(mean(x))), rep(2333.2333, ncol(m))))
  symbols <- suppressMessages(AnnotationDbi::mapIds(get(mapping), keys = rownames(m), keytype = "ENTREZID", column="SYMBOL", multiVals = 'first'))
  rownames(m) <- unname(symbols)
  rbind(m.go, m)
}
