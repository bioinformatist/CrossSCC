
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GMMClassifier

The goal of GMMClassifier is to classify single cells (as expression
matrix of scRNA-seq data) into subtypes, and to find out the mapping
relationship among subtypes across datasets (of different batch).

## Installation

You can install the latest version of GMMClassifier with:

``` r
install.packages("devtools")
devtools::install_github("GMMClassifier")
```

## About example data

The example data used in this project is part of
[GSE81861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81861)
of [this paper](https://www.nature.com/articles/ng.3818#accessions).
First, the FPKM matrix file of all cells was
downloaded:

``` bash
aria2c ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81861/suppl/GSE81861_Cell_Line_FPKM.csv.gz
pigz -d GSE81861_Cell_Line_FPKM.csv.gz
```

The matrix was read into R and the columns of two batch were selected as
subsets according to [the results part of the
paper](https://www.nature.com/articles/ng.3818#results) “To assess batch
effects, we performed scRNA–seq in two batches for GM12878
(lymphoblastoid) cells and also for H1 embryonic stem cells. Gene
expression was quantified as fragments per kilobase per million reads
(FPKM), and low-quality cells were discarded on the basis of multiple
metrics”:

``` r
library(data.table)
library(tidyverse)
library(org.Hs.eg.db)
library(usethis)
library(genefilter)
library(Biobase)
rowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
cl <- fread('GSE81861_Cell_Line_FPKM.csv')
names(cl)[1] <- 'Gene'
cl <- cl[, Ensembl := str_match(Gene, ".+_(.+)\\.\\d+$")[, 2]]
symbols <- mapIds(org.Hs.eg.db, keys = cl$Ensembl, keytype = "ENSEMBL", column="ENTREZID", multiVals = 'first')
cl <- cl[, Entrez := symbols[Ensembl]]
cl <- cl[, c('Gene', 'Ensembl'):= NULL]
cl.b1 <- cbind(cl[, .(Entrez)], cl[, .SD, .SDcols = names(cl) %like% "_B1_"])
cl.b2 <- cbind(cl[, .(Entrez)], cl[, .SD, .SDcols = names(cl) %like% "_B2_"])
cl <- lapply(list(cl.b1, cl.b2), function(x) x[, var := rowVars(.SD), .SDcols = -c('Entrez')])
cl <- lapply(cl, function(x) x[, max.var := max(var), by = 'Entrez'][var != 0 & max.var == var & !is.na(Entrez)])
cl <- lapply(cl, function(x) x[, grep("var", colnames(x)) := NULL])
cl <- lapply(cl, function(x) as.data.frame(x) %>% remove_rownames %>% column_to_rownames(var="Entrez"))
cl <- lapply(cl, function(x) new("ExpressionSet", exprs=as.matrix(x), annotation = 'org.Hs.eg.db'))
cl <- lapply(cl, function(x) nsFilter(x, require.entrez = TRUE, require.GOBP = TRUE, require.GOCC = TRUE, require.GOMF = TRUE, var.cutoff = 0.999))
cl.b1 <- exprs(cl[[1]][[1]])
cl.b2 <- exprs(cl[[2]][[1]])
use_data(cl.b1, compress = 'xz')
use_data(cl.b2, compress = 'xz')
```

Thus, the `matrix` object `cl.b1` contains *GM12878* and *H1* cell of
batch 1, and `cl.b2` contains those of batch 2.

## How to use?
