
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CrossSCC

The goal of CrossSCC is to classify **S**ingle-**C**ell data(as
expression matrix of scRNA-seq data) **Cross**ing batch into **C**lusters
using Gaussian Mixture Model, and to find out the mapping relationship
in clusters.

![](man/figures/readme.gif)

## Installation

You can install the latest version of CrossSCC with:

``` r
install.packages("remotes")
remotes::install_github("bioinformatist/CrossSCC")
```

## Example dataset

The example data used in this project is part of
[GSE81861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81861)
from [this paper](https://www.nature.com/articles/ng.3818#accessions).
First, the FPKM matrix file of all cells was
downloaded:

``` bash
aria2c ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81861/suppl/GSE81861_Cell_Line_FPKM.csv.gz
pigz -d GSE81861_Cell_Line_FPKM.csv.gz
```

### Two batches

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
cl[, Ensembl := str_match(Gene, ".+_(.+)\\.\\d+$")[, 2]]
symbols <- mapIds(org.Hs.eg.db, keys = cl$Ensembl, keytype = "ENSEMBL", column="ENTREZID", multiVals = 'first')
cl[, Entrez := symbols[Ensembl]]
cl[, c('Gene', 'Ensembl'):= NULL]
cl.b1 <- cbind(cl[, .(Entrez)], cl[, .SD, .SDcols = names(cl) %like% "_B1_"])
cl.b2 <- cbind(cl[, .(Entrez)], cl[, .SD, .SDcols = names(cl) %like% "_B2_"])
cl <- lapply(list(cl.b1, cl.b2), function(x) x[, var := rowVars(.SD), .SDcols = -c('Entrez')])
cl <- lapply(cl, function(x) x[, max.var := max(var), by = 'Entrez'][var != 0 & max.var == var & !is.na(Entrez), ])
cl <- lapply(cl, function(x) x[, grep("var", colnames(x)) := NULL])
cl <- lapply(cl, function(x) as.data.frame(x) %>% remove_rownames %>% column_to_rownames(var = "Entrez"))
cl <- lapply(cl, function(x) new("ExpressionSet", exprs=as.matrix(x), annotation = 'org.Hs.eg.db'))
cl.b1 <- cl[[1]]
cl.b2 <- cl[[2]]
# setwd() back to package root directory
use_data(cl.b1, compress = 'xz')
use_data(cl.b2, compress = 'xz')
```

Thus, the `matrix` object `cl.b1` contains *GM12878* and *H1* cell of
batch 1, and `cl.b2` contains those of batch 2.

### Five types in same batch

``` r
# Package library omitted
cl <- fread('GSE81861_Cell_Line_FPKM.csv')
names(cl)[1] <- 'Gene'
cl[, Ensembl := str_match(Gene, ".+_(.+)\\.\\d+$")[, 2]]
symbols <- mapIds(org.Hs.eg.db, keys = cl$Ensembl, keytype = "ENSEMBL", column="ENTREZID", multiVals = 'first')
cl[, Entrez := symbols[Ensembl]]
cl[, c('Gene', 'Ensembl'):= NULL]
cl <- cbind(cl[, .(Entrez)], cl[, .SD, .SDcols = !(names(cl) %like% "_B[[:digit:]]{1}_")])
cl[, var := rowVars(.SD), .SDcols = -c('Entrez')]
cl <- cl[, max.var := max(var), by = 'Entrez'][var != 0 & max.var == var & !is.na(Entrez), ]
cl[, grep("var", colnames(cl)) := NULL]
cl <- as.data.frame(cl) %>% remove_rownames %>% column_to_rownames(var = "Entrez")
cl <- cl[, !(names(cl) %in% 'Entrez')]
cl.no.batch <- new("ExpressionSet", exprs=as.matrix(cl), annotation = 'org.Hs.eg.db')
# setwd() back to package root directory
use_data(cl.no.batch, compress = 'xz')
```

`cl.no.batch` contains samples with five types in same batch.

## Usage && Assessment

### On [Batch 1](#two-batches)

``` r
library(CrossSCC)
data("cl.b1")
# Fixed seed means fixed mu and sigma during normalmixEM()
set.seed(920304)
handsome.zuo <- CrossSCC(cl.b1, ncores = 56, mean.posterior.cutoff = 0.31, verbose = FALSE)
```

To visualize the result (currently still too simple):

``` r
plot(handsome.zuo)
```

As we known, samples should be devided into two types, so just check
nodes at level 2 now:

``` r
library(data.tree)
#> 
#> Attaching package: 'data.tree'
#> The following object is masked from 'package:Biobase':
#> 
#>     Aggregate
handsome.zuo$Get('sampleNames', filterFun = isLeaf)
#> $`GO:0006958`
#> $`GO:0006958`[[1]]
#>  [1] "RHB1116__GM12878_B1__orange" "RHB1118__GM12878_B1__orange"
#>  [3] "RHB1140__GM12878_B1__orange" "RHB1146__GM12878_B1__orange"
#>  [5] "RHB1153__GM12878_B1__orange" "RHB1165__GM12878_B1__orange"
#>  [7] "RHB1167__GM12878_B1__orange" "RHB1169__GM12878_B1__orange"
#>  [9] "RHB1172__GM12878_B1__orange" "RHB1185__GM12878_B1__orange"
#> [11] "RHB1204__GM12878_B1__orange" "RHB1205__GM12878_B1__orange"
#> [13] "RHB1207__GM12878_B1__orange"
#> 
#> $`GO:0006958`[[2]]
#> [1] "RHB1129__GM12878_B1__orange" "RHB1141__GM12878_B1__orange"
#> [3] "RHB1189__GM12878_B1__orange" "RHB1200__GM12878_B1__orange"
#> [5] "RHB1201__GM12878_B1__orange" "RHB1206__GM12878_B1__orange"
#> [7] "RHB1208__GM12878_B1__orange"
#> 
#> 
#> $`GO:0046294`
#> $`GO:0046294`[[1]]
#> [1] "RHB1122__GM12878_B1__orange" "RHB1124__GM12878_B1__orange"
#> [3] "RHB1158__GM12878_B1__orange" "RHB1173__GM12878_B1__orange"
#> [5] "RHB1193__GM12878_B1__orange" "RHB1197__GM12878_B1__orange"
#> [7] "RHB1199__GM12878_B1__orange" "RHB1125__GM12878_B1__orange"
#> 
#> $`GO:0046294`[[2]]
#> [1] "RHB1147__GM12878_B1__orange" "RHB1162__GM12878_B1__orange"
#> [3] "RHB1171__GM12878_B1__orange" "RHB1188__GM12878_B1__orange"
#> [5] "RHC092__H1_B1__brown"       
#> 
#> 
#> $`GO:0007338`
#>  [1] "RHC069__H1_B1__brown" "RHC070__H1_B1__brown" "RHC071__H1_B1__brown"
#>  [4] "RHC072__H1_B1__brown" "RHC073__H1_B1__brown" "RHC074__H1_B1__brown"
#>  [7] "RHC075__H1_B1__brown" "RHC077__H1_B1__brown" "RHC078__H1_B1__brown"
#> [10] "RHC079__H1_B1__brown" "RHC080__H1_B1__brown" "RHC081__H1_B1__brown"
#> [13] "RHC082__H1_B1__brown" "RHC083__H1_B1__brown" "RHC084__H1_B1__brown"
#> [16] "RHC085__H1_B1__brown" "RHC087__H1_B1__brown" "RHC088__H1_B1__brown"
#> [19] "RHC089__H1_B1__brown" "RHC090__H1_B1__brown" "RHC093__H1_B1__brown"
#> [22] "RHC094__H1_B1__brown" "RHC095__H1_B1__brown" "RHC096__H1_B1__brown"
#> [25] "RHC097__H1_B1__brown" "RHC098__H1_B1__brown" "RHC100__H1_B1__brown"
#> [28] "RHC101__H1_B1__brown" "RHC102__H1_B1__brown" "RHC103__H1_B1__brown"
#> [31] "RHC104__H1_B1__brown" "RHC105__H1_B1__brown" "RHC107__H1_B1__brown"
#> [34] "RHC108__H1_B1__brown" "RHC109__H1_B1__brown" "RHC110__H1_B1__brown"
#> [37] "RHC111__H1_B1__brown" "RHC112__H1_B1__brown" "RHC113__H1_B1__brown"
#> [40] "RHC114__H1_B1__brown" "RHC115__H1_B1__brown" "RHC116__H1_B1__brown"
#> [43] "RHC117__H1_B1__brown" "RHC118__H1_B1__brown" "RHC120__H1_B1__brown"
#> [46] "RHC121__H1_B1__brown" "RHC122__H1_B1__brown" "RHC123__H1_B1__brown"
#> [49] "RHC124__H1_B1__brown" "RHC125__H1_B1__brown" "RHC127__H1_B1__brown"
#> [52] "RHC128__H1_B1__brown" "RHC129__H1_B1__brown" "RHC130__H1_B1__brown"
#> [55] "RHC131__H1_B1__brown" "RHC132__H1_B1__brown" "RHC134__H1_B1__brown"
#> [58] "RHC135__H1_B1__brown" "RHC136__H1_B1__brown" "RHC137__H1_B1__brown"
#> [61] "RHC138__H1_B1__brown" "RHC140__H1_B1__brown" "RHC141__H1_B1__brown"
#> [64] "RHC142__H1_B1__brown" "RHC143__H1_B1__brown" "RHC144__H1_B1__brown"
#> [67] "RHC145__H1_B1__brown" "RHC146__H1_B1__brown"
```

It seems only **one** sample `"RHC092__H1_B1__brown"` was incorrectly
classified.

To list full relationship between cluster and
samples:

``` r
(cl2 <- unname(lapply(rapply(handsome.zuo$Get('sampleNames', filterFun = isLeaf), enquote, how = 'unlist'), eval)))
#> [[1]]
#>  [1] "RHB1116__GM12878_B1__orange" "RHB1118__GM12878_B1__orange"
#>  [3] "RHB1140__GM12878_B1__orange" "RHB1146__GM12878_B1__orange"
#>  [5] "RHB1153__GM12878_B1__orange" "RHB1165__GM12878_B1__orange"
#>  [7] "RHB1167__GM12878_B1__orange" "RHB1169__GM12878_B1__orange"
#>  [9] "RHB1172__GM12878_B1__orange" "RHB1185__GM12878_B1__orange"
#> [11] "RHB1204__GM12878_B1__orange" "RHB1205__GM12878_B1__orange"
#> [13] "RHB1207__GM12878_B1__orange"
#> 
#> [[2]]
#> [1] "RHB1129__GM12878_B1__orange" "RHB1141__GM12878_B1__orange"
#> [3] "RHB1189__GM12878_B1__orange" "RHB1200__GM12878_B1__orange"
#> [5] "RHB1201__GM12878_B1__orange" "RHB1206__GM12878_B1__orange"
#> [7] "RHB1208__GM12878_B1__orange"
#> 
#> [[3]]
#> [1] "RHB1122__GM12878_B1__orange" "RHB1124__GM12878_B1__orange"
#> [3] "RHB1158__GM12878_B1__orange" "RHB1173__GM12878_B1__orange"
#> [5] "RHB1193__GM12878_B1__orange" "RHB1197__GM12878_B1__orange"
#> [7] "RHB1199__GM12878_B1__orange" "RHB1125__GM12878_B1__orange"
#> 
#> [[4]]
#> [1] "RHB1147__GM12878_B1__orange" "RHB1162__GM12878_B1__orange"
#> [3] "RHB1171__GM12878_B1__orange" "RHB1188__GM12878_B1__orange"
#> [5] "RHC092__H1_B1__brown"       
#> 
#> [[5]]
#>  [1] "RHC069__H1_B1__brown" "RHC070__H1_B1__brown" "RHC071__H1_B1__brown"
#>  [4] "RHC072__H1_B1__brown" "RHC073__H1_B1__brown" "RHC074__H1_B1__brown"
#>  [7] "RHC075__H1_B1__brown" "RHC077__H1_B1__brown" "RHC078__H1_B1__brown"
#> [10] "RHC079__H1_B1__brown" "RHC080__H1_B1__brown" "RHC081__H1_B1__brown"
#> [13] "RHC082__H1_B1__brown" "RHC083__H1_B1__brown" "RHC084__H1_B1__brown"
#> [16] "RHC085__H1_B1__brown" "RHC087__H1_B1__brown" "RHC088__H1_B1__brown"
#> [19] "RHC089__H1_B1__brown" "RHC090__H1_B1__brown" "RHC093__H1_B1__brown"
#> [22] "RHC094__H1_B1__brown" "RHC095__H1_B1__brown" "RHC096__H1_B1__brown"
#> [25] "RHC097__H1_B1__brown" "RHC098__H1_B1__brown" "RHC100__H1_B1__brown"
#> [28] "RHC101__H1_B1__brown" "RHC102__H1_B1__brown" "RHC103__H1_B1__brown"
#> [31] "RHC104__H1_B1__brown" "RHC105__H1_B1__brown" "RHC107__H1_B1__brown"
#> [34] "RHC108__H1_B1__brown" "RHC109__H1_B1__brown" "RHC110__H1_B1__brown"
#> [37] "RHC111__H1_B1__brown" "RHC112__H1_B1__brown" "RHC113__H1_B1__brown"
#> [40] "RHC114__H1_B1__brown" "RHC115__H1_B1__brown" "RHC116__H1_B1__brown"
#> [43] "RHC117__H1_B1__brown" "RHC118__H1_B1__brown" "RHC120__H1_B1__brown"
#> [46] "RHC121__H1_B1__brown" "RHC122__H1_B1__brown" "RHC123__H1_B1__brown"
#> [49] "RHC124__H1_B1__brown" "RHC125__H1_B1__brown" "RHC127__H1_B1__brown"
#> [52] "RHC128__H1_B1__brown" "RHC129__H1_B1__brown" "RHC130__H1_B1__brown"
#> [55] "RHC131__H1_B1__brown" "RHC132__H1_B1__brown" "RHC134__H1_B1__brown"
#> [58] "RHC135__H1_B1__brown" "RHC136__H1_B1__brown" "RHC137__H1_B1__brown"
#> [61] "RHC138__H1_B1__brown" "RHC140__H1_B1__brown" "RHC141__H1_B1__brown"
#> [64] "RHC142__H1_B1__brown" "RHC143__H1_B1__brown" "RHC144__H1_B1__brown"
#> [67] "RHC145__H1_B1__brown" "RHC146__H1_B1__brown"
```

To calculate ARI (Adjusted Rand Index):

``` r
library(clues)
# Convert sample names to binary values
cl1 <- ifelse(grepl('GM12878_B1', colnames(cl.b1)), 1, 2)
# Map samples to clusters in result
cl2 <- vapply(colnames(cl.b1), function(y) which(sapply(cl2, function(x) y %in% x)), 2333, USE.NAMES = FALSE)
# Calculate Adjusted Rand Index
adjustedRand(cl1, cl2, 'Rand')
#>      Rand 
#> 0.9138614
```

### On [Five Samples](#five-types-in-same-batch)

``` r
data('cl.no.batch')
handsome.zuo <- CrossSCC(cl.no.batch, ncores = 56, mean.posterior.cutoff = 0.31, verbose = FALSE)
```

``` r
handsome.zuo$Get('sampleNames', filterFun = isLeaf)
#> $`GO:0032496`
#> $`GO:0032496`[[1]]
#>  [1] "RHL1003__H1437__red"   "RHL1010__H1437__red"  
#>  [3] "RHL1019__H1437__red"   "RHL1033__H1437__red"  
#>  [5] "RHL1059__H1437__red"   "RHL1066__H1437__red"  
#>  [7] "RHL1067__H1437__red"   "RHL1075__H1437__red"  
#>  [9] "RHL1077__H1437__red"   "RHL1094__H1437__red"  
#> [11] "RHL1099__H1437__red"   "RHL1119__H1437__red"  
#> [13] "RHL1133__H1437__red"   "RHB1138__IMR90__black"
#> [15] "RHB1139__IMR90__black" "RHB1160__IMR90__black"
#> [17] "RHT021__K562__blue"    "RHT033__K562__blue"   
#> [19] "RHT055__K562__blue"    "RHT079__K562__blue"   
#> [21] "RHT085__K562__blue"    "RHT092__K562__blue"   
#> 
#> $`GO:0032496`[[2]]
#>  [1] "RHL995__H1437__red"     "RHL996__H1437__red"    
#>  [3] "RHL999__H1437__red"     "RHL1001__H1437__red"   
#>  [5] "RHL1005__H1437__red"    "RHL1011__H1437__red"   
#>  [7] "RHL1012__H1437__red"    "RHL1013__H1437__red"   
#>  [9] "RHL1015__H1437__red"    "RHL1022__H1437__red"   
#> [11] "RHL1050__H1437__red"    "RHL1057__H1437__red"   
#> [13] "RHL1062__H1437__red"    "RHL1071__H1437__red"   
#> [15] "RHL1073__H1437__red"    "RHL1088__H1437__red"   
#> [17] "RHL1095__H1437__red"    "RHL1100__H1437__red"   
#> [19] "RHL1102__H1437__red"    "RHL1104__H1437__red"   
#> [21] "RHL1105__H1437__red"    "RHL1106__H1437__red"   
#> [23] "RHL1109__H1437__red"    "RHL1110__H1437__red"   
#> [25] "RHH1143__HCT116__green" "RHH1148__HCT116__green"
#> [27] "RHH1149__HCT116__green" "RHB1133__IMR90__black" 
#> [29] "RHB1134__IMR90__black"  "RHB1136__IMR90__black" 
#> [31] "RHB1143__IMR90__black"  "RHB1152__IMR90__black" 
#> [33] "RHB1175__IMR90__black"  "RHB1176__IMR90__black" 
#> [35] "RHB1178__IMR90__black"  "RHB1180__IMR90__black" 
#> [37] "RHB1183__IMR90__black"  "RHB1187__IMR90__black" 
#> [39] "RHB1202__IMR90__black"  "RHT110__K562__blue"    
#> 
#> 
#> $`GO:0016233`
#> $`GO:0016233`[[1]]
#>  [1] "RHL1049__H1437__red"    "RHL1144__H1437__red"   
#>  [3] "RHH1111__HCT116__green" "RHH1112__HCT116__green"
#>  [5] "RHH1116__HCT116__green" "RHH1117__HCT116__green"
#>  [7] "RHH1118__HCT116__green" "RHH1120__HCT116__green"
#>  [9] "RHH1122__HCT116__green" "RHH1124__HCT116__green"
#> [11] "RHH1126__HCT116__green" "RHH1128__HCT116__green"
#> [13] "RHH1129__HCT116__green" "RHH1130__HCT116__green"
#> [15] "RHH1131__HCT116__green" "RHH1133__HCT116__green"
#> [17] "RHH1135__HCT116__green" "RHH1140__HCT116__green"
#> [19] "RHH1141__HCT116__green" "RHH1144__HCT116__green"
#> [21] "RHH1147__HCT116__green" "RHH1151__HCT116__green"
#> [23] "RHH1153__HCT116__green" "RHH1154__HCT116__green"
#> [25] "RHH1156__HCT116__green" "RHH1157__HCT116__green"
#> [27] "RHH1158__HCT116__green" "RHH1160__HCT116__green"
#> [29] "RHH1164__HCT116__green" "RHH1165__HCT116__green"
#> [31] "RHH1166__HCT116__green" "RHH1169__HCT116__green"
#> [33] "RHB1144__IMR90__black"  "RHB1168__IMR90__black" 
#> [35] "RHB1181__IMR90__black"  "RHB1203__IMR90__black" 
#> 
#> $`GO:0016233`[[2]]
#>  [1] "RHL1042__H1437__red"    "RHH1114__HCT116__green"
#>  [3] "RHH1119__HCT116__green" "RHH1121__HCT116__green"
#>  [5] "RHH1123__HCT116__green" "RHH1125__HCT116__green"
#>  [7] "RHH1127__HCT116__green" "RHH1132__HCT116__green"
#>  [9] "RHH1136__HCT116__green" "RHH1137__HCT116__green"
#> [11] "RHH1138__HCT116__green" "RHH1142__HCT116__green"
#> [13] "RHH1146__HCT116__green" "RHH1150__HCT116__green"
#> [15] "RHH1167__HCT116__green" "RHH1170__HCT116__green"
#> [17] "RHH1171__HCT116__green"
#> 
#> 
#> $`GO:1904431`
#>  [1] "RHL1016__H1437__red"    "RHL1078__H1437__red"   
#>  [3] "RHL1096__H1437__red"    "RHL1103__H1437__red"   
#>  [5] "RHL1135__H1437__red"    "RHL1136__H1437__red"   
#>  [7] "RHH1145__HCT116__green" "RHH1172__HCT116__green"
#>  [9] "RHB1130__IMR90__black"  "RHB1157__IMR90__black" 
#> [11] "RHB1174__IMR90__black"  "RHB1179__IMR90__black" 
#> [13] "RHT018__K562__blue"     "RHT019__K562__blue"    
#> [15] "RHT022__K562__blue"     "RHT023__K562__blue"    
#> [17] "RHT025__K562__blue"     "RHT029__K562__blue"    
#> [19] "RHT030__K562__blue"     "RHT035__K562__blue"    
#> [21] "RHT036__K562__blue"     "RHT038__K562__blue"    
#> [23] "RHT039__K562__blue"     "RHT040__K562__blue"    
#> [25] "RHT041__K562__blue"     "RHT044__K562__blue"    
#> [27] "RHT045__K562__blue"     "RHT046__K562__blue"    
#> [29] "RHT047__K562__blue"     "RHT048__K562__blue"    
#> [31] "RHT049__K562__blue"     "RHT051__K562__blue"    
#> [33] "RHT052__K562__blue"     "RHT053__K562__blue"    
#> [35] "RHT054__K562__blue"     "RHT056__K562__blue"    
#> [37] "RHT058__K562__blue"     "RHT059__K562__blue"    
#> [39] "RHT061__K562__blue"     "RHT062__K562__blue"    
#> [41] "RHT063__K562__blue"     "RHT067__K562__blue"    
#> [43] "RHT069__K562__blue"     "RHT071__K562__blue"    
#> [45] "RHT072__K562__blue"     "RHT074__K562__blue"    
#> [47] "RHT075__K562__blue"     "RHT080__K562__blue"    
#> [49] "RHT081__K562__blue"     "RHT082__K562__blue"    
#> [51] "RHT083__K562__blue"     "RHT086__K562__blue"    
#> [53] "RHT088__K562__blue"     "RHT089__K562__blue"    
#> [55] "RHT090__K562__blue"     "RHT091__K562__blue"    
#> [57] "RHT093__K562__blue"     "RHT094__K562__blue"    
#> [59] "RHT096__K562__blue"     "RHT098__K562__blue"    
#> [61] "RHT099__K562__blue"     "RHT100__K562__blue"    
#> [63] "RHT101__K562__blue"     "RHT102__K562__blue"    
#> [65] "RHT103__K562__blue"     "RHT105__K562__blue"    
#> [67] "RHT106__K562__blue"     "RHT109__K562__blue"    
#> [69] "RHT112__K562__blue"     "RHT113__K562__blue"    
#> 
#> $`GO:0048247`
#>  [1] "RHA015__A549__turquoise" "RHA016__A549__turquoise"
#>  [3] "RHA017__A549__turquoise" "RHA018__A549__turquoise"
#>  [5] "RHA028__A549__turquoise" "RHA029__A549__turquoise"
#>  [7] "RHA030__A549__turquoise" "RHA031__A549__turquoise"
#>  [9] "RHA032__A549__turquoise" "RHA033__A549__turquoise"
#> [11] "RHA034__A549__turquoise" "RHA035__A549__turquoise"
#> [13] "RHA036__A549__turquoise" "RHA037__A549__turquoise"
#> [15] "RHA038__A549__turquoise" "RHA039__A549__turquoise"
#> [17] "RHA040__A549__turquoise" "RHA041__A549__turquoise"
#> [19] "RHA042__A549__turquoise" "RHA043__A549__turquoise"
#> [21] "RHA044__A549__turquoise" "RHA045__A549__turquoise"
#> [23] "RHA046__A549__turquoise" "RHA047__A549__turquoise"
#> [25] "RHA048__A549__turquoise" "RHA049__A549__turquoise"
#> [27] "RHA050__A549__turquoise" "RHA051__A549__turquoise"
#> [29] "RHA052__A549__turquoise" "RHA053__A549__turquoise"
#> [31] "RHA054__A549__turquoise" "RHA055__A549__turquoise"
#> [33] "RHA056__A549__turquoise" "RHA057__A549__turquoise"
#> [35] "RHA058__A549__turquoise" "RHA059__A549__turquoise"
#> [37] "RHA060__A549__turquoise" "RHA061__A549__turquoise"
#> [39] "RHA062__A549__turquoise" "RHA063__A549__turquoise"
#> [41] "RHA064__A549__turquoise" "RHA065__A549__turquoise"
#> [43] "RHA066__A549__turquoise" "RHA067__A549__turquoise"
#> [45] "RHA068__A549__turquoise" "RHA069__A549__turquoise"
#> [47] "RHA070__A549__turquoise" "RHA071__A549__turquoise"
#> [49] "RHA072__A549__turquoise" "RHA073__A549__turquoise"
#> [51] "RHA074__A549__turquoise" "RHA075__A549__turquoise"
#> [53] "RHA076__A549__turquoise" "RHA077__A549__turquoise"
#> [55] "RHA078__A549__turquoise" "RHA079__A549__turquoise"
#> [57] "RHA080__A549__turquoise" "RHA081__A549__turquoise"
#> [59] "RHA082__A549__turquoise" "RHA084__A549__turquoise"
#> [61] "RHA085__A549__turquoise" "RHA086__A549__turquoise"
#> [63] "RHA087__A549__turquoise" "RHA088__A549__turquoise"
#> [65] "RHA089__A549__turquoise" "RHA091__A549__turquoise"
#> [67] "RHA092__A549__turquoise" "RHA093__A549__turquoise"
#> [69] "RHA094__A549__turquoise" "RHA095__A549__turquoise"
#> [71] "RHA096__A549__turquoise" "RHA097__A549__turquoise"
#> [73] "RHA098__A549__turquoise" "RHA099__A549__turquoise"
#> [75] "RHL1117__H1437__red"     "RHT042__K562__blue"     
#> [77] "RHT043__K562__blue"      "RHT050__K562__blue"     
#> [79] "RHT057__K562__blue"      "RHT060__K562__blue"     
#> [81] "RHT064__K562__blue"      "RHT097__K562__blue"     
#> [83] "RHT104__K562__blue"
cl2 <- unname(lapply(rapply(handsome.zuo$Get('sampleNames', filterFun = isLeaf), enquote, how = 'unlist'), eval))
cl2 <- vapply(colnames(cl.no.batch), function(y) which(sapply(cl2, function(x) y %in% x)), 2333, USE.NAMES = FALSE)
library(stringr)
cl1 <- factor(str_match(colnames(cl.no.batch), '__(.+)__')[, 2])
levels(cl1) <- seq_len(5)
cl1 <- as.character(cl1)
adjustedRand(cl1, cl2, 'Rand')
#>      Rand 
#> 0.8681033
```

## Interactive visualization

``` r
plot_CrossSCC(handsome.zuo)
```

![](man/figures/readme2.gif)

## Find best parameter combination by Bayesion Optimization (not finished yet)

``` r
test.CrossSCC <- function(mean.posterior.cutoff, ovl.cutoff, mean.posterior.weight, ovl.weight, lambda.cutoff) {
  handsome.zuo <- suppressMessages(CrossSCC(cl.b1[1:500], ncores = 10,
                                            mean.posterior.cutoff = mean.posterior.cutoff,
                           ovl.cutoff = ovl.cutoff, mean.posterior.weight = mean.posterior.weight,
                           ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff))
  cl2 <- unname(lapply(rapply(handsome.zuo$Get('sampleNames', filterFun = isLeaf),
                              enquote, how = 'unlist'), eval))
  cl2 <- vapply(colnames(cl.b1), function(y) which(vapply(cl2, function(x) y %in% x, logical(1))),
                2333, USE.NAMES = FALSE)
  cl1 <- ifelse(grepl('GM12878_B1', colnames(cl.b1)), 1, 2)
  list(Score = adjustedRand(cl1, cl2, 'Rand'), Pred = 0)
}

opt.res <- BayesianOptimization(test.CrossSCC,
                                bounds = list(mean.posterior.cutoff = c(0.15, 0.2), ovl.cutoff = c(0, 1),
                                              mean.posterior.weight = c(0, 1), ovl.weight = c(0, 1),
                                              lambda.cutoff = c(0, 1)),
                                init_points = 50, n_iter = 20)
```

## TODO list (sorted by priority in descending order)

  - \[x\] Finish the function for ranking features with rated
    components.
  - \[x\] Finish the main function with parameters.
  - \[x\] Design a method for visualizing the final pseudo decision
    tree.
  - \[x\] The result of fitting is somewhat random at the time. Try to
    make it reproducible when needed.
  - \[x\] Hide the plotting result of `boot.comp()`. Only keep it in
    test mode.
