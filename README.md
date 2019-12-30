
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CrossSCC

The goal of CrossSCC is to classify **S**ingle-**C**ell data(as
expression matrix of scRNA-seq data) **Cross**ing batch into
**C**lusters using Gaussian Mixture Model, and to find out the mapping
relationship in clusters.

![](man/figures/readme.gif)

## Installation

You can install the latest version of CrossSCC with:

``` r
install.packages("remotes")
remotes::install_github("bioinformatist/CrossSCC")
```

Besides, you also need a database
[org.HsSimple.eg.db](https://github.com/bioinformatist/org.HsSimple.eg.db)
provided by us: First download the [latest release
tarball](https://github.com/bioinformatist/org.HsSimple.eg.db/releases),
then run:

``` r
install.packages("./org.HsSimple.eg.db", repos=NULL)
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

### Five types as well as batch 1 of two types

``` r
cl <- fread('GSE81861_Cell_Line_FPKM.csv')
names(cl)[1] <- 'Gene'
cl[, Ensembl := str_match(Gene, ".+_(.+)\\.\\d+$")[, 2]]
symbols <- mapIds(org.Hs.eg.db, keys = cl$Ensembl, keytype = "ENSEMBL", column="ENTREZID", multiVals = 'first')
cl[, Entrez := symbols[Ensembl]]
cl[, c('Gene', 'Ensembl'):= NULL]
cl <- cbind(cl[, .(Entrez)], cl[, .SD, .SDcols = !(names(cl) %like% "_B2_")])
cl[, var := rowVars(.SD), .SDcols = -c('Entrez')]
cl <- cl[, max.var := max(var), by = 'Entrez'][var != 0 & max.var == var & !is.na(Entrez), ]
cl[, grep("var", colnames(cl)) := NULL]
cl <- as.data.frame(cl) %>% remove_rownames %>% column_to_rownames(var = "Entrez")
cl <- cl[, !(names(cl) %in% 'Entrez')]
cl.7 <- new("ExpressionSet", exprs=as.matrix(cl), annotation = 'org.Hs.eg.db')
# setwd() back to package root directory
use_data(cl.7, compress = 'xz')
```

## How to use & performance

``` r
library(CrossSCC)
#> Loading required package: R.utils
#> Loading required package: R.oo
#> Loading required package: R.methodsS3
#> R.methodsS3 v1.7.1 (2016-02-15) successfully loaded. See ?R.methodsS3 for help.
#> R.oo v1.23.0 successfully loaded. See ?R.oo for help.
#> 
#> Attaching package: 'R.oo'
#> The following object is masked from 'package:R.methodsS3':
#> 
#>     throw
#> The following objects are masked from 'package:methods':
#> 
#>     getClasses, getMethods
#> The following objects are masked from 'package:base':
#> 
#>     attach, detach, load, save
#> R.utils v2.9.2 successfully loaded. See ?R.utils for help.
#> 
#> Attaching package: 'R.utils'
#> The following object is masked from 'package:utils':
#> 
#>     timestamp
#> The following objects are masked from 'package:base':
#> 
#>     cat, commandArgs, getOption, inherits, isOpen, parse, warnings
#> Loading required package: snow
#> 
data("cl.7")
(handsome.zuo <- CrossSCC(cl.7, ncores = 16, mean.posterior.cutoff = 0.3475, ovl.cutoff = 0.2682, mean.posterior.weight = 0.0000, ovl.weight = 0.7094, lambda.cutoff = 0.8587, verbose = FALSE))
#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~

#> Warning in searchCommandline(parallel, cpus = cpus, type = type,
#> socketHosts = socketHosts, : Unknown option on commandline:
#> rmarkdown::render('/home/ysun/github.com/bioinformatist/GMMClassifier/
#> README.Rmd',~+~~+~encoding~+~
#>                                      levelName
#> 1  GO:0006826                                 
#> 2   ¦--LOC109864269                           
#> 3   ¦   ¦--GO:0016241                         
#> 4   ¦   ¦   ¦--GO:0032981                     
#> 5   ¦   ¦   ¦   ¦--GO:0032981\nsub-component 1
#> 6   ¦   ¦   ¦   °--GO:0032981\nsub-component 2
#> 7   ¦   ¦   °--RPL36A                         
#> 8   ¦   ¦       ¦--RPL36A\nsub-component 1    
#> 9   ¦   ¦       °--RPL36A\nsub-component 2    
#> 10  ¦   °--GO:0034080                         
#> 11  ¦       ¦--GO:0034080\nsub-component 1    
#> 12  ¦       °--GO:0034080\nsub-component 2    
#> 13  °--GO:0034341                             
#> 14      ¦--GO:0022900                         
#> 15      ¦   ¦--SNORD13                        
#> 16      ¦   ¦   ¦--SNORD13\nsub-component 1   
#> 17      ¦   ¦   °--SNORD13\nsub-component 2   
#> 18      ¦   °--GO:0048678                     
#> 19      ¦       ¦--GO:0048678\nsub-component 1
#> 20      ¦       °--GO:0048678\nsub-component 2
#> 21      °--GO:0000183                         
#> 22          ¦--GO:0060333                     
#> 23          ¦   ¦--GO:0060333\nsub-component 1
#> 24          ¦   °--GO:0060333\nsub-component 2
#> 25          °--TMSB10                         
#> 26              ¦--TMSB10\nsub-component 1    
#> 27              °--TMSB10\nsub-component 2
library(data.tree)
handsome.zuo$Get('sampleNames', filterFun = isLeaf, simplify = FALSE)
#> $`GO:0032981\nsub-component 1`
#>  [1] "RHB1118__GM12878_B1__orange" "RHB1169__GM12878_B1__orange"
#>  [3] "RHB1185__GM12878_B1__orange" "RHB1201__GM12878_B1__orange"
#>  [5] "RHL995__H1437__red"          "RHL1001__H1437__red"        
#>  [7] "RHL1003__H1437__red"         "RHL1010__H1437__red"        
#>  [9] "RHL1016__H1437__red"         "RHL1019__H1437__red"        
#> [11] "RHL1022__H1437__red"         "RHL1050__H1437__red"        
#> [13] "RHL1067__H1437__red"         "RHL1073__H1437__red"        
#> [15] "RHL1078__H1437__red"         "RHL1088__H1437__red"        
#> [17] "RHL1105__H1437__red"         "RHL1117__H1437__red"        
#> [19] "RHL1119__H1437__red"         "RHL1133__H1437__red"        
#> [21] "RHL1135__H1437__red"         "RHL1144__H1437__red"        
#> 
#> $`GO:0032981\nsub-component 2`
#>  [1] "RHL1033__H1437__red"    "RHL1059__H1437__red"   
#>  [3] "RHL1066__H1437__red"    "RHL1075__H1437__red"   
#>  [5] "RHL1077__H1437__red"    "RHL1094__H1437__red"   
#>  [7] "RHL1095__H1437__red"    "RHL1099__H1437__red"   
#>  [9] "RHL1104__H1437__red"    "RHH1143__HCT116__green"
#> [11] "RHB1138__IMR90__black"  "RHT021__K562__blue"    
#> [13] "RHT033__K562__blue"     "RHT035__K562__blue"    
#> [15] "RHT054__K562__blue"     "RHT055__K562__blue"    
#> [17] "RHT064__K562__blue"     "RHT079__K562__blue"    
#> [19] "RHT085__K562__blue"     "RHT092__K562__blue"    
#> [21] "RHT101__K562__blue"    
#> 
#> $`RPL36A\nsub-component 1`
#>  [1] "RHH1111__HCT116__green" "RHH1112__HCT116__green"
#>  [3] "RHH1114__HCT116__green" "RHH1116__HCT116__green"
#>  [5] "RHH1118__HCT116__green" "RHH1119__HCT116__green"
#>  [7] "RHH1120__HCT116__green" "RHH1121__HCT116__green"
#>  [9] "RHH1122__HCT116__green" "RHH1123__HCT116__green"
#> [11] "RHH1124__HCT116__green" "RHH1125__HCT116__green"
#> [13] "RHH1126__HCT116__green" "RHH1127__HCT116__green"
#> [15] "RHH1128__HCT116__green" "RHH1129__HCT116__green"
#> [17] "RHH1130__HCT116__green" "RHH1131__HCT116__green"
#> [19] "RHH1132__HCT116__green" "RHH1133__HCT116__green"
#> [21] "RHH1135__HCT116__green" "RHH1136__HCT116__green"
#> [23] "RHH1137__HCT116__green" "RHH1138__HCT116__green"
#> [25] "RHH1141__HCT116__green" "RHH1142__HCT116__green"
#> [27] "RHH1144__HCT116__green" "RHH1145__HCT116__green"
#> [29] "RHH1146__HCT116__green" "RHH1148__HCT116__green"
#> [31] "RHH1149__HCT116__green" "RHH1150__HCT116__green"
#> [33] "RHH1151__HCT116__green" "RHH1153__HCT116__green"
#> [35] "RHH1154__HCT116__green" "RHH1156__HCT116__green"
#> [37] "RHH1157__HCT116__green" "RHH1158__HCT116__green"
#> [39] "RHH1160__HCT116__green" "RHH1164__HCT116__green"
#> [41] "RHH1165__HCT116__green" "RHH1166__HCT116__green"
#> [43] "RHH1167__HCT116__green" "RHH1169__HCT116__green"
#> [45] "RHH1170__HCT116__green" "RHH1172__HCT116__green"
#> 
#> $`RPL36A\nsub-component 2`
#>  [1] "RHL996__H1437__red"     "RHL999__H1437__red"    
#>  [3] "RHL1005__H1437__red"    "RHL1011__H1437__red"   
#>  [5] "RHL1012__H1437__red"    "RHL1013__H1437__red"   
#>  [7] "RHL1015__H1437__red"    "RHL1042__H1437__red"   
#>  [9] "RHL1049__H1437__red"    "RHL1057__H1437__red"   
#> [11] "RHL1062__H1437__red"    "RHL1071__H1437__red"   
#> [13] "RHL1096__H1437__red"    "RHL1100__H1437__red"   
#> [15] "RHL1102__H1437__red"    "RHL1103__H1437__red"   
#> [17] "RHL1106__H1437__red"    "RHL1109__H1437__red"   
#> [19] "RHL1110__H1437__red"    "RHL1136__H1437__red"   
#> [21] "RHH1147__HCT116__green" "RHH1171__HCT116__green"
#> [23] "RHB1176__IMR90__black"  "RHB1181__IMR90__black" 
#> [25] "RHT047__K562__blue"     "RHT075__K562__blue"    
#> [27] "RHT110__K562__blue"    
#> 
#> $`GO:0034080\nsub-component 1`
#>  [1] "RHC071__H1_B1__brown" "RHC073__H1_B1__brown" "RHC075__H1_B1__brown"
#>  [4] "RHC077__H1_B1__brown" "RHC084__H1_B1__brown" "RHC087__H1_B1__brown"
#>  [7] "RHC089__H1_B1__brown" "RHC090__H1_B1__brown" "RHC093__H1_B1__brown"
#> [10] "RHC095__H1_B1__brown" "RHC097__H1_B1__brown" "RHC098__H1_B1__brown"
#> [13] "RHC100__H1_B1__brown" "RHC102__H1_B1__brown" "RHC103__H1_B1__brown"
#> [16] "RHC107__H1_B1__brown" "RHC108__H1_B1__brown" "RHC110__H1_B1__brown"
#> [19] "RHC112__H1_B1__brown" "RHC113__H1_B1__brown" "RHC114__H1_B1__brown"
#> [22] "RHC117__H1_B1__brown" "RHC121__H1_B1__brown" "RHC122__H1_B1__brown"
#> [25] "RHC125__H1_B1__brown" "RHC127__H1_B1__brown" "RHC128__H1_B1__brown"
#> [28] "RHC129__H1_B1__brown" "RHC130__H1_B1__brown" "RHC132__H1_B1__brown"
#> [31] "RHC134__H1_B1__brown" "RHC135__H1_B1__brown" "RHC137__H1_B1__brown"
#> [34] "RHC141__H1_B1__brown" "RHC142__H1_B1__brown" "RHC143__H1_B1__brown"
#> [37] "RHC145__H1_B1__brown"
#> 
#> $`GO:0034080\nsub-component 2`
#>  [1] "RHC069__H1_B1__brown" "RHC070__H1_B1__brown" "RHC072__H1_B1__brown"
#>  [4] "RHC074__H1_B1__brown" "RHC078__H1_B1__brown" "RHC079__H1_B1__brown"
#>  [7] "RHC080__H1_B1__brown" "RHC082__H1_B1__brown" "RHC083__H1_B1__brown"
#> [10] "RHC085__H1_B1__brown" "RHC088__H1_B1__brown" "RHC092__H1_B1__brown"
#> [13] "RHC094__H1_B1__brown" "RHC096__H1_B1__brown" "RHC101__H1_B1__brown"
#> [16] "RHC104__H1_B1__brown" "RHC105__H1_B1__brown" "RHC109__H1_B1__brown"
#> [19] "RHC111__H1_B1__brown" "RHC115__H1_B1__brown" "RHC116__H1_B1__brown"
#> [22] "RHC118__H1_B1__brown" "RHC123__H1_B1__brown" "RHC131__H1_B1__brown"
#> [25] "RHC136__H1_B1__brown" "RHC138__H1_B1__brown" "RHC144__H1_B1__brown"
#> 
#> $`SNORD13\nsub-component 1`
#>  [1] "RHA015__A549__turquoise"     "RHA018__A549__turquoise"    
#>  [3] "RHA029__A549__turquoise"     "RHA030__A549__turquoise"    
#>  [5] "RHA033__A549__turquoise"     "RHA035__A549__turquoise"    
#>  [7] "RHA037__A549__turquoise"     "RHA041__A549__turquoise"    
#>  [9] "RHA045__A549__turquoise"     "RHA057__A549__turquoise"    
#> [11] "RHA059__A549__turquoise"     "RHA062__A549__turquoise"    
#> [13] "RHA069__A549__turquoise"     "RHA070__A549__turquoise"    
#> [15] "RHA074__A549__turquoise"     "RHA075__A549__turquoise"    
#> [17] "RHA084__A549__turquoise"     "RHA092__A549__turquoise"    
#> [19] "RHA093__A549__turquoise"     "RHA096__A549__turquoise"    
#> [21] "RHB1116__GM12878_B1__orange" "RHB1146__GM12878_B1__orange"
#> [23] "RHB1158__GM12878_B1__orange" "RHB1172__GM12878_B1__orange"
#> [25] "RHB1205__GM12878_B1__orange" "RHB1207__GM12878_B1__orange"
#> [27] "RHB1130__IMR90__black"       "RHC081__H1_B1__brown"       
#> 
#> $`SNORD13\nsub-component 2`
#>  [1] "RHA016__A549__turquoise"     "RHA017__A549__turquoise"    
#>  [3] "RHA028__A549__turquoise"     "RHA032__A549__turquoise"    
#>  [5] "RHA039__A549__turquoise"     "RHA040__A549__turquoise"    
#>  [7] "RHA042__A549__turquoise"     "RHA043__A549__turquoise"    
#>  [9] "RHA046__A549__turquoise"     "RHA047__A549__turquoise"    
#> [11] "RHA048__A549__turquoise"     "RHA050__A549__turquoise"    
#> [13] "RHA053__A549__turquoise"     "RHA054__A549__turquoise"    
#> [15] "RHA055__A549__turquoise"     "RHA056__A549__turquoise"    
#> [17] "RHA058__A549__turquoise"     "RHA061__A549__turquoise"    
#> [19] "RHA064__A549__turquoise"     "RHA065__A549__turquoise"    
#> [21] "RHA066__A549__turquoise"     "RHA067__A549__turquoise"    
#> [23] "RHA068__A549__turquoise"     "RHA072__A549__turquoise"    
#> [25] "RHA073__A549__turquoise"     "RHA076__A549__turquoise"    
#> [27] "RHA078__A549__turquoise"     "RHA079__A549__turquoise"    
#> [29] "RHA080__A549__turquoise"     "RHA081__A549__turquoise"    
#> [31] "RHA082__A549__turquoise"     "RHA086__A549__turquoise"    
#> [33] "RHA087__A549__turquoise"     "RHA089__A549__turquoise"    
#> [35] "RHA091__A549__turquoise"     "RHA094__A549__turquoise"    
#> [37] "RHA095__A549__turquoise"     "RHA097__A549__turquoise"    
#> [39] "RHA099__A549__turquoise"     "RHB1206__GM12878_B1__orange"
#> 
#> $`GO:0048678\nsub-component 1`
#>  [1] "RHA034__A549__turquoise"     "RHA036__A549__turquoise"    
#>  [3] "RHA038__A549__turquoise"     "RHA044__A549__turquoise"    
#>  [5] "RHA049__A549__turquoise"     "RHA051__A549__turquoise"    
#>  [7] "RHA052__A549__turquoise"     "RHA060__A549__turquoise"    
#>  [9] "RHA063__A549__turquoise"     "RHA071__A549__turquoise"    
#> [11] "RHA077__A549__turquoise"     "RHA085__A549__turquoise"    
#> [13] "RHA098__A549__turquoise"     "RHB1140__GM12878_B1__orange"
#> [15] "RHB1147__GM12878_B1__orange" "RHB1165__GM12878_B1__orange"
#> [17] "RHB1200__GM12878_B1__orange" "RHH1117__HCT116__green"     
#> [19] "RHH1140__HCT116__green"     
#> 
#> $`GO:0048678\nsub-component 2`
#>  [1] "RHB1133__IMR90__black" "RHB1134__IMR90__black"
#>  [3] "RHB1136__IMR90__black" "RHB1139__IMR90__black"
#>  [5] "RHB1143__IMR90__black" "RHB1144__IMR90__black"
#>  [7] "RHB1152__IMR90__black" "RHB1157__IMR90__black"
#>  [9] "RHB1160__IMR90__black" "RHB1168__IMR90__black"
#> [11] "RHB1174__IMR90__black" "RHB1175__IMR90__black"
#> [13] "RHB1178__IMR90__black" "RHB1180__IMR90__black"
#> [15] "RHB1183__IMR90__black" "RHB1187__IMR90__black"
#> [17] "RHB1202__IMR90__black" "RHB1203__IMR90__black"
#> 
#> $`GO:0060333\nsub-component 1`
#>  [1] "RHA031__A549__turquoise" "RHA088__A549__turquoise"
#>  [3] "RHC124__H1_B1__brown"    "RHC140__H1_B1__brown"   
#>  [5] "RHT018__K562__blue"      "RHT019__K562__blue"     
#>  [7] "RHT022__K562__blue"      "RHT038__K562__blue"     
#>  [9] "RHT044__K562__blue"      "RHT050__K562__blue"     
#> [11] "RHT052__K562__blue"      "RHT053__K562__blue"     
#> [13] "RHT061__K562__blue"      "RHT062__K562__blue"     
#> [15] "RHT067__K562__blue"      "RHT069__K562__blue"     
#> [17] "RHT071__K562__blue"      "RHT072__K562__blue"     
#> [19] "RHT083__K562__blue"      "RHT088__K562__blue"     
#> [21] "RHT090__K562__blue"      "RHT091__K562__blue"     
#> [23] "RHT096__K562__blue"      "RHT098__K562__blue"     
#> [25] "RHT102__K562__blue"      "RHT109__K562__blue"     
#> 
#> $`GO:0060333\nsub-component 2`
#>  [1] "RHB1122__GM12878_B1__orange" "RHB1124__GM12878_B1__orange"
#>  [3] "RHB1129__GM12878_B1__orange" "RHB1141__GM12878_B1__orange"
#>  [5] "RHB1153__GM12878_B1__orange" "RHB1167__GM12878_B1__orange"
#>  [7] "RHB1171__GM12878_B1__orange" "RHB1173__GM12878_B1__orange"
#>  [9] "RHB1188__GM12878_B1__orange" "RHB1189__GM12878_B1__orange"
#> [11] "RHB1193__GM12878_B1__orange" "RHB1204__GM12878_B1__orange"
#> [13] "RHB1208__GM12878_B1__orange" "RHB1125__GM12878_B1__orange"
#> [15] "RHB1179__IMR90__black"      
#> 
#> $`TMSB10\nsub-component 1`
#>  [1] "RHT025__K562__blue" "RHT029__K562__blue" "RHT030__K562__blue"
#>  [4] "RHT036__K562__blue" "RHT039__K562__blue" "RHT040__K562__blue"
#>  [7] "RHT041__K562__blue" "RHT042__K562__blue" "RHT043__K562__blue"
#> [10] "RHT045__K562__blue" "RHT046__K562__blue" "RHT048__K562__blue"
#> [13] "RHT051__K562__blue" "RHT056__K562__blue" "RHT057__K562__blue"
#> [16] "RHT058__K562__blue" "RHT059__K562__blue" "RHT060__K562__blue"
#> [19] "RHT063__K562__blue" "RHT080__K562__blue" "RHT086__K562__blue"
#> [22] "RHT093__K562__blue" "RHT094__K562__blue" "RHT097__K562__blue"
#> [25] "RHT099__K562__blue" "RHT100__K562__blue" "RHT104__K562__blue"
#> [28] "RHT112__K562__blue" "RHT113__K562__blue"
#> 
#> $`TMSB10\nsub-component 2`
#>  [1] "RHB1162__GM12878_B1__orange" "RHB1197__GM12878_B1__orange"
#>  [3] "RHB1199__GM12878_B1__orange" "RHC120__H1_B1__brown"       
#>  [5] "RHC146__H1_B1__brown"        "RHT023__K562__blue"         
#>  [7] "RHT049__K562__blue"          "RHT074__K562__blue"         
#>  [9] "RHT081__K562__blue"          "RHT082__K562__blue"         
#> [11] "RHT089__K562__blue"          "RHT103__K562__blue"         
#> [13] "RHT105__K562__blue"          "RHT106__K562__blue"
library(stringr)
cl1 <- factor(str_match(colnames(cl.7), '__(.+?)_')[, 2])
levels(cl1) <- seq_len(7)
cl1 <- as.character(cl1)
cl2 <- unname(lapply(rapply(handsome.zuo$Get('sampleNames', filterFun = isLeaf), enquote, how = 'unlist'), eval))
cl2 <- vapply(colnames(cl.7), function(y) which(vapply(cl2, function(x) y %in% x,
                                                       logical(1))), 2333, USE.NAMES = FALSE)
library(clues)
adjustedRand(cl1, cl2, 'Rand')
#>      Rand 
#> 0.8946772
```

## Interactive visualization

``` r
plot_CrossSCC(handsome.zuo)
```

![](man/figures/readme2.gif)

## Find best parameter combination by Bayesion Optimization

``` r
library(rBayesianOptimization)
library(CrossSCC)
library(data.tree)
library(stringr)
library(clues)
data("cl.7")
test.CrossSCC <- function(mean.posterior.cutoff, ovl.cutoff, mean.posterior.weight, ovl.weight, lambda.cutoff) {
  tryCatch(
    {handsome.zuo <- CrossSCC(cl.7, ncores = 16, verbose = FALSE,
                                            mean.posterior.cutoff = mean.posterior.cutoff,
                           ovl.cutoff = ovl.cutoff, mean.posterior.weight = mean.posterior.weight,
                           ovl.weight = ovl.weight, lambda.cutoff = lambda.cutoff)
    cl1 <- factor(str_match(colnames(cl.7), '__(.+?)_')[, 2])
    levels(cl1) <- seq_len(7)
    cl1 <- as.character(cl1)
    cl2 <- unname(lapply(rapply(handsome.zuo$Get('sampleNames', filterFun = isLeaf), enquote, how = 'unlist'), eval))
    cl2 <- vapply(colnames(cl.7), function(y) which(vapply(cl2, function(x) y %in% x, 
                                                           logical(1))), 2333, USE.NAMES = FALSE)
    list(Score = adjustedRand(cl1, cl2, 'Rand'), Pred = 0)}, error=function(e) list(Score = 0, Pred = 0)
  )
}

opt.res <- BayesianOptimization(test.CrossSCC,
                                bounds = list(mean.posterior.cutoff = c(0, 0.5), ovl.cutoff = c(0, 0.5),
                                              mean.posterior.weight = c(0, 1), ovl.weight = c(0, 1),
                                              lambda.cutoff = c(0.5, 1)),
                                init_points = 20, n_iter = 100)
```

## Acknowledgement

  - Thanks to [ScreenToGif](https://github.com/NickeManarin/ScreenToGif)
    for producing gif images for this repo.
  - Thanks to [Dr. Qi Zhao](http://seqworld.com) at SYSUCC for
    suggestions on tree data structure and user experience.
