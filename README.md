
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
(handsome.zuo <- CrossSCC(cl.7, ncores = 16, mean.posterior.cutoff = 0.1036, ovl.cutoff = 0.1864, mean.posterior.weight = 0.2557,   ovl.weight = 0.6198, lambda.cutoff = 0.7750, verbose = FALSE))
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
#>                             levelName
#> 1 LOC109864269                       
#> 2  ¦--GO:0006119                     
#> 3  ¦   ¦--GO:0006119\nsub-component 1
#> 4  ¦   °--GO:0006119\nsub-component 2
#> 5  °--GO:1900116                     
#> 6      ¦--GO:1900116\nsub-component 1
#> 7      °--GO:1900116\nsub-component 2
library(data.tree)
handsome.zuo$Get('sampleNames', filterFun = isLeaf, simplify = FALSE)
#> $`GO:0006119\nsub-component 1`
#>   [1] "RHA015__A549__turquoise"     "RHA016__A549__turquoise"    
#>   [3] "RHA017__A549__turquoise"     "RHA018__A549__turquoise"    
#>   [5] "RHA028__A549__turquoise"     "RHA029__A549__turquoise"    
#>   [7] "RHA030__A549__turquoise"     "RHA031__A549__turquoise"    
#>   [9] "RHA032__A549__turquoise"     "RHA033__A549__turquoise"    
#>  [11] "RHA034__A549__turquoise"     "RHA035__A549__turquoise"    
#>  [13] "RHA036__A549__turquoise"     "RHA037__A549__turquoise"    
#>  [15] "RHA038__A549__turquoise"     "RHA039__A549__turquoise"    
#>  [17] "RHA041__A549__turquoise"     "RHA042__A549__turquoise"    
#>  [19] "RHA043__A549__turquoise"     "RHA044__A549__turquoise"    
#>  [21] "RHA045__A549__turquoise"     "RHA046__A549__turquoise"    
#>  [23] "RHA047__A549__turquoise"     "RHA048__A549__turquoise"    
#>  [25] "RHA049__A549__turquoise"     "RHA050__A549__turquoise"    
#>  [27] "RHA052__A549__turquoise"     "RHA053__A549__turquoise"    
#>  [29] "RHA054__A549__turquoise"     "RHA055__A549__turquoise"    
#>  [31] "RHA056__A549__turquoise"     "RHA057__A549__turquoise"    
#>  [33] "RHA058__A549__turquoise"     "RHA059__A549__turquoise"    
#>  [35] "RHA060__A549__turquoise"     "RHA061__A549__turquoise"    
#>  [37] "RHA062__A549__turquoise"     "RHA064__A549__turquoise"    
#>  [39] "RHA065__A549__turquoise"     "RHA066__A549__turquoise"    
#>  [41] "RHA067__A549__turquoise"     "RHA068__A549__turquoise"    
#>  [43] "RHA069__A549__turquoise"     "RHA070__A549__turquoise"    
#>  [45] "RHA071__A549__turquoise"     "RHA072__A549__turquoise"    
#>  [47] "RHA074__A549__turquoise"     "RHA075__A549__turquoise"    
#>  [49] "RHA077__A549__turquoise"     "RHA078__A549__turquoise"    
#>  [51] "RHA079__A549__turquoise"     "RHA080__A549__turquoise"    
#>  [53] "RHA081__A549__turquoise"     "RHA082__A549__turquoise"    
#>  [55] "RHA084__A549__turquoise"     "RHA085__A549__turquoise"    
#>  [57] "RHA086__A549__turquoise"     "RHA087__A549__turquoise"    
#>  [59] "RHA088__A549__turquoise"     "RHA089__A549__turquoise"    
#>  [61] "RHA091__A549__turquoise"     "RHA092__A549__turquoise"    
#>  [63] "RHA093__A549__turquoise"     "RHA094__A549__turquoise"    
#>  [65] "RHA095__A549__turquoise"     "RHA096__A549__turquoise"    
#>  [67] "RHA097__A549__turquoise"     "RHA098__A549__turquoise"    
#>  [69] "RHB1116__GM12878_B1__orange" "RHB1118__GM12878_B1__orange"
#>  [71] "RHB1122__GM12878_B1__orange" "RHB1129__GM12878_B1__orange"
#>  [73] "RHB1141__GM12878_B1__orange" "RHB1146__GM12878_B1__orange"
#>  [75] "RHB1147__GM12878_B1__orange" "RHB1153__GM12878_B1__orange"
#>  [77] "RHB1162__GM12878_B1__orange" "RHB1167__GM12878_B1__orange"
#>  [79] "RHB1169__GM12878_B1__orange" "RHB1172__GM12878_B1__orange"
#>  [81] "RHB1185__GM12878_B1__orange" "RHB1188__GM12878_B1__orange"
#>  [83] "RHB1197__GM12878_B1__orange" "RHB1199__GM12878_B1__orange"
#>  [85] "RHB1201__GM12878_B1__orange" "RHB1206__GM12878_B1__orange"
#>  [87] "RHB1207__GM12878_B1__orange" "RHB1208__GM12878_B1__orange"
#>  [89] "RHL995__H1437__red"          "RHL996__H1437__red"         
#>  [91] "RHL999__H1437__red"          "RHL1001__H1437__red"        
#>  [93] "RHL1003__H1437__red"         "RHL1005__H1437__red"        
#>  [95] "RHL1010__H1437__red"         "RHL1011__H1437__red"        
#>  [97] "RHL1012__H1437__red"         "RHL1015__H1437__red"        
#>  [99] "RHL1016__H1437__red"         "RHL1019__H1437__red"        
#> [101] "RHL1022__H1437__red"         "RHL1033__H1437__red"        
#> [103] "RHL1042__H1437__red"         "RHL1050__H1437__red"        
#> [105] "RHL1059__H1437__red"         "RHL1066__H1437__red"        
#> [107] "RHL1067__H1437__red"         "RHL1071__H1437__red"        
#> [109] "RHL1073__H1437__red"         "RHL1075__H1437__red"        
#> [111] "RHL1077__H1437__red"         "RHL1078__H1437__red"        
#> [113] "RHL1094__H1437__red"         "RHL1095__H1437__red"        
#> [115] "RHL1096__H1437__red"         "RHL1099__H1437__red"        
#> [117] "RHL1100__H1437__red"         "RHL1104__H1437__red"        
#> [119] "RHL1105__H1437__red"         "RHL1109__H1437__red"        
#> [121] "RHL1110__H1437__red"         "RHL1117__H1437__red"        
#> [123] "RHL1133__H1437__red"         "RHL1135__H1437__red"        
#> [125] "RHL1136__H1437__red"         "RHL1144__H1437__red"        
#> [127] "RHH1111__HCT116__green"      "RHH1112__HCT116__green"     
#> [129] "RHH1114__HCT116__green"      "RHH1116__HCT116__green"     
#> [131] "RHH1117__HCT116__green"      "RHH1118__HCT116__green"     
#> [133] "RHH1119__HCT116__green"      "RHH1120__HCT116__green"     
#> [135] "RHH1121__HCT116__green"      "RHH1122__HCT116__green"     
#> [137] "RHH1123__HCT116__green"      "RHH1124__HCT116__green"     
#> [139] "RHH1125__HCT116__green"      "RHH1126__HCT116__green"     
#> [141] "RHH1127__HCT116__green"      "RHH1128__HCT116__green"     
#> [143] "RHH1129__HCT116__green"      "RHH1130__HCT116__green"     
#> [145] "RHH1131__HCT116__green"      "RHH1132__HCT116__green"     
#> [147] "RHH1133__HCT116__green"      "RHH1135__HCT116__green"     
#> [149] "RHH1136__HCT116__green"      "RHH1137__HCT116__green"     
#> [151] "RHH1138__HCT116__green"      "RHH1140__HCT116__green"     
#> [153] "RHH1141__HCT116__green"      "RHH1142__HCT116__green"     
#> [155] "RHH1143__HCT116__green"      "RHH1144__HCT116__green"     
#> [157] "RHH1145__HCT116__green"      "RHH1146__HCT116__green"     
#> [159] "RHH1147__HCT116__green"      "RHH1148__HCT116__green"     
#> [161] "RHH1149__HCT116__green"      "RHH1150__HCT116__green"     
#> [163] "RHH1151__HCT116__green"      "RHH1153__HCT116__green"     
#> [165] "RHH1154__HCT116__green"      "RHH1156__HCT116__green"     
#> [167] "RHH1157__HCT116__green"      "RHH1158__HCT116__green"     
#> [169] "RHH1160__HCT116__green"      "RHH1164__HCT116__green"     
#> [171] "RHH1165__HCT116__green"      "RHH1166__HCT116__green"     
#> [173] "RHH1167__HCT116__green"      "RHH1169__HCT116__green"     
#> [175] "RHH1170__HCT116__green"      "RHH1171__HCT116__green"     
#> [177] "RHH1172__HCT116__green"      "RHB1125__GM12878_B1__orange"
#> [179] "RHB1130__IMR90__black"       "RHB1133__IMR90__black"      
#> [181] "RHB1134__IMR90__black"       "RHB1136__IMR90__black"      
#> [183] "RHB1138__IMR90__black"       "RHB1139__IMR90__black"      
#> [185] "RHB1143__IMR90__black"       "RHB1157__IMR90__black"      
#> [187] "RHB1168__IMR90__black"       "RHB1174__IMR90__black"      
#> [189] "RHB1176__IMR90__black"       "RHB1179__IMR90__black"      
#> [191] "RHB1180__IMR90__black"       "RHB1181__IMR90__black"      
#> [193] "RHB1183__IMR90__black"       "RHB1187__IMR90__black"      
#> [195] "RHB1202__IMR90__black"       "RHT021__K562__blue"         
#> [197] "RHT033__K562__blue"          "RHT035__K562__blue"         
#> [199] "RHT054__K562__blue"          "RHT055__K562__blue"         
#> [201] "RHT064__K562__blue"          "RHT079__K562__blue"         
#> [203] "RHT085__K562__blue"          "RHT092__K562__blue"         
#> [205] "RHT101__K562__blue"         
#> 
#> $`GO:0006119\nsub-component 2`
#>  [1] "RHA051__A549__turquoise"     "RHA073__A549__turquoise"    
#>  [3] "RHA076__A549__turquoise"     "RHA099__A549__turquoise"    
#>  [5] "RHB1140__GM12878_B1__orange" "RHB1158__GM12878_B1__orange"
#>  [7] "RHB1165__GM12878_B1__orange" "RHB1171__GM12878_B1__orange"
#>  [9] "RHB1173__GM12878_B1__orange" "RHB1189__GM12878_B1__orange"
#> [11] "RHB1193__GM12878_B1__orange" "RHB1200__GM12878_B1__orange"
#> [13] "RHB1205__GM12878_B1__orange" "RHL1013__H1437__red"        
#> [15] "RHL1049__H1437__red"         "RHL1088__H1437__red"        
#> [17] "RHL1102__H1437__red"         "RHL1103__H1437__red"        
#> [19] "RHL1106__H1437__red"         "RHL1119__H1437__red"        
#> [21] "RHB1160__IMR90__black"       "RHB1175__IMR90__black"      
#> [23] "RHB1178__IMR90__black"       "RHB1203__IMR90__black"      
#> [25] "RHT018__K562__blue"          "RHT019__K562__blue"         
#> [27] "RHT022__K562__blue"          "RHT023__K562__blue"         
#> [29] "RHT025__K562__blue"          "RHT029__K562__blue"         
#> [31] "RHT030__K562__blue"          "RHT036__K562__blue"         
#> [33] "RHT038__K562__blue"          "RHT039__K562__blue"         
#> [35] "RHT040__K562__blue"          "RHT041__K562__blue"         
#> [37] "RHT042__K562__blue"          "RHT043__K562__blue"         
#> [39] "RHT044__K562__blue"          "RHT045__K562__blue"         
#> [41] "RHT046__K562__blue"          "RHT047__K562__blue"         
#> [43] "RHT048__K562__blue"          "RHT049__K562__blue"         
#> [45] "RHT050__K562__blue"          "RHT051__K562__blue"         
#> [47] "RHT052__K562__blue"          "RHT053__K562__blue"         
#> [49] "RHT056__K562__blue"          "RHT057__K562__blue"         
#> [51] "RHT058__K562__blue"          "RHT059__K562__blue"         
#> [53] "RHT060__K562__blue"          "RHT061__K562__blue"         
#> [55] "RHT062__K562__blue"          "RHT063__K562__blue"         
#> [57] "RHT067__K562__blue"          "RHT069__K562__blue"         
#> [59] "RHT071__K562__blue"          "RHT072__K562__blue"         
#> [61] "RHT074__K562__blue"          "RHT075__K562__blue"         
#> [63] "RHT080__K562__blue"          "RHT081__K562__blue"         
#> [65] "RHT082__K562__blue"          "RHT083__K562__blue"         
#> [67] "RHT086__K562__blue"          "RHT088__K562__blue"         
#> [69] "RHT089__K562__blue"          "RHT090__K562__blue"         
#> [71] "RHT091__K562__blue"          "RHT093__K562__blue"         
#> [73] "RHT094__K562__blue"          "RHT097__K562__blue"         
#> [75] "RHT098__K562__blue"          "RHT099__K562__blue"         
#> [77] "RHT100__K562__blue"          "RHT102__K562__blue"         
#> [79] "RHT104__K562__blue"          "RHT105__K562__blue"         
#> [81] "RHT106__K562__blue"          "RHT109__K562__blue"         
#> [83] "RHT110__K562__blue"          "RHT112__K562__blue"         
#> [85] "RHT113__K562__blue"         
#> 
#> $`GO:1900116\nsub-component 1`
#>  [1] "RHB1144__IMR90__black" "RHB1152__IMR90__black"
#>  [3] "RHC069__H1_B1__brown"  "RHC070__H1_B1__brown" 
#>  [5] "RHC071__H1_B1__brown"  "RHC072__H1_B1__brown" 
#>  [7] "RHC073__H1_B1__brown"  "RHC074__H1_B1__brown" 
#>  [9] "RHC075__H1_B1__brown"  "RHC077__H1_B1__brown" 
#> [11] "RHC078__H1_B1__brown"  "RHC079__H1_B1__brown" 
#> [13] "RHC080__H1_B1__brown"  "RHC081__H1_B1__brown" 
#> [15] "RHC082__H1_B1__brown"  "RHC083__H1_B1__brown" 
#> [17] "RHC085__H1_B1__brown"  "RHC087__H1_B1__brown" 
#> [19] "RHC088__H1_B1__brown"  "RHC089__H1_B1__brown" 
#> [21] "RHC090__H1_B1__brown"  "RHC092__H1_B1__brown" 
#> [23] "RHC093__H1_B1__brown"  "RHC094__H1_B1__brown" 
#> [25] "RHC095__H1_B1__brown"  "RHC096__H1_B1__brown" 
#> [27] "RHC097__H1_B1__brown"  "RHC098__H1_B1__brown" 
#> [29] "RHC100__H1_B1__brown"  "RHC101__H1_B1__brown" 
#> [31] "RHC102__H1_B1__brown"  "RHC103__H1_B1__brown" 
#> [33] "RHC104__H1_B1__brown"  "RHC105__H1_B1__brown" 
#> [35] "RHC107__H1_B1__brown"  "RHC108__H1_B1__brown" 
#> [37] "RHC109__H1_B1__brown"  "RHC110__H1_B1__brown" 
#> [39] "RHC111__H1_B1__brown"  "RHC112__H1_B1__brown" 
#> [41] "RHC113__H1_B1__brown"  "RHC114__H1_B1__brown" 
#> [43] "RHC115__H1_B1__brown"  "RHC116__H1_B1__brown" 
#> [45] "RHC117__H1_B1__brown"  "RHC118__H1_B1__brown" 
#> [47] "RHC120__H1_B1__brown"  "RHC121__H1_B1__brown" 
#> [49] "RHC122__H1_B1__brown"  "RHC123__H1_B1__brown" 
#> [51] "RHC124__H1_B1__brown"  "RHC125__H1_B1__brown" 
#> [53] "RHC127__H1_B1__brown"  "RHC128__H1_B1__brown" 
#> [55] "RHC129__H1_B1__brown"  "RHC130__H1_B1__brown" 
#> [57] "RHC131__H1_B1__brown"  "RHC132__H1_B1__brown" 
#> [59] "RHC135__H1_B1__brown"  "RHC136__H1_B1__brown" 
#> [61] "RHC137__H1_B1__brown"  "RHC138__H1_B1__brown" 
#> [63] "RHC140__H1_B1__brown"  "RHC141__H1_B1__brown" 
#> [65] "RHC142__H1_B1__brown"  "RHC143__H1_B1__brown" 
#> [67] "RHC144__H1_B1__brown"  "RHC145__H1_B1__brown" 
#> [69] "RHC146__H1_B1__brown" 
#> 
#> $`GO:1900116\nsub-component 2`
#>  [1] "RHA040__A549__turquoise"     "RHA063__A549__turquoise"    
#>  [3] "RHB1124__GM12878_B1__orange" "RHB1204__GM12878_B1__orange"
#>  [5] "RHL1057__H1437__red"         "RHL1062__H1437__red"        
#>  [7] "RHC084__H1_B1__brown"        "RHC134__H1_B1__brown"       
#>  [9] "RHT096__K562__blue"          "RHT103__K562__blue"
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
#> 0.7035908
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
