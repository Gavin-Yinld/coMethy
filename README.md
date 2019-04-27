<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/co-methy.gif"/></div>

# coMethly
# Introduction
coMethy is an R package for grouping the the genomic loci with similar methylation pattern. According to the methylation profiles of genomic loci in target methylomes, co-methylation analysis was performed in combination of k-means clustering and WGCNA analysis. For each co-methylation module, PCA analysis was performed to select a subset of loci as the eigen loci representing methylation trend.

# Current Features
* Co-methylation analysis to cluster pCSM loci with similar methylation pattern into co-methylated modules
* PCA analysis to extract eigen-pCSM loci representing methylation trend of ecah co-methylation module

# Installation
In R console,
```R
source("http://bioconductor.org/biocLite.R")
biocLite("WGCNA")
library("devtools")
devtools::install_github("Gavin-Yinld/coMethly")
```
# How to Use

## Step 1. K-means clustering to divide the genomic with distinct methylation level
`csmFinder` takes the methylation profile of input genomic loci in target methylomes. A numeric matrix is needed as input with each row denotes a genomic loci and each column denotes a sample.
```R
library("coMethy")
# get the demo dataset
file=paste(system.file(package="coMethy"),"extdata/co_methy.test.data.txt",sep='/')
meth_data <- read.table(file,sep='\t',header=T,row.names=1)

# a typical input data looks like this:
head(meth_data)
                         glia_7wk neu_7wk glia_12mo_F neu_12mo_F exc_neu pv_neu
chr7_82243987_82244107       0.39    0.77        0.21       0.73    0.81   0.89
chr2_166158817_166158941     0.64    0.90        0.59       0.92    0.90   0.46
chr5_30670219_30670424       0.69    0.84        0.64       0.88    0.85   0.22
chr7_66421727_66421857       0.66    0.76        0.68       0.80    0.74   0.80
chr4_131771412_131771552     0.47    0.91        0.33       0.93    0.94   0.93
chr3_135529338_135529388     0.92    1.00        1.00       0.97    0.76   0.73
                         vip_neu astrocyte oligodendrocyte
chr7_82243987_82244107      0.65      0.72            0.73
chr2_166158817_166158941    0.62      0.81            0.59
chr5_30670219_30670424      0.90      0.93            0.83
chr7_66421727_66421857      0.84      0.64            0.71
chr4_131771412_131771552    0.10      0.36            0.30
chr3_135529338_135529388    0.90      0.91            0.92

```
```R
kmeans_cluster <- co_methylation_step1(meth_data)
module <- co_methylation_step2(data=meth_data,
                               kmeans_result=kmeans_cluster,
                               softPower_list=c(16,20,16),plot=T)
eigen_loci <- extract_eigen(module$profile,module$module_id,100,plot=T)
```
<div align=center><img width="700" height="525" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/power.png"/></div>
<div align=center><img width="700" height="525" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/wgcna.png"/></div>
<div align=center><img width="700" height="525" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/eigen_loci.png"/></div>
