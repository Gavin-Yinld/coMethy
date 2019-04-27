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
