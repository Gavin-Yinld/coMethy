<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/co-methy.gif"/></div>

# coMethy
# Introduction
coMethy is an R package for grouping the the genomic loci with similar methylation pattern. According to the methylation profiles of genomic loci in target methylomes, co-methylation analysis was performed in combination of k-means clustering and WGCNA analysis. For each co-methylation module, PCA analysis was performed to select a subset of loci as the eigen loci representing methylation trend.

# Installation
In R console,
```R
source("http://bioconductor.org/biocLite.R")
biocLite("WGCNA")
library("devtools")
devtools::install_github("Gavin-Yinld/coMethy")
```
# How to Use

## Step 1. K-means clustering to divide the genomic with distinct methylation level
`coMethy` takes the methylation profile of input genomic loci in target methylomes. A numeric matrix is needed as input with each row denotes a genomic loci ,each column denotes a sample and each value denotes the methylation level.
```R
library("coMethy")

# get the demo dataset
file=paste(system.file(package="coMethy"),"extdata/co_methy.test.data.txt",sep='/')
meth_data <- read.table(file,sep='\t',header=T,row.names=1)

# a typical input data looks like this:
head(meth_data)
                         sample1 sample2 sample3 sample4 sample5 sample6
chr7_82243987_82244107      0.39    0.77    0.21    0.73    0.81    0.89
chr2_166158817_166158941    0.64    0.90    0.59    0.92    0.90    0.46
chr5_30670219_30670424      0.69    0.84    0.64    0.88    0.85    0.22
chr7_66421727_66421857      0.66    0.76    0.68    0.80    0.74    0.80
chr4_131771412_131771552    0.47    0.91    0.33    0.93    0.94    0.93
chr3_135529338_135529388    0.92    1.00    1.00    0.97    0.76    0.73
                         sample7 sample8 sample9
chr7_82243987_82244107      0.65    0.72    0.73
chr2_166158817_166158941    0.62    0.81    0.59
chr5_30670219_30670424      0.90    0.93    0.83
chr7_66421727_66421857      0.84    0.64    0.71
chr4_131771412_131771552    0.10    0.36    0.30
chr3_135529338_135529388    0.90    0.91    0.92
```
Firstly, K-means clustering analysis is called to divide pCSM loci into hypo/mid/hyper-methylated groups. In addinion, `pickSoftThreshold` in `WGCNA` package is called to show the topological properties of the network in each group. 

```R
kmeans_cluster <- co_methylation_step1(meth_data)
# A file named "parameter.pdf" will be generated to show the the topological properties.
```
<div align=center><img width="700" height="525" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/power.png"/></div>

## Step 2. WGCNA analysis to detect the co-methylation module
Network construction was performed using the `blockwiseModules` function in the `WGCNA` package, which allows the network construction for the entire data set. For each of kmeans-group, a pair-wise correlation matrix is computed, and an adjacency matrix is calculated by raising the correlation matrix to a power. The proper power need to be chosen using the scale-free topology criterion in step 1. For example, the power of 16,20,16 are chosen for the networks built from each kmeans-group.
```R
module <- co_methylation_step2(data=meth_data,
                               kmeans_result=kmeans_cluster,
                               softPower_list=c(16,20,16),plot=T)
# By setting "plot=T", a file named "wgcna.module.pdf" will be generated to show the methylation level of the loci in each group.
```
<div align=center><img width="700" height="525" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/wgcna.png"/></div>

## Step 3. PCA analysis to extract eigen-loci from each co-methylation module
PCA analysis is adopted to pick a set of pCSM loci with the largest loading in PC1 as eigen-loci for the corresponding module to represent methylation trend.
```R
eigen_loci <- extract_eigen(methy_data=module$profile,
                            all_label=module$module_id,
                            number_of_eig=100,plot=T)
#By setting "plot=T", a file named "eigen_loci.pdf" will be generated to show the methylation level of the eigen-loci in each group.
```
<div align=center><img width="700" height="525" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/eigen_loci.png"/></div>
