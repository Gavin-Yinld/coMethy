<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/co-methy.gif"/></div>

# coMethly
Co-methylation analysis to cluster the genomic loci with similar methylation pattern
# Installation
csmFinder needs the following tools to be installed and available in the `PATH` environment:
1.  [R](https://www.r-project.org/)(>=3.4.4)
2.  [python2](https://www.python.org/downloads/) (>=2.7.10), to process the bismark extractor results
3.  [bedtools2](https://github.com/arq5x/bedtools2) (>=2.28.0), to merge the overlapped pCSM segments into pCSM loci

In R console,
```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges","WGCNA"))
library("devtools")
devtools::install_github("Gavin-Yinld/coMethy")
```



```R
library("coMethy")
file=paste(system.file(package="coMethy"),"extdata/co_methy.test.data.txt",sep='/')
meth_data <- read.table(file,sep='\t',header=T,row.names=1)
kmeans_cluster <- co_methylation_step1(meth_data)
module <- co_methylation_step2(data=meth_data,
                               kmeans_result=kmeans_cluster,
                               softPower_list=c(16,20,16),plot=T)
eigen_loci <- extract_eigen(module$profile,module$module_id,100,plot=T)
```
