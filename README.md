<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/co-methy.gif"/></div>

# coMethly
Co-methylation analysis to cluster the genomic loci with similar methylation pattern
# Installation

In R console,
```R
source("http://bioconductor.org/biocLite.R")
biocLite("WGCNA")
library("devtools")
devtools::install_github("Gavin-Yinld/coMethly")
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
<div align=center><img width="800" height="600" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/power.png"/></div>
<div align=center><img width="600" height="800" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/wgcna.png"/></div>
<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/eigen_loci.pdf"/></div>
