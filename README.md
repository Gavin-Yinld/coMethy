<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/coMethly/blob/master/figures/co-methy.gif"/></div>

# coMethly
Co-methylation analysis to cluster the genomic loci with similar methylation pattern

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
