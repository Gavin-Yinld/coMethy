co_methylation_step1 <- function(csm_ml_matrix){
options(stringsAsFactors=F)
set.seed(123)
data <- as.matrix(csm_ml_matrix)
km.res <- kmeans(data, 3, nstart = 25)
par(mfrow = c(3,2),mar=c(3,5,2,2));
for(i in 1:3)
{
  pdf("parameter.pdf")
  require(WGCNA)
	wgcna.data <- t(data[which(km.res$cluster==i),])
	powers = c(c(1:10), seq(from = 12, to=30, by=2))
	###Call the network topology analysis function
	sft = pickSoftThreshold(wgcna.data, powerVector = powers, verbose = 5,networkType="signed")
	####Plot the results:

	cex1 = 0.9;
	#######Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology \n Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
	############this line corresponds to using an R^2 cut-off of h
	abline(h=0.80,col="red")
	#######Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	#save(sft,file=paste0("kmeans_",i,"sft.Rdata"))
 dev.off()

}
return(km.res)
}
	###########################################################################














