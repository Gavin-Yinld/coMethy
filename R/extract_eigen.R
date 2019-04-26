extract_eigen <- function(methy_data ,all_label , number_of_eig,plot=F )
{	
	eig_num_per_module <- function(label,total_number=10000)
	{
	freq <- table(label)
	base_line <- floor(total_number/length(freq))
	min_cluster_num <- min(freq)
	pick_number <- 0
	result <- data.frame(module_id=names(freq),number=0)
	rownames(result) <- names(freq)
	picked_num <- sum(result$number)
	while(picked_num < total_number)
	{	#print(picked_num)
		if(base_line>min_cluster_num)
		{
			index <- which(freq<base_line)
			result[names(freq)[index],2] <- freq[index]
			picked_num <- sum(result$number)
			freq <- freq[-index]
			base_line <- floor((total_number-picked_num)/length(freq))
			min_cluster_num <- min(freq)
		} else {
			result[names(freq),2] <- base_line
			picked_num <- sum(result$number)
			if(picked_num<total_number)
			{
				rest <- total_number - picked_num
				#index <- which.max(freq)
				index <- order(freq,decreasing=TRUE)[1:rest]
				result[names(freq)[index],2] <- result[names(freq)[index],2] + 1
				picked_num <- sum(result$number)
			}
		}
	}
	if(picked_num > total_number)
	{
		rest <- picked_num - total_number
		index <- order(freq,decreasing=TRUE)[1:rest]
		result[index,2] <- result[index,2] - 1
		picked_num <- sum(result$number)
	}
	return(result)
	}


	number_per_module <- eig_num_per_module(label=all_label,total_number=number_of_eig)
	eigen_loci=NULL
	eigen_loci.label <- NULL
	for(label in unique(all_label))
	{
		#print(label)
		temp.data <- methy_data[which(all_label==label),]
		p <-  prcomp(t(temp.data), scale=FALSE)
		temp.res <- temp.data[order(abs(p$rotation[,1]),decreasing=T)[1:number_per_module[label,2]],]
		eigen_loci <- rbind(eigen_loci,temp.res)
		eigen_loci.label <- c(eigen_loci.label,rep(label,nrow(temp.res)))
	}
	if(plot==T)
	{ 
	  label.list <- sort(unique(all_label))
	  pdf("eigen_loci.pdf")
	  par(mfrow=c(3,3))
	  for(label in sort(unique(eigen_loci.label)))
	  {
	    temp.data <- eigen_loci[which(eigen_loci.label==label),]
	    plot(1:ncol(temp.data),seq(0,1,length.out=ncol(temp.data)),type='n',xlab='Samples',
	         ylab='Methylation level',main=paste0("Module ",label,',','\n',"n=",nrow(temp.data)))
	    for(j in 1:nrow(temp.data))
	    {
	      lines(x=1:ncol(temp.data),y=temp.data[j,],col='red',lwd=0.1)
	    }
	    #lines(x=1:ncol(temp.data),y=colMeans(temp.data),col='red',lwd=0.1)
	  }
	  dev.off()
	}
	  
	
	
	return(list(methy_prof=eigen_loci,eigen_loci.label=eigen_loci.label))
}
