#' @title Split clusters based on temperatures 
#' 
#' @description Split clusters into new clusters based on temperature for the deterministic annealing algorithm
#' @param pca.components Expects a matrix with n rows as the cells and m columns as the principal components.
#' @param parameter.estimates Expects a list of parameters (mu, covariance, and weights) from initialize_gmm or maximization_step
#' @param clusters Expect the most likely cluster assignments
#' @param temperature Current temperature to use to split clusters
#' @param decay.step Data.frame with the first column as the step number and the second column as the temperature for that step
#' @param max.clusters Maximum number of clusters for the dataset
#' @export

cluster_split <- function(pca.components,parameter.estimates,clusters,temperature,decay.step,max.clusters,num.cores) {
	
	initial.clusters <- as.numeric(as.character(levels(clusters)))
	
	critical.temps <- vector()
	
	for (i in 1:length(levels(clusters))) {
		
	suppressWarnings(critical.temps[i] <- irlba::prcomp_irlba(pca.components[clusters==i,],n=1)$sdev^2*2
	)
	
	}
	
	if (any(temperature>critical.temps)) {
		temperature <- decay.step[which(min(critical.temps) > decay.step[,2])[1],"temperature"]
		
	}
	
	split.clusters <- data.frame("initial"=NULL,"split"=NULL)
	
	if((sum(critical.temps>temperature)*2)>max.clusters) {
		num.to.split <- max.clusters-length(initial.clusters)
		clusters.to.split <- order(critical.temps,decreasing=T)[1:num.to.split]
		
	for (i in clusters.to.split) {
		
		clusters <- factor(clusters,levels=c(seq(1,length(levels(clusters))+1,1)))
		
	if (nrow(split.clusters)==0) {
		
				split.clusters <- data.frame("initial"=levels(clusters)[i],"split"=levels(clusters)[length(levels(clusters))])
				
	} else {
		
		split.clusters <- rbind(split.clusters,data.frame("initial"=levels(clusters)[i],"split"=levels(clusters)[length(levels(clusters))]))
		
	}
	
	half.cluster.length <- round(length(which(clusters==i))/2)
			clusters[which(clusters==i)[(half.cluster.length+1):(2*half.cluster.length)]] <- levels(clusters)[length(levels(clusters))]
			
			}
			
	} else {
	
	for (i in 1:length(critical.temps)) {
		
		if (critical.temps[i]>temperature) {
			
			clusters <- factor(clusters,levels=c(seq(1,length(levels(clusters))+1,1)))
			
			if (nrow(split.clusters)==0) {
				
				split.clusters <- data.frame("initial"=levels(clusters)[i],"split"=levels(clusters)[length(levels(clusters))])
				
			} else {
				
				split.clusters <- rbind(split.clusters,data.frame("initial"=levels(clusters)[i],"split"=levels(clusters)[length(levels(clusters))]))
				
			}
			
			half.cluster.length <- round(length(which(clusters==i))/2)
			clusters[which(clusters==i)[(half.cluster.length+1):(2*half.cluster.length)]] <- levels(clusters)[length(levels(clusters))]
			
		}
		
	}
	
	}
	
	new.clusters <- setdiff(as.numeric(as.character(levels(clusters))),initial.clusters)
	
	if (length(new.clusters)<1) {
		return(list("clusters"=clusters,"maximization"=parameter.estimates))
		
	} else {
	
	new.cluster.params <- vector("list",length(levels(clusters)))
	
	for (i in new.clusters) {
		subset <- pca.components[clusters==i,]
		new.cluster.params[[i]]$mu <- apply(subset,2,mean)
		new.cluster.params[[i]]$covariance <- cov(subset)
		new.cluster.params[[i]]$weight <- nrow(subset)/nrow(pca.components)
	
	}
	
	parameter.estimates <- c(parameter.estimates[initial.clusters],new.cluster.params[new.clusters])
	
	for (i in 1:nrow(split.clusters)) {
		
	parameter.estimates[[as.numeric(as.character(split.clusters$initial[i]))]]$weight <- parameter.estimates[[as.numeric(as.character(split.clusters$split[i]))]]$weight
		
	}
	
	}
	
	expect_results <- expectation_step(pca.components,clusters,parameter.estimates,temperature,num.cores)
	
return(list("clusters"=clusters,"maximization"=parameter.estimates,"conditional.prob"=expect_results[["conditional.prob"]],"posterior.probability"=expect_results[["posterior.probability"]],temperature=temperature))

	
}