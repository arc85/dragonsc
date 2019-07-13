#' @title Initialization of the annealing model
#'
#' @description This function uses PCA components as input, and splits the data into two clusters at a given temperature. It also provides initial parameters for the expectation and maximization.
#' @param pca.components Expects a matrix with n rows as the cells and m columns as the principal components.
#' @param temperature Initial temperature for annealing.
#' @export 

initialize.gmm <- function(pca.components,temperature,num.cores) {
		
	clusters <- as.vector(rep(NA,nrow(pca.components)))
	clusters[1:round(nrow(pca.components)/2)] <- 1
	clusters[(round(nrow(pca.components)/2)+1):nrow(pca.components)] <- 2
	clusters <- as.factor(clusters)
	
	cluster.summaries <- vector("list",length=length(levels(clusters)))

	for (i in 1:length(levels(clusters))) {
		subset <- pca.components[clusters==i,]
		cluster.summaries[[i]]$mu <- apply(subset,2,mean)
		cluster.summaries[[i]]$covariance <- cov(subset)
		cluster.summaries[[i]]$weight <- nrow(subset)/nrow(pca.components)
	
	}	

em_results <- expectation_step(pca.components,clusters,cluster.summaries,temperature,num.cores)
	return(list("clusters"=clusters,"maximization"=cluster.summaries,"likelihood.probability"=em_results[["likelihood.probability"]],"posterior.probability"=em_results[["posterior.probability"]]))
	
}