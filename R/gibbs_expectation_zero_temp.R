#' @title Deterministic annealing expectation
#'
#' @description Uses a Euclidean distance metric to determine the distance of each cell from cluster centers, and uses a Gibbs sampler to evaluate the likelihood probability and posterior probability
#' @param pca.components Expects a matrix with n rows as the cells and m columns as the principal components.
#' @param clusters Expect the most likely cluster assignments
#' @param parameter.estimates Expects a list of parameters (mu, covariance, and weights) from initialize_gmm or maximization_step
#' @importFrom foreach %dopar%
#' @export

expectation_step_zero <- function(pca.components,clusters,parameter.estimates,num.cores) {
	
	if (num.cores==1) {
	
	dist.measure <- vector("list",length=length(levels(clusters)))
	llh.prob <- vector("list",length=length(levels(clusters)))
	post.prob <- vector("list",length=length(levels(clusters)))
	
	for (i in 1:length(levels(clusters))) {
		
dist.measure[[i]] <- mapply(function(x,y) mahalanobis(pca.components[x,],parameter.estimates[[y]]$mu,parameter.estimates[[y]]$cov),1:nrow(pca.components),i)

llh.prob[[i]] <- exp(-dist.measure[[i]])
	
	}
	
	llh.prob.sum <- Reduce("+",llh.prob)
	
	for (i in 1:length(levels(clusters))) {
				
	post.prob[[i]] <- llh.prob[[i]]/llh.prob.sum
		
	}
	
	names(llh.prob) <- paste("Cluster_",seq(1,length(levels(clusters)),1),sep="")
	llh.prob.mat <- as.matrix(dplyr::bind_rows(llh.prob))
	rownames(llh.prob.mat) <- rownames(pca.components)
	
	names(post.prob) <- paste("Cluster_",seq(1,length(levels(clusters)),1),sep="")
	post.prob.mat <- as.matrix(dplyr::bind_rows(post.prob))
	rownames(post.prob.mat) <- rownames(pca.components)
	
	results.list <- vector("list",length=2)
	names(results.list) <- c("likelihood.probability","posterior.probability")

	results.list[[1]] <- llh.prob.mat
	results.list[[2]] <- post.prob.mat
	
	return(results.list)
	
	} else {
		
	dist.measure <- vector()
	llh.prob <- vector()
	post.prob <- vector("list",length=length(levels(clusters)))
	
llh.prob <- foreach::foreach(i=1:length(levels(clusters))) %dopar% {
		
dist.measure <- mapply(function(x,y) mahalanobis(pca.components[x,],parameter.estimates[[y]]$mu,parameter.estimates[[y]]$cov),1:nrow(pca.components),i)

exp(-dist.measure)
	
	}
	
	llh.prob.sum <- Reduce("+",llh.prob)
	
	for (i in 1:length(levels(clusters))) {
				
	post.prob[[i]] <- llh.prob[[i]]/llh.prob.sum
		
	}
	
	names(llh.prob) <- paste("Cluster_",seq(1,length(levels(clusters)),1),sep="")
	llh.prob.mat <- as.matrix(dplyr::bind_rows(llh.prob))
	rownames(llh.prob.mat) <- rownames(pca.components)
	
	names(post.prob) <- paste("Cluster_",seq(1,length(levels(clusters)),1),sep="")
	post.prob.mat <- as.matrix(dplyr::bind_rows(post.prob))
	rownames(post.prob.mat) <- rownames(pca.components)
	
	results.list <- vector("list",length=2)
	names(results.list) <- c("likelihood.probability","posterior.probability")

	results.list[[1]] <- llh.prob.mat
	results.list[[2]] <- post.prob.mat
	
	return(results.list)

	}

} 