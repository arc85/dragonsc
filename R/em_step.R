#' @title EM step of the algorithm
#'
#' @description Wrapper to go between the expectation and maximization steps until the maximum iterations or the convergence of log likelihood
#' @param pca.components Expects a matrix with n rows as the cells and m columns as the principal components.
#' @param clusters Expect the most likely cluster assignments
#' @param parameter.estimates Expects a list of parameters (mu, covariance, and weights) from initialize_gmm or maximization_step
#' @param max.iterations The maximum number of iterations of the EM
#' @param delta.log.li Change in the log likelihood to satisfy convergence
#' @export

expect_max <- function(pca.components,tsne.coordinates,clusters,parameter.estimates,max.iterations,delta.log.li,temp,num.cores) {
	
for (i in 1:max.iterations) {
	
if (i==1) {
	
	log.li <- vector()
	cur.log.li <- vector()
	
	expectation <-expectation_step(pca.components,clusters,parameter.estimates,temp,num.cores)
	
	post.sums <- rowSums(expectation[["likelihood.probability"]])
	
	if (any(post.sums==0)) {
		
	zero.index <- which(post.sums==0)
	
	pca.components <- pca.components[-zero.index,]
	tsne.coordinates <- tsne.coordinates[-zero.index,]
	clusters <- clusters[-zero.index]
	expectation[[1]] <- expectation[[1]][-zero.index,]
	expectation[[2]] <- expectation[[2]][-zero.index,]
	
}

	clusters <- as.factor(apply(unlist(expectation[["posterior.probability"]]),1,function(x) which.max(x)))
	
	parameter.estimates <- maximization_step(pca.components,clusters,parameter.estimates,posteriors=expectation[["posterior.probability"]],num.cores)
	
	weights <- sapply(parameter.estimates,"[[",3)
	
	if (any(weights==0)) {
		zero.index <- which(weights==0)
		parameter.estimates <- parameter.estimates[-zero.index]
	}
	
	cur.log.li <- sum(log(apply(expectation[["likelihood.probability"]],1,sum),base=exp(1)))
	
	if (cur.log.li=="-Inf") {
		zero.index <- which(log(apply(expectation[["likelihood.probability"]]*weights,2,sum),base=exp(1))=="-Inf")
		pca.components <- pca.components[-zero.index,]
		tsne.coordinates <- tsne.coordinates[-zero.index,]
		clusters <- clusters[-zero.index]
		expectation[[1]] <- expectation[[1]][-zero.index,]
		expectation[[2]] <- expectation[[2]][-zero.index,]
		cur.log.li <- sum(log(apply(expectation[["likelihood.probability"]],1,sum),base=exp(1)))
	}
	
	log.li[i] <- cur.log.li
	
	print(paste("Current log likelihood=",cur.log.li,sep=""))
	
} else {

	expectation <-expectation_step(pca.components,clusters,parameter.estimates,temp,num.cores)
	
	post.sums <- rowSums(expectation[["likelihood.probability"]])
	
	if (any(post.sums==0)) {
		
	zero.index <- which(post.sums==0)
	
	pca.components <- pca.components[-zero.index,]
	tsne.coordinates <- tsne.coordinates[-zero.index,]
	clusters <- clusters[-zero.index]
	expectation[[1]] <- expectation[[1]][-zero.index,]
	expectation[[2]] <- expectation[[2]][-zero.index,]
	
}

	clusters <- as.factor(apply(unlist(expectation[["posterior.probability"]]),1,function(x) which.max(x)))
	
	parameter.estimates <- maximization_step(pca.components,clusters,parameter.estimates,posteriors=expectation[["posterior.probability"]],num.cores)
	
	weights <- sapply(parameter.estimates,"[[",3)
	
	if (any(weights==0)) {
		zero.index <- which(weights==0)
		parameter.estimates <- parameter.estimates[-zero.index]
	}
	
	cur.log.li <- sum(log(apply(expectation[["likelihood.probability"]],1,sum),base=exp(1)))
	
		if (cur.log.li=="-Inf") {
		zero.index <- which(log(apply(expectation[["likelihood.probability"]]*weights,2,sum),base=exp(1))=="-Inf")
		pca.components <- pca.components[-zero.index,]
		tsne.coordinates <- tsne.coordinates[-zero.index,]
		clusters <- clusters[-zero.index]
		expectation[[1]] <- expectation[[1]][-zero.index,]
		expectation[[2]] <- expectation[[2]][-zero.index,]
		cur.log.li <- sum(log(apply(expectation[["likelihood.probability"]],1,sum),base=exp(1)))
	}
	
	log.li[i] <- cur.log.li
	
	print(paste("Current log likelihood=",cur.log.li,sep=""))
	
		if (abs(cur.log.li-log.li[i-1])<delta.log.li) {
	
return(list("pca.components"=pca.components,"tsne.coordinates"=tsne.coordinates,"expectation"=expectation,"maximization"=parameter.estimates,"clusters"=clusters,"MLE"=log.li))
	
	}
		
}

}	
	
return(list("pca.components"=pca.components,"tsne.coordinates"=tsne.coordinates,"expectation"=expectation,"maximization"=parameter.estimates,"clusters"=clusters,"MLE"=log.li))

}