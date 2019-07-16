#' @title Maximization step of the algorithm
#'
#' @description Use the expectation step to generate parameters of the mixture model
#'
#' @param pca.components Expects a matrix with n rows as the cells and m columns as the principal components.
#' @param clusters Expect the most likely cluster assignments.
#' @param parameter.estimates Expects a list of parameters (mu, covariance, and weights) from initialize_gmm or maximization_step.
#' @param posteriors Data.frame of priors from the expectation step.
#' @param num.cores Number of cores used for parallel calculations
#' @importFrom foreach %dopar%
#' @export

maximization_step <- function(pca.components,clusters,parameter.estimates,posteriors,num.cores) {

	if (num.cores==1) {

	maximization.parameters <- vector("list",length=length(levels(clusters)))

	post.sum <- apply(posteriors,2,sum)

	for (i in 1:length(levels(clusters))) {

		maximization.parameters[[i]]$mu <- mapply(function(x,y) sum(posteriors[,x]*pca.components[,y])/post.sum[x],i,1:ncol(pca.components))
		names(maximization.parameters[[i]]$mu) <- colnames(pca.components)

		trans.pca.mu <- t(apply(pca.components,1,function(x) x-maximization.parameters[[i]]$mu))
		maximization.parameters[[i]]$covariance <- t(posteriors[,i]*trans.pca.mu) %*% (posteriors[,i]*trans.pca.mu)/post.sum[i]

		maximization.parameters[[i]]$weight <- post.sum[i]/nrow(pca.components)

	}

	return(maximization.parameters)

	} else {

	mu <- vector()
	covariance <- vector()
	weight <- vector()
	maximization.parameters <- vector()

	post.sum <- apply(posteriors,2,sum)

	combine.func <- function(x) {unlist(x,recursive=F)}

	maximization.parameters <- foreach::foreach(i=1:length(levels(clusters))) %dopar% {

		mu <- mapply(function(x,y) sum(posteriors[,x]*pca.components[,y])/post.sum[x],i,1:ncol(pca.components))
		names(mu) <- colnames(pca.components)

		trans.pca.mu <- t(apply(pca.components,1,function(x) x-mu))
		covariance <- t(posteriors[,i]*trans.pca.mu) %*% (posteriors[,i]*trans.pca.mu)/post.sum[i]

		weight <- post.sum[i]/nrow(pca.components)

		list("mu"=mu,"covariance"=covariance,"weight"=weight)

	}

	return(maximization.parameters)

	}

	}
