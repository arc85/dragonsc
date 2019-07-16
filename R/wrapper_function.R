#' @title DRAGON clustering algorithm
#'
#' @description Takes input from, PCA components, temperature decay steps, maximum iterations, delta log likelihood for convergence at each step, maximum number of resulting clusters, and number of cores for parallelization. Has options for saving intermediate files and for verbose output.
#' @param pca.components Expects a matrix with n rows as the cells and m columns as the principal components.
#' @param temp.decay.steps Temperature progressed used for the algorithm. Expects a two column data frame, with columns named "step" and "temperature"
#' @param max.iterations Maximum number of iterations of the expectation maximization at each temperature
#' @param delta.log.likelihood Change in the log likelihood between steps of the expectation maximization for convergence
#' @param max.clusters Maximum number of clusters in the final step of the algorithm
#' @param save.intermediate.files Whether to save files at each temperature of the algorithm. Default is FALSE, but can be very informative to view intermediate results.
#' @param num.cores Number of cores to use for the algorithm. Uses doParallel as a backend.
#' @param verbose Whether to display the current step of the algorithm and the log likelihood for each iteration of the expectation maximization algorithm. Default is FALSE.
#' @export

dragon <- function(pca.components,temp.decay.steps,max.iterations,delta.log.likelihood,max.clusters,save.intermediate.files=FALSE,num.cores,verbose=FALSE) {

	if (num.cores>1) {
	cl <- parallel::makeCluster(num.cores)
	doParallel::registerDoParallel(cl,cores=num.cores)
	}

	for (i in 1:nrow(temp.decay.steps)) {

		if (i==1) {

		if (verbose==TRUE) {
			print(paste("Initializing DRAGON clustering at temperature=",signif(temp.decay.steps[i,2],digits=2),sep=""))
		}

		em_results <- initialize_gmm(pca.components,temp.decay.steps[i,2],num.cores)

		if (verbose==TRUE) {
		print(paste("Starting loop ",i," at temperature=",signif(temp.decay.steps[i,2],digits=2),sep=""))
	}

em_results <- expect_max(pca.components,clusters=em_results[["clusters"]],em_results[["maximization"]],max.iterations,delta.log.likelihood,temp.decay.steps[i,2],num.cores,verbose)

if (verbose==TRUE) {
print(paste("Completed loop ",i," at temperature=",signif(temp.decay.steps[i,2],digits=2),sep=""))
}

if (save.intermediate.files==TRUE){
cur.dir <- getwd()
if (any(list.files()=="dragon_intermediate_files")) {
	setwd(paste(cur.dir,"/dragon_intermediate_files/",sep=""))
} else {
dir.create(paste(cur.dir,"/dragon_intermediate_files/",sep=""))
new.dir <- paste(cur.dir,"/dragon_intermediate_files/",sep="")
setwd(new.dir)
}
saveRDS(em_results,file=paste("em_results_temp_",signif(temp.decay.steps[i,2],digits=2),".RDS",sep=""))
}

} else {

if (length(levels(em_results[["clusters"]]))<max.clusters) {

	if (i>2) {
				if (temp.decay.steps[i,2]==temp.use) {
			temp.use <- temp.decay.steps[i+1,2]
		} else { if (temp.decay.steps[i,2]>temp.use) {
			temp.use <- temp.decay.steps[min(which(temp.decay.steps[,2]<temp.use)),2]
		} else {temp.use <- temp.decay.steps[i,2]}
		}
	} else {temp.use <- temp.decay.steps[i,2]}

	split.clusters <- cluster_split(pca.components,parameter.estimates=em_results[["maximization"]],clusters=em_results[["clusters"]],temperature=temp.use,temp.decay.steps,max.clusters,num.cores)

	temp.use <- split.clusters[["temperature"]]

	if (verbose==TRUE) {
	print(paste("Starting loop ",i," at temperature=",signif(split.clusters[["temperature"]],digits=2),sep=""))
}

	em_results <- expect_max(pca.components,split.clusters[["clusters"]],split.clusters[["maximization"]],max.iterations,delta.log.likelihood,temp=split.clusters[["temperature"]],num.cores,verbose)

	if (verbose==TRUE) {
	print(paste("Completed loop ",i," at temperature=", signif(split.clusters[["temperature"]],digits=2),sep=""))
}

if (save.intermediate.files==TRUE) {
saveRDS(em_results,file=paste("em_results_temp_",signif(temp.use,digits=2),".RDS",sep=""))
}

} else {

		if (temp.decay.steps[i,2]==temp.use) {
			temp.use <- temp.decay.steps[i+1,2]
		} else { if (temp.decay.steps[i,2]>temp.use) {
			temp.use <- temp.decay.steps[min(which(temp.decay.steps[,2]<temp.use)),2]
		} else {temp.use <- temp.decay.steps[i,2]}
		}

	if (temp.use>0) {

if (verbose==TRUE) {
		print(paste("Starting loop ",i," at temperature=",signif(temp.use,digits=2),sep=""))
	}

	em_results <- expect_max(pca.components,em_results[["clusters"]],em_results[["maximization"]],max.iterations,delta.log.likelihood,temp=temp.use,num.cores,verbose)

if (verbose==TRUE) {
	print(paste("Completed loop ",i," at temperature=", signif(temp.use,digits=2),sep=""))
}

if (save.intermediate.files==TRUE) {
saveRDS(em_results,file=paste("em_results_temp_",signif(temp.use,digits=2),".RDS",sep=""))
}

} else {

	if (verbose==TRUE) {
	print(paste("Starting loop ",i," at temperature=0",sep=""))
}

	em_results <- expect_max_zero(pca.components,em_results[["clusters"]],em_results[["maximization"]],max.iterations,delta.log.likelihood,num.cores,verbose)

if (verbose==TRUE) {
	print(paste("Completed loop ",i," at temperature=0",sep=""))
}

if (save.intermediate.files==TRUE) {
saveRDS(em_results,file="em_results_temp_0.RDS")
}

	return(em_results)

				} #close if for temp>0

			} #close else for max.clusters

		} #close i > 1 else step

	} #close nrow(temp.decay.steps) loop

	parallel::stopCluster(cl)

} #close function
