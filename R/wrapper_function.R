#' @title Master wrapper for daGMM
#'
#' @description Takes input from the TSNE cluster coordinates, PCA components, decay step, maximum iterations, delta log likelihood for convergence at each step, and the maximum number of resulting clusters
#' @param tsne.coordinates Expects a data.frame with two columns (tSNE_1 and tSNE_2, respectively) and rows of cells
#' @param pca.components Expects a matrix with n rows as the cells and m columns as the principal components.
#' @param decay.step Expects a data.frame of two columns, with the first as an index of temperature and the second as the temperatures for the determinisitic annealing
#' @param max.iterations Maximum number of iterations of the EM for each step
#' @param delta.log.likelihood Change in the log likelihood between steps of the EM for convergence
#' @param max.clusters Maximum number of clusters in the final step of the algorithm
#' @export

daGMM <- function(tsne.coordinates,pca.components,decay.step,max.iterations,delta.log.likelihood,max.clusters,tsne.plot.x.lims,tsne.plot.y.lims,num.cores) {
	
	if (num.cores>1) {
	cl <- parallel::makeCluster(num.cores)
	doParallel::registerDoParallel(cl,cores=num.cores)
	}
	
	initial.pca.components <- pca.components
	initial.tsne.coordinates <- tsne.coordinates
	
	for (i in 1:nrow(decay.step)) {
		
		if (i==1) {
			
		em_results <- initialize.gmm(pca.components,decay.step[i,2],num.cores)
		
em_results <- expect_max(initial.pca.components,initial.tsne.coordinates,clusters=em_results[["clusters"]],em_results[["maximization"]],max.iterations,delta.log.likelihood,decay.step[i,2],num.cores)
	
print(paste("Completed loop ",i," at temperature=",signif(decay.step[i,2],digits=2),sep=""))

overall.data <- data.frame(em_results[["tsne.coordinates"]],cluster=em_results[["clusters"]])

px <- ggplot2::ggplot(data=overall.data,aes_string(colnames(tsne.coordinates)[1],colnames(tsne.coordinates)[2])) +
stat_density_2d(geom="polygon",aes(alpha=..level..,fill=cluster),bins=5,size=2) + scale_x_continuous(limits=tsne.plot.x.lims,expand=c(0,0)) + scale_y_continuous(limits=tsne.plot.y.lims,expand=c(0,0)) + labs(alpha="Level",fill="Cluster") + theme(panel.background=element_blank(),axis.line=element_line(colour="black"),legend.position="bottom") + guides(alpha=guide_legend(nrow=2),fill=guide_legend(nrow=2)) + coord_equal()

pdf(paste("density_plot_temp_",decay.step[i,2],".pdf",sep=""))
print(px)
dev.off()

saveRDS(em_results,file=paste("em_results_temp_",decay.step[i,2],".RDS",sep=""))

} else {

if (length(levels(em_results[["clusters"]]))<max.clusters) {
	
	if (i>2) {
				if (decay.step[i,2]==temp.use) {
			temp.use <- decay.step[i+1,2]
		} else { if (decay.step[i,2]>temp.use) {
			temp.use <- decay.step[min(which(decay.step[,2]<temp.use)),2]
			} else {temp.use <- decay.step[i,2]}
		}
	} else {temp.use <- decay.step[i,2]}
	
	split.clusters <- cluster_split(initial.pca.components,parameter.estimates=em_results[["maximization"]],clusters=em_results[["clusters"]],temperature=temp.use,decay.step,max.clusters,num.cores)
	
	temp.use <- split.clusters[["temperature"]]
	
	print(paste("Starting loop ",i," at temperature=",signif(split.clusters[["temperature"]],digits=2),sep=""))
	
	em_results <- expect_max(initial.pca.components,initial.tsne.coordinates,split.clusters[["clusters"]],split.clusters[["maximization"]],max.iterations,delta.log.likelihood,temp=split.clusters[["temperature"]],num.cores)
	
	print(paste("Completed loop ",i," at temperature=", signif(split.clusters[["temperature"]],digits=2),sep=""))

overall.data <- data.frame(em_results[["tsne.coordinates"]],cluster=em_results[["clusters"]])

px <- ggplot2::ggplot(data=overall.data,aes_string(colnames(tsne.coordinates)[1],colnames(tsne.coordinates)[2])) +
stat_density_2d(geom="polygon",aes(alpha=..level..,fill=cluster),bins=5,size=2) + scale_x_continuous(limits=as.numeric(tsne.plot.x.lims),expand=c(0,0)) + scale_y_continuous(limits=as.numeric(tsne.plot.x.lims),expand=c(0,0)) + labs(alpha="Level",fill="Cluster") + theme(panel.background=element_blank(),axis.line=element_line(colour="black"),legend.position="bottom") + guides(alpha=guide_legend(nrow=2),fill=guide_legend(nrow=2)) + coord_equal()

hc <- hclust(dist(t(sapply(em_results[["maximization"]],function(x) x[["mu"]]))))

pdf(paste("density_plot_temp_",signif(split.clusters[["temperature"]],digits=2),".pdf",sep=""))

vp.left <- grid::viewport(height=unit(1,"npc"),width=unit(0.66,"npc"),just=c("right"),y=0.5,x=0.66)
par(mfrow=c(1,3))
plot.new()
plot.new()
print(px,vp=vp.left)
plot(hc,main=NULL,axes=F,xlab=NULL,ylab=NULL,ann=F)

dev.off()

saveRDS(em_results,file=paste("em_results_temp_",signif(split.clusters[["temperature"]],digits=2),".RDS",sep=""))
	
} else {
			
		if (decay.step[i,2]==temp.use) {
			temp.use <- decay.step[i+1,2]
		} else { if (decay.step[i,2]>temp.use) {
			temp.use <- decay.step[min(which(decay.step[,2]<temp.use)),2]
			} else {temp.use <- decay.step[i,2]}
		}
				
	if (temp.use>0) {
	
		print(paste("Starting loop ",i," at temperature=",signif(temp.use,digits=2),sep=""))
	
	em_results <- expect_max(initial.pca.components,initial.tsne.coordinates,em_results[["clusters"]],em_results[["maximization"]],max.iterations,delta.log.likelihood,temp=temp.use,num.cores)

	print(paste("Completed loop ",i," at temperature=", signif(temp.use,digits=2),sep=""))

overall.data <- data.frame(em_results[["tsne.coordinates"]],cluster=em_results[["clusters"]])

px <- ggplot2::ggplot(data=overall.data,aes_string(colnames(tsne.coordinates)[1],colnames(tsne.coordinates)[2])) +
stat_density_2d(geom="polygon",aes(alpha=..level..,fill=cluster),bins=5,size=2) + scale_x_continuous(limits=tsne.plot.x.lims,expand=c(0,0)) + scale_y_continuous(limits=tsne.plot.y.lims,expand=c(0,0)) + labs(alpha="Level",fill="Cluster") + theme(panel.background=element_blank(),axis.line=element_line(colour="black"),legend.position="bottom") + guides(alpha=guide_legend(nrow=2),fill=guide_legend(nrow=2)) + coord_equal()

hc <- hclust(dist(t(sapply(em_results[["maximization"]],function(x) x[["mu"]]))))

pdf(paste("density_plot_temp_",signif(temp.use,digits=2),".pdf",sep=""))

vp.left <- grid::viewport(height=unit(1,"npc"),width=unit(0.66,"npc"),just=c("right"),y=0.5,x=0.66)
par(mfrow=c(1,3))
plot.new()
plot.new()
print(px,vp=vp.left)
plot(hc,main=NULL,axes=F,xlab=NULL,ylab=NULL,ann=F)

dev.off()

saveRDS(em_results,file=paste("em_results_temp_",signif(temp.use,digits=2),".RDS",sep=""))

} else {
	
	print(paste("Starting loop ",i," at temperature=0",sep=""))
	
	em_results <- expect_max_zero(initial.pca.components,initial.tsne.coordinates,em_results[["clusters"]],em_results[["maximization"]],max.iterations,delta.log.likelihood,num.cores)
		
	print(paste("Completed loop ",i," at temperature=0",sep=""))

overall.data <- data.frame(em_results[["tsne.coordinates"]],cluster=em_results[["clusters"]])

px <- ggplot2::ggplot(data=overall.data,aes_string(colnames(tsne.coordinates)[1],colnames(tsne.coordinates)[2])) +
stat_density_2d(geom="polygon",aes(alpha=..level..,fill=cluster),bins=5,size=2) + scale_x_continuous(limits=tsne.plot.x.lims,expand=c(0,0)) + scale_y_continuous(limits=tsne.plot.y.lims,expand=c(0,0)) + labs(alpha="Level",fill="Cluster") + theme(panel.background=element_blank(),axis.line=element_line(colour="black"),legend.position="bottom") + guides(alpha=guide_legend(nrow=2),fill=guide_legend(nrow=2)) + coord_equal()

hc <- hclust(dist(t(sapply(em_results[["maximization"]],function(x) x[["mu"]]))))

pdf(paste("density_plot_temp_0.pdf",sep=""))

vp.left <- grid::viewport(height=unit(1,"npc"),width=unit(0.66,"npc"),just=c("right"),y=0.5,x=0.66)
par(mfrow=c(1,3))
plot.new()
plot.new()
print(px,vp=vp.left)
plot(hc,main=NULL,axes=F,xlab=NULL,ylab=NULL,ann=F)

dev.off()

p.final <- ggplot2::ggplot(data=overall.data,aes_string(colnames(tsne.coordinates)[1],colnames(tsne.coordinates)[2])) +
geom_point(aes(colour=cluster)) +
theme(panel.background=element_blank(),axis.line=element_line(colour="black"),legend.position="bottom") + labs(colour="Cluster")

pdf("final_gmm_clusters_plot.pdf")
print(p.final)
dev.off()

saveRDS(em_results,file="em_results_temp_0.RDS")

	return(em_results)

				} #close if for temp>0
				
			} #close else for max.clusters 
			
		} #close i > 1 else step
		
	} #close nrow(decay.step) loop
	
	parallel::stopCluster(cl)
	
} #close function