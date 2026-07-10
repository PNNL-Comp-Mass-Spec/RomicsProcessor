#' FeatureCorrelation()
#' @description calculate the correlation matrix for the features of a romics_object
#' @param romics_object  an object of type romics_object
#' @param corr_type has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.
#' @param use has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'
#' @details calculate a correlation matrix for the features of an Romics_object
#' @return This function will return a correlation matrix for the features of an romics_object
#' @author Geremy Clair
#' @export
FeatureCorrelation<-function(romics_object,corr_type="pearson",use='everything'){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(corr_type)){corr_type="pearson"}
  if(!corr_type %in% c("pearson", "kendall","spearman")){stop("<corr_type> has to be either 'pearson', 'kendall',or 'spearman'")}
  if(missing(use)){use="everything"}
  if(!use %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop("<use> has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'")}

  d<-romics_object$data
  feat.cor <- cor(t(d), use = use, method = corr_type)

  return(feat.cor)
  }

#' SampleCorrelation()
#' @description calculate the correlation matrix for the samples of a romics_object
#' @param romics_object  an object of type romics_object
#' @param corr_type has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.
#' @param use has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'
#' @details calculate a correlation matrix for the samples of an Romics_object
#' @return This function will return a correlation matrix for the samples of an romics_object
#' @author Geremy Clair
#' @export
SampleCorrelation<-function(romics_object,corr_type="pearson",use='everything'){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(corr_type)){corr_type="pearson"}
  if(!corr_type %in% c("pearson", "kendall","spearman")){stop("<corr_type> has to be either 'pearson', 'kendall',or 'spearman'")}
  if(missing(use)){use="everything"}

  d<-romics_object$data
  NCores<-parallel::detectCores()

  if(ncol(d)>250){
    library(parallel)
    library(doParallel)
    cores<- NCores-3
    print(paste0("Due to the computing power required for such a large dataset, the computation was paralellized using",cores, "of your computer."))
    options('mc.cores' = cores)
    doParallel::registerDoParallel(cores)
    pb <- txtProgressBar(max=100, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)

    sample.cor <- foreach(i = seq_len(ncol(d)),
                 .combine = rbind,
                 .multicombine = TRUE,
                 .inorder = FALSE,
                 .packages = c('data.table', 'doParallel')) %dopar% {
                   cor(d[,i], d, method = corr_type,use=use)
                 }
    close(pb)}else{sample.cor <- cor(d, use = use, method = corr_type)}
  return(sample.cor)
}

#' FeatureCorrelationHclust()
#' @description calculate the feature correlation matrix for an romics_object and generate a hclust object.
#' @param romics_object  an object of type romics_object
#' @param corr_type has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.
#' @param use has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'
#' @param hclust_method has to be in 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'
#' @param dist_method has to be in 'correlation.dissimilarity','euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' @details calculate the feature correlation matrix for the samples of an romics object and generate a hclust object.
#' @return This function will return a hclust result that can be used to generate hclust plots.
#' @author Geremy Clair
#' @export
FeatureCorrelationHclust<-function(romics_object,corr_type="pearson",use='everything',hclust_method="ward.D",dist_method="euclidean"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(corr_type)){corr_type="pearson"}
  if(!corr_type %in% c("pearson", "kendall","spearman")){stop("<corr_type> has to be either 'pearson', 'kendall',or 'spearman'")}
  if(missing(use)){use="everything"}
  if(!use %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop("<use> has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'")}
  if(missing(hclust_method)){hclust_method<-"ward.D"}
  if(missing(dist_method)){dist_method<-"euclidean"}
  if(!dist_method %in% c("correlation.dissimilarity","euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")){
    stop("<dist_method> has to be in 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'")
  }
  if(!hclust_method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")){
    stop("<hclust_method> has to be in 'correlation.dissimilarity','ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'")
  }

  cor<-FeatureCorrelation(romics_object,corr_type = corr_type,use=use)
  if(sum(is.na(cor))>0){
    warning("Your data contains missing values, it is not recommended to use this function on non full data.
            The <NA> correlations were replaced by 0.")
    cor[is.na(cor)]<-0}

  if(dist_method=="correlation.dissimilarity"){hclust_Feature<-hclust(as.dist(1-cor),method = hclust_method)}else{
    hclust_Feature<-hclust(dist(cor,method=dist_method),method =  hclust_method)
  }

  return(hclust_Feature)
}

#' SampleCorrelationHclust()
#' @description calculate the sample correlation matrix for an romics_object and generate a hclust object.
#' @param romics_object  an object of type romics_object
#' @param corr_type has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.
#' @param use has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'
#' @param hclust_method has to be in 'correlation.dissimilarity','ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'
#' @param dist_method has to be in 'correlation.dissimilarity','euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' @details calculate the sample correlation matrix for the samples of an romics object and generate a hclust object.
#' @return This function will return a hclust result that can be used to generate hclust plots.
#' @author Geremy Clair
#' @export
SampleCorrelationHclust<-function(romics_object,corr_type="pearson",use='everything',hclust_method="ward.D",dist_method="euclidean"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(corr_type)){corr_type="pearson"}
  if(!corr_type %in% c("pearson", "kendall","spearman")){stop("<corr_type> has to be either 'pearson', 'kendall',or 'spearman'")}
  if(missing(use)){use="everything"}
  if(!use %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop("<use> has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'")}
  if(missing(hclust_method)){hclust_method<-"ward.D"}
  if(missing(dist_method)){dist_method<-"euclidean"}
  if(!dist_method %in% c("correlation.dissimilarity","euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")){
    stop("<dist_method> has to be in 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'")
  }
  if(!hclust_method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")){
    stop("<hclust_method> has to be in 'correlation.dissimilarity','ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'")
  }

  cor<-SampleCorrelation(romics_object,corr_type = corr_type,use=use)
  if(sum(is.na(cor))>0){
    warning("Your data contains missing values, it is not recommended to use this function on non full data.
            The <NA> correlations were replaced by 0.")
    cor[is.na(cor)]<-0}

  if(dist_method=="correlation.dissimilarity"){hclust_Sample<-hclust(as.dist(1-cor))}else{
    hclust_Sample<-hclust(dist(cor,method=dist_method),method =  hclust_method)
  }

  return(hclust_Sample)
}

#' FeatureCorrelationHclustPlot()
#' @description calculate the Feature correlation matrix for the samples of an romics_object and generate a hclust plot.
#' @param romics_object  an object of type romics_object
#' @param corr_type has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.
#' @param use has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'
#' @param hclust_method has to be in 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'
#' @param dist_method has to be in 'correlation.dissimilarity','euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' @param plot_type has to be either 'dendrogram', 'unrooted', or 'fan'
#' @param number_of_clusters has to be a numeric value to indicate the number of clusters
#' @param cex has to be a numeric value to indicate the label font size
#' @param cluster_colors has to be a character vector containing the list of colors used for the clusters.
#' @details calculate the Feature correlation matrix for the samples of an romics object and generate a hclust plot
#' @return This function will return a hclust plot generated with the base R plot() function.
#' @author Geremy Clair
#' @export
FeatureCorrelationHclustPlot<-function(romics_object,corr_type="pearson",use="everything",hclust_method="ward.D",dist_method="euclidean",plot_type = "dendrogram" ,number_of_clusters=0,cex=1,cluster_colors= c("c1","c2","c3","c4")){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(corr_type)){corr_type="pearson"}
  if(!corr_type %in% c("pearson", "kendall","spearman")){stop("<corr_type> has to be either 'pearson', 'kendall',or 'spearman'")}
  if(missing(use)){use="everything"}
  if(!use %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop("<use> has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'")}
  if(missing(number_of_clusters)){number_of_clusters=0}
  if(length(cex)!=1|!is.numeric(cex)){stop("<cex> has to be a numerical value.")}
  if(missing(cex)){cex<-1}
  if(missing(cluster_colors)){if(number_of_clusters==0){cluster_colors<-"gray50"}else{cluster_colors<- ROP_colors }}
  if(missing(hclust_method)){hclust_method<-"ward.D"}
  if(missing(dist_method)){dist_method<-"euclidean"}
  if(!dist_method %in% c("correlation.dissimilarity","euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")){
    stop("<dist_method> has to be in 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'")
  }
  if(!hclust_method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")){
    stop("<hclust_method> has to be in 'correlation.dissimilarity','ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'")
  }

  colors<-colors[1:number_of_clusters]

  cor<-FeatureCorrelation(romics_object,corr_type=corr_type,use=use)
  if(sum(is.na(cor))>0){
    warning("Your data contains missing values, it is not recommended to use this function on non full data.
            The <NA> correlations were replaced by 0.")
    cor[is.na(cor)]<-0}

  if(dist_method=="correlation.dissimilarity"){hclust_Feature<-hclust(as.dist(1-cor))}else{
    hclust_Feature<-hclust(dist(cor,method=dist_method),method =  hclust_method)
  }

  if(number_of_clusters>0){
    clust<-cutree(hclust_Feature,number_of_clusters)
    colors<-colors[clust]
    }else{colors<-rep("gray50",nrow(romics_object$data))}

  if(plot_type=="dendrogram"){
    dend <- hclust_Feature %>%
            as.dendrogram() %>%
            set("labels_cex",cex)

    if(number_of_clusters>0){
      dend<- dendextend::color_branches(dend, k=number_of_clusters,col=ROP_colors[1:number_of_clusters])
      dend<-dendextend::color_labels(dend, k=number_of_clusters,col=ROP_colors[1:number_of_clusters])
    }
  plot(dend)
  }

  if(plot_type=="unrooted"){
  plot(as.phylo(hclust_Feature),type="unrooted",tip.color=colors,cex=cex)
  }

  if(plot_type=="fan"){
  plot(as.phylo(hclust_Feature),type="fan",tip.color=colors,cex=cex)
  }
}

#' SampleCorrelationHclustPlot()
#' @description calculate the Sample correlation matrix for the samples of an romics_object and generate a hclust plot.
#' @param romics_object  an object of type romics_object
#' @param corr_type has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.
#' @param use has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'
#' @param hclust_method has to be in 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'
#' @param dist_method has to be in 'correlation.dissimilarity','euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' @param plot_type has to be either 'dendrogram', 'unrooted', or 'fan'
#' @param number_of_clusters has to be a numeric value to indicate the number of clusters
#' @param cex has to be a numeric value to indicate the label font size
#' @param cluster_colors has to be a character vector containing the list of colors used for the clusters.
#' @details calculate the Sample correlation matrix for the samples of an romics object and generate a hclust plot
#' @return This function will return a hclust plot generated with the base R plot() function.
#' @author Geremy Clair
#' @export
SampleCorrelationHclustPlot<-function(romics_object,corr_type="pearson",use="everything",hclust_method="ward.D",dist_method="euclidean",plot_type = "dendrogram" ,number_of_clusters=0,cex=1,cluster_colors= c("c1","c2","c3","c4")){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(corr_type)){corr_type="pearson"}
  if(!corr_type %in% c("pearson", "kendall","spearman")){stop("<corr_type> has to be either 'pearson', 'kendall',or 'spearman'")}
  if(missing(use)){use="everything"}
  if(!use %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop("<use> has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'")}
  if(missing(number_of_clusters)){number_of_clusters=0}
  if(length(cex)!=1|!is.numeric(cex)){stop("<cex> has to be a numerical value.")}
  if(missing(cex)){cex<-1}
  if(missing(cluster_colors)){if(number_of_clusters==0){cluster_colors<-"gray50"}else{cluster_colors<- ROP_colors }}
  if(missing(hclust_method)){hclust_method<-"ward.D"}
  if(missing(dist_method)){dist_method<-"euclidean"}
  if(!dist_method %in% c("correlation.dissimilarity","euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")){
    stop("<dist_method> has to be in 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'")
  }
  if(!hclust_method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")){
    stop("<hclust_method> has to be in 'correlation.dissimilarity','ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'")
  }

  colors<-colors[1:number_of_clusters]

  cor<-SampleCorrelation(romics_object,corr_type=corr_type,use=use)
  if(sum(is.na(cor))>0){
    warning("Your data contains missing values, it is not recommended to use this function on non full data.
            The <NA> correlations were replaced by 0.")
    cor[is.na(cor)]<-0}

  if(dist_method=="correlation.dissimilarity"){hclust_Sample<-hclust(as.dist(1-cor))}else{
    hclust_Sample<-hclust(dist(cor,method=dist_method),method =  hclust_method)
  }

  if(number_of_clusters>0){
    clust<-cutree(hclust_Sample,number_of_clusters)
    colors<-colors[clust]
  }else{colors<-rep("gray50",nrow(romics_object$data))}

  if(plot_type=="dendrogram"){
    dend <- hclust_Sample %>%
      as.dendrogram() %>%
      set("labels_cex",cex)

    if(number_of_clusters>0){
      dend<- dendextend::color_branches(dend, k=number_of_clusters,col=ROP_colors[1:number_of_clusters])
      dend<-dendextend::color_labels(dend, k=number_of_clusters,col=ROP_colors[1:number_of_clusters])
    }
    plot(dend)
  }

  if(plot_type=="unrooted"){
    plot(as.phylo(hclust_Sample),type="unrooted",tip.color=colors,cex=cex)
  }

  if(plot_type=="fan"){
    plot(as.phylo(hclust_Sample),type="fan",tip.color=colors,cex=cex)
  }
}

#' FeatureCorrelationHeatmap()
#' @description calculate the feature correlation matrix of an romics_object and hclust the results are displayed as an heatmap.
#' @param romics_object  an object of type romics_object
#' @param corr_type has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.
#' @param use has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'
#' @param hclust_method has to be in 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'
#' @param dist_method has to be in 'correlation.dissimilarity','euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' @param number_of_clusters has to be a numeric value to indicate the number of clusters
#' @param cex has to be a numeric value to indicate the label font size
#' @param cluster_colors has to be a character vector containing the list of colors used for the clusters.
#' @param heatmap_col has to be a color gradient (e.g.viridis(20), greenred(50))
#' @param cellnote has to be TRUE or FALSE to indicate if the correlation numbers should be displayed onto the heatmap.
#' @return This function will return a an heatmap generated with gplots::heatmap.2().
#' @author Geremy Clair
#' @export
FeatureCorrelationHeatmap<-function(romics_object,corr_type="pearson",use="everything",number_of_clusters=0,cex=1,cluster_colors= c("c1","c2","c3","c4"), heatmap_col=viridis(20),dist_method="euclidean",hclust_method="ward.D",cellnote=F){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(corr_type)){corr_type="pearson"}
  if(!corr_type %in% c("pearson", "kendall","spearman")){stop("<corr_type> has to be either 'pearson', 'kendall',or 'spearman'")}
  if(missing(use)){use="everything"}
  if(!use %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop("<use> has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'")}
  if(missing(number_of_clusters)){number_of_clusters=0}
  if(length(cex)!=1|!is.numeric(cex)){stop("<cex> has to be a numerical value.")}
  if(missing(cex)){cex<-1}
  if(missing(cluster_colors)){if(number_of_clusters==0){cluster_colors<-"gray50"}else{cluster_colors<- ROP_colors }}
  if(missing(hclust_method)){hclust_method<-"ward.D"}
  if(missing(dist_method)){dist_method<-"euclidean"}
  if(!dist_method %in% c("correlation.dissimilarity","euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")){
    stop("<dist_method> has to be in 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'")
  }
  if(!hclust_method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")){
    stop("<hclust_method> has to be in 'correlation.dissimilarity','ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'")
  }
  if(missing(heatmap_col)){heatmap_col=viridis(20)}
  if(missing(cellnote)){cellnote=F}

  cor<-FeatureCorrelation(romics_object,corr_type=corr_type,use=use)
  if(sum(is.na(cor))>0){
    warning("Your data contains missing values, it is not recommended to use this function on non full data.
            The <NA> correlations were replaced by 0.")
    cor[is.na(cor)]<-0}

  if(dist_method=="correlation.dissimilarity"){hclust_Feature<-hclust(as.dist(1-cor))}else{
    hclust_Feature<-hclust(dist(cor,method=dist_method),method =  hclust_method)
  }
  if(missing(cellnote)){cellnote=FALSE}
  if(cellnote!=TRUE&cellnote!=FALSE){stop("cellnote should either be TRUE or FALSE")}

  #create colors
  colors<-cluster_colors[1:number_of_clusters]
  if(number_of_clusters>0){
    clust<-cutree(hclust_Feature,number_of_clusters)
    colors<-colors[clust]
  }else{colors<-rep("gray50",nrow(romics_object$data))}

  #generate dend
  dend <- hclust_Feature %>%
    as.dendrogram() %>%
    dendextend::set("labels_cex",cex)
  if(number_of_clusters>0){
    dend<- color_branches(dend, k=number_of_clusters,col=ROP_colors[1:number_of_clusters])
    dend<-color_labels(dend, k=number_of_clusters,col=ROP_colors[1:number_of_clusters])
  }

  if(cellnote==T){gplots::heatmap.2(cor,Rowv = dend,Colv = dend,trace="none",cellnote = round(cor,2),notecol="gray50",col =heatmap_col)}else{
    gplots::heatmap.2(cor,Rowv = dend,Colv = dend,col =heatmap_col,trace="none")
  }
}

#' SampleCorrelationHeatmap()
#' @description calculate the sample correlation matrix of an romics_object and hclust the results are displayed as an heatmap.
#' @param romics_object  an object of type romics_object
#' @param corr_type has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.
#' @param use has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'
#' @param hclust_method has to be in 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'
#' @param dist_method has to be in 'correlation.dissimilarity','euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' @param number_of_clusters has to be a numeric value to indicate the number of clusters
#' @param cex has to be a numeric value to indicate the label font size
#' @param cluster_colors has to be a character vector containing the list of colors used for the clusters.
#' @param heatmap_col has to be a color gradient (e.g.viridis(20), greenred(50))
#' @param cellnote has to be TRUE or FALSE to indicate if the correlation numbers should be displayed onto the heatmap.
#' @return This function will return a an heatmap generated with gplots::heatmap.2().
#' @author Geremy Clair
#' @export
SampleCorrelationHeatmap<-function(romics_object,corr_type="pearson",use="everything",number_of_clusters=0,cex=1,cluster_colors= c("c1","c2","c3","c4"), heatmap_col=viridis(20),dist_method="euclidean",hclust_method="ward.D",cellnote=F){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(corr_type)){corr_type="pearson"}
  if(!corr_type %in% c("pearson", "kendall","spearman")){stop("<corr_type> has to be either 'pearson', 'kendall',or 'spearman'")}
  if(missing(use)){use="everything"}
  if(!use %in% c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")){
    stop("<use> has to be either 'everything', 'all.obs', 'complete.obs', 'na.or.complete', or 'pairwise.complete.obs'")}
  if(missing(number_of_clusters)){number_of_clusters=0}
  if(length(cex)!=1|!is.numeric(cex)){stop("<cex> has to be a numerical value.")}
  if(missing(cex)){cex<-1}
  if(missing(cluster_colors)){if(number_of_clusters==0){cluster_colors<-"gray50"}else{cluster_colors<- ROP_colors }}
  if(missing(hclust_method)){hclust_method<-"ward.D"}
  if(missing(dist_method)){dist_method<-"euclidean"}
  if(!dist_method %in% c("correlation.dissimilarity","euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")){
    stop("<dist_method> has to be in 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'")
  }
  if(!hclust_method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")){
    stop("<hclust_method> has to be in 'correlation.dissimilarity','ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'")
  }
  if(missing(heatmap_col)){heatmap_col=viridis(20)}
  if(missing(cellnote)){cellnote=F}

  cor<-SampleCorrelation(romics_object,corr_type=corr_type,use=use)
  if(sum(is.na(cor))>0){
    warning("Your data contains missing values, it is not recommended to use this function on non full data.
            The <NA> correlations were replaced by 0.")
    cor[is.na(cor)]<-0}

  if(dist_method=="correlation.dissimilarity"){hclust_Sample<-hclust(as.dist(1-cor))}else{
    hclust_Sample<-hclust(dist(cor,method=dist_method),method =  hclust_method)
  }
  if(missing(cellnote)){cellnote=FALSE}
  if(cellnote!=TRUE&cellnote!=FALSE){stop("cellnote should either be TRUE or FALSE")}

  #create colors
  colors<-cluster_colors[1:number_of_clusters]
  if(number_of_clusters>0){
    clust<-cutree(hclust_Sample,number_of_clusters)
    colors<-colors[clust]
  }else{colors<-rep("gray50",nrow(romics_object$data))}

  #generate dend
  dend <- hclust_Sample %>%
    as.dendrogram() %>%
    dendextend::set("labels_cex",cex)
  if(number_of_clusters>0){
    dend<- color_branches(dend, k=number_of_clusters,col=ROP_colors[1:number_of_clusters])
    dend<-color_labels(dend, k=number_of_clusters,col=ROP_colors[1:number_of_clusters])
  }

  if(cellnote==T){gplots::heatmap.2(cor,Rowv = dend,Colv = dend,trace="none",cellnote = round(cor,2),notecol="gray50",col =heatmap_col)}else{
    gplots::heatmap.2(cor,Rowv = dend,Colv = dend,col =heatmap_col,trace="none")
  }
}
