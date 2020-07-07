#' romicsHclust()
#' @description Plots the hierarchical clustering of the samples contaned in the romics_object calculated using the functions dist() and hclust() (see documentation of those R functions for more details). The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param method_dist the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param method_hclust the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @details This function uses the dist() and hclust() functions to calculate the hierachical clustering and then plots the hclust with colors based on the current main_factor of the romics_object.
#' @return a hierarchical clustering tree plot with its branches colored by factor.
#' @author Geremy Clair
#' @export
romicsHclust<-function(romics_object, method_dist="euclidean", method_hclust="ward.D" ){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if (missing(method_dist)){method_dist<-"euclidean"}
  if(!method_dist %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){stop("The method_dist used is not appropriate")}
  if(missing(method_hclust)){method_hclust<-"ward.D"}
  if(!method_hclust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")){stop("The method_hclust used is not appropriate")}

  #extract the data
  data<-romics_object$data
  data<-t(data)

  #extract the colors
  colors_dend <- as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))

  #define the clustering method
  distance_samples<-dist(data, method=method_dist)

  #define the hclust agglomeration method
  hc<- hclust(distance_samples, method=method_hclust)

  #convert the hclust into a dendrogram
  dd <- as.dendrogram(hc)

  #order the colors as the hclust
  colors_dend<-colors_dend[order=hc$order]

  #color the branches
  dd <- color_branches(dd,k = NULL, h = NULL,col=as.vector(colors_dend))

  #plot the dendrogram to choose the number of clusters
  par(cex=0.7, mar=c(8, 8, 4, 1))
  plot(dd,main="Hierarchical clustering samples", sub="")
}

#' romicsPCA()
#' @description Calculate the PCA of the data layer of the romics_object using the package FactoMineR. If the data layer contains some missing values those will be imputed using the missMDA::imputePCA() method (see the documentation of this function for more details). This function will return the PCA results and not a romics_object
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param ncp inherited from missmda::imputePCA().
#' @param scale  inherited from missmda::imputePCA(). boolean. TRUE implies a same weight for each variable
#' @param method  inherited from missmda::imputePCA(). "Regularized" by default or "EM". TRUE implies a same weight for each variable
#' @param row.w inherited from missmda::imputePCA(). row weights (by default, a vector of 1 for uniform row weights)
#' @param ncp.min used only if ncp is not set. inherited from missmda::estim_ncpPCA().integer corresponding to the minimum number of components to test
#' @param ncp.max used only if ncp is not set. inherited from missmda::estim_ncpPCA().integer corresponding to the minimum number of components to test
#' @param ... further arguments passed to or from other methods
#' @details This function uses the dist() and hclust() functions to calculate the hierachical clustering and then plots the hclust with colors based on the current main_factor of the romics_object.
#' @return Return the results of the PCA performed on the current version of the romics_object
#' @author Geremy Clair
#' @export
romicsPCA<-function(romics_object, ncp=5, scale = TRUE, method=c("Regularized","EM"), ncp.min = 0, ncp.max = 5, method.cv = c("gcv","loo","Kfold"),...){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  data<-t(romics_object$data)

  if(sum(is.na(data))>0){
    warning("Your romics_object data layer was containing some missing values. The missmda::imputePCA() method was used to impute these values")
    if(missing("scale")){scale=TRUE}
    if(missing("method")){method="Regularized"}
    if(missing(ncp)){
      if(missing(ncp.min)){ncp.min=0}
      if(missing(ncp.max)){ncp.max=5}
      if(missing(method.cv)){method.cv="Kfold"}
      warning("ncp was not defined, it will be automatically decided using the function missMDA::estim_ncpPCA() function, this might take a few minuts.")
      ncp <- missMDA::estim_ncpPCA(data,ncp.min=ncp.min,ncp.max=ncp.max,method.cv = method.cv, ...)
      print("Below is the result of the estim_ncpPCA()")
      print(ncp)
      ncp <- ncp$ncp
      }
    data<-missMDA::imputePCA(data, ncp = ncp, scale = scale, method = method, ...)
    }

  #run the PCA without showing the graphics
  pdf(file = NULL)
  pca_results<-FactoMineR::PCA(data)
  dev.off()

  return(pca_results)
}

#' indPCAplot()
#' @description Plots the samples on the PCA sample plot, the percentage of explained variance, or both on demand using ggplots2. The colors used for the plotting will correspond to the main_factor of the romics_object. The axis to be plotted can be chosen. The PCA results are calculated using the function romicsPCA (see the documentation of this function for more details).
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param Xcomp numerical/double. Indicate the cp to plot on the X axis
#' @param Ycomp numerical/double. Indicate the cp to plot on the Y axis
#' @param label boolean. Indicate if the sample name label should be plotted
#' @param plotType should be one of the following options to indicate the type of plot to be returned : 'dual'(for both), 'individual', or 'percentage'
#' @param ... further arguments passed to or from other methods
#' @details This function will plot the results of a PCA calculated on the romics_object data layer using the function romicsPCA() (see documentation for more details). The function can plot different plots based on user input. The first type of plot is a classical sample PCA plot the principal components to be plotted on each axis depend on the user input. The second type of plot is the calculation of the percentage of variance explained by each component. By using dual both plot will be returned. the plots are generated using ggplot2 and are subsequently adjustable using ggplot2 commands.
#' @return Returns either one ggplot2 or a combination plot generated with grid.arrange. On the sample PCA plot the colors plotted correspond to the main_factor utilized.
#' @author Geremy Clair
#' @export
indPCAplot <- function(romics_object, Xcomp=1, Ycomp=2, label=TRUE, plotType="dual", ... ){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(label)){label=TRUE}
  if(missing(plotType)){plotType="dual"}
  if(!plotType %in% c("dual","individual","percentage")){stop("Your <plotType> should be either dual, individual, or percentage")}

  #load data, factor and colors
  data<-t(romics_object$data)
  factor<-as.character(t(romics_object$metadata[rownames(romics_object$metadata)==romics_object$main_factor,]))
  colors<-as.character(t(romics_object$metadata[nrow(romics_object$metadata),]))

  #run the PCA and extract the results explained
  pca_results<-romicsPCA(romics_object, ...)
  pca_results_explained <- data.frame(pca_results$eig[,2])
  pca_results_explained <- round(as.vector(pca_results_explained[,1]),2)

  #verify the Xcomp and Ycomp
  if(missing(Xcomp) | !(is.numeric(Xcomp)|is.double(Xcomp))){
    warning("The component to be plotted on the first axis was not define PC1 will be used")
    Xcomp=1}
  if(missing(Ycomp)){
    warning("The component to be plotted on the first axis was not define PC2 will be used")
    Ycomp=2}
  if(Xcomp>nrow(pca_results$eig)){stop("The Xcomp selected is too large for the data you are using")}
  if(Ycomp>nrow(pca_results$eig)){stop("The Ycomp selected is too large for the data you are using")}

  #plots the component percentage of variance explained plot
  explained_plot<- fviz_screeplot(pca_results,barcolor= "gray20", barfill = "gray20")+theme_ROP()+theme(axis.text.x = element_text(angle = 0))+ggtitle("Percentage of exp. variance")

  #select the appropriate scale
  pca_plot_scale<- max(c(ceiling(max(pca_results$ind$coord[,c(Xcomp,Ycomp)])/10)+1)*10,abs(min(floor(min(pca_results$ind$coord[,c(Xcomp,Ycomp)])/10)-1)*10))

  #pca coordinate object
  pca_coordinates<-as.data.frame(pca_results$ind$coord[,c(Xcomp,Ycomp)])
  pca_coordinates$factor<-factor

  #create the indiv plot
  pca_indiv_plot<- ggplot(pca_coordinates, aes(x=pca_coordinates[,1], y=pca_coordinates[,2]))+
    geom_point(aes(colour= factor),size = 3,alpha=I(.8)) +
    scale_y_continuous(limits=c(-pca_plot_scale, pca_plot_scale))+
    scale_x_continuous(limits=c(-pca_plot_scale, pca_plot_scale))+
    scale_color_manual(values = unique (colors))+
    xlab(paste0("PC",Xcomp,"(",pca_results_explained[Xcomp],"%)"))+ylab(paste0("PC",Ycomp,"(",pca_results_explained[Ycomp],"%)"))+
    ggtitle("Principal component analysis")+
    labs(colour = romics_object$main_factor)+
    theme_ROP()

  #if label is TRUE add labels
  if(label==TRUE){pca_indiv_plot<-pca_indiv_plot+geom_text(aes(colour=factor), size = 3,label=row.names(data))}

  #plot the results
  if(plotType=="dual"){grid.arrange(explained_plot,pca_indiv_plot,ncol=2)}
  if(plotType=="individual"){plot(pca_indiv_plot)}
  if(plotType=="percentage"){plot(explained_plot)}

}

#' indPCA3D()
#' @description Plots the coordinate of the samples on a PCA plot in 3D. The colors used for the plotting will correspond to the main_factor of the romics_object. The axis to be plotted can be chosen. The PCA results are calculated using the function romicsPCA (see the documentation of this function for more details).
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param Xcomp numerical/double. Indicate the cp to plot on the X axis
#' @param Ycomp numerical/double. Indicate the cp to plot on the Y axis
#' @param Zcomp numerical/double. Indicate the cp to plot on the Z axis
#' @param ... further arguments passed to or from other methods
#' @details This function will plot the results of a PCA calculated on the romics_object data layer using the function romicsPCA() (see documentation for more details).
#' @return Returns a 3D plot generated with the plot3D::scatter3D() function.
#' @author Geremy Clair
#' @export
indPCA3D <- function(romics_object, Xcomp=1, Ycomp=2, Zcomp=3, ... ){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  #load data, factor and colors
  data<-t(romics_object$data)
  factor<-as.character(t(romics_object$metadata[rownames(romics_object$metadata)==romics_object$main_factor,]))
  colors<-as.character(t(romics_object$metadata[nrow(romics_object$metadata),]))

  #run the PCA and extract the results explained
  pca_results<-romicsPCA(romics_object, ...)
  pca_results_explained <- data.frame(pca_results$eig[,2])
  pca_results_explained <- round(as.vector(pca_results_explained[,1]),2)

  #verify the Xcomp and Ycomp
  if(missing(Xcomp) | !(is.numeric(Xcomp)|is.double(Xcomp))){
    warning("The component to be plotted on the first axis was not define PC1 will be used")
    Xcomp=1}
  if(missing(Ycomp)){
    warning("The component to be plotted on the first axis was not define PC2 will be used")
    Ycomp=2}
  if(missing(Zcomp)){
    warning("The component to be plotted on the first axis was not define PC3 will be used")
    Zcomp=3}
  if(Xcomp>nrow(pca_results$eig)){stop("The Xcomp selected is too large for the data you are using")}
  if(Ycomp>nrow(pca_results$eig)){stop("The Ycomp selected is too large for the data you are using")}
  if(Zcomp>nrow(pca_results$eig)){stop("The Zcomp selected is too large for the data you are using")}

  #select the appropriate scale
  pca_plot_scale<- max(c(ceiling(max(pca_results$ind$coord[,c(Xcomp,Ycomp, Zcomp)])/10)+1)*10,abs(min(floor(min(pca_results$ind$coord[,c(Xcomp,Ycomp, Zcomp)])/10)-1)*10))

  #pca coordinate object
  pca_coordinates<-as.data.frame(pca_results$ind$coord[,c(Xcomp,Ycomp)])

  #create the 3D plot
  scatter3D(pca_results$ind$coord[,Xcomp],
            pca_results$ind$coord[,Ycomp],
            pca_results$ind$coord[,Zcomp],
            colvar= as.integer(as.factor(factor)),col=colors,
            pch = 16,
            bty = "b2",
            phi=60,
            ticktype = "detailed",
            main ="Principal component analysis",
            col.panel ="gray",
            col.grid = "gray95",
            type = "h",
            colkey = list(plot = FALSE),
            xlim = c(-pca_plot_scale,pca_plot_scale),
            ylim = c(-pca_plot_scale,pca_plot_scale),
            zlim = c(-pca_plot_scale,pca_plot_scale),
            xlab = paste0("PC",Xcomp,"(",pca_results_explained[Xcomp],"%)"),
            ylab =paste0("PC",Ycomp,"(",pca_results_explained[Ycomp],"%)"),
            zlab = paste0("PC",Zcomp,"(",pca_results_explained[Zcomp],"%)"))
}
