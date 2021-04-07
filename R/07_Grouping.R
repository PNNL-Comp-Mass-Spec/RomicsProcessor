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

  pca_coord<-data.frame(pca_results$ind$coord[,c(Xcomp,Ycomp, Zcomp)])
  colnames(pca_coord)<-c("PCA","PCB","PCC")
  pca_coord<-cbind(pca_coord,factor=factor)
  pca_coord$factor<-as.factor(pca_coord$factor)

  fig <- plot_ly(pca_coord, x = ~PCA, y = ~PCB, z = ~PCC, color = ~factor, colors = unique(colors))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = paste0("PC",Xcomp,"(",pca_results_explained[Xcomp],"%)")),
                                     yaxis = list(title = paste0("PC",Ycomp,"(",pca_results_explained[Ycomp],"%)")),
                                     zaxis = list(title = paste0("PC",Zcomp,"(",pca_results_explained[Zcomp],"%)"))))

  fig
}

#' romicsUmap()
#' @description Calculate the umap of the data layer of the romics_object using the package umap. This function require the data to be full in the layer data to be run.
#' Plots the samples on the PCA sample plot, the percentage of explained variance, or both on demand using ggplots2. The colors used for the plotting will correspond to the main_factor of the romics_object. The axis to be plotted can be chosen. The PCA results are calculated using the function romicsPCA (see the documentation of this function for more details).
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param umap_config object of class umap.config (see umap() function documentation)
#' @param method character, implementation. Available methods are 'naive' (an implementation written in pure R) and 'umap-learn' (requires python package 'umap-learn')
#' @param lock_seed has to be TRUE or FALSE to indicate if the seed is locked to enable reproducing plotting
#' @param seed numeric value indicating what seed to use when random seed are used
#' @param ... Arguments passed to umap function.
#' @details This function will
#' @return Returns a umap result
#' @author Geremy Clair
#' @export
romicsUmap<- function(romics_object, umap_config=umap.defaults, method = c("naive", "umap-learn"), lock_seed=TRUE, seed = 42){
  if(missing(romics_object)){stop("romics_object is missing")}
  if(class(romics_object)!="romics_object"){stop("your romics_object was not created using the function romicsCreateObject")}
  if(romics_object$steps[1]!="romics_object"){stop("romics_object is not in the appropriate format")}
  if(missing(umap_config)){umap_config=umap.defaults}
  if(missing(method)){method="naive"}
  if(missing(lock_seed)){lock_seed=TRUE}
  if(!is.logical(lock_seed)){stop("lock_seed has to be either TRUE or FALSE")}
  if(missing(seed)){seed<-42}

  #set seed so the plot will always look the sameif lock seed was TRUE else use a random seed.
  if (lock_seed==TRUE){
    set.seed(seed)}else{
      set.seed(Sys.time())}

  #load data
  data<-t(romics_object$data)

  #run the umap
  umap <- umap(data, umap_config, method=method)

  #return result
  return(umap)
}

#' romicsUmapPlot()
#' @description Plots the samples on the umap. Can be colored either by a factor or by the values of a given variable present in the dataset. This function require the data to be full in the layer data to be run. Calculate the umap of the data layer of the romics_object using the package umap.
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param umap_config object of class umap.config (see umap() function documentation)
#' @param method character, implementation. Available methods are 'naive' (an implementation written in pure R) and 'umap-learn' (requires python package 'umap-learn')
#' @param lock_seed has to be TRUE or FALSE to indicate if the seed is locked to enable reproducing plotting
#' @param color_by has to be either 'factor' or 'variable to indicate if the coloring of the sample points is done using the factor set with the parameter 'factor' or with the numerical values of a variable feature set with the factor 'variable'
#' @param factor has to be a factor of the romics_object, the existing factors names can be found using the function romicsFactorNames()
#' @param color_palette has to be a R gradient color palette (viridis(n=20) by default) this gradient color palette will be used to color the points when the 'color_by' parameter is set to 'variable'
#' @param scale has to be TRUE or FALSE to indicate if the palette defined in color_palette has to be displayed or not.
#' @param seed numeric value indicating what seed to use when random seed are used
#' @param ... further arguments passed to or from other methods
#' @return Returns either a ggplot2 showing the sample umap plot
#' @author Geremy Clair
#' @export

romicsUmapPlot<-function(romics_object, umap_config=umap.defaults, label=TRUE, method = c("naive", "umap-learn"), lock_seed=TRUE, seed = 42, color_by=c("factor","variable"), factor = "main", variable="variable_name",color_palette=viridis(n=20) , scale=TRUE){
  if(missing(romics_object)){stop("romics_object is missing")}
  if(class(romics_object)!="romics_object"){stop("your romics_object was not created using the function romicsCreateObject")}
  if(romics_object$steps[1]!="romics_object"){stop("romics_object is not in the appropriate format")}
  if(missing(umap_config)){umap_config=umap.defaults}
  if(missing(method)){method="naive"}
  if(missing(lock_seed)){lock_seed=TRUE}
  if(!is.logical(lock_seed)){stop("lock_seed as to be either TRUE or FALSE")}
  if(missing(seed)){seed<-42}
  if(missing(color_by)){color_by="factor"}
  if(!color_by %in% c("factor","variable")){stop("'color_by' has to be either 'variable' or 'factor'.")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){stop("factor has to be a either 'main' or a factor_name of the romics_object, please use the function romicsFactorNames() to identify the usable factors.")}
  if(missing(variable)){variable=rownames(romics_object$data)[1]}
  if(!is.character(variable)&length(variable!=1)){stop("<variable> should be a character vector of lenght 1")}
  if(lock_seed==TRUE){
    set.seed(seed)}else{
      set.seed(Sys.time())}
  if(missing(color_palette)){color_palette =viridis::viridis(20)}
  if(missing(scale)){scale=T}
  #load data
  data<-t(romics_object$data)

  #run the umap
  umap <- umap(data, umap_config, method=method)

  #coordinates object creation
  umap_coordinates<-data.frame(umap$layout)

  if(color_by=="factor"){
    if(factor=="main"){
      labels<-romics_object$metadata[rownames(romics_object$metadata)==romics_object$main_factor,]
      colors<-romics_object$metadata[rownames(romics_object$metadata)=="colors_romics",]
    }else{
      changed_factor<-romicsChangeFactor(romics_object,main_factor = factor)
      labels<-changed_factor$metadata[rownames(changed_factor$metadata)==factor,]
      colors<-changed_factor$metadata[rownames(changed_factor$metadata)=="colors_romics",]
    }

    umap_coordinates<-merge(umap_coordinates,t(labels),by=0)

    rownames(umap_coordinates)<-umap_coordinates$Row.names
    umap_coordinates<-umap_coordinates[-1]
    umap_coordinates<-merge(umap_coordinates,t(colors),by=0)
    rownames(umap_coordinates)<-umap_coordinates$Row.names
    umap_coordinates<-umap_coordinates[-1]
    colnames(umap_coordinates)<-c("umap1", "umap2","labels","colors")

    umap_plot<- ggplot(umap_coordinates, aes(x=umap1, y=umap2, color=labels)) +
      geom_point(size=4, alpha=I(.8))+
      theme_ROP()+
      ggtitle("umap plot")+
      scale_color_manual(values = unique(umap_coordinates$colors)) +labs(color = romics_object$main_factor)

  }else{
    variable<-colnames(data)[grepl(variable,colnames(data))]

    if(length(variable)>1){
      warning("More than one variable contained the variable choosen:")
      warning(variable)
      stop()
    }
    if(length(variable)==0){stop("your variable was not present in the romics_object")}
    variable_val<-data[,colnames(data)==variable]

    if(scale==T){
      scaled_variable <-as.numeric(scale(variable_val))
      names(scaled_variable)<-names(variable_val)
      variable_val<-scaled_variable
    }

    umap_coordinates<-cbind(umap_coordinates,variable_val)
    colnames(umap_coordinates)<-c("umap1", "umap2","variable_val")

    umap_plot<- ggplot(umap_coordinates, aes(x=umap1, y=umap2, color=variable_val)) +
      geom_point(size=4, alpha=I(.8))+ theme_ROP()+ ggtitle(paste0("umap plot colored by ",variable))+
      scale_color_gradientn(colours = color_palette)
  }

  if(label==TRUE){
    umap_plot<-umap_plot+geom_text(colour=umap_coordinates$colors, size = 3,label=rownames(data))
  }

  plot(umap_plot)
}

