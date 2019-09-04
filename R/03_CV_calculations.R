#' romicsCVs()
#' Calculate the CVs (works properly for non logged romics_objects only) for each level of a metadata factor and eventually plot the CV boxplot and/or barplot (on demand)
#'
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param factor A string indicating a factor from the romics_object. The factors usable for a given romics_object names can be obtain using the function romicsFactorNames(romics_object)
#' @param plot can be either FALSE, 'boxplot', 'barplot', or 'all'
#'
#' @details This function will calculate the CVs, the percentage CVs, and on demand will plot the barplot and/or the boxplot of the CVs per sampletype. The CVs are calculated on unlogged romics_object only.
#' @return
#'
#' @author Geremy Clair
#' @export
#'

#this function needs to be modified to be able to act by factor
RomicsCVs<-function(romics_object, factor="main", plot="all"){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor)|factor=="main"){factor=romics_object$main_factor}
  if(!factor %in% romicsFactorNames(romics_object)){
    print("Your factor does not exist in your data, here is the list of the existing factors:")
    print(romicsFactorNames(romics_object))
    stop()
  }
  if(missing(plot)){
    print("<plot> was not defined, by defaults the CVs will not be plotted.")
    plot=FALSE}
  if(!plot %in% c(FALSE,"boxplot","barplot","all")){
    print("<plot> was not properly defined, by defaults the CVs will not be plotted.")
    plot=FALSE}

  #obtain the levels_factor of the selected factor and the colors from the romics_object
  levels_factor<-unique(as.character(t(romics_object$metadata[rownames(romics_object$metadata)==factor,])))
  colors<-ROP_colors[1:length(levels_factor)]

  #Create the CVs table
  CVs<-data.frame(matrix(nrow=nrow(romics_object$data), ncol=length(levels_factor)))
  rownames(CVs)<-rownames(romics_object$data)
  colnames(CVs)<-levels_factor

  #Create a table that will contain the mean and median of the CVs
  CVs_mean_and_median <- data.frame(matrix(nrow = length(levels_factor), ncol = 2))
  colnames(CVs_mean_and_median)<- c("CV_Mean", "CV_Median")
  rownames(CVs_mean_and_median)<- levels_factor

  Stdev<-as.numeric()
  for (i in 1:length(levels_factor)){
    for(j in 1:nrow(romics_object$data)){
      Stdev[j]<-sd(romics_object$data[j,romics_object$metadata[rownames(romics_object$metadata)==factor,]==levels_factor[i]])
    }
    Means<-rowMeans(romics_object$data[,romics_object$metadata[rownames(romics_object$metadata)==factor,]==levels_factor[i]] )

    CVs[,i]<-(Stdev/Means)*100

    CVs_mean_and_median[i,2]<-median(CVs[,i],na.rm=T)
  }

  CVs_mean_and_median$CV_Mean<-colMeans(CVs, na.rm=T)

  if(plot=="boxplot"|plot=="all"){boxplot(CVs,ylims=c(0,100),col=colors)}

  if(plot=="barplot"|plot=="all"){
    CVs_freq<-data.frame(matrix(ncol=ncol(CVs),nrow=10))
    rownames(CVs_freq) <- c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")
    colnames(CVs_freq) <- colnames(CVs)
    for(i in 1:ncol(CVs_freq)){
      CVs_freq[,i]<- c( sum(0<=CVs[,i] &CVs[,i]<10,na.rm=T),
                        sum(10<=CVs[,i] &CVs[,i]<20,na.rm=T),
                        sum(20<=CVs[,i] &CVs[,i]<30,na.rm=T),
                        sum(30<=CVs[,i] &CVs[,i]<40,na.rm=T),
                        sum(40<=CVs[,i] &CVs[,i]<50,na.rm=T),
                        sum(50<=CVs[,i] &CVs[,i]<60,na.rm=T),
                        sum(60<=CVs[,i] &CVs[,i]<70,na.rm=T),
                        sum(70<=CVs[,i] &CVs[,i]<80,na.rm=T),
                        sum(80<=CVs[,i] &CVs[,i]<90,na.rm=T),
                        sum(90<=CVs[,i] &CVs[,i]<10,na.rm=T))
    }
    barplot(t(CVs_freq),ylab="Number of peptides", beside=T,col=colors )
    legend("topright", levels_factor, cex=1.3, bty="n", fill=colors)
  }

  results<-list(CVs_table=CVs, CVs_means_and_medians=CVs_mean_and_median)

  romics_object<-romicsUpdateSteps(romics_object,arguments)

  if(plot==FALSE){return(results)}
}




