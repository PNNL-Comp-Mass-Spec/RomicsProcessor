#' romicsMean()
#' @description Calculates the means of each variable within each level of the selected factor and add the generated columns in the statistics layer of the romics_object.
#' @param romics_object A object created using the function romicsCreateObject().
#' @param main_factor Either 'main' OR any factor from the romics_object, the list of factors from a romics object can be obtained using the function romicsFactorNames().
#' @details Adds the means columns to the statistics layer. enable to choose a different factor than the main one to do those calculation
#' @return This function returns a modified romics object, containing mean columns in the statistics layer.
#' @author Geremy Clair
#' @export
romicsMean<-function(romics_object, factor="main"){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){
    warning("The selected factor is not in the list of factors of the romics_object")
    warning(romicsFactorNames(romics_object))
    }

  #load data
  data<-romics_object$data

  #check if the $statistic object already is a part of the Romics_object (if not create it)
  if(is.null(romics_object$statistics)){
    print("The Statistics layer was added to your object")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if the statistics object does not have the same number of rows as the data replace it by a null statistics object
  if(!is.null(romics_object$statistics)&&nrow(romics_object$statistics)!=nrow(data)){
    warning("The Statistics layer was not containing the same number of rows as your data, it was replaced by an empty statistics layer")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if factor is main extract the factor from the romics_object$main_factor
  if(factor=="main"){factor<-romics_object$main_factor}
  #extract the factor from the metadata
  factor<-as.factor(as.character(t(romics_object$metadata[romicsFactorNames(romics_object)==factor,])))

  #find the levels for this factor
  levels_factor<-levels(factor)
  #create a list columns
  means_names <- character()
  #create a means object
  means<-data.frame(matrix(nrow=nrow(data),ncol=0))
  rownames(means)<-rownames(romics_object$statistics)
  #fill this table
  for(i in 1:length(levels_factor)){
    means_names[i] <- paste(levels_factor[i],"_mean",sep="")
    for (j in 1:nrow(romics_object$statistics)){
      means[j,i]<-mean(as.numeric(data[j,factor==levels_factor[i]]))}}

  colnames(means) <- means_names

  #remove old columns if means were previously existing
  if(sum(means_names %in% colnames(romics_object$statistics))>0){
    warning("The means were previously calculated, the old mean columns were removed.")
    romics_object$statistics<- romics_object$statistics[,!colnames(romics_object$statistics) %in% means_names]
  }

  #add the new columns
  romics_object$statistics <- cbind(romics_object$statistics,means)

  #print info
  print("Means columns (*_mean) were added to the statistics")

  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return the new object
  return(romics_object)
}

#' romicsSd()
#' @description Calculates the standard deviation of each variable within each level of the selected factor and add the generated columns in the statistics layer of the romics_object.
#' @param romics_object A object created using the function romicsCreateObject().
#' @param main_factor Either 'main' OR any factor from the romics_object, the list of factors from an romics object can be obtained using the function romicsFactorNames().
#' @details Adds the sd columns to the statistics layer. enable to choose a different factor than the main one to do those calculation
#' @return This function returns a modified romics object, containing sd columns in the statistics layer.
#' @author Geremy Clair
#' @export
romicsSd<-function(romics_object, factor="main"){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){
    warning("The selected factor is not in the list of factors of the romics_object")
    warning(romicsFactorNames(romics_object))
  }

  #load data
  data<-romics_object$data

  #check if the $statistic object already is a part of the Romics_object (if not create it)
  if(is.null(romics_object$statistics)){
    print("The Statistics layer was added to your object")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if the statistics object does not have the same number of rows as the data replace it by a null statistics object
  if(!is.null(romics_object$statistics)&&nrow(romics_object$statistics)!=nrow(data)){
    warning("The Statistics layer was not containing the same number of rows as your data, it was replaced by an empty statistics layer")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if factor is main extract the factor from the romics_object$main_factor
  if(factor=="main"){factor<-romics_object$main_factor}
  #extract the factor from the metadata
  factor<-as.factor(as.character(t(romics_object$metadata[romicsFactorNames(romics_object)==factor,])))

  #find the levels for this factor
  levels_factor<-levels(factor)
  #create a list columns
  sd_names <- character()
  #create a means object
  sd<-data.frame(matrix(nrow=nrow(data),ncol=0))
  rownames(sd)<-rownames(romics_object$statistics)
  #fill this table
  for(i in 1:length(levels_factor)){
    sd_names[i] <- paste(levels_factor[i],"_sd",sep="")
    for (j in 1:nrow(romics_object$statistics)){
      sd[j,i] <- sd(as.numeric(data[j,factor==levels_factor[i]]))}}
  colnames(sd) <- sd_names

  #remove old columns if sd were previously existing
  if(sum(sd_names %in% colnames(romics_object$statistics))>0){
    warning("The standard deviations were previously calculated, the old mean columns were removed.")
    romics_object$statistics<- romics_object$statistics[,!colnames(romics_object$statistics) %in% sd_names]
  }

  #add the new columns
  romics_object$statistics <- cbind(romics_object$statistics,sd)

  #print info
  print("The standard deviation columns (*_sd) were added to the statistics")

  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return romics_object
  return(romics_object)
}

#' romicsTtest()
#' @description Performs all possible paired-T.tests for each variable using the levels of the selected factor of the romics_object. The results are added as new columns in the statistics layers.
#' @param romics_object has to be an romics_object created with the function romicsCreateObject(),
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param paired a logical indicating whether you want a paired t-test.
#' @param pairing_factor name of a factor contained in an r_object, the list of the available factor can be obtained by using the function romics_factors()
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param padj a logical variable indincating wheter to perform or not adjustment of pvalues
#' @param padj_method correction method. Must be in  {"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr","none"}
#' @param mode must be in {"vs", "enrichments"} indicates if the groups should be compaired by pair or against all other groups
#' @param factor a character string indicating the factor to use for the test, the list of the available factor can be obtained by using the function romics_factors(), if missing the function will use the main factor of the object
#' @param percentage_completeness a numerical value comprised between 0 and 100 to indicate the minimum completeness required in at least one group calculate the T.test (if set completeness is not met, p and fold change will be NA)
#' @param reverse_order a boolean to indicate if the order of the factors needs to be reversed (this will make the calculated fold changes values A/B become B/A or log2(A/B) become log2(B/A))
#' @details When paired T.test are performed it is possible to include a second factor to generate the pairs. This function will also calculate the fold-changes or log2(fold-change). Please, note that the test will automatically determine if a log tranformation was performed to the object, subsequently we recommend to import not pre-logged data.frames when creating the object. For paired T.tests, it is possible to set a second factor containing the pairs, if missing the function will consider the pairs based on the column order in the romics_object.
#' @return an romics_object with the statistical layer containing the newly generated t.tests and fold-changes
#' @author Geremy Clair
#' @export
romicsTtest<-function(romics_object, alternative="two.sided", paired = FALSE, pairing_factor="none", var.equal=FALSE, factor = "main", padj=TRUE, padj_method="BH",mode="vs",percentage_completeness=0,reverse_order=FALSE, ...){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){
    warning("The selected factor is not in the list of factors of the romics_object")
    warning(romicsFactorNames(romics_object))
  }
  if(missing(paired)){paired<-FALSE}
  if(missing(reverse_order)){reverse_order<-FALSE}
  if(missing(var.equal)){var.equal<-FALSE}
  if(missing(padj)){padj<-TRUE}
  if(missing(padj_method)){padj_method<-"BH"}
  if(!padj_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){stop("padj_method has to be in the following list: holm, hochberg, hommel, bonferroni, BH, BY,fdr, none")}
  if(missing(pairing_factor)){pairing_factor ="none"}
  if(missing(mode)){mode="vs"}
  if(!mode %in% c("vs", "enrichment")){stop("'mode' has to be either 'vs' or 'enrichment' to indicate if the comparison should be done between all possible groups or if the features within a given group should be compared to all the other groups (enrichment).")}
  if(missing(percentage_completeness)){percentage_completeness<-0}
  if(percentage_completeness<0 & percentage_completeness>100){stop("the completeness has to be comprised between 0 and 100 %")}

  #load data
  data<-romics_object$data

  #load missing
  missingdata<-romics_object$missingdata

  #check if the $statistic object already is a part of the Romics_object (if not create it)
  if(is.null(romics_object$statistics)){
    print("The Statistics layer was added to your object")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if the statistics object does not have the same number of rows as the data replace it by a null statistics object
  if(!is.null(romics_object$statistics)&&nrow(romics_object$statistics)!=nrow(data)){
    warning("The Statistics layer was not containing the same number of rows as your data, it was replaced by an empty statistics layer")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if factor is main extract the factor from the romics_object$main_factor
  if(factor=="main"){factor<-romics_object$main_factor}
  #extract the factor from the metadata
  factor<-as.factor(as.character(t(romics_object$metadata[romicsFactorNames(romics_object)==factor,])))
  factor<-as.factor(t(factor))

  if(pairing_factor=="none"){
    pairing_factor<-romics_object$main_factor
  }else{
    if(!pairing_factor %in% rownames(romics_object$metadata)){stop("Your paired_factor is not present in the metadata")}}

  #extract paired_factor from the romics_object
  pairing_factor<-romics_object$metadata[pairing_factor==rownames(romics_object$metadata),]
  pairing_factor<-as.factor(t(pairing_factor))

  #order the data, the factor AND the paired_factor based on the factor and the paired_factor
  data<-data[,order(factor,pairing_factor)]
  missingdata<-missingdata[,order(factor,pairing_factor)]

  #create an object containing both factor to sort the two at the time
  both_factor<-paste0(factor,"@",pairing_factor)
  both_factor<-both_factor[order(factor,pairing_factor)]

  #re-extract factors
  factor<-as.factor(sapply(strsplit(both_factor,"@"), function(x) x[1]))
  pairing_factor<-as.factor(sapply(strsplit(both_factor,"@"), function(x) x[2]))

  #extract the levels to be considered
  levels_factor<-levels(factor)

  #create T_results and fold_change
  t_result <- vector(mode="numeric",length=nrow(data))
  fold_change <- vector(mode="numeric",length=nrow(data))
  t_padj<- vector(mode="numeric",length=nrow(data))

  #Create a T_table
  T_table<- data.frame(matrix(nrow=nrow(data),ncol=0))

  #determine if the data was log transformed
  if(sum(grepl("log10transform",romics_object$steps))+sum(grepl("log2transform",romics_object$steps))>0){log_transformed<-TRUE}else{log_transformed<-FALSE}

  if(mode=="vs"){
    #determine the list of combinations to consider
    by2combinations<- t(combn(levels_factor,2))
    if(reverse_order==TRUE){by2combinations<-by2combinations[,2:1]}

    #loop calculating pval,  fold changes
    for(i in 1:nrow(by2combinations)){
      for(j in 1:nrow(data)){
        if(log_transformed==TRUE){
          fold_change[j] <- mean(as.numeric(data[j,factor==by2combinations[i,2]]),na.rm = T)-mean(as.numeric(data[j,factor==by2combinations[i,1]]),na.rm = T)
        }else{
          fold_change[j] <- mean(as.numeric(data[j,factor==by2combinations[i,2]]),na.rm = T)/mean(as.numeric(data[j,factor==by2combinations[i,1]]),na.rm = T)
        }
        #calculate pvalues column
        if(log_transformed==TRUE & fold_change[j]==0 || log_transformed==FALSE && fold_change[j]==1){
          t_result[j] <- NA
          t_padj[j]<-NA
        }else{
          t_result[j] <- t.test(as.numeric(data[j,factor==by2combinations[i,1]]),as.numeric(data[j,factor==by2combinations[i,2]]),alternative=alternative, paired = paired, var.equal=var.equal)$p.value
        }
      }

      if(percentage_completeness>0){
        replicates_factor <- as.double(table(factor))
        names(replicates_factor) <- levels_factor
        replicates_factor <- replicates_factor[names(replicates_factor)%in%by2combinations[i,]]
        max_empty <- floor((replicates_factor)*(1-percentage_completeness/100))
        m1<-rowSums(missingdata[,factor==by2combinations[i,1]])>max_empty[names(max_empty)==by2combinations[i,1]]
        m2<-rowSums(missingdata[,factor==by2combinations[i,2]])>max_empty[names(max_empty)==by2combinations[i,2]]
        fold_change[m1*m2==1]<-NA
        t_result[m1*m2==1]<-NA
        }

      #add T.test p to T_table
      T_table[,paste(by2combinations[i,2],"_vs_",by2combinations[i,1],"_Ttest_p",sep="")]<- t_result

      ##add the adjusted p if padj=TRUE
      if(padj==TRUE){
        t_padj<-p.adjust(t_result, method=padj_method)
        T_table[,paste(by2combinations[i,2],"_vs_",by2combinations[i,1],"_Ttest_padj",sep="")]<- t_padj
      }

      #add (log(fold-change)) to the T_table
      if(log_transformed==TRUE){
        T_table[,paste("log(",by2combinations[i,2],"/",by2combinations[i,1],")",sep="")]<- fold_change
      }

      if(log_transformed==FALSE){
        T_table[,paste("(",by2combinations[i,2],"/",by2combinations[i,1],")",sep="")]<- fold_change
      }
    }
  }else{
    for(i in 1:length(levels_factor)){
      for(j in 1:nrow(data)){
        #calculate fold change (or log2(foldchange)) column
        if(log_transformed==TRUE){
          fold_change[j] <- mean(as.numeric(data[j,factor==levels_factor[i]]),na.rm = T)-mean(as.numeric(data[j,factor!=levels_factor[i]]),na.rm = T)
        }else{
          fold_change[j] <- mean(as.numeric(data[j,factor==levels_factor[i]]),na.rm = T)/mean(as.numeric(data[j,factor!=levels_factor[i]]),na.rm = T)
        }
        #calculate pvalues column
        if(log_transformed==TRUE & fold_change[j]==0 || log_transformed==FALSE && fold_change[j]==1){
          t_result[j] <- NA
          t_padj[j]<-NA
        }else{
          t_result[j] <- t.test(x = as.numeric(data[j,factor==levels_factor[i]]),y=as.numeric(data[j,factor!=levels_factor[i]]),alternative=alternative, paired = paired,...)$p.value
        }
      }

      if(percentage_completeness>0){
        replicates_factor <- as.double(table(factor))
        names(replicates_factor) <- levels_factor
        replicates_factor <- replicates_factor[names(replicates_factor)]
        other_factors <- sum(replicates_factor[names(replicates_factor)!=levels_factor[i]])
        replicates_factor<-c(replicates_factor[names(replicates_factor)==levels_factor[i]],other_factors=other_factors)
        max_empty <- floor((replicates_factor)*(1-percentage_completeness/100))
        m1<-rowSums(missingdata[,factor==levels_factor[i]])>max_empty[1]
        m2<-rowSums(missingdata[,factor!=levels_factor[i]])>max_empty[2]
        fold_change[m1*m2==1]<-NA
        t_result[m1*m2==1]<-NA
      }

      #add t.test p to T_table
      T_table[,paste0(levels_factor[i],"_vs_others_Ttest_p")]<- t_result

      ##add the adjusted p if padj=TRUE
      if(padj==TRUE){
        t_padj<-p.adjust(t_result, method=padj_method)
        T_table[,paste0(levels_factor[i],"_vs_others_Ttest_padj")]<- t_padj
      }

      #add (log(fold-change)) to the T_table
      if(log_transformed==TRUE){
        T_table[,paste0("log(",levels_factor[i],"/others)")]<- fold_change
      }

      if(log_transformed==FALSE){
        T_table[,paste0("(",levels_factor[i],"/others)")]<- fold_change
      }
    }}

  romics_object$statistics <- cbind(romics_object$statistics,T_table)

  #print info
  print("T_test columns were added to the statistics")

  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return romics_object
  return(romics_object)

}

#' romicsWilcoxTest()
#' @description Performs all possible paired Wilcoxon signed-rank test for each variable using the levels of the selected factor of the romics_object. The results are added as new columns in the statistics layers.
#' @param romics_object has to be an romics_object created with the function romicsCreateObject(),
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param paired a logical indicating whether you want a paired t-test.
#' @param pairing_factor name of a factor contained in an r_object, the list of the available factor can be obtained by using the function romics_factors()
#' @param padj a logical variable indincating wheter to perform or not adjustment of pvalues
#' @param padj_method correction method. Must be in  {"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr","none"}
#' @param factor a character string indicating the factor to use for the test, the list of the available factor can be obtained by using the function romics_factors(), if missing the function will use the main factor of the object
#' @param mode 'vs' or 'enrichment' indicating if all the between group comparisons should be performed (vs) or if the features in a given group should be compared to all the other groups (enrichment).
#' @param ... other parameters can be passed down to the wilox.test() function from the 'stat' package.
#' @details When paired Wilcox.test are performed it is possible to include a second factor to generate the pairs. This function will also calculate the fold-changes or log2(fold-change). Please, note that the test will automatically determine if a log tranformation was performed to the object, subsequently we recommend to import not pre-logged data.frames when creating the object. For paired T.tests, it is possible to set a second factor containing the pairs, if missing the function will consider the pairs based on the column order in the romics_object.
#' @return an romics_object with the statistical layer containing the newly generated t.tests and fold-changes
#' @author Geremy Clair
#' @export
romicsWilcoxTest<-function(romics_object, alternative="two.sided", paired = FALSE, pairing_factor="none", factor = "main", padj=TRUE, padj_method="BH",mode="vs",...){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){
    warning("The selected factor is not in the list of factors of the romics_object")
    warning(romicsFactorNames(romics_object))
  }
  if(missing(paired)){paired<-FALSE}
  if(missing(padj)){padj<-TRUE}
  if(missing(padj_method)){padj_method<-"BH"}
  if(!padj_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){stop("padj_method has to be in the following list: holm, hochberg, hommel, bonferroni, BH, BY,fdr, none")}
  if(missing(pairing_factor)){pairing_factor ="none"}
  if(missing(mode)){mode="vs"}
  if(!mode %in% c("vs", "enrichment")){stop("'mode' has to be either 'vs' or 'enrichment' to indicate if the comparison should be done between all possible groups or if the features within a given group should be compared to all the other groups (enrichment).")}
  #load data
  data<-romics_object$data

  #check if the $statistic object already is a part of the Romics_object (if not create it)
  if(is.null(romics_object$statistics)){
    print("The Statistics layer was added to your object")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if the statistics object does not have the same number of rows as the data replace it by a null statistics object
  if(!is.null(romics_object$statistics)&&nrow(romics_object$statistics)!=nrow(data)){
    warning("The Statistics layer was not containing the same number of rows as your data, it was replaced by an empty statistics layer")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if factor is main extract the factor from the romics_object$main_factor
  if(factor=="main"){factor<-romics_object$main_factor}
  #extract the factor from the metadata
  factor<-as.factor(as.character(t(romics_object$metadata[romicsFactorNames(romics_object)==factor,])))
  factor<-as.factor(t(factor))

  if(pairing_factor=="none"){
    pairing_factor<-romics_object$main_factor
  }else{
    if(!pairing_factor %in% rownames(romics_object$metadata)){stop("Your paired_factor is not present in the metadata")}}

  #extract paired_factor from the romics_object
  pairing_factor<-romics_object$metadata[pairing_factor==rownames(romics_object$metadata),]
  pairing_factor<-as.factor(t(pairing_factor))

  #order the data, the factor AND the paired_factor based on the factor and the paired_factor
  data<-data[,order(factor,pairing_factor)]

  #create an object containing both factor to sort the two at the time
  both_factor<-paste0(factor,"@",pairing_factor)
  both_factor<-both_factor[order(factor,pairing_factor)]

  #re-extract factors
  factor<-as.factor(sapply(strsplit(both_factor,"@"), function(x) x[1]))
  pairing_factor<-as.factor(sapply(strsplit(both_factor,"@"), function(x) x[2]))

  #extract the levels to be considered
  levels_factor<-levels(factor)

  #determine if the data was log transformed
  if(sum(grepl("log10transform",romics_object$steps))+sum(grepl("log2transform",romics_object$steps))>0){log_transformed<-TRUE}else{log_transformed<-FALSE}

  #create Wilcox_results and fold_change
  wilcox_result <- vector(mode="numeric",length=nrow(data))
  fold_change <- vector(mode="numeric",length=nrow(data))
  wilcox_padj<- vector(mode="numeric",length=nrow(data))

  #Create a T_table
  wilcox_table<- data.frame(matrix(nrow=nrow(data),ncol=0))

  if(mode=="vs"){
  #determine the list of combinations to consider
    by2combinations<- t(combn(levels_factor,2))

  #loop calculating pval,  fold changes
    for(i in 1:nrow(by2combinations)){
      for(j in 1:nrow(data)){
          #calculate fold change (or log2(foldchange)) column
          if(log_transformed==TRUE){
            fold_change[j] <- mean(as.numeric(data[j,factor==by2combinations[i,2]]))-mean(as.numeric(data[j,factor==by2combinations[i,1]]))
          }else{
            fold_change[j] <- mean(as.numeric(data[j,factor==by2combinations[i,2]]))/mean(as.numeric(data[j,factor==by2combinations[i,1]]))
          }
          #calculate pvalues column
          if(log_transformed==TRUE & fold_change[j]==0 || log_transformed==FALSE && fold_change[j]==1){
            t_result[j] <- NA
            t_padj[j]<-NA
          }else{
            wilcox_result[j] <- wilcox.test(x = as.numeric(data[j,factor==by2combinations[i,1]]),y= as.numeric(data[j,factor==by2combinations[i,2]]),alternative=alternative, paired = paired,...)$p.value
          }
        }

    #add Wilcox.test p to wilcox_table
    wilcox_table[,paste(by2combinations[i,2],"_vs_",by2combinations[i,1],"_Wilcox_test_p",sep="")]<- wilcox_result

    ##add the adjusted p if padj=TRUE
    if(padj==TRUE){
      wilcox_padj<-p.adjust(wilcox_result, method=padj_method)
      wilcox_table[,paste(by2combinations[i,2],"_vs_",by2combinations[i,1],"_Wilcox_test_padj",sep="")]<- wilcox_padj
      }

    #add (log(fold-change)) to the T_table
    if(log_transformed==TRUE){
      wilcox_table[,paste("log(",by2combinations[i,2],"/",by2combinations[i,1],")",sep="")]<- fold_change
      }

    if(log_transformed==FALSE){
      wilcox_table[,paste("(",by2combinations[i,2],"/",by2combinations[i,1],")",sep="")]<- fold_change
      }
    }
  }else{
    for(i in 1:length(levels_factor)){
        for(j in 1:nrow(data)){
          #calculate fold change (or log2(foldchange)) column
          if(log_transformed==TRUE){
            fold_change[j] <- mean(as.numeric(data[j,factor==levels_factor[i]]))-mean(as.numeric(data[j,factor!=levels_factor[i]]))
          }else{
            fold_change[j] <- mean(as.numeric(data[j,factor==levels_factor[i]]))/mean(as.numeric(data[j,factor!=levels_factor[i]]))
          }
          #calculate pvalues column
          if(log_transformed==TRUE & fold_change[j]==0 || log_transformed==FALSE && fold_change[j]==1){
            t_result[j] <- NA
            t_padj[j]<-NA
          }else{
            wilcox_result[j] <- wilcox.test(x = as.numeric(data[j,factor==levels_factor[i]]),y=as.numeric(data[j,factor!=levels_factor[i]]),alternative=alternative, paired = paired,...)$p.value
          }
        }

        #add Wilcox.test p to wilcox_table
        wilcox_table[,paste0(levels_factor[i],"_vs_others_Wilcox_test_p")]<- wilcox_result

        ##add the adjusted p if padj=TRUE
        if(padj==TRUE){
          wilcox_padj<-p.adjust(wilcox_result, method=padj_method)
          wilcox_table[,paste0(levels_factor[i],"_vs_others_Wilcox_test_padj")]<- wilcox_padj
          }

        #add (log(fold-change)) to the T_table
        if(log_transformed==TRUE){
          wilcox_table[,paste0("log(",levels_factor[i],"/others)")]<- fold_change
          }

        if(log_transformed==FALSE){
          wilcox_table[,paste0("(",levels_factor[i],"/others)")]<- fold_change
          }
      }}

  romics_object$statistics <- cbind(romics_object$statistics,wilcox_table)

  #print info
  print("wilcox_test columns were added to the statistics")

  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return romics_object
  return(romics_object)

}

#' romicsANOVA()
#' @description Performs the ANOVA for each variable contained in the data layer of the romics_object. The factor of the romics_object to be used for the analysis can be selected.
#' @param romics_object A romics_object created with the function romicsCreateObject(),
#' @param padj Boolean indincating wheter to perform or not adjustment of pvalues
#' @param padj_method correction method. Must be in  {"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr","none"}
#' @param factor A character string indicating the factor to use for the test, the list of the available factor can be obtained by using the function romics_factors(), if missing the function will use the main factor of the object.
#' @details perform an ANOVA for each variable of the romics_object the factor used will be the main_factor of the romics_object unless specified differently.
#' @return an romics_object with the statistical layer containing the newly generated ANOVA columns
#' @author Geremy Clair
#' @export
romicsANOVA<-function(romics_object, padj=TRUE, padj_method="BH", factor="main"){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){
    warning("The selected factor is not in the list of factors of the romics_object")
    warning(romicsFactorNames(romics_object))
  }
  if(missing(padj)){padj<-TRUE}
  if(missing(padj_method)){padj_method<-"BH"}
  if(!padj_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){stop("padj_method has to be in the following list: holm, hochberg, hommel, bonferroni, BH, BY,fdr, none")}

  #extract data from the romics_object
  data<-romics_object$data
  t_data<-t(data)

  #check if the $statistic object already is a part of the Romics_object (if not create it)
  if(is.null(romics_object$statistics)){
    print("The Statistics layer was added to your object")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if the statistics object does not have the same number of rows as the data replace it by a null statistics object
  if(!is.null(romics_object$statistics)&&nrow(romics_object$statistics)!=nrow(data)){
    warning("The Statistics layer was not containing the same number of rows as your data, it was replaced by an empty statistics layer")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if factor is main extract the factor from the romics_object$main_factor
  if(factor=="main"){factor<-romics_object$main_factor}

  #extract the factor from the metadata
  factor<-as.factor(as.character(t(romics_object$metadata[romicsFactorNames(romics_object)==factor,])))
  factor<-as.factor(t(as.character(factor)))

  #extrapolate the levels
  levels_factor<-levels(factor)

  #run the ANOVA
  ANOVA_results <- lapply(1:length(colnames(t_data)), function(i) summary(aov(t_data[,i] ~ factor)))

  # Put the names to the variables
  ANOVA_results <- setNames(ANOVA_results,rownames(data))
  ANOVA_results <- unlist(ANOVA_results)
  ANOVA_results <- ANOVA_results[grep(".Pr(>F)1",names(ANOVA_results), fixed = TRUE)]
  names(ANOVA_results) <- gsub(".Pr(>F)1","",names(ANOVA_results), fixed =TRUE)
  ANOVA_results<- data.frame(ANOVA_p=ANOVA_results)

  # if padj demanded calculate padj
  if(padj==TRUE){
    padj<-p.adjust(ANOVA_results$ANOVA_p, method=padj_method)
    ANOVA_results[,"ANOVA_padj"]<-padj
    print("The ANOVA columns (ANOVA_p and ANOVA_padj) were added to the statistics")
    }else{print("The ANOVA column (ANOVA_p) was added to the statistics")}

  romics_object$statistics <- cbind(romics_object$statistics,ANOVA_results)

  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return romics_object
  return(romics_object)

}

#' romicsZscores()
#' @description Calculates the Zscores for each cell of the data layer of the romics_object. the generated results are added to the statistics Layer as new columns.
#' @param romics_object A romics_object created with the function romicsCreateObject(),
#' @details adds the Zscores columns to the statistics Layer of a romics_object
#' @return an romics_object with the statistical layer containing the newly generated Zscores columns
#' @author Geremy Clair
#' @export
romicsZscores<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  data<-romics_object$data
  Means<-rowMeans(data,na.rm = T)
  Stdev<- apply(data,1, sd, na.rm = TRUE)


  Z_scores<-data
  colnames(Z_scores)<-paste0("Z_scores_",colnames(Z_scores))

  for (i in 1:nrow(data)){
    Z_scores[i,]<-(data[i,]-Means[i])/Stdev[i]
  }


  #check if the $statistic object already is a part of the Romics_object (if not create it)
  if(is.null(romics_object$statistics)){
    print("The Statistics layer was added to your object")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }


  #if the statistics object does not have the same number of rows as the data replace it by a null statistics object
  if(!is.null(romics_object$statistics)&&nrow(romics_object$statistics)!=nrow(data)){
    warning("The Statistics layer was not containing the same number of rows as your data, it was replaced by an empty statistics layer")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(data),ncol=0))
    rownames(romics_object$statistics)<-rownames(data)
  }

  #if Z_scores column previously exist remove them to replace by the new ones
  if(sum(grepl("Z_score_",colnames(romics_object$statistics)))>0){
    warning("The romics object contained previously some <Z_score_> columns those were removed and replace by the newly calculated ones")
    romics_object$statistics<-romics_object$statistics[,!grepl("Z_score_",colnames(romics_object$statistics))]
  }

  romics_object$statistics<-cbind(romics_object$statistics,Z_scores)

  print("Z_score_ columns were added to the statistics")

  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  return(romics_object)
  }

#' pFrequencyPlot()
#' @description makes a frequency plot of the pvalues and adjustedpvalues (columns of the statistics layer ending by '_p' and '_padj')
#' @param romics_object A romics_object created with the function romicsCreateObject()
#' @param p_columns 'all' by default, otherwise it can be a text vector containing the columns to be plotted.
#' @param p indicate the target pvalue to be plotted with a red dotted bar.
#' @param bin_width numeric vector, by default 0.01 indicate the width of the frequency bins
#' @details plot all or a specified list of pvalue and adjusted pvalues frequency plots
#' @return returns one or multiple plots
#' @author Geremy Clair
#' @export

pFrequencyPlot<-function(romics_object,p_columns="all",p=0.05,bin_width=0.01){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(p_columns)){p_columns="all"}
  if(!is.character(p_columns)){stop("p_columns has to be a character vector.")}
  if(missing(p)|!is.numeric(p)){p=0.05}

  pcol<-romics_object$statistics[, grepl(".*_p$",colnames(romics_object$statistics))|grepl(".*_padj$",colnames(romics_object$statistics))]

  if(p_columns!="all"){pcol<-pcol[,colnames(pcol) %in% p_columns]}

  for (i in 1:ncol(pcol)){
  p1<-p
  pval<-data.frame(ids=rownames(pcol), p=as.numeric(t(pcol[i])))
  print(paste0(sum(pval<p)," with ",colnames(pcol)[i],"<",p))
  print(ggplot(pval, aes(p)) +
    geom_histogram(binwidth = bin_width)+
    ggtitle(paste0("Frequency plot: ",colnames(pcol[i])))+geom_vline(xintercept=p,linetype="dashed", color = "red")+
    geom_text(aes(x=p1,y=max(hist(pval$p, seq(0,1,by=bin_width), plot = FALSE)$counts)/2), label=paste0("\n",colnames(pcol)[i],"<",p),colour="red",angle=90)+
    theme_ROP())
  }
}
