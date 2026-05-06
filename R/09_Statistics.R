#' romicsStatLayer()
#' @description This function will add a statistic layer if it does not exist, it will also verify if the existing stat layer has the right number of rows (same as the data), if not it will create a blank stat layer.
#' @param romics_object An romics_object containing a statistic layer
#' @details Adds the statistics layer.
#' @return This function returns a modified romics object, containing statistics layer
#' @author Geremy Clair
#' @export
romicsStatLayer<-function(romics_object=romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format.")}
  if(is.null(romics_object$statistics)){
    print("The Statistics layer was added to your object")
    romics_object$statistics<-data.frame(matrix(nrow=nrow(romics_object$data),ncol=0))
    rownames(romics_object$statistics)<-rownames(romics_object$data)
  }else{
    if(!is.null(romics_object$statistics) && nrow(romics_object$statistics)!=nrow(romics_object$data)){
      warning("The <romics_object$statistics> layer was not containing the same number of rows as your data, it was replaced by an empty <romics_object$statistics>.")
      romics_object$statistics<-data.frame(matrix(nrow=nrow(romics_object$data),ncol=0))
      rownames(romics_object$statistics)<-rownames(romics_object$data)
    }
  }
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
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(p_columns)){p_columns="all"}
  if(!is.character(p_columns)){stop("p_columns has to be a character vector.")}
  if(missing(p)|!is.numeric(p)){p=0.05}

  pcol<-romics_object$statistics[, grepl(".*_p$",colnames(romics_object$statistics))|grepl(".*_padj$",colnames(romics_object$statistics))]

  if(p_columns!="all"){pcol<-pcol[,colnames(pcol) %in% p_columns]}

  for (i in 1:ncol(pcol)){
    p1<-p
    pval<-data.frame(ids=rownames(pcol), p=as.numeric(t(pcol[i])))
    print(paste0(sum(pval<p,na.rm = T)," with ",colnames(pcol)[i],"<",p))
    print(ggplot(pval, aes(p)) +
            geom_histogram(binwidth = bin_width)+
            ggtitle(paste0("Frequency plot: ",colnames(pcol[i])))+geom_vline(xintercept=p,linetype="dashed", color = "red")+
            geom_text(aes(x=p1,y=max(hist(pval$p, seq(0,1,by=bin_width), plot = FALSE)$counts)/2), label=paste0("p=",p),colour="red",angle=90)+
            theme_ROP())
  }
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
#' @param ... other parameters can be passed down to the t.test() function from the 'stat' package.
#' @details When paired T.test are performed it is possible to include a second factor to generate the pairs. This function will also calculate the fold-changes or log2(fold-change). Please, note that the test will automatically determine if a log tranformation was performed to the object, subsequently we recommend to import not pre-logged data.frames when creating the object. For paired T.tests, it is possible to set a second factor containing the pairs, if missing the function will consider the pairs based on the column order in the romics_object.
#' @return an romics_object with the statistical layer containing the newly generated t.tests and fold-changes
#' @author Geremy Clair
#' @export
romicsTtest<-function(romics_object, alternative="two.sided", paired = FALSE, pairing_factor="none", var.equal=FALSE, factor = "main", padj=TRUE, padj_method="BH",mode="vs",percentage_completeness=0,reverse_order=FALSE, ...){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
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
  if(percentage_completeness<0 && percentage_completeness>100){stop("the completeness has to be comprised between 0 and 100 %")}

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
      # Pre-calculate which rows meet the percentage_completeness threshold
      if(percentage_completeness>0){
        replicates_factor <- as.double(table(factor))
        names(replicates_factor) <- levels_factor
        replicates_factor <- replicates_factor[names(replicates_factor)%in%by2combinations[i,]]
        max_empty <- floor((replicates_factor)*(1-percentage_completeness/100))
        m1<-rowSums(missingdata[,factor==by2combinations[i,1]])>max_empty[names(max_empty)==by2combinations[i,1]]
        m2<-rowSums(missingdata[,factor==by2combinations[i,2]])>max_empty[names(max_empty)==by2combinations[i,2]]
        meets_threshold <- !(m1 | m2)
      }else{
        meets_threshold <- rep(TRUE, nrow(data))
      }

      for(j in 1:nrow(data)){
        if(log_transformed==TRUE){
          fold_change[j] <- mean(as.numeric(data[j,factor==by2combinations[i,2]]),na.rm = T)-mean(as.numeric(data[j,factor==by2combinations[i,1]]),na.rm = T)
        }else{
          fold_change[j] <- mean(as.numeric(data[j,factor==by2combinations[i,2]]),na.rm = T)/mean(as.numeric(data[j,factor==by2combinations[i,1]]),na.rm = T)
        }
        #calculate pvalues column
        if(!meets_threshold[j] || !is.finite(fold_change[j]) || (log_transformed==TRUE && fold_change[j]==0) || (log_transformed==FALSE && fold_change[j]==1)){
          t_result[j] <- NA
          t_padj[j]<-NA
        }else{
          t_result[j] <- t.test(as.numeric(data[j,factor==by2combinations[i,1]]),as.numeric(data[j,factor==by2combinations[i,2]]),alternative=alternative, paired = paired, var.equal=var.equal)$p.value
        }
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
      # Pre-calculate which rows meet the percentage_completeness threshold
      if(percentage_completeness>0){
        replicates_factor <- as.double(table(factor))
        names(replicates_factor) <- levels_factor
        replicates_factor <- replicates_factor[names(replicates_factor)]
        other_factors <- sum(replicates_factor[names(replicates_factor)!=levels_factor[i]])
        replicates_factor<-c(replicates_factor[names(replicates_factor)==levels_factor[i]],other_factors=other_factors)
        max_empty <- floor((replicates_factor)*(1-percentage_completeness/100))
        m1<-rowSums(missingdata[,factor==levels_factor[i]])>max_empty[1]
        m2<-rowSums(missingdata[,factor!=levels_factor[i]])>max_empty[2]
        meets_threshold <- !(m1 | m2)
      }else{
        meets_threshold <- rep(TRUE, nrow(data))
      }

      for(j in 1:nrow(data)){
        #calculate fold change (or log2(foldchange)) column
        if(log_transformed==TRUE){
          fold_change[j] <- mean(as.numeric(data[j,factor==levels_factor[i]]),na.rm = T)-mean(as.numeric(data[j,factor!=levels_factor[i]]),na.rm = T)
        }else{
          fold_change[j] <- mean(as.numeric(data[j,factor==levels_factor[i]]),na.rm = T)/mean(as.numeric(data[j,factor!=levels_factor[i]]),na.rm = T)
        }
        #calculate pvalues column
        if(!meets_threshold[j] || !is.finite(fold_change[j]) || (log_transformed==TRUE && fold_change[j]==0) || (log_transformed==FALSE && fold_change[j]==1)){
          t_result[j] <- NA
          t_padj[j]<-NA
        }else{
          t_result[j] <- t.test(x = as.numeric(data[j,factor==levels_factor[i]]),y=as.numeric(data[j,factor!=levels_factor[i]]),alternative=alternative, paired = paired, var.equal=var.equal,...)$p.value
        }
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
#' @param ... other parameters can be passed down to the wilcox.test() function from the 'stat' package.
#' @details When paired Wilcox.test are performed it is possible to include a second factor to generate the pairs. This function will also calculate the fold-changes or log2(fold-change). Please, note that the test will automatically determine if a log tranformation was performed to the object, subsequently we recommend to import not pre-logged data.frames when creating the object. For paired T.tests, it is possible to set a second factor containing the pairs, if missing the function will consider the pairs based on the column order in the romics_object.
#' @return an romics_object with the statistical layer containing the newly generated Wilcox tests and fold-changes
#' @author Geremy Clair
#' @export
romicsWilcoxTest<-function(romics_object, alternative="two.sided", paired = FALSE, pairing_factor="none", factor = "main", padj=TRUE, padj_method="BH",mode="vs",...){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
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

  #Create a wilcox_table
  wilcox_table<- data.frame(matrix(nrow=nrow(data),ncol=0))

  if(mode=="vs"){
    #determine the list of combinations to consider
    by2combinations<- t(combn(levels_factor,2))

    #loop calculating pval,  fold changes
    for(i in 1:nrow(by2combinations)){
      # Pre-calculate which rows meet the percentage_completeness threshold
      if(percentage_completeness>0){
        replicates_factor <- as.double(table(factor))
        names(replicates_factor) <- levels_factor
        replicates_factor <- replicates_factor[names(replicates_factor)%in%by2combinations[i,]]
        max_empty <- floor((replicates_factor)*(1-percentage_completeness/100))
        m1<-rowSums(missingdata[,factor==by2combinations[i,1]])>max_empty[names(max_empty)==by2combinations[i,1]]
        m2<-rowSums(missingdata[,factor==by2combinations[i,2]])>max_empty[names(max_empty)==by2combinations[i,2]]
        meets_threshold <- !(m1 | m2)
      }else{
        meets_threshold <- rep(TRUE, nrow(data))
      }

      for(j in 1:nrow(data)){
        #calculate fold change (or log2(foldchange)) column
        if(log_transformed==TRUE){
          fold_change[j] <- mean(as.numeric(data[j,factor==by2combinations[i,2]]),na.rm=TRUE)-mean(as.numeric(data[j,factor==by2combinations[i,1]]),na.rm=TRUE)
        }else{
          fold_change[j] <- mean(as.numeric(data[j,factor==by2combinations[i,2]]),na.rm=TRUE)/mean(as.numeric(data[j,factor==by2combinations[i,1]]),na.rm=TRUE)
        }
        #calculate pvalues column
        if(!meets_threshold[j] || !is.finite(fold_change[j]) || (log_transformed==TRUE && fold_change[j]==0) || (log_transformed==FALSE && fold_change[j]==1)){
          wilcox_result[j] <- NA
          wilcox_padj[j]<-NA
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

      #add (log(fold-change)) to the wilcox_table
      if(log_transformed==TRUE){
        wilcox_table[,paste("log(",by2combinations[i,2],"/",by2combinations[i,1],")",sep="")]<- fold_change
      }

      if(log_transformed==FALSE){
        wilcox_table[,paste("(",by2combinations[i,2],"/",by2combinations[i,1],")",sep="")]<- fold_change
      }
    }
  }else{
    for(i in 1:length(levels_factor)){
      # Pre-calculate which rows meet the percentage_completeness threshold
      if(percentage_completeness>0){
        replicates_factor <- as.double(table(factor))
        names(replicates_factor) <- levels_factor
        replicates_factor <- replicates_factor[names(replicates_factor)]
        other_factors <- sum(replicates_factor[names(replicates_factor)!=levels_factor[i]])
        replicates_factor<-c(replicates_factor[names(replicates_factor)==levels_factor[i]],other_factors=other_factors)
        max_empty <- floor((replicates_factor)*(1-percentage_completeness/100))
        m1<-rowSums(missingdata[,factor==levels_factor[i]])>max_empty[1]
        m2<-rowSums(missingdata[,factor!=levels_factor[i]])>max_empty[2]
        meets_threshold <- !(m1 | m2)
      }else{
        meets_threshold <- rep(TRUE, nrow(data))
      }

      for(j in 1:nrow(data)){
        #calculate fold change (or log2(foldchange)) column
        if(log_transformed==TRUE){
          fold_change[j] <- mean(as.numeric(data[j,factor==levels_factor[i]]),na.rm=TRUE)-mean(as.numeric(data[j,factor!=levels_factor[i]]),na.rm=TRUE)
        }else{
          fold_change[j] <- mean(as.numeric(data[j,factor==levels_factor[i]]),na.rm=TRUE)/mean(as.numeric(data[j,factor!=levels_factor[i]]),na.rm=TRUE)
        }
        #calculate pvalues column
        if(!meets_threshold[j] || !is.finite(fold_change[j]) || (log_transformed==TRUE && fold_change[j]==0) || (log_transformed==FALSE && fold_change[j]==1)){
          wilcox_result[j] <- NA
          wilcox_padj[j]<-NA
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

      #add (log(fold-change)) to the wilcox_table
      if(log_transformed==TRUE){
        wilcox_table[,paste0("log(",levels_factor[i],"/others)")]<- fold_change
      }

      if(log_transformed==FALSE){
        wilcox_table[,paste0("(",levels_factor[i],"/others)")]<- fold_change
      }
    }}

  romics_object$statistics <- cbind(romics_object$statistics,wilcox_table)

  #print info
  print("Wilcox_test columns were added to the statistics")

  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return romics_object
  return(romics_object)

}


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
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format.")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){stop("The selected <factor> is not in the list of factors of the <romics_object>.")}
  #statistics layer checks
  romics_object<-romicsStatLayer(romics_object)
  #load data
  data<-romics_object$data
  #extract the factor from the metadata
  factor<-romicsExtractFactor(romics_object,factor = factor)
  #calculate the means for each factor
  means <- sapply(levels(factor), function(level) {
    rowMeans(data[, factor == level, drop = FALSE], na.rm = TRUE)
  })
  #add "_means" to the columns names
  colnames(means)<-paste0(colnames(means),"_mean")
  #convert the NaNs into NAs
  means[is.nan(means)]<-NA
  #add the new columns to the statistics layer
  romics_object$statistics <- cbind(romics_object$statistics,means)
  #print info
  print("The 'means' columns were added to the <romics_object$statistics> layer.")
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
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format.")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){stop("The selected <factor> is not in the list of factors of the <romics_object>.")}
  #statistics layer checks
  romics_object<-romicsStatLayer(romics_object)
  #load data
  data<-romics_object$data
  #extract the factor from the metadata
  factor<-romicsExtractFactor(romics_object,factor = factor)
  #calculate the sd
  sd <- sapply(levels(factor), function(level) {
    apply(data[, factor == level, drop = FALSE], 1, function(row) sd(as.numeric(row), na.rm = TRUE))
  })
  #add "_sd" to the columns names
  colnames(sd)<-paste0(colnames(sd),"_sd")
  #add the new columns to the statistics layer
  romics_object$statistics <- cbind(romics_object$statistics,sd)
  #print info
  print("The 'sd' columns were added to the <romics_object$statistics> layer.")
  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  #return the new object
  return(romics_object)
}

#' romicsPercentComplete()
#' @description Calculates the percentage of completeness for each level of the selected factor and add the generated columns in the statistics layer of the romics_object.
#' @param romics_object A object created using the function romicsCreateObject().
#' @param main_factor Either 'main' OR any factor from the romics_object, the list of factors from a romics object can be obtained using the function romicsFactorNames().
#' @details Adds the percentage_completeness columns to the statistics layer. Enable to choose a different factor than the main one to do those calculation.
#' @return This function returns a modified romics object, containing percentage_completeness columns in the statistics layer.
#' @author Geremy Clair
#' @export
romicsPercentComplete<-function(romics_object, factor="main"){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format.")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){stop("The selected <factor> is not in the list of factors of the <romics_object>.")}
  #statistics layer checks
  romics_object<-romicsStatLayer(romics_object)
  #load data
  data<- !romics_object$missingdata
  #extract the factor from the metadata
  factor<-romicsExtractFactor(romics_object,factor = factor)

  #calculate the means for each factor
  percentage_completeness <- sapply(levels(factor), function(level) {
    subset_data<- data[, factor == level, drop = FALSE]
    completeness<-rowSums(subset_data, na.rm = TRUE)/ncol(subset_data)*100
    return(completeness)
  })
  #add "_percentage_completeness" to the columns names
  colnames(percentage_completeness)<-paste0(colnames(percentage_completeness),"_percentage_completeness")

  #convert the NaNs into NAs
  percentage_completeness[is.nan(percentage_completeness)]<-NA
  #add the new columns to the statistics layer
  romics_object$statistics <- cbind(romics_object$statistics,percentage_completeness)
  #print info
  print("The 'percentage_completeness' columns were added to the <romics_object$statistics> layer.")
  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  class(romics_object)<-"romics_object"
  #return the new object
  return(romics_object)
}

#' romicsPercentMissing()
#' @description Calculates the percentage of Missingness for each level of the selected factor and add the generated columns in the statistics layer of the romics_object.
#' @param romics_object A object created using the function romicsCreateObject().
#' @param main_factor Either 'main' OR any factor from the romics_object, the list of factors from a romics object can be obtained using the function romicsFactorNames().
#' @details Adds the percentage_missingness columns to the statistics layer. Enable to choose a different factor than the main one to do those calculation.
#' @return This function returns a modified romics object, containing percentage_missingness columns in the statistics layer.
#' @author Geremy Clair
#' @export
romicsPercentMissing<-function(romics_object, factor="main"){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format.")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){stop("The selected <factor> is not in the list of factors of the <romics_object>.")}
  #statistics layer checks
  romics_object<-romicsStatLayer(romics_object)
  #load data
  data<- romics_object$missingdata
    #extract the factor from the metadata
    factor<-romicsExtractFactor(romics_object,factor = factor)

  #calculate the means for each factor
  percentage_missingness <- sapply(levels(factor), function(level) {
    subset_data<- data[, factor == level, drop = FALSE]
    missingness<-rowSums(subset_data, na.rm = TRUE)/ncol(subset_data)*100
    return(missingness)
  })

  #add "_percentage_missingness" to the columns names
  colnames(percentage_missingness)<-paste0(colnames(percentage_missingness),"_percentage_missingness")
  #convert the NaNs into NAs
  percentage_missingness[is.nan(percentage_missingness)]<-NA
  #add the new columns to the statistics layer
  romics_object$statistics <- cbind(romics_object$statistics,percentage_missingness)
  #print info
  print("The 'percentage_missingness' columns were added to the <romics_object$statistics> layer.")
  #update steps
  class(romics_object)<-"romics_object"
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  #return the new object
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
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format.")}
  #statistics layer checks
  romics_object<-romicsStatLayer(romics_object)
  #load data
  data<-romics_object$data
  #calculate the means and sd
  means <- rowMeans(data, na.rm = TRUE)
  sd <- apply(data, 1, sd, na.rm = TRUE)
  #if Z_scores column previously exist remove them to replace by the new ones
  if(sum(grepl("Z_scores_",colnames(romics_object$statistics)))>0){
    warning("The romics object contained previously some 'Z_scores_' columns those were removed and replace by the newly calculated ones.")
    romics_object$statistics<-romics_object$statistics[,!grepl("Z_scores_",colnames(romics_object$statistics))]
  }
  #calculate the Zscores
  Z <- (data - means) / sd
  #modify the column names to begin by "Z_scores_"
  colnames(Z)<-paste0("Z_scores_",colnames(Z))
  romics_object$statistics<-cbind(romics_object$statistics,Z)
  #message
  print("'Z_scores_' columns were added to the <romics_object$statistics> layer.")
  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  #return the object
  return(romics_object)
}

#' romicsANOVA()
#' @description Performs the ANOVA for each variable contained in the data layer of the romics_object. The factor of the romics_object to be used for the analysis can be selected.
#' @param romics_object A romics_object created with the function romicsCreateObject(),
#' @param padj Boolean indincating wheter to perform or not adjustment of pvalues
#' @param padj_method correction method. Must be in  c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
#' @param factor A character string indicating the factor to use for the test, the list of the available factor can be obtained by using the function romics_factors(), if missing the function will use the main factor of the object.
#' @details perform an ANOVA for each variable of the romics_object the factor used will be the main_factor of the romics_object unless specified differently.
#' @return an romics_object with the statistical layer containing the newly generated ANOVA columns
#' @author Geremy Clair
#' @export
romicsANOVA<-function(romics_object, padj=TRUE, padj_method="BH", factor="main"){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format.")}
  if(missing(factor)){factor="main"}
  if(!factor %in% c("main",romicsFactorNames(romics_object))){
    warning("The selected <factor> is not in the list of factors of the <romics_object>.")
    warning(romicsFactorNames(romics_object))
  }
  if(missing(padj)){padj<-TRUE}
  if(missing(padj_method)){padj_method<-"BH"}
  if(!padj_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr")){stop("<padj_method> has to be in the following list: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY',
                                                                                                     fdr'.")}
  #extract data from the romics_object
  data<-data.frame(t(romics_object$data))
  #if factor is main extract the factor from the romics_object$main_factor
  if(factor=="main"){factor<-romics_object$main_factor}
  factor_name<-factor
  #extract the factor from the metadata
  factor<-romicsExtractFactor(romics_object,factor)
  #create data_frame for tests
  df<-cbind(factor=factor,data)
  #test function
  anova_feature <- function(feature, factor) {
    data_frame <- data.frame(value = feature, factor = factor)
    data_frame <- data_frame[!is.na(data_frame$value),]
    if(length(unique(data_frame$factor))<2){p<-NA}else{
      model <- aov(value ~ factor, data = data_frame)
      p <- summary(model)[[1]]$`Pr(>F)`[1]}
    return(p)
  }
  #run the tests
  ANOVA_results <-data.frame(p=apply(df[2:ncol(df)], 2, anova_feature, factor = df$factor))
  #adjust_the_colnames
  colnames(ANOVA_results)<-paste0("ANOVA_",factor_name,"_p")
  #if padj demanded calculate padj
  if(padj==TRUE){
    ANOVA_results[,2]<-p.adjust(ANOVA_results[,1], method=padj_method)
    colnames(ANOVA_results)[2]<-paste0(colnames(ANOVA_results)[1],"adj")
  }
  #message
  print("The 'ANOVA columns' were added to the <romics_object$statistics> layer:")
  print(colnames(ANOVA_results))
  #add to the statistics layer
  romics_object$statistics <- cbind(romics_object$statistics,ANOVA_results)
  #update steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  #return romics_object
  return(romics_object)
}

#' romicsGlmBinomial()
#' @description Performs generalized linear model (GLM) binomial tests for statistical significance testing with support for pairwise and enrichment modes
#' @param romics_object A romics_object with data and metadata layers
#' @param factor Character string specifying which factor to use for grouping. Default: "main" uses the main factor
#' @param cluster_factor Character string specifying an optional clustering factor for subset analysis. Default: "none"
#' @param padj Logical indicating whether to calculate adjusted p-values. Default: TRUE
#' @param padj_method Character string specifying p-value adjustment method ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"). Default: "BH"
#' @param mode Character string either "vs" for pairwise comparisons or "enrichment" for each group vs all others. Default: "vs"
#' @param reverse_order Logical indicating whether to reverse comparison order. Default: FALSE
#' @param suppress_warnings Logical indicating whether to suppress warning messages. Default: TRUE
#' @details GLM binomial test compares proportions of significant features between groups. Results include p-values, adjusted p-values, and directionality (1=up, -1=down, 0=no change)
#' @return romics_object with statistical test results added to the $statistics layer
#' @author Geremy Clair
#' @export
romicsGlmBinomial <- function(romics_object,
                              factor = "main",
                              cluster_factor = "none",
                              padj = TRUE,
                              padj_method = "BH",
                              mode = "vs",
                              reverse_order = FALSE,
                              suppress_warnings = TRUE) {
  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("<romics_object> is missing or is not in the appropriate format.")
  }
  if (missing(factor)) {
    factor <- "main"
  }
  if (!factor %in% c("main", romicsFactorNames(romics_object))) {
    stop("The selected <factor> is not in the list of factors of the <romics_object>.")
  }
  if (missing(padj)) {
    padj <- TRUE
  }
  if (missing(padj_method)) {
    padj_method <- "BH"
  }
  if (!padj_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
    stop("The <padj_method> has to be in the following list: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'.")
  }
  if (missing(mode)) {
    mode <- "vs"
  }
  if (!mode %in% c("vs", "enrichment")) {
    stop("The <mode> has to be either 'vs' or 'enrichment' to indicate if the comparison should be done between all possible groups or if the features within a given group should be compared to all the other groups (enrichment).")
  }
  if (missing(reverse_order)) {
    reverse_order <- FALSE
  }
  if (missing(suppress_warnings)) {
    suppress_warnings <- TRUE
  }
  if (missing(cluster_factor)) {
    cluster_factor <- "none"
  }

  # Validate cluster_factor
  if (cluster_factor != "none" && !cluster_factor %in% romicsFactorNames(romics_object)) {
    stop(paste("Cluster factor", cluster_factor, "not found in available factors."))
  }

  # Check if the $statistics object already is a part of the romics_object (if not create it)
  romics_object <- romicsStatLayer(romics_object)

  # Load data
  data <- romics_object$data
  metadata <- data.frame(t(romics_object$metadata), stringsAsFactors = TRUE)

  # Create missing data matrix if it doesn't exist
  if (is.null(romics_object$missingdata)) {
    romics_object$missingdata <- is.na(romics_object$data)
    romics_object$missingdata <- romics_object$missingdata * 1  # Convert to 0/1
  }
  missingdata <- romics_object$missingdata

  # If factor is main, extract the factor from the romics_object$main_factor
  if (factor == "main") {
    factor <- romics_object$main_factor
  }

  # Extract the factor from the metadata
  factor_values <- romicsExtractFactor(romics_object, factor)

  # Determine cluster levels
  if (cluster_factor != "none") {
    cluster_levels <- levels(as.factor(metadata[[cluster_factor]]))
    message(paste("GLM Binomial test will be performed for", length(cluster_levels),
                  "cluster levels:", paste(cluster_levels, collapse = ", ")))
  } else {
    cluster_levels <- "all"
  }

  # Function to calculate the GLM test results
  GLM_feature <- function(row, f) {
    row <- as.numeric(row)
    g <- data.frame(pop = as.character(f), missing = row)
    tryCatch({
      if (suppress_warnings) {
        model_result <- suppressWarnings(anova(glm(missing ~ pop, family = "binomial", data = g), test = "Chisq"))
      } else {
        model_result <- anova(glm(missing ~ pop, family = "binomial", data = g), test = "Chisq")
      }
      p <- model_result[2, "Pr(>Chi)"]
      return(p)
    }, error = function(e) {
      return(NA)
    })
  }

  # Initialize results data frame
  all_results <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
  rownames(all_results) <- rownames(data)

  # Process each cluster level
  for (cluster_level in cluster_levels) {
    if (cluster_factor != "none") {
      message(paste("Processing cluster level:", cluster_level))
      # Subset data and metadata for this cluster
      cluster_indices <- metadata[[cluster_factor]] == cluster_level
      subset_data <- missingdata[, cluster_indices, drop = FALSE]
      subset_metadata <- metadata[cluster_indices, , drop = FALSE]

      # Extract factor from subset metadata
      subset_factor <- subset_metadata[[factor]]
      if (!is.factor(subset_factor)) {
        subset_factor <- as.factor(subset_factor)
      }

      cluster_suffix <- paste0("_within_", cluster_level)
    } else {
      subset_data <- missingdata
      subset_factor <- factor_values
      cluster_suffix <- ""
    }

    # Convert factor to factor type if needed
    if (!is.factor(subset_factor)) {
      subset_factor <- as.factor(subset_factor)
    }

    # Check if we have enough samples
    if (ncol(subset_data) < 3) {
      warning(paste("Insufficient samples in cluster", cluster_level, ". Skipping."))
      next
    }

    # Determine the list of combinations to consider
    if (mode == "vs") {
      by2combinations <- data.frame(t(combn(levels(subset_factor), 2)))
      if (reverse_order == TRUE) {
        by2combinations <- by2combinations[, 2:1]
      }
    } else {
      by2combinations <- data.frame(X1 = levels(subset_factor), X2 = "other")
      if (reverse_order == TRUE) {
        by2combinations <- by2combinations[, 2:1]
      }
    }

    # Create the output table for this cluster
    cluster_results <- data.frame(matrix(ncol = 0, nrow = nrow(subset_data)))
    rownames(cluster_results) <- rownames(subset_data)

    # Calculate the stats for each factor combination
    for (i in 1:nrow(by2combinations)) {
      # Extract data for first group (now by2combinations[i,2] since we renamed to i,2 vs i,1)
      if (as.character(by2combinations[i, 2]) != "other") {
        df1 <- as.data.frame(subset_data[, as.character(subset_factor) == as.character(by2combinations[i, 2])])
        f1 <- as.character(subset_factor)[subset_factor == as.character(by2combinations[i, 2])]
      } else {
        df1 <- subset_data[, as.character(subset_factor) != as.character(by2combinations[i, 1])]
        f1 <- rep("other", ncol(df1))
      }

      # Extract data for second group (now by2combinations[i,1])
      if (as.character(by2combinations[i, 1]) != "other") {
        df2 <- as.data.frame(subset_data[, as.character(subset_factor) == as.character(by2combinations[i, 1])])
        f2 <- as.character(subset_factor)[subset_factor == as.character(by2combinations[i, 1])]
      } else {
        df2 <- subset_data[, as.character(subset_factor) != as.character(by2combinations[i, 2])]
        f2 <- rep("other", ncol(df2))
      }

      # Calculate directionality
      # 1 = more present in second group (i,1), -1 = more present in first group (i,2), 0 = equal
      dir <- rep(0, nrow(subset_data))
      dir[rowSums(df2) / ncol(df2) > rowSums(df1) / ncol(df1)] <- 1
      dir[rowSums(df1) / ncol(df1) > rowSums(df2) / ncol(df2)] <- -1

      # Check if we have enough samples
      if (ncol(df1) < 2 || ncol(df2) < 2) {
        p <- rep(NA, nrow(subset_data))
        if (ncol(df1) < 2) {
          message(paste0("The level <", as.character(by2combinations[i, 2]),
                         "> has less than 2 samples in cluster ", cluster_level,
                         ", the GLM Binomial test will not be calculated for the comparison: ",
                         as.character(by2combinations[i, 2]), "_vs_",
                         as.character(by2combinations[i, 1])))
        } else {
          message(paste0("The level <", as.character(by2combinations[i, 1]),
                         "> has less than 2 samples in cluster ", cluster_level,
                         ", the GLM Binomial test will not be calculated for the comparison: ",
                         as.character(by2combinations[i, 2]), "_vs_",
                         as.character(by2combinations[i, 1])))
        }
      } else {
        # Combine data and factor
        f <- as.factor(c(f1, f2))
        df <- cbind(df1, df2)
        # Calculate p-values
        p <- as.numeric(apply(df, 1, GLM_feature, f = f))
      }

      # Apply p-value adjustment if requested
      if (padj == TRUE) {
        adj <- p.adjust(p, method = padj_method)
        r <- data.frame(p = p, adj = adj, dir = dir)
      } else {
        r <- data.frame(p = p, dir = dir)
      }

      # Create column names with cluster suffix
      # In enrichment mode, use the non-"other" level (level being tested)
      if (mode == "enrichment" && as.character(by2combinations[i, 2]) == "other") {
        level_name <- as.character(by2combinations[i, 1])
      } else {
        level_name <- as.character(by2combinations[i, 2])
      }

      colnames(r)[colnames(r) == "p"] <- paste0(
        level_name, "_vs_others",
        cluster_suffix, "_glmBinomialTest_p"
      )
      colnames(r)[colnames(r) == "adj"] <- paste0(
        level_name, "_vs_others",
        cluster_suffix, "_glmBinomialTest_padj"
      )
      colnames(r)[colnames(r) == "dir"] <- paste0(
        level_name, "_vs_others",
        cluster_suffix, "_directionality"
      )

      cluster_results <- cbind(cluster_results, r)
    }

    # Add cluster results to overall results
    all_results <- cbind(all_results, cluster_results)
  }

  # Add results to statistics layer
  if (ncol(all_results) > 0) {
    romics_object$statistics <- cbind(romics_object$statistics, all_results)
  }

  # Message
  n_comparisons <- ncol(all_results)
  if (n_comparisons > 0) {
    if (cluster_factor != "none") {
      message(paste("GLM Binomial test complete.", n_comparisons,
                    "statistical columns were added to the statistics layer for",
                    length(cluster_levels), "cluster(s)."))
    } else {
      message(paste("GLM Binomial test complete.", n_comparisons,
                    "statistical columns were added to the statistics layer."))
    }
  } else {
    warning("No results were generated. Check your data and factor specifications.")
  }

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  # Return romics_object
  return(romics_object)
}

#' romicsLMM()
#' @description Performs linear mixed model analysis for each variable using specified fixed and random effects. The results are added as new columns in the statistics layer.
#' @param romics_object An romics_object created with the function romicsCreateObject().
#' @param fixed_effects A character string specifying the fixed effects formula (e.g., "condition" or "condition + treatment").
#' @param random_effects A character string specifying the random effects formula (e.g., "(1|donor)" or "(1|donor) + (1|batch)").
#' @param cluster_factor A character string indicating a clustering factor to subset data by (e.g., "Leiden_clust_names"). If "none", analysis is performed on the full dataset. Default: "none".
#' @param method A character string specifying the method for p-value calculation: "satterthwaite" (using lmerTest) or "likelihood_ratio" (comparing nested models). Default: "satterthwaite".
#' @param padj A logical variable indicating whether to perform p-value adjustment. Default: TRUE.
#' @param padj_method Correction method. Must be in c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"). Default: "BH".
#' @param percentage_completeness A numerical value between 0 and 100 indicating the minimum completeness required to calculate the LMM (if not met, p-value will be NA). Default: 0.
#' @param REML A logical indicating whether to use REML estimation. Default: TRUE for satterthwaite method, FALSE for likelihood_ratio method.
#' @param handle_singular A character string specifying how to handle singular fits: "keep" (keep results with warning), "simplify" (try simpler random effects), "drop" (return NA). Default: "simplify".
#' @param suppress_warnings A logical indicating whether to suppress convergence warnings. Default: TRUE.
#' @param pairwise_comparisons A logical indicating whether to perform all pairwise comparisons for factors with 3+ levels. Default: TRUE.
#' @details This function fits linear mixed models to test fixed effects while accounting for random effects structure. When pairwise_comparisons is TRUE and a factor has 3 or more levels, all pairwise comparisons are computed. When models fail to converge or produce singular fits, the function can automatically simplify the random effects structure or fall back to fixed effects only models. When cluster_factor is specified, separate models are fitted for each cluster level.
#' @return An romics_object with the statistics layer containing the newly generated LMM results including effect sizes (log fold changes or differences), p-values, and adjusted p-values.
#' @author Geremy Clair
#' @export
romicsLMM <- function(romics_object,
                      fixed_effects,
                      random_effects,
                      cluster_factor = "none",
                      method = "satterthwaite",
                      padj = TRUE,
                      padj_method = "BH",
                      percentage_completeness = 0,
                      REML = NULL,
                      handle_singular = "simplify",
                      suppress_warnings = TRUE,
                      pairwise_comparisons = TRUE) {
  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("The <romics_object> is missing or is not in the appropriate format.")
  }
  if (missing(fixed_effects) || !is.character(fixed_effects)) {
    stop("<fixed_effects> must be a character string specifying the fixed effects formula.")
  }
  if (missing(random_effects) || !is.character(random_effects)) {
    stop("<random_effects> must be a character string specifying the random effects formula.")
  }
  if (!method %in% c("satterthwaite", "likelihood_ratio")) {
    stop("<method> must be either 'satterthwaite' or 'likelihood_ratio'.")
  }
  if (!padj_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
    stop("The <padj_method> has to be in the following list: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'.")
  }
  if (percentage_completeness < 0 || percentage_completeness > 100) {
    stop("The <percentage_completeness> has to be between 0 and 100%.")
  }
  if (!handle_singular %in% c("keep", "simplify", "drop")) {
    stop("<handle_singular> must be 'keep', 'simplify', or 'drop'.")
  }

  # Set default REML based on method
  if (is.null(REML)) {
    REML <- ifelse(method == "satterthwaite", TRUE, FALSE)
  }

  # Check required packages
  required_packages <- c("lme4", "lmerTest", "emmeans")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required. Please install it with: install.packages('", pkg, "')", sep = ""))
    }
  }

  # Configure emmeans to suppress degree of freedom adjustment warnings
  if (suppress_warnings) {
    emmeans::emm_options(
      pbkrtest.limit = Inf,
      lmerTest.limit = Inf,
      disable.pbkrtest = TRUE
    )
  }

  # Check statistics layer
  romics_object <- romicsStatLayer(romics_object)

  # Load data and metadata
  data <- romics_object$data
  metadata <- data.frame(t(romics_object$metadata), stringsAsFactors = TRUE)

  # Check if data is log-transformed
  is_log_transformed <- romicsLogCheck(romics_object)

  # Extract and validate factors from fixed and random effects
  fixed_factors <- trimws(unlist(strsplit(gsub("\\*|\\:", "+", fixed_effects), "\\+")))

  # Better parsing for random effects that handles multiple formats
  random_effects_clean <- random_effects
  random_effects_clean <- gsub("\\(1\\|", "", random_effects_clean)
  random_effects_clean <- gsub("[\\(\\)]", "", random_effects_clean)
  random_factors <- unlist(strsplit(random_effects_clean, "[+/]"))
  random_factors <- trimws(random_factors)
  random_factors <- random_factors[random_factors != ""]

  all_factors <- colnames(metadata)

  if (!all(fixed_factors %in% all_factors)) {
    missing_fixed <- fixed_factors[!fixed_factors %in% all_factors]
    stop(paste("Fixed effects not found in metadata:", paste(missing_fixed, collapse = ", ")))
  }
  if (!all(random_factors %in% all_factors)) {
    missing_random <- random_factors[!random_factors %in% all_factors]
    stop(paste("Random effects not found in metadata:", paste(missing_random, collapse = ", ")))
  }

  # Validate cluster factor
  if (cluster_factor != "none" && !cluster_factor %in% all_factors) {
    stop(paste("Cluster factor", cluster_factor, "not found in metadata."))
  }

  # Determine cluster levels
  if (cluster_factor != "none") {
    cluster_levels <- levels(as.factor(metadata[[cluster_factor]]))
    message(paste("Analysis will be performed for", length(cluster_levels), "cluster levels:", paste(cluster_levels, collapse = ", ")))
  } else {
    cluster_levels <- "all"
  }

  # Determine if we should do pairwise comparisons
  do_pairwise <- FALSE
  pairwise_contrasts <- NULL

  if (pairwise_comparisons && length(fixed_factors) == 1) {
    factor_levels <- levels(as.factor(metadata[[fixed_factors[1]]]))
    if (length(factor_levels) >= 3) {
      do_pairwise <- TRUE
      # Generate all pairwise combinations
      pairwise_contrasts <- combn(factor_levels, 2, simplify = FALSE)
      message(paste("Pairwise comparisons enabled for factor", fixed_factors[1],
                    "with", length(factor_levels), "levels."))
      message(paste("Will perform", length(pairwise_contrasts), "pairwise comparisons."))
    }
  }

  # Create fallback random effects formulas
  fallback_formulas <- list()
  if (length(random_factors) > 1 && handle_singular == "simplify") {
    for (rf in random_factors) {
      fallback_formulas[[length(fallback_formulas) + 1]] <- paste0("(1|", rf, ")")
    }
  }

  # Function to fit LMM for a single feature with fallback
  fit_lmm_feature <- function(response_data, metadata_subset, feature_name) {
    # Check completeness
    completeness <- sum(!is.na(response_data)) / length(response_data) * 100
    if (completeness < percentage_completeness) {
      return(list(estimates = NA, p_values = NA, contrasts = NULL, model_type = "insufficient_data"))
    }

    # Prepare model data
    model_data <- data.frame(response = as.numeric(response_data))
    model_data <- cbind(model_data, metadata_subset)
    model_data <- model_data[!is.na(model_data$response), ]

    # Check if we have enough data
    if (nrow(model_data) < 3) {
      return(list(estimates = NA, p_values = NA, contrasts = NULL, model_type = "insufficient_data"))
    }

    # Check if random effects have enough levels
    random_effect_valid <- TRUE
    for (rf in random_factors) {
      if (length(unique(model_data[[rf]])) < 2) {
        random_effect_valid <- FALSE
        break
      }
    }

    # Try full model first, then fallbacks, then fixed effects only
    formulas_to_try <- list()
    model_types <- list()

    if (random_effect_valid) {
      formulas_to_try[[1]] <- random_effects
      model_types[[1]] <- "full_mixed"

      if (length(fallback_formulas) > 0 && handle_singular == "simplify") {
        for (i in seq_along(fallback_formulas)) {
          formulas_to_try[[length(formulas_to_try) + 1]] <- fallback_formulas[[i]]
          model_types[[length(model_types) + 1]] <- "simplified_mixed"
        }
      }
    }

    if (handle_singular == "simplify") {
      formulas_to_try[[length(formulas_to_try) + 1]] <- NULL
      model_types[[length(model_types) + 1]] <- "fixed_only"
    }

    for (idx in seq_along(formulas_to_try)) {
      re_formula <- formulas_to_try[[idx]]
      model_type <- model_types[[idx]]

      tryCatch({
        if (method == "satterthwaite") {
          # Decide whether to use mixed or fixed model
          if (is.null(re_formula)) {
            formula_full <- as.formula(paste("response ~", fixed_effects))
            if (suppress_warnings) {
              model <- suppressMessages(suppressWarnings(lm(formula_full, data = model_data)))
            } else {
              model <- lm(formula_full, data = model_data)
            }
            is_singular <- FALSE
          } else {
            formula_full <- as.formula(paste("response ~", fixed_effects, "+", re_formula))
            if (suppress_warnings) {
              model <- suppressMessages(suppressWarnings(lmerTest::lmer(formula_full, data = model_data, REML = REML)))
            } else {
              model <- lmerTest::lmer(formula_full, data = model_data, REML = REML)
            }
            is_singular <- lme4::isSingular(model)
          }

          # Handle singular fits
          if (is_singular) {
            if (handle_singular == "drop") {
              next
            }
            if (handle_singular == "simplify" && idx < length(formulas_to_try)) {
              next
            }
          }

          # Extract standard coefficients
          coeffs <- summary(model)$coefficients

          if (nrow(coeffs) < 2) {
            next
          }

          fixed_terms <- rownames(coeffs)[-1]
          estimates <- coeffs[-1, "Estimate", drop = FALSE]
          p_values <- coeffs[-1, "Pr(>|t|)", drop = FALSE]
          names(estimates) <- fixed_terms
          names(p_values) <- fixed_terms

          # If pairwise comparisons requested, compute them using emmeans
          pairwise_results <- NULL
          if (do_pairwise && !is.null(pairwise_contrasts)) {
            tryCatch({
              emm <- suppressMessages(emmeans::emmeans(model, specs = fixed_factors[1]))
              pairwise_results <- list()

              for (contrast_pair in pairwise_contrasts) {
                contrast_name <- paste(contrast_pair[1], "vs", contrast_pair[2], sep = "_")

                # Create contrast vector
                contrast_vector <- rep(0, length(levels(model_data[[fixed_factors[1]]])))
                names(contrast_vector) <- levels(model_data[[fixed_factors[1]]])
                contrast_vector[contrast_pair[1]] <- 1
                contrast_vector[contrast_pair[2]] <- -1

                # Compute contrast
                contr_result <- suppressMessages(emmeans::contrast(emm, method = list(contrast_vector)))
                contr_summary <- summary(contr_result)

                # FIXED: Ensure we get scalar values
                pairwise_results[[contrast_name]] <- list(
                  estimate = as.numeric(contr_summary$estimate[1]),
                  p_value = as.numeric(contr_summary$p.value[1])
                )
              }
            }, error = function(e) {
              NULL
            })
          }

          return(list(estimates = estimates, p_values = p_values,
                      contrasts = pairwise_results, model_type = model_type))

        } else if (method == "likelihood_ratio") {
          # Likelihood ratio method
          if (is.null(re_formula)) {
            formula_null <- as.formula("response ~ 1")
            formula_full <- as.formula(paste("response ~", fixed_effects))

            if (suppress_warnings) {
              model_null <- suppressMessages(suppressWarnings(lm(formula_null, data = model_data)))
              model_full <- suppressMessages(suppressWarnings(lm(formula_full, data = model_data)))
            } else {
              model_null <- lm(formula_null, data = model_data)
              model_full <- lm(formula_full, data = model_data)
            }

            coeffs <- summary(model_full)$coefficients
            if (nrow(coeffs) < 2) {
              next
            }

            fixed_terms <- rownames(coeffs)[-1]
            estimates <- coeffs[-1, "Estimate", drop = FALSE]
            names(estimates) <- fixed_terms

            f_test <- anova(model_null, model_full)
            p_value <- f_test[2, "Pr(>F)"]
            test_name <- paste("F_test", gsub("\\+", "_", gsub(" ", "", fixed_effects)), sep = "_")
            p_values <- setNames(p_value, test_name)

            # Pairwise comparisons for fixed effects model
            pairwise_results <- NULL
            if (do_pairwise && !is.null(pairwise_contrasts)) {
              tryCatch({
                emm <- suppressMessages(emmeans::emmeans(model_full, specs = fixed_factors[1]))
                pairwise_results <- list()

                for (contrast_pair in pairwise_contrasts) {
                  contrast_name <- paste(contrast_pair[1], "vs", contrast_pair[2], sep = "_")

                  contrast_vector <- rep(0, length(levels(model_data[[fixed_factors[1]]])))
                  names(contrast_vector) <- levels(model_data[[fixed_factors[1]]])
                  contrast_vector[contrast_pair[1]] <- 1
                  contrast_vector[contrast_pair[2]] <- -1

                  contr_result <- suppressMessages(emmeans::contrast(emm, method = list(contrast_vector)))
                  contr_summary <- summary(contr_result)

                  # FIXED: Ensure we get scalar values
                  pairwise_results[[contrast_name]] <- list(
                    estimate = as.numeric(contr_summary$estimate[1]),
                    p_value = as.numeric(contr_summary$p.value[1])
                  )
                }
              }, error = function(e) {
                NULL
              })
            }

            return(list(estimates = estimates, p_values = p_values,
                        contrasts = pairwise_results, model_type = model_type))
          } else {
            # Mixed model with likelihood ratio
            formula_null <- as.formula(paste("response ~", re_formula))
            formula_full <- as.formula(paste("response ~", fixed_effects, "+", re_formula))

            if (suppress_warnings) {
              model_null <- suppressMessages(suppressWarnings(lme4::lmer(formula_null, data = model_data, REML = REML)))
              model_full <- suppressMessages(suppressWarnings(lme4::lmer(formula_full, data = model_data, REML = REML)))
            } else {
              model_null <- lme4::lmer(formula_null, data = model_data, REML = REML)
              model_full <- lme4::lmer(formula_full, data = model_data, REML = REML)
            }

            is_singular <- lme4::isSingular(model_full)

            if (is_singular) {
              if (handle_singular == "drop") {
                next
              }
              if (handle_singular == "simplify" && idx < length(formulas_to_try)) {
                next
              }
            }

            coeffs <- summary(model_full)$coefficients
            if (nrow(coeffs) < 2) {
              next
            }

            fixed_terms <- rownames(coeffs)[-1]
            estimates <- coeffs[-1, "Estimate", drop = FALSE]
            names(estimates) <- fixed_terms

            lrt_result <- anova(model_null, model_full)
            p_value <- lrt_result[2, "Pr(>Chisq)"]
            lrt_name <- paste("LRT", gsub("\\+", "_", gsub(" ", "", fixed_effects)), sep = "_")
            p_values <- setNames(p_value, lrt_name)

            # Pairwise comparisons for mixed model
            pairwise_results <- NULL
            if (do_pairwise && !is.null(pairwise_contrasts)) {
              tryCatch({
                emm <- suppressMessages(emmeans::emmeans(model_full, specs = fixed_factors[1]))
                pairwise_results <- list()

                for (contrast_pair in pairwise_contrasts) {
                  contrast_name <- paste(contrast_pair[1], "vs", contrast_pair[2], sep = "_")

                  contrast_vector <- rep(0, length(levels(model_data[[fixed_factors[1]]])))
                  names(contrast_vector) <- levels(model_data[[fixed_factors[1]]])
                  contrast_vector[contrast_pair[1]] <- 1
                  contrast_vector[contrast_pair[2]] <- -1

                  contr_result <- suppressMessages(emmeans::contrast(emm, method = list(contrast_vector)))
                  contr_summary <- summary(contr_result)

                  # FIXED: Ensure we get scalar values
                  pairwise_results[[contrast_name]] <- list(
                    estimate = as.numeric(contr_summary$estimate[1]),
                    p_value = as.numeric(contr_summary$p.value[1])
                  )
                }
              }, error = function(e) {
                NULL
              })
            }

            return(list(estimates = estimates, p_values = p_values,
                        contrasts = pairwise_results, model_type = model_type))
          }
        }
      }, error = function(e) {
        NULL
      })
    }

    return(list(estimates = NA, p_values = NA, contrasts = NULL, model_type = "failed"))
  }

  # Process each cluster level
  all_results <- NULL
  model_type_counts <- list()
  singular_count <- 0
  failed_count <- 0

  for (cluster_level in cluster_levels) {
    message(paste("Processing cluster level:", cluster_level))

    # Subset data and metadata if cluster factor is specified
    if (cluster_factor != "none") {
      cluster_indices <- metadata[[cluster_factor]] == cluster_level
      subset_data <- data[, cluster_indices, drop = FALSE]
      subset_metadata <- metadata[cluster_indices, , drop = FALSE]
      cluster_suffix <- paste0("_within_", cluster_level)
    } else {
      subset_data <- data
      subset_metadata <- metadata
      cluster_suffix <- ""
    }

    if (ncol(subset_data) < 3) {
      warning(paste("Insufficient samples in cluster", cluster_level, ". Skipping."))
      next
    }

    # Process each feature
    p_values_list <- list()
    estimates_list <- list()
    contrasts_list <- list()
    model_types_list <- list()

    for (i in seq_len(nrow(subset_data))) {
      feature_name <- rownames(subset_data)[i]
      response_data <- subset_data[i, ]
      result <- fit_lmm_feature(response_data, subset_metadata, feature_name)

      p_values_list[[i]] <- result$p_values
      estimates_list[[i]] <- result$estimates
      contrasts_list[[i]] <- result$contrasts
      model_types_list[[i]] <- result$model_type

      if (length(result$p_values) == 1 && is.na(result$p_values)) {
        failed_count <- failed_count + 1
      }
    }

    # Track model types
    model_type_table <- table(unlist(model_types_list))
    model_type_counts[[cluster_level]] <- model_type_table

    # Create results dataframe for this cluster
    cluster_results <- data.frame(matrix(ncol = 0, nrow = nrow(subset_data)))
    rownames(cluster_results) <- rownames(subset_data)

    # Add pairwise comparison results if available
    if (do_pairwise && !is.null(pairwise_contrasts)) {
      for (contrast_pair in pairwise_contrasts) {
        contrast_name <- paste(contrast_pair[1], "vs", contrast_pair[2], sep = "_")

        # FIXED: Initialize with NA for all rows to ensure correct length
        contrast_estimates <- rep(NA_real_, nrow(subset_data))
        contrast_pvals <- rep(NA_real_, nrow(subset_data))

        # Fill in values where available with EXTRA robust checking
        for (i in seq_along(contrasts_list)) {
          tryCatch({
            contr <- contrasts_list[[i]]

            # Check if contrasts exist and are valid
            if (!is.null(contr) &&
                is.list(contr) &&
                length(contr) > 0 &&
                contrast_name %in% names(contr)) {

              # Get the contrast element
              contr_element <- contr[[contrast_name]]

              # Check that the contrast element itself is valid
              if (!is.null(contr_element) &&
                  is.list(contr_element) &&
                  "estimate" %in% names(contr_element) &&
                  "p_value" %in% names(contr_element)) {

                # Extract estimate - force to single scalar
                est_val <- contr_element$estimate
                if (!is.null(est_val) && length(est_val) > 0) {
                  est_val <- as.numeric(est_val)[1]  # Take only first element
                  if (is.finite(est_val)) {  # Check it's a valid number
                    contrast_estimates[i] <- est_val
                  }
                }

                # Extract p-value - force to single scalar
                p_val <- contr_element$p_value
                if (!is.null(p_val) && length(p_val) > 0) {
                  p_val <- as.numeric(p_val)[1]  # Take only first element
                  if (is.finite(p_val) && p_val >= 0 && p_val <= 1) {  # Valid p-value
                    contrast_pvals[i] <- p_val
                  }
                }
              }
            }
          }, error = function(e) {
            # Silently skip this feature if there's any error
            NULL
          })
        }

        # Final safety check: verify the vectors have the correct length
        if (length(contrast_estimates) != nrow(subset_data)) {
          warning(paste("Length mismatch for", contrast_name, "estimates: expected",
                        nrow(subset_data), "got", length(contrast_estimates)))
          # Force correct length
          contrast_estimates <- rep(NA_real_, nrow(subset_data))
        }

        if (length(contrast_pvals) != nrow(subset_data)) {
          warning(paste("Length mismatch for", contrast_name, "p-values: expected",
                        nrow(subset_data), "got", length(contrast_pvals)))
          # Force correct length
          contrast_pvals <- rep(NA_real_, nrow(subset_data))
        }

        # Create column names with new format
        if (is_log_transformed) {
          effect_col_name <- paste0("log(", contrast_pair[1], "/", contrast_pair[2], ")", cluster_suffix)
        } else {
          effect_col_name <- paste0(contrast_name, "_diff", cluster_suffix)
        }
        p_col_name <- paste0(contrast_name, cluster_suffix, "_LMMtest_p")

        # Final verification before assignment
        stopifnot(length(contrast_estimates) == nrow(subset_data))
        stopifnot(length(contrast_pvals) == nrow(subset_data))

        # Add to results - now absolutely guaranteed to have correct length
        cluster_results[[effect_col_name]] <- contrast_estimates
        cluster_results[[p_col_name]] <- contrast_pvals

        # Add adjusted p-values if requested
        if (padj) {
          if (sum(!is.na(contrast_pvals)) > 0) {
            adj <- p.adjust(contrast_pvals, method = padj_method)
            padj_col_name <- paste0(contrast_name, cluster_suffix, "_LMMtest_padj")

            # Verify adjusted p-values also have correct length
            stopifnot(length(adj) == nrow(subset_data))
            cluster_results[[padj_col_name]] <- adj
          }
        }
      }
    } else {
      # Original behavior for non-pairwise comparisons
      # Get all unique term names from standard model output
      all_terms <- unique(unlist(lapply(p_values_list, names)))
      all_terms <- all_terms[!is.na(all_terms)]

      if (length(all_terms) > 0) {
        # Create matrices for p-values and estimates
        p_matrix <- matrix(NA, nrow = nrow(subset_data), ncol = length(all_terms))
        est_matrix <- matrix(NA, nrow = nrow(subset_data), ncol = length(all_terms))
        colnames(p_matrix) <- all_terms
        colnames(est_matrix) <- all_terms
        rownames(p_matrix) <- rownames(subset_data)
        rownames(est_matrix) <- rownames(subset_data)

        # Fill in values
        for (i in seq_along(p_values_list)) {
          p_vals <- p_values_list[[i]]
          est_vals <- estimates_list[[i]]
          if (length(p_vals) > 0 && !all(is.na(p_vals))) {
            valid_names <- names(p_vals)[!is.na(names(p_vals))]
            if (length(valid_names) > 0) {
              p_matrix[i, valid_names] <- p_vals[valid_names]
              est_matrix[i, valid_names] <- est_vals[valid_names]
            }
          }
        }

        # Add results for each comparison with new naming convention
        for (j in seq_len(ncol(est_matrix))) {
          term_name <- colnames(est_matrix)[j]

          # Parse the term name to create cleaner comparison names
          # Remove the factor name prefix (e.g., "condition" from "conditionCTL")
          comparison_name <- term_name
          if (length(fixed_factors) == 1) {
            comparison_name <- gsub(paste0("^", fixed_factors[1]), "", term_name)
          }

          # Get factor levels to determine reference
          if (length(fixed_factors) == 1) {
            if (cluster_factor != "none") {
              factor_levels <- levels(subset_metadata[[fixed_factors[1]]])
            } else {
              factor_levels <- levels(metadata[[fixed_factors[1]]])
            }
            reference_level <- factor_levels[1]

            # Create pairwise name (e.g., "CTL_vs_CHR")
            pairwise_name <- paste(comparison_name, "vs", reference_level, sep = "_")
          } else {
            pairwise_name <- comparison_name
            reference_level <- NULL
          }

          # Create column names with new format
          if (is_log_transformed) {
            if (!is.null(reference_level)) {
              effect_col_name <- paste0("log(", comparison_name, "/", reference_level, ")", cluster_suffix)
            } else {
              effect_col_name <- paste0("log(", comparison_name, ")", cluster_suffix)
            }
          } else {
            effect_col_name <- paste0(pairwise_name, "_diff", cluster_suffix)
          }
          p_col_name <- paste0(pairwise_name, cluster_suffix, "_LMMtest_p")

          # Add effect sizes and p-values
          cluster_results[[effect_col_name]] <- est_matrix[, j]
          cluster_results[[p_col_name]] <- p_matrix[, j]

          # Add adjusted p-values if requested
          if (padj) {
            p_col <- p_matrix[, j]
            if (sum(!is.na(p_col)) > 0) {
              adj <- p.adjust(p_col, method = padj_method)
              padj_col_name <- paste0(pairwise_name, cluster_suffix, "_LMMtest_padj")
              cluster_results[[padj_col_name]] <- adj
            }
          }
        }
      }
    }

    # Combine results
    if (is.null(all_results)) {
      all_results <- cluster_results
    } else {
      all_results <- cbind(all_results, cluster_results)
    }
  }

  # Add results to statistics layer
  if (!is.null(all_results) && ncol(all_results) > 0) {
    romics_object$statistics <- cbind(romics_object$statistics, all_results)
  }

  # Print summary message
  n_comparisons <- if (!is.null(all_results)) ncol(all_results) else 0
  if (n_comparisons > 0) {
    message(paste("Linear Mixed Model analysis complete.", n_comparisons, "statistical columns were added to the statistics layer."))
    if (is_log_transformed) {
      message("Data detected as log-transformed. Effect sizes represent log fold changes.")
    } else {
      message("Data detected as non-log-transformed. Effect sizes represent arithmetic differences.")
    }

    if (do_pairwise) {
      message(paste("Pairwise comparisons performed:", length(pairwise_contrasts), "comparisons per cluster."))
    }

    message("Results include: effect sizes, p-values, and adjusted p-values.")

    if (failed_count > 0) {
      message(paste("Note:", failed_count, "features had NA results due to model fitting issues."))
    }

    if (length(model_type_counts) > 0) {
      message("\nModel fitting summary by cluster:")
      for (cluster_name in names(model_type_counts)) {
        message(paste("  Cluster:", cluster_name))
        counts <- model_type_counts[[cluster_name]]
        for (model_type in names(counts)) {
          message(paste("    -", model_type, ":", counts[model_type]))
        }
      }
    }

    if (handle_singular == "simplify") {
      message("\nNote: Singular fits were handled by simplifying random effects structure or using fixed effects only models.")
    }
  } else {
    warning("No results were generated. Check your data and model specifications.")
  }

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)
  return(romics_object)
}

#' romicsExtractSignificantFeatures()
#' @description Extract the names of features that are significant based on any statistical test in the statistics layer
#' @param romics_object A romics_object with calculated statistics
#' @param p_threshold Numeric value for significance threshold. Default: 0.05
#' @param padj_type Character string indicating which p-value type to use: "p" for raw p-values or "padj" for adjusted p-values. Default: "p"
#' @param specific_comparisons Character vector of specific comparisons to consider (e.g., "Alport_vs_Ctrl"). If NULL, considers all comparisons. Default: NULL
#' @param test_types Character vector of test types to consider (e.g., "Ttest", "Wilcox_test", "glmBinomialTest"). If NULL, considers all tests. Default: NULL
#' @param directionality Character string indicating which direction of change to consider: "all" (both up and down), "up" (upregulated only), "down" (downregulated only). Default: "all"
#' @param fc_threshold Numeric value for fold change threshold (used with directionality). For log-transformed data, this is the log fold change threshold. Default: 0 (no fold change filtering)
#' @details This function searches through statistical test columns in the statistics layer and identifies
#'          features that meet the significance criteria and directionality requirements in any test.
#'          For t-test and Wilcox test: up = log(FC) > 0 (or FC > 1), down = log(FC) < 0 (or FC < 1).
#'          For GLM binomial test: up = directionality = 1, down = directionality = -1.
#'          It automatically detects p-value columns based on the naming convention and the specified padj_type.
#' @return Character vector of significant feature names meeting the directionality criteria
#' @author Geremy Clair
#' @export
romicsExtractSignificantFeatures <- function(romics_object,
                                             p_threshold = 0.05,
                                             padj_type = "p",
                                             specific_comparisons = NULL,
                                             test_types = NULL,
                                             directionality = "all",
                                             fc_threshold = 0) {
  # Input validation
  if (!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if (!padj_type %in% c("p", "padj")) {
    stop("padj_type must be either 'p' or 'padj'")
  }
  if (!is.numeric(p_threshold) || p_threshold < 0 || p_threshold > 1) {
    stop("p_threshold must be a numeric value between 0 and 1")
  }
  if (!directionality %in% c("all", "up", "down")) {
    stop("directionality must be one of: 'all', 'up', 'down'")
  }
  if (!is.numeric(fc_threshold) || fc_threshold < 0) {
    stop("fc_threshold must be a non-negative numeric value")
  }
  if (is.null(romics_object$statistics) || ncol(romics_object$statistics) == 0) {
    warning("No statistics found in romics_object")
    return(character(0))
  }

  # Get all column names from statistics
  stat_columns <- colnames(romics_object$statistics)

  # Identify p-value columns based on padj_type
  if (padj_type == "padj") {
    target_columns <- grep("_padj$", stat_columns, value = TRUE)
  } else {
    target_columns <- grep("_p$", stat_columns, value = TRUE)
  }

  if (length(target_columns) == 0) {
    warning(paste("No", padj_type, "columns found in statistics"))
    return(character(0))
  }

  # Filter by test types if specified
  if (!is.null(test_types)) {
    test_pattern <- paste(test_types, collapse = "|")
    target_columns <- target_columns[grepl(paste0("_(", test_pattern, ")_", padj_type, "$"), target_columns)]
    if (length(target_columns) == 0) {
      warning(paste("No", padj_type, "columns found for the specified test types:", paste(test_types, collapse = ", ")))
      return(character(0))
    }
  }

  # Filter by specific comparisons if provided
  if (!is.null(specific_comparisons)) {
    comparison_pattern <- paste(specific_comparisons, collapse = "|")
    target_columns <- target_columns[grepl(comparison_pattern, target_columns)]
    if (length(target_columns) == 0) {
      warning(paste("No", padj_type, "columns found for the specified comparisons:", paste(specific_comparisons, collapse = ", ")))
      return(character(0))
    }
  }

  # Extract the statistics data for target columns
  stats_data <- romics_object$statistics[, target_columns, drop = FALSE]

  # Find significant features based on p-value threshold
  significant_matrix <- stats_data < p_threshold & !is.na(stats_data)

  # Apply directionality filter if not "all"
  if (directionality != "all") {

    # Initialize directionality matrix
    direction_matrix <- matrix(TRUE, nrow = nrow(stats_data), ncol = ncol(stats_data))
    rownames(direction_matrix) <- rownames(stats_data)
    colnames(direction_matrix) <- colnames(stats_data)

    for (i in 1:ncol(stats_data)) {
      col_name <- colnames(stats_data)[i]

      # Determine test type from column name
      if (grepl("_Ttest_", col_name) || grepl("_Wilcox_test_", col_name)) {
        # Use fold change columns for t-test and Wilcox test

        # Extract comparison name (e.g., "Alport_vs_Ctrl" from "Alport_vs_Ctrl_Ttest_p")
        comparison_name <- sub("_Ttest.*", "", col_name)
        comparison_name <- sub("_Wilcox_test.*", "", comparison_name)

        # Look for corresponding fold change column: log(A/B) format
        fc_col_pattern <- paste0("log\\(", gsub("_vs_", "/", comparison_name), "\\)")
        fc_col <- grep(fc_col_pattern, stat_columns, value = TRUE)

        if (length(fc_col) > 0) {
          fc_values <- romics_object$statistics[[fc_col[1]]]

          if (directionality == "up") {
            # Up: log(FC) > fc_threshold (default fc_threshold = 0, so log(FC) > 0)
            direction_matrix[, i] <- fc_values > fc_threshold & !is.na(fc_values)
          } else if (directionality == "down") {
            # Down: log(FC) < -fc_threshold (default fc_threshold = 0, so log(FC) < 0)
            direction_matrix[, i] <- fc_values < -fc_threshold & !is.na(fc_values)
          }
        } else {
          warning(paste("Could not find fold change column for", col_name, "- looking for pattern:", fc_col_pattern))
          direction_matrix[, i] <- FALSE
        }

      } else if (grepl("_glmBinomialTest_", col_name)) {
        # Use directionality column for GLM binomial test

        # Extract comparison name (e.g., "Alport_vs_Ctrl" from "Alport_vs_Ctrl_glmBinomialTest_p")
        comparison_name <- sub("_glmBinomialTest.*", "", col_name)
        dir_col_name <- paste0(comparison_name, "_directionality")

        if (dir_col_name %in% stat_columns) {
          dir_values <- romics_object$statistics[[dir_col_name]]

          if (directionality == "up") {
            # Up: directionality = 1
            direction_matrix[, i] <- dir_values == 1 & !is.na(dir_values)
          } else if (directionality == "down") {
            # Down: directionality = -1
            direction_matrix[, i] <- dir_values == -1 & !is.na(dir_values)
          }
        } else {
          warning(paste("Could not find directionality column for", col_name, "- looking for:", dir_col_name))
          direction_matrix[, i] <- FALSE
        }
      } else {
        # Unknown test type, assume all pass directionality filter
        warning(paste("Unknown test type for column", col_name, "- ignoring directionality filter"))
        direction_matrix[, i] <- TRUE
      }
    }

    # Combine significance and directionality filters
    final_matrix <- significant_matrix & direction_matrix
  } else {
    # No directionality filter
    final_matrix <- significant_matrix
  }

  # Get features that meet criteria in any test
  significant_features <- rownames(stats_data)[rowSums(final_matrix, na.rm = TRUE) > 0]

  # Provide informative message
  message(paste("Found", length(significant_features), "significant features out of",
                nrow(stats_data), "total features"))
  message(paste("Using", padj_type, "with threshold:", p_threshold))
  message(paste("Directionality filter:", directionality))
  if (directionality != "all") {
    if (fc_threshold > 0) {
      message(paste("Fold change threshold:", fc_threshold))
      message(paste("  Up: log(FC) >", fc_threshold, "or directionality = 1"))
      message(paste("  Down: log(FC) < -", fc_threshold, "or directionality = -1"))
    } else {
      message("  Up: log(FC) > 0 or directionality = 1")
      message("  Down: log(FC) < 0 or directionality = -1")
    }
  }
  if (!is.null(test_types)) {
    message(paste("Test types considered:", paste(test_types, collapse = ", ")))
  }
  if (!is.null(specific_comparisons)) {
    message(paste("Comparisons considered:", paste(specific_comparisons, collapse = ", ")))
  }
  message(paste("P-value columns used:", paste(target_columns, collapse = ", ")))

  return(significant_features)
}

#' featureSignificanceToTests()
#' @description Evaluate which statistical tests find a specific feature significant
#' @param romics_object A romics_object with calculated statistics
#' @param feature Character string of feature name (must match a rowname in romics_object$data)
#' @param ptype Character string indicating which p-value type to use: "padj" for adjusted p-values or "p" for raw p-values. Default: "padj"
#' @param p_threshold Numeric value for significance threshold. Default: 0.05
#' @param test_types Character vector of test types to consider (e.g., c("Ttest", "Wilcox_test", "glmBinomialTest")). If NULL, considers all tests. Default: NULL
#' @param specific_comparisons Character vector of specific comparisons to consider (e.g., "Alport_vs_Ctrl"). If NULL, considers all comparisons. Default: NULL
#' @details Returns a data frame with one row per test, showing whether the feature is significant.
#'          Directionality (up/down) is included for pairwise comparison tests (t-test, Wilcox, GLM binomial)
#'          but not for multi-group tests (ANOVA). Direction is NA for non-significant results.
#' @return Data frame with columns: test_name (character), significant (logical), direction (character: "up"/"down"/NA)
#' @examples
#' \dontrun{
#'   # Get significance across all tests for a specific feature
#'   results <- featureSignificanceToTests(romics_proteins, feature = "ProteinA")
#'
#'   # Get significance for t-tests only
#'   results <- featureSignificanceToTests(romics_proteins, feature = "ProteinA", test_types = "Ttest")
#'
#'   # Use raw p-values instead of adjusted
#'   results <- featureSignificanceToTests(romics_proteins, feature = "ProteinA", ptype = "p")
#' }
#' @author Geremy Clair
#' @export
featureSignificanceToTests <- function(
  romics_object,
  feature,
  ptype = c("padj", "p"),
  p_threshold = 0.05,
  test_types = NULL,
  specific_comparisons = NULL
) {

  # Validate romics_object
  if (!is.romicsObject(romics_object)) {
    stop("Input must be a valid romics_object")
  }

  # Match ptype argument
  ptype <- match.arg(ptype)

  # Check feature exists
  if (!feature %in% rownames(romics_object$data)) {
    stop("Feature '", feature, "' not found in romics_object. Available features: ",
         paste(head(rownames(romics_object$data), 5), collapse = ", "), "...")
  }

  # Validate p_threshold
  if (!is.numeric(p_threshold) || p_threshold < 0 || p_threshold > 1) {
    stop("p_threshold must be a numeric value between 0 and 1")
  }

  # Check if statistics layer exists
  if (is.null(romics_object$statistics) || nrow(romics_object$statistics) == 0) {
    message("No statistical tests have been run on this romics_object")
    return(data.frame(test_name = character(), significant = logical(), direction = character(),
                      stringsAsFactors = FALSE))
  }

  stat_columns <- colnames(romics_object$statistics)

  # Identify target p-value columns
  if (ptype == "padj") {
    target_columns <- grep("_padj$", stat_columns, value = TRUE)
    if (length(target_columns) == 0) {
      warning("No adjusted p-values found. Using raw p-values instead.")
      target_columns <- grep("_p$", stat_columns, value = TRUE)
      # Exclude _padj columns if they exist (they end in 'j' not 'p')
      target_columns <- grep("_padj$", target_columns, value = TRUE, invert = TRUE)
    }
  } else {
    target_columns <- grep("_p$", stat_columns, value = TRUE)
    # Exclude _padj columns
    target_columns <- grep("_padj$", target_columns, value = TRUE, invert = TRUE)
  }

  if (length(target_columns) == 0) {
    message("No p-value columns found in statistics layer")
    return(data.frame(test_name = character(), significant = logical(), direction = character(),
                      stringsAsFactors = FALSE))
  }

  # Filter by test_types if specified
  if (!is.null(test_types)) {
    test_patterns <- paste0("_(", paste(test_types, collapse = "|"), ")_")
    target_columns <- target_columns[grepl(test_patterns, target_columns)]
  }

  # Filter by specific_comparisons if specified
  if (!is.null(specific_comparisons)) {
    target_columns <- target_columns[grepl(paste(specific_comparisons, collapse = "|"), target_columns)]
  }

  if (length(target_columns) == 0) {
    message("No tests match the specified filters")
    return(data.frame(test_name = character(), significant = logical(), direction = character(),
                      stringsAsFactors = FALSE))
  }

  # Extract p-values for the feature as a named vector
  feature_pvals <- as.numeric(romics_object$statistics[feature, target_columns])
  names(feature_pvals) <- target_columns

  # Initialize result data frame
  results <- data.frame(
    test_name = character(),
    significant = logical(),
    direction = character(),
    stringsAsFactors = FALSE
  )

  # Process each test column
  for (col in target_columns) {
    pval <- feature_pvals[col]

    # Determine test type
    test_type <- if (grepl("_Ttest_", col)) {
      "Ttest"
    } else if (grepl("_Wilcox_test_", col)) {
      "Wilcox_test"
    } else if (grepl("_glmBinomialTest_", col)) {
      "glmBinomialTest"
    } else if (grepl("_ANOVA_", col)) {
      "ANOVA"
    } else {
      "Other"
    }

    # Extract comparison name (remove test type and p/padj suffix)
    comparison_name <- col
    if (ptype == "padj") {
      comparison_name <- sub("_padj$", "", comparison_name)
    } else {
      comparison_name <- sub("_p$", "", comparison_name)
    }

    # Determine if significant
    is_significant <- !is.na(pval) && pval < p_threshold

    # Determine direction (only for pairwise comparison tests and if significant)
    direction <- NA_character_
    if (is_significant && test_type %in% c("Ttest", "Wilcox_test")) {
      # Extract comparison to find fold-change column
      comparison_base <- sub(paste0("_", test_type, "$"), "", comparison_name)
      fc_pattern <- gsub("_vs_", "/", comparison_base)

      # Look for fold-change column in format: log(A/B) or (A/B)
      fc_cols <- grep(paste0("(log)?\\(", fc_pattern, "\\)"), stat_columns, value = TRUE)

      if (length(fc_cols) > 0) {
        fc_col_name <- fc_cols[1]
        fc_val <- romics_object$statistics[feature, fc_col_name]
        if (!is.na(fc_val)) {
          # Determine if this is log-transformed based on column name
          is_log_transformed <- grepl("^log\\(", fc_col_name)

          if (is_log_transformed) {
            # For log-transformed: > 0 is up, < 0 is down
            direction <- if (fc_val > 0) "up" else if (fc_val < 0) "down" else NA_character_
          } else {
            # For raw ratio: > 1 is up, < 1 is down
            direction <- if (fc_val > 1) "up" else if (fc_val < 1) "down" else NA_character_
          }
        }
      }
    } else if (is_significant && test_type == "glmBinomialTest") {
      # Extract directionality from {comparison}_directionality column
      # Note: 1 = more in others (down), -1 = more in level (up), so signs are reversed
      comparison_base <- sub("_glmBinomialTest", "", comparison_name)
      dir_col <- paste0(comparison_base, "_directionality")

      if (dir_col %in% stat_columns) {
        dir_val <- romics_object$statistics[feature, dir_col]
        if (!is.na(dir_val)) {
          # Reverse the interpretation: -1 (more in level) = "up", 1 (more in others) = "down"
          direction <- if (dir_val < 0) "up" else if (dir_val > 0) "down" else NA_character_
        }
      }
    }
    # For ANOVA and other multi-group tests, direction remains NA

    # Add row to results
    results <- rbind(results, data.frame(
      test_name = comparison_name,
      significant = is_significant,
      direction = direction,
      stringsAsFactors = FALSE
    ))
  }

  # Sort by test_name for readability
  results <- results[order(results$test_name), ]
  rownames(results) <- NULL

  return(results)
}

#' romicsPlotSignificantFeatures()
#' @description Create a barplot showing the number of significant features for each statistical test
#' @param romics_object A romics_object with calculated statistics
#' @param p_threshold Numeric value for significance threshold. Default: 0.05
#' @param padj_type Character string indicating which p-value type to use: "p" for raw p-values or "padj" for adjusted p-values. Default: "p"
#' @param specific_comparisons Character vector of specific comparisons to consider (e.g., "Alport_vs_Ctrl"). If NULL, considers all comparisons. Default: NULL
#' @param test_types Character vector of test types to consider (e.g., "Ttest", "Wilcox_test", "glmBinomialTest"). If NULL, considers all tests. Default: NULL
#' @param split_directions Logical indicating whether to split up and downregulated features into separate bars. Default: FALSE
#' @param fc_threshold Numeric value for fold change threshold (used with directionality). Default: 0
#' @param colors Character vector of colors for bars. If split_directions=TRUE, provide colors for c("Up", "Down"). If FALSE, provide single color. Default: c("#E31A1C", "#1F78B4") for up/down or "#2cbcb2" for total
#' @param title Character string for plot title. If NULL, generates automatic title. Default: NULL
#' @param rotate_x_labels Logical indicating whether to rotate x-axis labels 45 degrees. Default: TRUE
#' @details This function uses romicsExtractSignificantFeatures() to count significant features for each test type
#'          and creates a barplot using ggplot2 with theme_ROP(). Can show total counts or split by up/down regulation.
#' @return ggplot object
#' @author [Your Name]
#' @export
romicsPlotSignificantFeatures <- function(romics_object,
                                          p_threshold = 0.05,
                                          padj_type = "p",
                                          specific_comparisons = NULL,
                                          test_types = NULL,
                                          split_directions = FALSE,
                                          fc_threshold = 0,
                                          colors = NULL,
                                          title = NULL,
                                          rotate_x_labels = TRUE) {
  # Input validation
  if (!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }

  library(ggplot2)

  # Set default colors
  if (is.null(colors)) {
    if (split_directions) {
      colors <- c("Up" = "#E31A1C", "Down" = "#1F78B4")
    } else {
      colors <- "#2cbcb2"
    }
  }

  # Get statistics columns
  stat_columns <- colnames(romics_object$statistics)

  # Identify p-value columns
  if (padj_type == "padj") {
    p_columns <- grep("_padj$", stat_columns, value = TRUE)
  } else {
    p_columns <- grep("_p$", stat_columns, value = TRUE)
  }

  if (length(p_columns) == 0) {
    stop(paste("No", padj_type, "columns found in statistics"))
  }

  # Filter by test types if specified
  if (!is.null(test_types)) {
    test_pattern <- paste(test_types, collapse = "|")
    p_columns <- p_columns[grepl(paste0("_(", test_pattern, ")_", padj_type, "$"), p_columns)]
  }

  # Filter by specific comparisons if provided
  if (!is.null(specific_comparisons)) {
    comparison_pattern <- paste(specific_comparisons, collapse = "|")
    p_columns <- p_columns[grepl(comparison_pattern, p_columns)]
  }

  if (length(p_columns) == 0) {
    stop("No p-value columns found matching the specified criteria")
  }

  # Create data frame to store results
  if (split_directions) {
    plot_data <- data.frame(
      Test = character(0),
      Direction = character(0),
      Count = numeric(0),
      stringsAsFactors = FALSE
    )

    # For each p-value column, get up and down counts
    for (col in p_columns) {
      # Extract test name for display
      test_name <- col
      test_name <- gsub("_padj$|_p$", "", test_name)
      test_name <- gsub("_", " ", test_name)

      # Get p-values for this test
      p_values <- romics_object$statistics[[col]]

      # Find significant features for this test
      sig_mask <- p_values < p_threshold & !is.na(p_values)

      # Extract comparison name to find fold change column
      # Handle Ttest, Wilcox_test, and other test types
      comparison_name <- sub("_Ttest.*", "", col)
      comparison_name <- sub("_Wilcox_test.*", "", comparison_name)
      comparison_name <- sub("_glmBinomialTest.*", "", comparison_name)

      # Initialize counts
      up_count <- 0
      down_count <- 0

      # Determine directionality based on test type
      if (grepl("_Ttest_", col) || grepl("_Wilcox_test_", col)) {
        # For t-test and Wilcox test, use fold change
        fc_col_pattern <- paste0("log\\(", gsub("_vs_", "/", comparison_name), "\\)")
        fc_cols <- grep(fc_col_pattern, colnames(romics_object$statistics), value = TRUE)

        if (length(fc_cols) > 0) {
          fc_values <- romics_object$statistics[[fc_cols[1]]]
          up_mask <- sig_mask & fc_values > fc_threshold & !is.na(fc_values)
          down_mask <- sig_mask & fc_values < -fc_threshold & !is.na(fc_values)
          up_count <- sum(up_mask, na.rm = TRUE)
          down_count <- sum(down_mask, na.rm = TRUE)
        }
      } else if (grepl("_glmBinomialTest_", col)) {
        # For GLM binomial test, use directionality column
        dir_col_name <- paste0(comparison_name, "_directionality")
        if (dir_col_name %in% colnames(romics_object$statistics)) {
          dir_values <- romics_object$statistics[[dir_col_name]]
          up_mask <- sig_mask & dir_values == 1 & !is.na(dir_values)
          down_mask <- sig_mask & dir_values == -1 & !is.na(dir_values)
          up_count <- sum(up_mask, na.rm = TRUE)
          down_count <- sum(down_mask, na.rm = TRUE)
        }
      }

      # Add to plot data
      plot_data <- rbind(plot_data,
                         data.frame(Test = test_name, Direction = "Up", Count = up_count),
                         data.frame(Test = test_name, Direction = "Down", Count = down_count))
    }

    # Create plot with split directions
    p <- ggplot(plot_data, aes(x = Test, y = Count, fill = Direction)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = colors) +
      labs(y = "Number of Significant Features",
           x = "Statistical Test",
           fill = "Direction") +
      theme_ROP()

  } else {
    # Total counts only
    plot_data <- data.frame(
      Test = character(0),
      Count = numeric(0),
      stringsAsFactors = FALSE
    )

    # For each p-value column, get total count
    for (col in p_columns) {
      # Extract test name for display
      test_name <- col
      test_name <- gsub("_padj$|_p$", "", test_name)
      test_name <- gsub("_", " ", test_name)

      # Get total significant features for this test
      p_values <- romics_object$statistics[[col]]
      sig_count <- sum(p_values < p_threshold & !is.na(p_values))

      # Add to plot data
      plot_data <- rbind(plot_data,
                         data.frame(Test = test_name, Count = sig_count))
    }

    # Create plot with total counts
    p <- ggplot(plot_data, aes(x = Test, y = Count)) +
      geom_bar(stat = "identity", fill = colors, alpha = 0.8) +
      labs(y = "Number of Significant Features",
           x = "Statistical Test") +
      theme_ROP()
  }

  # Add text labels on bars
  if (split_directions) {
    p <- p + geom_text(aes(label = Count),
                       position = position_dodge(width = 0.9),
                       vjust = -0.5, size = 3)
  } else {
    p <- p + geom_text(aes(label = Count),
                       vjust = -0.5, size = 3)
  }

  # Rotate x-axis labels if requested (add after theme_ROP)
  if (rotate_x_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  # Add title
  if (is.null(title)) {
    if (split_directions) {
      title <- paste0("Significant Features by Test Type and Direction\n(",
                      padj_type, " < ", p_threshold, ")")
    } else {
      title <- paste0("Significant Features by Test Type\n(",
                      padj_type, " < ", p_threshold, ")")
    }

    if (fc_threshold > 0) {
      title <- paste0(title, ", FC threshold: ", fc_threshold)
    }
  }

  p <- p + ggtitle(title)

  # Print summary message
  message("Barplot created with the following parameters:")
  message(paste("  P-value type:", padj_type))
  message(paste("  P-value threshold:", p_threshold))
  message(paste("  Split directions:", split_directions))
  message(paste("  Fold change threshold:", fc_threshold))
  message(paste("  Tests included:", paste(gsub("_padj$|_p$", "", gsub("_", " ", p_columns)), collapse = ", ")))

  return(p)
}

#' romicsTransferStatistics()
#' @description Transfers statistics columns from a source romics_object to a target romics_object, matching features by row name
#' @param source_romics_object A romics_object containing the statistics to transfer
#' @param target_romics_object A romics_object that will receive the statistics
#' @param statistics Character vector specifying which statistics columns to transfer (default: NULL, transfers all)
#' If NULL, all statistics columns from source will be transferred
#' @details This function matches features between source and target romics_objects by their row names.
#' Statistics columns are transferred for matching features. Features present in target but not in source will have NA values.
#' If target_romics_object does not have a statistics layer, one will be created.
#' @return Returns the target_romics_object with transferred statistics added/updated
#' @author Geremy Clair
#' @export
#' @examples
#' # Transfer all statistics from source to target
#' target_obj <- romicsTransferStatistics(
#'   source_romics_object = clustered_obj,
#'   target_romics_object = original_obj
#' )
#'
#' # Transfer specific statistics columns
#' target_obj <- romicsTransferStatistics(
#'   source_romics_object = clustered_obj,
#'   target_romics_object = original_obj,
#'   statistics = c("Cmeans_cluster", "Cmeans_membership_C1")
#' )
romicsTransferStatistics <- function(source_romics_object,
                                   target_romics_object,
                                   statistics = NULL) {
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(source_romics_object) || missing(source_romics_object)) {
    stop("source_romics_object is missing or is not in the appropriate format")
  }

  if(!is.romicsObject(target_romics_object) || missing(target_romics_object)) {
    stop("target_romics_object is missing or is not in the appropriate format")
  }

  if(is.null(source_romics_object$statistics)) {
    stop("source_romics_object does not have a statistics layer")
  }

  message("Transferring statistics from source to target romics_object...")
  message("  Source features: ", nrow(source_romics_object$data))
  message("  Target features: ", nrow(target_romics_object$data))

  # Determine which statistics to transfer
  if(is.null(statistics)) {
    stats_to_transfer <- colnames(source_romics_object$statistics)
    message("  Transferring all statistics (", length(stats_to_transfer), " columns)")
  } else {
    if(!is.character(statistics)) {
      stop("statistics must be NULL or a character vector")
    }
    missing_stats <- setdiff(statistics, colnames(source_romics_object$statistics))
    if(length(missing_stats) > 0) {
      stop("The following statistics columns were not found in source_romics_object: ",
           paste(missing_stats, collapse = ", "))
    }
    stats_to_transfer <- statistics
    message("  Transferring ", length(stats_to_transfer), " statistics columns")
  }

  # Initialize statistics layer in target if needed
  if(is.null(target_romics_object$statistics)) {
    target_romics_object$statistics <- data.frame(matrix(nrow = nrow(target_romics_object$data), ncol = 0))
    rownames(target_romics_object$statistics) <- rownames(target_romics_object$data)
  }

  # Get feature names
  source_features <- rownames(source_romics_object$statistics)
  target_features <- rownames(target_romics_object$data)

  # Transfer statistics
  for(stat_col in stats_to_transfer) {
    target_romics_object$statistics[[stat_col]] <- NA

    # Match features and transfer values
    for(feat in target_features) {
      source_idx <- which(source_features == feat)
      if(length(source_idx) == 1) {
        target_romics_object$statistics[feat, stat_col] <- source_romics_object$statistics[source_idx, stat_col]
      }
    }

    # Count matches
    n_matched <- sum(!is.na(target_romics_object$statistics[[stat_col]]))
    message("  Transferred '", stat_col, "': ", n_matched, " of ", length(target_features), " features matched")
  }

  # Update processing steps
  target_romics_object <- romicsUpdateSteps(target_romics_object, arguments)

  message("Statistics transfer complete!")
  return(target_romics_object)
}
