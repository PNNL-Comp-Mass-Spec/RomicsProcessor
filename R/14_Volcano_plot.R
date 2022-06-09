#' romicsVolcano()
#' @description generate one or multiple volcano plots from t.test or wilcox.test run using the functions romicsTtest() or romicsWilcoxTest() respectively
#' @param p_type 'p' or 'padj' to indicate the type of tpvalue to consider for the volcano plots.
#' @param p numeric value to indicate the maximum pvalue to use for the coloring of the volcano plot
#' @param colors numeric vector of lenght 3 used to color the features lower, not significantly changing and higher in the considered comparisons.
#' @param stat_type 't.test' or 'wilcox.test' to indicate what statistics to use for the volcano plots generation.
#' @param plot either 'all' if all paired comparisons have to be displayed OR a vector of numeric values comprised between 1 and the maximum number of possible plots generated.
#' @param plot_type 'plotly' or 'ggplot' to indicate the type of plot to be returned.
#' @details generate one or multiple volcano plots from t.test or wilcox.test run using the functions romicsTtest() or romicsWilcoxTest() respectively
#' @return This function will print different plot requested
#' @author Geremy Clair
#' @export
romicsVolcano<-function(romics_object, p_type= "p", p= 0.05, min_fold_change=0.6,colors = c("#2cbcb2", "#242021", "#d44e28"),stat_type="t.test",plot="all", plot_type="ggplot"){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(p_type)){p_type="p"}
  if(!p_type %in% c("p","padj")){stop("'p_type' has to be either 'p' or 'padj'")}
  if(missing(colors)){colors=c("#2cbcb2", "#242021", "#d44e28")}
  if(!is.character(colors)| length(colors)!=3){
    warning("'colors' should be a color character vector of lenght 3. The defaults colors were used")
    colors=c("#2cbcb2", "#242021", "#d44e28")
  }
  if(missing(p)){p=0.05}
  if(!is.numeric(p)|p>1|p<0){stop("'p' should be numeric and comprised between 0 and 1.")}
  if(missing(min_fold_change)){min_fold_change=0.6}
  if(!is.numeric(min_fold_change)|min_fold_change<0){stop("'min_fold_change' should be numeric and comprised higher than 0.")}
  if(missing(plot_type)){plot_type="ggplot"}
  if(!plot_type %in% c("plotly","ggplot")){
    warning("'plot_type' was not either 'plotly' or 'ggplot' the 'plotly' type was used by default.")
    plot_type="plotly"
    }
  if(missing(plot)){plot="all"}
  #extract stats
  stats<-romics_object$statistics
  #remove duplicate columns
  stats<-stats[,unique(colnames(stats))]

  #verify if Ttest and or WilcoxTest exist
  if(sum(grepl("_Ttest_p",colnames(stats)))+sum(grepl("_Wilcox_test_p",colnames(stats)))<=0){stop("Prior to plot the volcano plot(s) either T.tests or Wilcox.tests have to be run, to run these test use the functions romicsTtest() and/or romicsWilcoxTest()")}
  #if stat_type missing use t.test by default unless does not exist (then use wilcox.test)
  if(missing(stat_type)){
    if(sum(grepl("_Ttest_p",colnames(stats)))>0){
      stat_type="t.test"
      print("'stat_type' was missing 't.test' were used by default")}else{
      stat_type="wilcox.test"
      print("'stat_type' was missing 'wilcox.test' were used by default")}
    }
  if(!stat_type %in% c("t.test","wilcox.test")){stop("'stat_type' has to be either 't.test' or 'wilcox.test'.")}

  #remove columns with p or padj (depending on the p_type)
  if (p_type == "p") {
    stats <- stats[, !grepl("_padj$", colnames(stats))]
  }
  if (p_type == "padj") {
    stats <- stats[, !grepl("_p$", colnames(stats))]
  }

  #identify Ttest or WilcoxTests columns
  if(stat_type=="t.test"){
    test_col<-stats[grepl("_Ttest_",colnames(stats))]
    colnames(test_col)<-sub("_Ttest.*","",colnames(test_col))
    colnames(test_col)<-sub("_vs_","\\/",colnames(test_col))
  }else{
    test_col<-stats[grepl("_Wilcox_test_",colnames(stats))]
    colnames(test_col)<-sub("_Wilcox_test_.*","",colnames(test_col))
    colnames(test_col)<-sub("_vs_","\\/",colnames(test_col))
    }

  #transform the test columns to calculate the -log10(p(test))
    test_col<- log10(test_col)*-1
    fc_col <- stats[grepl("\\/", colnames(stats))] #collected earlier in function (following creation of 'stats' object)

  #identify the fold_change columns
  if(sum(grepl("log\\(",colnames(fc_col)))==ncol(fc_col)){fc_log=TRUE}else{
  if(sum(grepl("log\\(",colnames(fc_col)))==0){fc_log=FALSE}else{
    warning("some of the fold-changes in the statistics layer of the 'romics_object' were calculated both prior and after log_transform.")
    warning("Only the log transformed will be used to generate the Volcano plots")
    fc_col<-fc_col[,grepl("log\\(",colnames(fc_col))]
    }
  }

  #if not logged then log transform the fold change columns
  if(fc_log==FALSE){
    fc_col=log2(fc_col)
    log_type=2
    min_fold_change<-log2(min_fold_change)
    }else{
      if(romicsLogCheck(romics_object)&grepl("fun\\|log2",romics_object$steps[grepl("fun\\|log",romics_object$steps)])){log_type=2}else{log_type=10}
    }

  #format the colnames so they are identical to the pvalues ones
  colnames(fc_col)<-sub("log\\(","",colnames(fc_col))
  colnames(fc_col)<-sub("\\)$","",colnames(fc_col))

  minus_log_p<-log10(p)*-1

  if(sum(colnames(fc_col) %in% colnames(test_col))!=ncol(fc_col)){warning("Some of the fold-change columns were not having a equivalent statistical test to generate a Volcano plot.")}

  if(plot=="all"){plot<-1:ncol(fc_col)}

  if(!is.numeric(plot) & sum(!plot %in% 1:ncol(fc_col))!=0){stop(paste0("'plot' as to be either 'all' or a numeric vector with values comprised between 1 and ",ncol(fc_col),"."))}else{
    for(i in plot){
      df<-cbind(rownames(fc_col),fc_col[i], test_col[colnames(test_col)==colnames(fc_col)[i]])
      colnames(df)<-c("ID","fc","p")

      class<- rep("non_significant",nrow(df))
      class[df$p>minus_log_p&df$fc<(min_fold_change*-1)]<-"down"
      class[df$p>minus_log_p&df$fc>(min_fold_change)]<-"up"
      class<-paste0(class,"_in_",sub("\\/.*","",colnames(fc_col[i])))
      df$class<-class

      if(plot_type=="ggplot"){
      fig<-ggplot(df,aes(x=fc,y=p,colour=class))+geom_point(alpha=0.5)+
      theme_ROP()+ggtitle(paste0("Volcano plot for ",colnames(fc_col[i])))+
      xlab(paste0("log",log_type,"(",colnames(fc_col[i]),")"))+
      ylab(paste0("-log10(",p_type,"_",stat_type,"_",colnames(fc_col[i]),")"))+
      scale_colour_manual(values=colors)
      plot(fig)}else{
        title=paste0("Volcano plot for ",colnames(fc_col[i]))
      print(plot_ly(x = df$fc,
                   y = df$p,
                   color =df$class,
                   colors=colors,
                   type = "scatter",mode="markers",
                   text=paste("ID=",df$ID)) %>% layout(title=paste0("Volcano plot for ",colnames(fc_col[i])),
                            xaxis=list(title=paste0("log",log_type,"(",colnames(fc_col[i]),")")),
                            yaxis=list(title=paste0("-log10(",p_type,"_",stat_type,"_",colnames(fc_col[i]),")"))))
      }
    }}
  }

