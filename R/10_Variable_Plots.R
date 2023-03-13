#' singleVariablePlot()
#' @description Plots the abundances for a given variable based on a given factor contained in the metadata of the romics_object. This function will use the factor selected to group the samples and will color the plots based on the level of this factor. This function will function with partial matches.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param variable A string, the function will look for any string containing this name (partial match work, however if multiple rows contain the same partial name, the function will stop and indicate the multiple options)
#' @param type Must be in the following list: 'jitter','boxplot','violin','jb',or 'jv'. Indicates the type of plot to be returned. 'jb' is for both a jitter and a boxplot. 'jv' is for both a jitter and a violin
#' @param factor has to be either 'main" or a factor of the romics_object, the list of factor can be retrieved using the function romicsFactorNames()
#' @param limits has to be a numerical vector of length 2 corresponding to the lower and higher y_limits of the plot, by default ggplot will automatically define these, has to be a numeric vector of length 2 c(limit_min, limit_max)
#' @param title String. Indicate the title of the plot. automatically it will use the full feature name found for the chosen variable (in the case of partial match).
#' @param pval character vector of length one. Has to be either "none", "p", or "padj". "none" by default , "p" or "padj" is employed it indicates that pairwise ttest pvalues should be displayed above the groups of the the main factor (these have to be calculated prior running this function).
#' @param y_bracket_pos has to be a numerical vector of length one indicating the position and position increments for the pval bracket (only necessary when pval in c('p', 'padj')).
#' @details This function allows to quickly plot a the values for a given variable across multiple levels of a given factor
#' @return a ggplot figure
#' @author Geremy Clair
#' @export
singleVariablePlot<-function(romics_object, variable="variable", type = "jb", factor="main", limits=c(-10,10), title="auto", pval = "p",y_bracket_pos=0.5){
  #general checkings
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
    if(!type %in% c("jitter", "boxplot", "violin", "jb", "jv")){stop("type has to be in 'jitter','boxplot','violin','jb',or 'jv')")}
  if(missing(factor)){factor="main"}
  if(factor=="main"){factor=romics_object$main_factor}
  if(!factor %in% romicsFactorNames(romics_object)){
    print("the chosen <factor> is not in the list of factor for the selected romics object, this list comprises:")
    print(romicsFactorNames(romics_object))
    stop()
  }
  if(!is.character(variable)&length(variable!=1)){stop("<variable> should be a character vector of lenght 1")}
  if(missing(title)){title="auto"}
  if(missing(pval)){pval<-"p"}
  if(length(pval)!=1 | !pval %in% c("none", "p","padj")){stop("<pval> has to be 'none','p',or 'padj'.")}
  if(factor != romics_object$main_factor){romicsChangeFactor(romics_object,main_factor = factor)}
  if(missing(y_bracket_pos)){y_bracket_pos=0.5}
  if(!is.numeric(y_bracket_pos) & y_bracket_pos<=0 & length(y_bracket_pos)!=1){stop("<y_bracket_pos> has to be a positive numerical vector of lenght 1.")}
  #find the variable indicated using a grepl function
  variable<-rownames(romics_object$data)[variable==rownames(romics_object$data)]

  if(length(variable)==0){stop("your variable was not present in the romics_object")}

  #gather the data and the groups and place this in a data.frame named data
  data<- t(as.character(romics_object$data[rownames(romics_object$data)==variable,]))
  group<-as.character(t(romics_object$metadata[rownames(romics_object$metadata)==factor,]))
  fill<-as.character(t(romics_object$metadata[rownames(romics_object$metadata)=="colors_romics",]))
  data<-rbind(intensity=data, group=group,fill=fill)
  data<-data.frame(t(data))
  colnames(data)<-c("intensity","group","fill")
  data$sample<-colnames(romics_object$data)

  data$intensity<-as.numeric(data$intensity)
  data<-data[!is.na(data$intensity),]
  fill<-unique(data$fill)

  #create the plots
  plot<-ggplot(data, aes(x = group, y = intensity))+ theme_ROP()

  if(!missing(limits)){plot<-plot+scale_y_continuous(limits=limits)}

  if(type %in% c("boxplot","jb")){
    plot<- plot +
      geom_boxplot(alpha=0.5)
  }

  if(type %in% c("violin","jv")){
    plot<-plot + geom_violin(aes(fill=group, alpha=0.25))+scale_fill_manual(values=unique(data$fill))
  }

  if(type %in% c("jitter","jb","jv")){
    plot<- plot +
      geom_jitter(aes(color=group),position=position_jitter(0.25),size=3) +
      scale_color_manual(values=fill)
  }

  if(type %in% c("jitter","violin","jv")){
    plot<- plot + stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1),
                                geom="pointrange", color = "black",shape=3)
  }

  if(title=="auto"){plot<-plot+ggtitle(variable)}else{plot<-plot+ggtitle(title)}

  if(pval!="none"){
    f<-factor
    factor<-as.factor(as.character(t(romics_object$metadata[romicsFactorNames(romics_object)==factor,])))
    factor<-as.factor(t(factor))
    levels_factor<-levels(factor)
    by2combinations<- data.frame(t(combn(levels_factor,2)))
    if(!paste0(by2combinations[1,1],"_vs_",by2combinations[1,2],"_Ttest_",pval) %in% colnames(romics_object$statistics)){
      if(!paste0(by2combinations[1,2],"_vs_",by2combinations[1,1],"_Ttest_",pval) %in% colnames(romics_object$statistics)){
        warning(paste0("The ttests for the selected for the factor <",f,"> were not previously calculated. The pvalues won't be displayed on the plot."))
      }
    }else{
    if(!paste0(by2combinations[1,1],"_vs_",by2combinations[1,2],"_Ttest_",pval) %in% colnames(romics_object$statistics)){
      by2combinations<- by2combinations[,2:1]}
    y_increment<-y_bracket_pos
    y_bracket_pos<-max(data$intensity)+y_increment
    for(i in 1:nrow(by2combinations)){
      plot<- plot + geom_bracket(xmin = by2combinations[i,1],
                                 xmax = by2combinations[i,2],
                                 y.position = y_bracket_pos,
                                 label =formatC(romics_object$statistics[variable==rownames(romics_object$statistics),
                                                                         colnames(Figs$statistics)==paste0(by2combinations[1,1],"_vs_",by2combinations[1,2],"_Ttest_",pval)],
                                                format = "f"))
    y_bracket_pos<-y_increment+y_bracket_pos


    }
    }
  }
  return(plot)
  }

