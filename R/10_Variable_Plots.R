#' singleVariablePlot()
#' @description Plots the abundances for a given variable based on a given factor contained in the metadata of the romics_object. This function will use the factor selected to group the samples and will color the plots based on the level of this factor. This function will function with partial matches.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param variable A string, the function will look for any string containing this name (partial match work, however if multiple rows contain the same partial name, the function will stop and indicate the multiple options)
#' @param type Must be in the following list: 'jitter','boxplot','violin','jb',or 'jv'. Indicates the type of plot to be returned. 'jb' is for both a jitter and a boxplot. 'jv' is for both a jitter and a violin
#' @param factor has to be either 'main" or a factor of the romics_object, the list of factor can be retrieved using the function romicsFactorNames()
#' @param limits is corresponding to the y_limits of the plot, by default ggplot will automatically define these, has to be a numeric vector of length 2 c(limit_min, limit_max)
#' @param title String. Indicate the title of the plot. automatically it will use the full feature name found for the chosen variable (in the case of partial match).
#' @details This function allows to quickly plot a the values for a given variable accross multiple levels of a given factor
#' @return a ggplot figure
#' @author Geremy Clair
#' @export
singleVariablePlot<-function(romics_object, variable="variable", type = "jb", factor="main", limits=c(-10,10), title="auto"){
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

  #find the variable indicated using a grepl function
  variable<-rownames(romics_object$data)[grepl(variable,rownames(romics_object$data))]
  if(length(variable)>1){
    warning("More than one variable contained the variable choosen:")
    warning(variable)
    stop()
  }
  if(length(variable)==0){stop("your variable was not present in the romics_object")}

  #gather the data and the groups and place this in a data.frame named data
  data<- t(as.character(romics_object$data[rownames(romics_object$data)==variable,]))
  group<-as.character(t(romics_object$metadata[rownames(romics_object$metadata)==factor,]))
  fill<-as.character(t(romics_object$metadata[rownames(romics_object$metadata)=="colors_romics",]))
  data<-rbind(intensity=data, group=group,fill=fill)
  data<-data.frame(t(data))
  colnames(data)<-c("intensity","group","fill")

  data$intensity<-as.numeric(data$intensity)

  fill<-unique(as.character(t(fill)))

  #create the plots
  plot<-ggplot(data, aes(x = group, y = intensity))+ theme_ROP()

  if(!missing(limits)){plot<-plot+scale_y_continuous(limits=limits)}

  if(type %in% c("jitter","jb","jv")){
    plot<- plot +
    geom_jitter(aes(color=group),position=position_jitter(0.25),size=3) +
    scale_color_manual(values=fill)
    }

  if(type %in% c("boxplot","jb")){
    plot<- plot +
      geom_boxplot(alpha=0.5)
  }

  if(type %in% c("violin","jv")){
    plot<-plot + geom_violin(aes(fill=group, alpha=0.25))+scale_fill_manual(values=unique(data$fill))
  }

  if(type %in% c("jitter","violin","jv")){
    plot<- plot + stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1),
                                geom="pointrange", color = "black",shape=3)
  }

  if(title=="auto"){plot<-plot+ggtitle(variable)}else{plot<-plot+ggtitle(title)}

  return(plot)
  }
