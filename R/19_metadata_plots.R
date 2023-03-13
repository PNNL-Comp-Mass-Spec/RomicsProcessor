#' metadata_relations()
#' @description display relationships between different factors contained in the metadata layers of an Romics_object
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param factors A character vector (lenght>2) containing the list of factors of the romics_object. the list of factors can be identified by using the function romicsFactorNames()
#'
#' @details generates a plot displaying relationships between different factors
#' @return a ggplot
#' @author Geremy Clair
#' @export
#'
metadataRelations<-function(romics_object=romics_object,factors=c("factor1", "factor2", "factor3"),textsize = 3){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(!all(factors %in% romicsFactorNames(romics_object)) & length(factors)<2){stop("'factors' has to be a character vector of lenght >2 containing a list of factors from the romics_object, factors can be identified by using the function romicsFactorNames()")}
  if(!is.numeric(textsize) | length(textsize)!=1){stop("'textsize' has to be a single numerical value.")}
  if(missing(textsize)){textsize=3}

  meta<-as.data.frame(t(romics_object$metadata))
  m<-data.frame(matrix(ncol=length(factors),nrow=nrow(meta)))
  for(i in 1:length(factors)){
    m[,i]<-meta[colnames(meta)==factors[i]]
    colnames(m)[i]<-factors[i]
  }

  ml<-to_lodes_form(m)

  ggplot(ml, aes(x = x, stratum = stratum,
                               alluvium = alluvium, fill = stratum, label = stratum)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "NA") +
    geom_stratum(alpha = 0.5) +
    geom_text(stat = "stratum",size=textsize) +
    theme_ROP() +
    xlab("Metadata") +
    ylab("Frequency")+
    scale_fill_manual(values=ROP_colors[1:length(unique(ml$stratum))])
}


#' factorCountLevels()
#' @description display relationships between different factors contained in the metadata layers of an Romics_object
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param factors A character vector (lenght=1) containing a factors of the romics_object. the list of factors can be identified by using the function romicsFactorNames()
#'
#' @details generates a bar graph indicating how many times each level is seen in a factor
#' @return a ggplot
#' @author Geremy Clair
#' @export
#'
factorCountLevels<-function(romics_object=romics_object,factor="factor",textsize = 3){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(!all(factor %in% romicsFactorNames(romics_object)) & length(factor)!=1){stop("'factor' has to be a character vector of lenght=1 contained in the romics_object, factors can be identified by using the function romicsFactorNames()")}
  if(missing(factor)){factor=romics_object$main_factor}
  if(!is.numeric(textsize) | length(textsize)!=1){stop("'textsize' has to be a single numerical value.")}
  if(missing(textsize)){textsize=3}

  meta<-as.data.frame(t(romics_object$metadata))
  meta<-meta[,colnames(meta)==factor]
  f<-as.factor(meta)
  t<-table(f)
  t<-data.frame(levels=names(t),frequency=as.numeric(t))

  ggplot(t, aes(x = levels, y= frequency)) +
    geom_bar(stat = "identity",fill=ROP_colors[2])+
    geom_text(aes(label=frequency),vjust=-0.3, color="black", size=textsize)+
    theme_ROP() + xlab("Levels") +
    ylab("Frequency")+
    ggtitle(factor)+
    theme(axis.text.x = element_text(angle = 0, hjust=0))
}



