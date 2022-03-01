#' lmp()
#' @description This function extract the overall p-value of a linear model
#' @param modelobject as to be an object of the class 'lm'
#' @details This function extract the p-value of an object of class 'lm'
#' @return This function returns the p-value of an object of class 'lm'
#' @author Geremy Clair and Feng Song
#' @export
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("The function lmp has been design to work on objects of class 'lm', your object was not of this class")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

#' getTrend()
#' @description This function extract the linear and quadratic trend pvalues and direction based on the pvalue it attributes the best fitted trend
#' @param sData is a serie of data. It should be a vector of numerical values corresponding to the different points of the serie (e.g. time serie)
#' @param m is the metadata associated to the serie (e.g. for a time-serie, it would be the time) if m is null a numeric sequence from 1:n will be used.
#' @param p is the minimum fitting pvalue (by default p = 0.05)
#' @param type indicates the type of functions to be fitted for the data serie. Has to be either "quadratic" (quadratic only), "linear" (linear only), or "both" (linear+quardatic)
#' @details this function enables to calculate and extract the linear or quardratic trend p values and the best fitted trend (e.g. linear_increasing, linear_decreasing, quadratic_concave,quadratic_convex)
#' @return This function returns a list
#' @author Geremy Clair and Feng Song
#' @export
getTrend <- function(sData, m = NULL, p=0.05, type = "both" ) {
    len = length(sData)
    if (is.null(m)) {
        m <- seq(len)
    } else {
        if (length(m) != length(sData)) {
            warning("The sData and m are not the same size! Using a numerical sequence as m!")
            m <- seq(len)
        }
    }
    if (len < 100) {
        modelT <- seq(min(m), max(m), (max(m)-min(m))/100)
    } else {
        modelT <- m
    }

    if(type %in% c("both", "linear")){
        liModel <- lm(sData ~ m)
        linear_p <- lmp(liModel)
        if(linear_p<0.05){
            if(liModel$coefficients[2]>0){linear_direction = "linear_increasing"}else{linear_direction = "linear_decreasing"}
        }else{linear_direction ="other_trend"}
    }

    if(type %in% c("both", "quadratic")){
        qModel <- lm(sData ~ poly(m,2))
        quadratic_p <- lmp(qModel)
        quadratic_direction <- qModel$coefficients[3]
        if(quadratic_p<0.05){
            if(qModel$coefficients[3]>0){quadratic_type = "quadratic_convex"}else{quadratic_type = "quadratic_concave"}
        }else{quadratic_type="other_trend"}
    }

    if(type == "both"){
        if(linear_p < 0.05 || quadratic_p < 0.05){
            if (quadratic_p > linear_p){best_fitted_trend = linear_direction}else{best_fitted_trend = quadratic_type}
            }else{best_fitted_trend="other_trend"}
        results <- list(linear_p=linear_p,quadratic_p=quadratic_p,best_fitted_trend=best_fitted_trend)
        }

    if(type == "linear"){results <- list(linear_p=linear_p,best_fitted_trend=linear_direction)}

    if(type == "quadratic"){results <- list(quadratic_p=quadratic_p,best_fitted_trend=quadratic_type)}

    return(results)
    }

#' romicsTrend()
#' @description This function extract the linear and quadratic trend pvalues and direction based on the pvalue it attributes the best fitted trend
#' @param romics_object has to be a romics_object created using the function romicsCreateObject() (see dedicated documentation) and has to be numeric.
#' @param factor has to be a factor contained in the romics_object created, the list of factor can be extracted using the function romicsFactorNames() (see dedicated documentation)
#' @param log_factor Boolean, indicate if the factor needs to be log transformed (TRUE), or not (FALSE). By default the factor is not transformed
#' @param type indicates if linear and/or quadratic trends analysis have to be performed, has to be 'both', 'linear',or 'quadratic'.
#' @param p indicates the minimal fitting p-value to be considered.
#' @details This function allows to calculate the linear and quadratic p value and estimate the best fitting model
#' @return This function returns a romics_object with new columns in the stat layer corresponding to the best fitting model, and the linear and/or quadratic p-values
#' @export
romicsTrend <- function(romics_object, factor="main",  log_factor = FALSE, type = "both", p=0.05) {
    arguments<-as.list(match.call())
    if(!is.romics_object(romics_object) || missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
    if(missing(factor)){factor="main"}
    if(!factor %in% c("main",romicsFactorNames(romics_object))){
        warning("The selected factor is not in the list of factors of the romics_object")
        warning(romicsFactorNames(romics_object))
    }
    #extract_data
    data <- romics_object$data

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
    f<-romics_object$metadata[romicsFactorNames(romics_object)==factor,]
        options(warn=-1)
    f<-as.numeric(f)
    if(length(f)==1){if(is.na(f)){stop("The factor selected could not be converted to a list of numbers.")}}
    options(warn=0)


    #order the data and f
    data <- data[, order(f)]

    if (log_factor == TRUE) {
        f <- log(sort(f))
    } else {
        f <- sort(f)
    }

    for(i in 1:nrow(data)){
        if(i==1){r=t(getTrend(as.numeric(data[i,]),f,p,type))}else{r=rbind(r,t(getTrend(as.numeric(data[i,]),f,p,type)))}
        }

    r<-data.frame(r)
    r[,1]<-as.numeric(r[,1])
    if(ncol(r)==2){
    r[,2]<-as.character(r[,2])
    }else{
    r[,2]<-as.numeric(r[,2])
    r[,3]<-as.character(r[,3])
    }

    r<-data.frame(r)
    romics_object$statistics <- cbind(romics_object$statistics, r)

    #print info
    print("The trend analysis columns were added to the statistics")

    #update steps
    romics_object<-romicsUpdateSteps(romics_object,arguments)

    return(romics_object)
}

#' singleVariableTrend()
#' @description This function plots the best fitted trend and a dotplot of the analytes intensities.
#' @param romics_object has to be a romics_object created using the function romicsCreateObject() (see dedicated documentation) and has to be numeric.
#' @param variable has to be a character vector. Romics will try to find an exact match or a partial match of the variable in the analytes names (colnames(romics_object$data))
#' @param factor has to be a factor contained in the romics_object created, the list of factor can be extracted using the function romicsFactorNames() (see dedicated documentation)
#' @param log_factor Boolean, indicate if the factor needs to be log transformed (TRUE), or not (FALSE). By default the factor is not transformed
#' @param title this has to be a character vector, either "auto" (in this case the rowname containing the variable will be displayed as title) or a title of your choice.
#' @details this function enables to calculate and extract the linear or quardratic trend p values and the best fitted trend (e.g. linear_increasing, linear_decreasing, quadratic_concave,quadratic_convex)
#' @return This function returns a list
#' @author Geremy Clair and Feng Song
#' @export
singleVariableTrend<-function(romics_object, variable="variable", factor="main", log_factor = FALSE, title="auto"){
    #general checkings
    trendArgs <- romics_object$steps[max(which(grepl("fun\\|romicsTrend",romics_object$steps)))]
    if (grepl("factor=", trendArgs)) {
        factor <- regmatches(trendArgs, regexec("factor='(.*?)'", trendArgs))[[1]][2]
    }
    if (grepl("log_factor=", trendArgs)) {
        factorscale <- regmatches(trendArgs, regexec("log_factor='(.*?)'", trendArgs))[[1]][2]
    }
    if(!is.romics_object(romics_object) || missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
    if(missing(factor)){factor="main"}
    if(!factor %in% c("main",romicsFactorNames(romics_object))){
        warning("The selected factor is not in the list of factors of the romics_object")
        warning(romicsFactorNames(romics_object))
    }

    #if factor is main extract the factor from the romics_object$main_factor
    if(factor=="main"){factor<-romics_object$main_factor}
    #extract the factor from the metadata
    f<-romics_object$metadata[romicsFactorNames(romics_object)==factor,]

    options(warn=-1)
    f<-as.numeric(f)
    if(length(f)==1){if(is.na(f)){stop("The factor selected could not be converted to a list of numbers.")}}
    options(warn=0)

    if(!is.character(variable)&&length(variable!=1)){stop("<variable> should be a character vector of lenght 1")}
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
    #order the data and f
    data <- data[, order(f)]

    if (log_factor==TRUE) {
        f <- log(sort(f))
    } else {
        f <- sort(f)
    }

    trendType <- romics_object$statistics[[variable, "best_fitted_trend"]]

    if (trendType == "linear_increasing" || trendType == "linear_decreasing") {
        model <- lm(data ~ f)
        tm <- seq(min(f), max(f), (max(f) - min(f))/100)
        ym <- predict(model, list(f=tm))
        trendPvalue<-romics_object$statistics[[variable, "linear_p"]]
    } else if (trendType == "quadratic_concave" || trendType == "quadratic_convex") {
        model <- lm(data ~ poly(f,2))
        tm <- seq(min(f), max(f), (max(f) - min(f))/100)
        ym <- predict(model, list(f=tm))
        trendPvalue<-romics_object$statistics[[variable, "quadratic_p"]]
    } else {
        ma <- function(x, n = 3){as.vector(filter(x, rep(1 / n, n), sides = 2))}
        tm <- f
        ym <- ma(data)
    }

    data<-rbind(t=f, y=data)
    data<-data.frame(t(data))
    colnames(data)<-c("t","y")
    data$y<-as.numeric(data$y)
    data$t<-as.numeric(data$t)

    ms<-rbind(t=tm, y=ym)
    ms<-data.frame(t(ms))
    colnames(ms)<-c("t","y")
    ms$y<-as.numeric(ms$y)
    ms$t<-as.numeric(ms$t)

    #create the plots
    plot <- ggplot() + theme_ROP() + geom_point(aes(x = t, y = y), data) +
        geom_line(aes(x = t, y = y), ms) +
        xlab(paste0(factor, " (", factorscale, "-scale)")) +
        ylab("data")

    if(is.na(trendPvalue)){
        plot<- plot +
            annotate("text", x=(max(data$t) - min(data$t))/2, y=max(data$y), label= trendType)
    } else {
        plot<- plot +
            annotate("text", x=(max(data$t) - min(data$t))/2, y=max(data$y), label= paste0(trendType, ", p=", trendPvalue))
    }

    if(title=="auto"){plot<-plot+ggtitle(variable)}else{plot<-plot+ggtitle(title)}

    return(plot)
}

#' romicsTrendHeatmap()
#' @description This function plots a heatmap of all the proteins by trend type??
#' @param romics_object has to be a romics_object created using the function romicsCreateObject() (see dedicated documentation) and has to be numeric.
#' @param factor has to be a factor contained in the romics_object created, the list of factor can be extracted using the function romicsFactorNames() (see dedicated documentation)
#' @param log_factor Boolean, indicate if the factor needs to be log transformed (TRUE), or not (FALSE). By default the factor is not transformed
#' @param ... arguments to be passed to heatmap.2() see this function documentation for more details
#' @details this function enables to calculate and extract the linear or quardratic trend p values and the best fitted trend (e.g. linear_increasing, linear_decreasing, quadratic_concave,quadratic_convex)
#' @return This function returns a list
#' @author Geremy Clair and Feng Song
#' @export
romicsTrendHeatmap<-function(romics_object,factor="main", log_factor = FALSE, ...){
    trendArgs <- romics_object$steps[max(which(grepl("fun\\|romicsTrend",romics_object$steps)))]
    if(!is.romics_object(romics_object) || missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
    if(is.null(romics_object$statistics$best_fitted_trend)){print("The trend analysis was not perform on the romics_object selected, please run the trend analysis using the function romicsTrend()")}
    if (grepl("factor=", trendArgs)) {
        factor <- regmatches(trendArgs, regexec("factor='(.*?)'", trendArgs))[[1]][2]
    }
    if (grepl("log_factor=", trendArgs)) {
        factorscale <- regmatches(trendArgs, regexec("factorscale='(.*?)'", trendArgs))[[1]][2]
    }
    if(missing(factor)){factor="main"}
    if(!factor %in% c("main",romicsFactorNames(romics_object))){
        warning("The selected factor is not in the list of factors of the romics_object")
        warning(romicsFactorNames(romics_object))
    }
    if(factor=="main"){factor<-romics_object$main_factor}

    data<-romics_object$data
    trend<-romics_object$statistics$best_fitted_trend

    #extract the factor from the metadata
    f<-romics_object$metadata[romicsFactorNames(romics_object)==factor,]
    options(warn=-1)
    f<-as.numeric(f)
    if(length(f)==1){if(is.na(f)){stop("The factor selected could not be converted to a list of numbers.")}}
    options(warn=0)

    data<-data[trend!="other_trend",]
    trend<-trend[trend!="other_trend"]

    data<-data[order(trend),order(f)]
    trend<-trend[order(trend)]

    rowCols <- trend
    rowCols[rowCols=="linear_decreasing"]<-"#ff0066"
    rowCols[rowCols=="linear_increasing"]<-"#0099ff"
    rowCols[rowCols=="quadratic_concave"]<-"#ffcc00"
    rowCols[rowCols=="quadratic_convex"]<-"#9900ff"

    if (log_factor == TRUE) {
        f <- log(f)
        f <- sort(f)
    } else {
        f <- sort(f)
    }

    colnames(data)<-paste0(colnames(data),"(",round(f,6),")")

    heatmap.2(as.matrix(data),
              dendrogram='none',
              Rowv = FALSE,
              Colv = FALSE,
              scale = "row",
              trace = "none",
              col=viridis(100),
              breaks=seq(-2, 2, length.out=101),
              RowSideColors=rowCols,
              ...
    )
}


