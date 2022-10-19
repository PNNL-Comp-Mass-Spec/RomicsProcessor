#' romicsCombineObjects()
#' @description Combine 2 Romics objects, remove the calculated stats ensure the log transformation are identical.
#' @param romics_object1 A romics_object created using romicsCreateObject().
#' @param romics_object2 A romics_object created using romicsCreateObject().
#' @param suffix1 A character string to append to the rownames of the romics_object1
#' @param suffix2 A character string to append to the rownames of the romics_object2
#'
#' @details Combine 2 romics objects, only works with romics object that have identical metadata (in the same order), it is therefore recommended to create the object with the same column order if you want to combine the two object later.
#' @return A combined romics_object containing the new layer romics_object$steps_old.
#' @author Geremy Clair
#' @export
romicsCombineObjects<-function(romics_object1,romics_object2,suffix1="_neg",suffix2="_pos"){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object1) | missing(romics_object1)) {stop("<romics_object1> is missing or is not in the appropriate format")}
  if(!is.romics_object(romics_object2) | missing(romics_object2)) {stop("<romics_object2> is missing or is not in the appropriate format")}
  if(missing(suffix1)){suffix1="romics_object1"}
  if(missing(suffix2)){suffix2="romics_object2"}
  if(!(is.character(suffix1)|length(suffix1)==1)){stop("<suffix1> is in the wrong format")}
  if(!(is.character(suffix2)|length(suffix2)==1)){stop("<suffix2> is in the wrong format")}
  romics_object<-romics_object1

  #the steps of the previous romics object will be recorded in the steps_old layer
  steps1<- paste0(romics_object1$steps,"|",suffix1)
  steps2<- paste0(romics_object2$steps,"|",suffix2)
  steps_old<-c(steps1, steps2)

  romics_object$steps<-romics_object1$steps[1:3]
  romicsUpdateSteps(romics_object,arguments)

  #if data was logged, unlog it prior combination
    log1<-0
    log2<-0

    if((sum(romicsSteps(romics_object1,show_dates = F,show_details = F)=="log2transform")-
       sum(romicsSteps(romics_object1,show_dates = F,show_details = F)=="unlog2data")) %% 2 == 1 ){
      log1<-log1+2
    }

    if((sum(romicsSteps(romics_object2,show_dates = F,show_details = F)=="log2transform")-
        sum(romicsSteps(romics_object2,show_dates = F,show_details = F)=="unlog2data")) %% 2 == 1 ){
      log2<-log2+2
    }

    if((sum(romicsSteps(romics_object1,show_dates = F,show_details = F)=="log10transform")-
        sum(romicsSteps(romics_object1,show_dates = F,show_details = F)=="unlog10data")) %% 2 == 1 ){
      log1<-log1+10
    }

    if((sum(romicsSteps(romics_object2,show_dates = F,show_details = F)=="log10transform")-
        sum(romicsSteps(romics_object2,show_dates = F,show_details = F)=="unlog10data")) %% 2 == 1 ){
      log2<-log2+10
    }

if(log1!=log2){stop("The <romics_object1> and <romics_object2> were not logged the same way.")}
if(log1==log2& log1>10){stop("The <romics_object1> and <romics_object2> were not logged more than once.")}
if(log1==log2& log1==2){
  print("both romics_object were log2 transformed the combination conserved the log2 transfomration")
  romics_object$steps<-c(romics_object$steps,paste0("date|",gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),"|romics_log2"))
  romics_object$steps<-c(romics_object$steps, "fun|log2transform(romics_object=romics_object_old)")
  }

if(log1==log2& log1==10){
  print("both romics_object were log10 transformed the combination conserved the log10 transfomration")
  romics_object$steps<-c(romics_object$steps,paste0("date|",gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),"|romics_log10"))
  romics_object$steps<-c(romics_object$steps, "fun|log10transform(romics_object=romics_object_old)")
  }

    #check if the $statistic object already is a part of the Romics_objects (if so remove stats)
    if(!is.null(romics_object1$statistics)){
      print("The Statistics layer was present in the <romics_object1> it was removed")
      romics_object1<-romics_object1[!names(romics_object1)=="statistics"]
    }
    if(!is.null(romics_object2$statistics)){
      print("The Statistics layer was present in the <romics_object2> it was removed")
      romics_object2<-romics_object2[!names(romics_object2)=="statistics"]
    }

    #if needed change the main factor of romics_object2
    if(romics_object1$main_factor!=romics_object2$main_factor){
      print("The main factor of <romics_object1> and <romics_object2> were different, the main factor of <romics_object1> was conserved")
    }
    romics_object2<-romicsChangeFactor(romics_object2,main_factor = romics_object1$main_factor)
    romics_object1<-romicsChangeFactor(romics_object1,main_factor = romics_object1$main_factor)


    #check if metadata of the two romics_objects are the same
    meta1<-romics_object1$metadata
    meta2<-romics_object2$metadata

    if(!setequal(meta1,meta2)){stop("The metadata of the two romics objects were different.")}

    odata1<-romics_object1$original_data
    odata2<-romics_object2$original_data
    rownames(odata1)<-paste0(rownames(odata1),suffix1)
    rownames(odata2)<-paste0(rownames(odata2),suffix2)
    romics_object$original_data<-rbind(odata1,odata2)

    data1<-romics_object1$data
    data2<-romics_object2$data
    rownames(data1)<-paste0(rownames(data1),suffix1)
    rownames(data2)<-paste0(rownames(data2),suffix2)
    romics_object$data<-rbind(data1,data2)

    mdata1<-romics_object1$missingdata
    mdata2<-romics_object2$missingdata
    rownames(mdata1)<-paste0(rownames(mdata1),suffix1)
    rownames(mdata2)<-paste0(rownames(mdata2),suffix2)
    romics_object$missingdata<-rbind(mdata1,mdata2)

    romics_object$steps_old<-steps_old

    return(romics_object)
  }
