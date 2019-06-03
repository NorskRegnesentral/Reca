#' Estimate catch-at-age model
#'
#' 
#' @param stoxdata list with input data
#' @param common list with common parameters
#' @param win indicator variable, =1 if windows platform, 0 otherwise
#' @keywords c-program mcmc
#' @export
## Created 30.10.2017
## Modified
##
run.fit <- function(stoxdata,common,win)
{
  stoxdata$AgeLength$DataMatrix<-bootstrap.subsamples(stoxdata$AgeLength$DataMatrix)
  ageList<-make.ageList(stoxdata$AgeLength,common)
  lgaList<-make.lgaList(stoxdata$AgeLength,common)
  stoxdata$WeightLength$DataMatrix<-bootstrap.subsamples(stoxdata$WeightLength$DataMatrix)
  wglList<-make.wglList(stoxdata$WeightLength,common)
  
  numpar<-make.numpar(ageList,wglList,common)
  
  fixed.par<-make.fixed(common$nSamples,lganew=NULL,lganew.cc=NULL,wglnew=NULL,wglnew.cc=NULL)

  inc.cts.var<-sum(ageList$cell$info$continuous*(1-ageList$cell$info$in.landings))>0
  if(inc.cts.var)hszList<-make.hszList(ageList)
  
  priorList<-make.priorList(ageList,wglList,common)
  m.i.n<-make.length.intervals(stoxdata$AgeLength$DataMatrix,common$resolution)
  a.vec<-common$minage:common$maxage
  nAges<-length(a.vec)

  common.par<-list(seed=as.integer(common$seed),
                   mcmc.par=as.integer(c(common$burnin,common$thin,common$nSamples)),
                   num.par.model1=as.integer(numpar$numpar1),
                   num.par.model2=as.integer(numpar$wgl.numpar),
                   usedebug=common$usedebug,
                   sim.ar=as.integer(common$sim.ar),
                   print.boat=as.integer(common$print.boat),
                   inputfolder=common$inputfolder,
                   filename.mcmc1=common$filename.mcmc1,
                   filename.mcmc2=common$filename.mcmc2,
                   filename.hsz.mcmc2=common$filename.hsz.mcmc2,
                   filename.hsz.it=common$filename.hsz.it,
                   filename.hsz.hauleff=common$filename.hsz.hauleff,
                   coastal.cod=as.integer(common$cc), # Moved from ageList in order to appear in both age/lga and wgl model
                   CCerror=make.CCerror(common$CCerror,nAges,stoxdata$AgeLength$CCerrorList),
                   m.i.n=m.i.n,
                   inc.haulsize=as.integer(inc.cts.var))

  if(is.null(common$print.format))
    common.par$print.format <- as.integer(0)  #print.format =0 (binary), 1 (ascii - only for testing)
  else
    common.par$print.format <- as.integer(common$print.format)
  if(is.null(common$old.version))
    common.par$old.version <- as.integer(0) # only for testing!!
  else
    common.par$old.version <- as.integer(common$old.version)     
  
  priorList<<-priorList
  common.par<<-common.par
  ageList<<-ageList
  lgaList<<-lgaList
  wglList<<-wglList
  if(inc.cts.var)hszList<<-hszList
  
  write.ascii.C.list(common.par,common$inputfolder,"common_par_fit_ascii")  
  #write.bin.C.list(common.par,common$inputfolder,"common_par_fit")
  write.bin.C.list(priorList,common$inputfolder,"priorList")
  write.bin.C.list(ageList,common$inputfolder,"stoxdata_age")
  write.bin.C.list(lgaList,common$inputfolder,"stoxdata_lga")
  write.bin.C.list(wglList,common$inputfolder,"stoxdata_wgl")
  if(inc.cts.var)write.bin.C.list(hszList,common$inputfolder,"stoxdata_hsz")
  if(win){
    caa.call <- system.file("bin/caa_main_model1.exe",package="Reca")
  } else {
    caa.call <- system.file("bin/caa_main_model1",package="Reca")
  }
  s <- system(paste(shQuote(caa.call)," ",common$inputfolder,sep=""), intern=TRUE)
  filename <- paste0(common$inputfolder,"log.txt")
  if(!file.exists(filename)) stop("Error in fitting age and lga model")
  fp <- file(filename,"r") 
  log <- readLines(fp)
  close(fp)
  if(log[length(log)]=="OK")cat("Success fitting age and lga model\n")
  else stop("Error fitting age and lga model")
  file.remove(filename)  # To make sure that you don't read an old file
  
  if(1){
    # common.par$inc.haulsize<-0  Not needed here?
    if(win){
      caa.call <- system.file("bin/caa_main_model2.exe",package="Reca")
    } else {
      caa.call <- system.file("bin/caa_main_model2",package="Reca")
    }
    s <- system(paste(shQuote(caa.call)," ",common$inputfolder,sep=""), intern=TRUE)
    filename <- paste0(common$inputfolder,"log.txt")
    if(!file.exists(filename)) stop("Error in fitting wgl model")
    fp <- file(filename,"r")
    log <- readLines(fp)
    close(fp)
    if(log[length(log)]=="OK")cat("Success fitting wgl model\n")
    else stop("Error fitting wgl model")
    file.remove(filename)
  }
}

