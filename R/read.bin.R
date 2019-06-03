#' Read binary file from prediction with mcmc samples
#'
#' @param filename filename (with correct path) of binary file to read mcmc results from estimation
#' @export
## Created 01.11.2017
## Modified
## 11.12.2017    HWR   Changed names in result lists
read.predict.bin<-function(filename.predict){
  fp<-file(filename.predict,"rb")
  nMCMC <- readBin(fp,integer(),n=1,endian="little")
  nAges <- readBin(fp,integer(),n=1,endian="little")
  AgeCategories <- readBin(fp,integer(),n=nAges,endian="little")
  nlint <- readBin(fp,integer(),n=1,endian="little")
  LengthIntervalsLog <- readBin(fp,double(),n=nlint,endian="little")
  TotalCount <- array(dim=c(nlint,nAges,nMCMC))
  MeanLength <- MeanWeight <- array(dim=c(nAges,nMCMC))                
  for(i in 1:nMCMC){for(j in 1:nAges){
    TotalCount[,j,i]<-readBin(fp,double(),n=nlint,endian="little")}
    MeanLength[,i]<-readBin(fp,double(),n=nAges,endian="little")
    MeanWeight[,i]<-readBin(fp,double(),n=nAges,endian="little")
  }
  close(fp)
  list(TotalCount=TotalCount, MeanLength=MeanLength, MeanWeight=MeanWeight,
       AgeCategories=AgeCategories, LengthIntervalsLog=LengthIntervalsLog)
}

#########################################################################################
#' Read binary file from estimation with mcmc samples
#'
#' @param filename filename (with correct path) of binary file to read mcmc results from estimation
#' @param stoxdata list with input data
#' @export
## Created 01.11.2017
## Modified
## 11.12.2017    HWR   Changed names in result lists
read.fit.bin<-function(filename.mcmc1,filename.mcmc2,filename.hsz,stoxdata){
  info1<-stoxdata$AgeLength$info
  info2<-stoxdata$WeightLength$info
  fp1 <- file(filename.mcmc1,"rb")
  par1<-read.params.1(fp1)
  age<-make.fit.structure.age(par1$age.cov,info1,par1$nMCMC,par1$print.boat,stoxdata$GlobalParameters$CCerror)
  lga<-make.fit.structure.lgawgl(par1$lga.cov,info1,par1$nMCMC)
  if(par1$ga.model){
    lga$nonlin <- vector(length=par1$nMCMC,mode="numeric")
  }
  if(par1$cc){
    lga.cc<-make.fit.structure.lgawgl(par1$lga.cov.cc,info1,par1$nMCMC)
    if(par1$ga.model){
      lga.cc$nonlin <- vector(length=par1$nMCMC,mode="numeric")
    }
  }
  for(i in 1:par1$nMCMC){
    age<-(read.mcmc.it(fp1,i,age,par1$age.cov,par1$print.boat))$structure
    lga<-(read.mcmc.it(fp1,i,lga,par1$lga.cov,par1$print.boat))$structure
    if(par1$ga.model){
      tmp<-readBin(fp1,double(),n=3,endian="little")
      lga$nonlin[i] <- tmp[3]
    }
    if(par1$cc){
      lga.cc<-(read.mcmc.it(fp1,i,lga.cc,par1$lga.cov.cc,par1$print.boat))$structure
      if(stoxdata$GlobalParameters$CCerror){
        tmp<-readBin(fp1,double(),n=2,endian="little")
        age$CCerror$kC[i] <- tmp[1]
        age$CCerror$kA[i] <- tmp[2]
      }
      if(par1$ga.model){
        tmp<-readBin(fp1,double(),n=3,endian="little")
        lga.cc$nonlin[i] <- tmp[3]
      }
    }
  }
  # remove renormalisation when changing to continuous g-function
  #amin.log <- log(min(par1$age.vec)+1/12)
  #amax.log <- log(max(par1$age.vec)+1)
  #lga <- renorm.lga(lga,amin.log,amax.log)
  #if(par1$cc)
  #  lga.cc <- renorm.lga(lga.cc,amin.log,amax.log)
  close(fp1)
 
  fp2 <- file(filename.mcmc2,"rb")
  par2<-read.params.2(fp2)
  wgl<-make.fit.structure.lgawgl(par2$wgl.cov,info2,par2$nMCMC)
  if(par2$cc){
    wgl.cc<-make.fit.structure.lgawgl(par2$wgl.cov.cc,info2,par2$nMCMC)
  }
  for(i in 1:par2$nMCMC){
    wgl<-(read.mcmc.it(fp2,i,wgl,par2$wgl.cov,par2$print.boat))$structure
    if(par2$cc){
      wgl.cc<-(read.mcmc.it(fp2,i,wgl.cc,par2$wgl.cov.cc,par2$print.boat))$structure
    }
  }
  close(fp2)

  res <- list(ProportionAtAge=age, LengthGivenAge=lga, WeightGivenLength=wgl)
  
  if(par1$cc)
    res$LengthGivenAgeCC <- lga.cc
  if(par2$cc)
    res$WeightGivenLengthCC <- wgl.cc

  inc.cts.var<-sum(stoxdata$AgeLength$info[,"continuous"]*(1-stoxdata$AgeLength$info[,"in.landings"]))>0
  if(inc.cts.var){
    fp3 <- file(filename.hsz,"rb")
    par3 <- read.params.2(fp3)
    hsz <- make.fit.structure.lgawgl(par3$wgl.cov,info1,par3$nMCMC)
    for(i in 1:par3$nMCMC){
      hsz <-(read.mcmc.it(fp3,i,hsz,par3$wgl.cov,0,1))$structure
    }
    res$HaulsizeModel <- hsz
    close(fp3)
  }
  
  return(res)
}

#########################################################################################
#' @export
read.mcmc.it<-function(fp,it,structure,cov,print.boat,hsz=0){
  ## Read linear effects
  for(i in 1:cov$nxcov){
    nm <- names(structure[[i]]$cov)
    ncov<-cov$ncov[i]
    if(i==1){
      nlev<-cov$int.nlev
      if(cov$nxcov==1){ # Haul effects are not printed for the age model in the mcmc file
        ncov<-ncov-1
        nlev<-nlev[-length(nlev)]
      }
    } else {
      nlev<-cov$slp.nlev
    }    
    for(a in 1:cov$ncat){
      for(j in 1:ncov){
        if(nm[j]=="boat" & print.boat==1){
          structure[[i]]$cov[[j]][a,,it]<-readBin(fp,double(),n=nlev[j],endian="little") 
        } else if(nm[j]=="boat" & print.boat==0){ #No boat effect printed in the mcmc file
          structure[[i]]$cov[[j]][a,,it]<-NA
        } else {
          structure[[i]]$cov[[j]][a,,it]<-readBin(fp,double(),n=nlev[j],endian="little")
        }
      }
    }
    if(hsz & (i==1)){ 
      eff_hsz<-readBin(fp,double(),n=cov$nhaul,endian="little")
    }
  }
  ## Read ar-coef
  if(cov$int.ispat==1){
    structure[[1]]$CAR[[1]][it]<-readBin(fp,double(),n=1,endian="little")
  }
  ## Read precisions for random effects
  ncov<-cov$ncov[1]
  if(cov$nxcov==1) ncov <- ncov-1
  ind <- 1
  for(j in 1:ncov){
    if(cov$int.fix[j]==0){
      structure[[1]]$tau[[ind]][it]<-readBin(fp,double(),n=1,endian="little")
      ind <- ind+1
    }
  }
  ## Read observation precision 
  if(cov$nxcov==1)
    structure[[1]]$tau$catchSample[it]<-readBin(fp,double(),n=1,endian="little")
  else
    structure[[1]]$tau$fish[it]<-readBin(fp,double(),n=1,endian="little")
  ## Read loglikelihood 
  structure$LogLikelihood[it]<-readBin(fp,double(),n=1,endian="little")

  list(structure=structure)
}

#########################################################################################
#' @export
make.names <- function(info,age=T,slope=F,hsz=F,print.boat=F){
  use <- T
  if(slope) use<-use&info[,"in.slopeModel"]==1
  if(!age) use<-use&info[,"continuous"]!=1
  if(age&!print.boat) use<-use&rownames(info)!="boat"
  names.cov <- rownames(info)[use]
  names.tau <- NULL
  use.tau <- use&info[,"random"]>0
  if(any(use.tau==TRUE))
    names.tau <- rownames(info)[use.tau]
  if(sum(info[use,"interaction"])>1){
    names.cov <- c(names.cov,"cell")
    names.tau <- c(names.tau,"cell")
  }
  if(!age&!slope){
    names.cov <- c(names.cov,"catchSample")
    names.tau <- c(names.tau,"catchSample")
  }
  #nloop <- length(names.cov)
  names.CAR <- NULL
  if(sum(info[use,"CAR"])>0) names.CAR <- rownames(info)[info[use,"CAR"]==1]
  #nloop<-c(nloop,length(names.other))
  #if(length(random.names)>0)names.other<-c(names.other,paste0("tau.",random.names))
  #nloop<-c(nloop,length(random.names))
  list(names.cov=names.cov,names.tau=names.tau,names.CAR=names.CAR)
}

#########################################################################################
#' @export
make.fit.structure<-function(names,nlev,nage,nMCMC){
  cov <- vector("list",length=length(names$names.cov))
  names(cov) <- names$names.cov
  for(i in 1:length(nlev))
    cov[[i]] <- array(dim=c(nage,nlev[i],nMCMC))
  tau <- NULL
  if(!is.null(names$names.tau)){
    tau <- vector("list",length=length(names$names.tau))
    names(tau) <- names$names.tau
    for(i in 1:length(names$names.tau))
      tau[[i]] <- vector(length=nMCMC,mode="numeric")
  }
  CAR <- NULL
  if(!is.null(names$names.CAR)){
    spatial <- vector(length=nMCMC,mode="numeric")
    CAR <- list(spatial=spatial)
  }
  structure <- list(cov=cov,tau=tau,CAR=CAR)
  structure
}

#########################################################################################
#' @export
make.fit.structure.age <- function(params,info,nMCMC,print.boat,CCerror){
  nm <- make.names(info,age=T,slope=F,hsz=F,print.boat=print.boat)
  nm$names.tau <- c(nm$names.tau,"catchSample")
  structure <- make.fit.structure(nm,params$int.nlev[-length(params$int.nlev)],params$ncat,nMCMC)
  LogLikelihood <- vector(length=nMCMC,mode="numeric")
  if(CCerror){
    kC <- vector(length=nMCMC,mode="numeric")
    kA <- vector(length=nMCMC,mode="numeric")
    CCerror <- list(kC=kC,kA=kA)
    res <- list(Intercept=structure,LogLikelihood=LogLikelihood,CCerror=CCerror)
  } else {
    res <- list(Intercept=structure,LogLikelihood=LogLikelihood)
  }
  res
}

#########################################################################################
#' @export
make.fit.structure.lgawgl <- function(params,info,nMCMC){
  nm.int<-make.names(info,age=F)
  nm.slp<-make.names(info,age=F,slope=T)
  nm.int$names.tau <- c(nm.int$names.tau,"fish")
  structure.int <- make.fit.structure(nm.int,params$int.nlev,1,nMCMC)
  structure.slp <- make.fit.structure(nm.slp,params$slp.nlev,1,nMCMC)
  list(Intercept=structure.int, Slope=structure.slp)
}

#########################################################################################
#' @export
read.params.1<-function(fp){
  nMCMC <- readBin(fp,integer(),n=1,endian="little")
  num.par <- readBin(fp,integer(),n=5,endian="little")
  print.boat <- readBin(fp,integer(),n=1,endian="little")
  age.cov <- read.cov(fp)
  nage <- age.cov$ncat
  age.vec<-readBin(fp,integer(),n=nage,endian="little")
  delta.age<-readBin(fp,double(),n=1,endian="little")
  lga.cov <- read.cov(fp)
  ga.model<-readBin(fp,integer(),n=1,endian="little")
  cc<-readBin(fp,integer(),n=1,endian="little")
  lga.cov.cc<-NULL
  if(cc){
    lga.cov.cc <- read.cov(fp)
    tmp<-readBin(fp,integer(),n=1,endian="little") # ga.model(coastal cod)
  }
  list(nMCMC=nMCMC,num.par=num.par,print.boat=print.boat,age.cov=age.cov,age.vec=age.vec,delta.age=delta.age,
       lga.cov=lga.cov,ga.model=ga.model,cc=cc,lga.cov.cc=lga.cov.cc)
}

#########################################################################################
#' @export
read.params.2<-function(fp){
  nMCMC <- readBin(fp,integer(),n=1,endian="little")
  num.par <- readBin(fp,integer(),n=2,endian="little")
  print.boat <- readBin(fp,integer(),n=1,endian="little")
  wgl.cov <- read.cov(fp)
  cc<-readBin(fp,integer(),n=1,endian="little")
  wgl.cov.cc<-NULL
  if(cc){
    wgl.cov.cc <- read.cov(fp)
  }
  list(nMCMC=nMCMC,num.par=num.par,print.boat=print.boat,wgl.cov=wgl.cov,cc=cc,wgl.cov.cc=wgl.cov.cc)
}

#########################################################################################
#' @export
read.cov<-function(fp){
  ncat<-readBin(fp,integer(),n=1,endian="little")
  nhaul<-readBin(fp,integer(),n=1,endian="little")
  nxcov<-readBin(fp,integer(),n=1,endian="little")
  ncov<-rep(NA,nxcov)
  slp.nlev<-slp.fix<-slp.interaction<-slp.in.landings<-NULL
  for(i in 1:nxcov){
    ncov[i]<-readBin(fp,integer(),n=1,endian="little")
    if(i==1){
      int.nlev<-readBin(fp,integer(),n=ncov[i],endian="little")
      int.fix<-readBin(fp,integer(),n=ncov[i],endian="little")
      int.interaction<-readBin(fp,integer(),n=ncov[i],endian="little")
      int.in.landings<-readBin(fp,integer(),n=ncov[i],endian="little")
    } else {
      slp.nlev<-readBin(fp,integer(),n=ncov[i],endian="little")
      slp.fix<-readBin(fp,integer(),n=ncov[i],endian="little")
      slp.interaction<-readBin(fp,integer(),n=ncov[i],endian="little")
      slp.in.landings<-readBin(fp,integer(),n=ncov[i],endian="little")
    }
  }
  tmp<-readBin(fp,integer(),n=5,endian="little") # internal c-code information (ispat,ihaul,icell,iboat,ihaulsize)
  if(tmp[1]>0) # spatial covariate included
    ispat<-1
  else
    ispat<-(-1)
  list(ncat=ncat,nhaul=nhaul,nxcov=nxcov,ncov=ncov,int.nlev=int.nlev,int.fix=int.fix,
       int.interaction=int.interaction,int.in.landings=int.in.landings,int.ispat=ispat,
       slp.nlev=slp.nlev,slp.fix=slp.fix,slp.interaction=slp.interaction,slp.in.landings=slp.in.landings) 
}
#########################################################################################
#' @export
renorm.lga<-function(lga,amin.log,amax.log){
  lga <- lga
  lga$Intercept$cov$constant <- lga$Intercept$cov$constant-lga$Slope$cov$constant*amin.log/(amax.log-amin.log)
  lga$Slope$cov$constant <- lga$Slope$cov$constant/(amax.log-amin.log)
  lga
}



