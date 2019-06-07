#######################################################################
#' Make input list for age data
#'
#' @param data list with input data
#' @param common list with common parameters
#' @export   
## Created 
make.ageList<-function(data,common){
  if(any(data$DataMatrix$part.year<=0 | data$DataMatrix$part.year>1))
    stop("Error: part.year must be >0 and <=1")
  a.vec<-common$minage:common$maxage
  nAges<-length(a.vec)
  if(common$cc)a.vec<-c(a.vec,a.vec)
  ncat<-length(a.vec)
  cell<-make.cell(data,ncat) 
  sortmat<-sort.mat(data$DataMatrix)
  if(common$age.error)
    A2A<-make.A2A(common$A2A,common$minage,common$maxage)
  else
    A2A<-NULL
  # Test the neighbourhood structure
  if(any(data$info[,"CAR"]==1)){
    narea<-length(data$CARNeighbours$numNeighbours)
    amat<-matrix(0,nrow=narea,ncol=narea)
    nn<-cumsum(c(0,data$CARNeighbours$numNeighbours))
    id<-data$CARNeighbours$idNeighbours
    for(i in 1:narea){
      use<-(nn[i]+1):nn[i+1]
      a<-id[use]
      amat[i,a]<-1
    }
    if(sum(amat!=t(amat))!=0){
      stop("Error: The area structure is not symmetric")
    }
  }
  if(any(!is.finite(data$CovariateMatrix$haulcount)))
    stop("Error: Missing values in haulcount is not allowed")
  if(common$CC){
    ind <- data$DataMatrix$otolithtype%in%c(1,2,4,5,NA)
    if(any(ind==FALSE))
       stop("Error: Otolithtype different from 1, 2, 4 or 5 is not allowed when running coastal cod model")
  }
  list(sortmat=sortmat,cell=cell,delta.age=as.numeric(common$delta.age),
       a.vec=a.vec,A2A=A2A,
       age.errors=as.integer(common$age.error),
       num.adj.area=as.integer(data$CARNeighbours$numNeighbours),
       adj.area=as.integer(data$CARNeighbours$idNeighbours)
  )
}
#######################################################################
#' Make input list for lga data
#'
#' @param data list with input data
#' @param common list with common parameters
#' @export   
make.lgaList<-function(data,common){
  if(any(is.na(data$DataMatrix$lengthCM))){
    stop("Error: No missing lengths are allowed in DataMatrix")
  }
  if(common$lgamodel=="log-linear")
  {
    g.a.model = 0
    g.a.par.init = 0
  }
  else if(common$lgamodel=="non-linear")
  {
    g.a.model = 1
    g.a.par.init<-c(1,1,0.05)
  }
  cell<-make.cell(data,1) 
  list(g.a.model=as.integer(g.a.model),g.a.par.init=g.a.par.init,res=cell$res)
}
#######################################################################
#' Make input list for wgl data
#'
#' @param data list with input data
#' @param common list with common parameters
#' @export   
make.wglList<-function(data,common){
  if(any(is.na(data$DataMatrix$weightKG))){
    stop("Error: No missing weights are allowed in DataMatrix")
  }
  data$DataMatrix$samplingID<-as.integer(data$DataMatrix$samplingID)
  data$DataMatrix$otolithtype<-as.integer(data$DataMatrix$otolithtype)
  cell<-make.cell(data,1) 
  haul<-data$DataMatrix$samplingID
  nFishBoat<-table(haul)
  DataMatrix<-data$DataMatrix
  for (icol in 1:ncol(DataMatrix))DataMatrix[,icol]<-code.na(DataMatrix[,icol],-99999)
  list(cell=cell,nFishBoat=as.integer(nFishBoat),
       ncol.data=as.integer(ncol(data$DataMatrix)),DataMatrix=DataMatrix,
       num.adj.area=as.integer(data$CARNeighbours$numNeighbours),
       adj.area=as.integer(data$CARNeighbours$idNeighbours))
}
#######################################################################
#' Make input list for haulsize data
#'
#' @param ageList list with age data
#' @export   
make.hszList<-function(ageList){
  info<-ageList$cell$info
  cov<-ageList$cell$covmat
  ihsz<-ageList$cell$info$continuous==1
  hsz<-cov[,ihsz]
  cov<-cov[,!ihsz,drop=F]
  info$random<-as.integer(info$random[!ihsz])
  info$CAR<-as.integer(info$CAR[!ihsz])
  info$continuous<-as.integer(info$continuous[!ihsz])
  info$in.landings<-as.integer(info$in.landings[!ihsz])
  info$nlev<-as.integer(info$nlev[!ihsz])
  info$interaction<-as.integer(rep(0,sum(!ihsz)))
  info$in.slopeModel<-as.integer(info$in.slopeModel[!ihsz])
  list(hsz=hsz,n.col.cov=length(cov),covmat=cov,info=info)
}
#######################################################################
#' Sort matrix
#'
#' @param DataMatrix matrix with input data
#' @export   
# 07.02.2018 - new function due to table() is changed in new R-version
sort.mat<-function(DataMatrix){
  DataMatrix<-DataMatrix
  DataMatrix$age<-as.integer(DataMatrix$age)
  DataMatrix$otolithtype<-as.integer(DataMatrix$otolithtype)
  haul<-DataMatrix$samplingID
  no.age<-is.na(DataMatrix$age)
  order<-order(haul,no.age)
  DataMatrix<-DataMatrix[order,]
  haul<-haul[order]
  tab<-table(haul,DataMatrix$age,useNA="always")
  nc<-ncol(tab)
  num.noAge<-tab[-nrow(tab),nc]
  nFishBoat<-table(haul)
  haulstart<-match(unique(haul),haul)
  start.noAge<-haulstart+nFishBoat-num.noAge-1
  start.noAge[num.noAge==0]<-0
  for (icol in 1:ncol(DataMatrix))DataMatrix[,icol]<-code.na(DataMatrix[,icol],-99999)
  list(n.col.data=as.integer(ncol(DataMatrix)),DataMatrix=DataMatrix,num.noAge=as.integer(num.noAge),
       nFishBoat=as.integer(nFishBoat),start.noAge=as.integer(start.noAge))
}
#######################################################################
#' Create length interval information
#'
#' @export   
make.length.intervals<-function(data,resolution){
  r.len = range(data$lengthCM,na.rm=TRUE)
  n.int.len = 1+as.integer((r.len[2]-r.len[1])/resolution)
  int.len.vec<-seq(r.len[1],r.len[2],resolution)
  int.len.lim<-int.len.vec-0.5*resolution
  int.len.vec<-log(int.len.vec)
  int.len.lim<-c(log(int.len.lim),max(log(r.len[2]+100.0),99999.9))
  n.int.len<-length(int.len.vec)
  list(n.int.len=n.int.len,int.len.lim=int.len.lim,int.len.vec=int.len.vec)
}
#######################################################################
#' @export   
bootstrap.subsamples<-function(data){
  data<-data
  if(!is.null(data$partnumber)){
    tab<-table(data$partnumber,data$samplingID)
    tab<-matrix(tab,nrow=length(unique(data$partnumber)),dimnames=dimnames(tab))
    nsub<-colSums(tab>0)
    multi.part<-colnames(tab)[nsub>1]
    if(length(multi.part)>0){
      for(name in multi.part){
        if(sum(is.na(data$partcount[data$samplingID==name]))==0)
          data[data$samplingID==name,]<-resample.ss(data[data$samplingID==name,])
        else data<-data[data$samplingID!=name,]}
    }
  }
  data
}
#######################################################################
#' @export   
resample.ss<-function(data){
  ## 23.05.2018 - Changed catchnumbers to partcount, and subsample to partnumber
  ## catchweight - corresponds to old FVEKT?
  if(is.matrix(data)|is.data.frame(data)){
    data$partcount[is.na(data$partcount)]<-0
    new.index<-sample(sum(data$partcount>0),replace=T,prob=data$partcount)
    newdata<-data[new.index,]
    if(!is.null(newdata$age)){## 23.05.2018 - Changed to sort only on length if age not available (in WeighLength)
      order<-order(newdata$age,newdata$lengthCM)
    } else {
      order<-order(newdata$lengthCM)
    }
    newdata<-newdata[order,]
  }
  else newdata<-data
  newdata$partcount<-sum(data$partcount[match(unique(newdata$partnumber),data$partnumber)])
#  newdata$catchweight<-sum(data$catchweight[match(unique(newdata$subsample),data$subsample)])
  newdata
}

#######################################################################
#' @export   
make.continuous.age<-function (d,m,y,age){
  dom<-c(0,31,28,31,30,31,30,31,31,30,31,30,31)
  cumd<-cumsum(dom)
  div<-365
  leapyear<-floor(y/4)==y/4
  div<-365+leapyear
  doy<-cumd[m]
  doy<-doy[leapyear&m>2]<-doy+1
  cont.age<-age+doy/div
  cont.age
}

#######################################################################
#' @export   
tab.order<-function(x){
  table(x)[match(unique(as.character(x)),names(table(x)))]
}
#######################################################################
#' @export   
tapply.order<-function(y,x){
  tapply(y,x,FUN=sum)[match(unique(as.character(x)),names(table(x)))]
}
#######################################################################
#' @export   
add.cell.info<-function(info,ncell){
  mat<-info[,c("random","CAR","continuous","in.landings","nlev","interaction","in.slopeModel")]
  rbind(mat,cell=c(1,0,0,0,ncell,0,0))
}
#######################################################################
#' @export   
make.info.list<-function(info){
  list(random=as.integer(info[,"random"]),CAR=as.integer(info[,"CAR"]),
       continuous=as.integer(info[,"continuous"]),
       in.landings=as.integer(info[,"in.landings"]),nlev=as.integer(info[,"nlev"]),
       interaction=as.integer(info[,"interaction"]),in.slopeModel=as.integer(info[,"in.slopeModel"]))
}
#######################################################################
#' @export   
make.cell<-function(data,ncat){
  info<-data$info
  covmat<-data$CovariateMatrix
  cellvar<-info[,"interaction"]
  if(sum(cellvar)>1) {
    nlev<-info[,"nlev"][cellvar==1]
    real.cell<-enumerate.cell(data$CovariateMatrix[,cellvar==1],nlev)
    cell<-as.integer(as.factor(real.cell))
    res = make.constr.fit(ncat,nlev,sort(unique(cell)))
  ### CHANGED 12/1-18
  ### Use all combination of cells - real.cell
  ##covmat$cell<-cell 
  ##info<-add.cell.info(info,max(cell))
  ## res not used anymore?
    covmat$cell <- as.integer(real.cell)
    info<-add.cell.info(info,prod(nlev))
  } else {
    res = list(cell=as.integer(0),Sigma=matrix(1),constr=matrix(1),n.constr.cell=as.integer(0))
    real.cell<-NULL
  }
  info<-make.info.list(info)  
  covmat<-format.covmat(covmat,info)
  list(n.col.cov=length(covmat),covmat=covmat,res=res,info=info,real.cell=real.cell)
}
######################################################################
#' @export   
format.covmat<-function(covmat,info){
  covmat<-covmat
  for (i in 1:length(covmat)){
    if(info$continuous[i]==1)covmat[[i]]<-as.double(covmat[[i]])
    else covmat[[i]]<-as.integer(covmat[[i]])
  }
  covmat
}
######################################################################
#' @export   
make.CCerror<-function(CCerror,nAges,CCerrorList){
  newCCerror<-vector("list")
  newCCerror$CCerror<-CCerror
  if(CCerror==0)
  {
    newCCerror$ptype1.CC<-newCCerror$ptype2.CC<-newCCerror$ptype4.S<-newCCerror$ptype5.S<-rep(1,nAges)
    newCCerror$ptype1.S<-newCCerror$ptype2.S<-newCCerror$ptype4.CC<-newCCerror$ptype5.CC<-rep(0,nAges)
  }
  else
  {
    newCCerror$ptype1.CC = CCerrorList$ptype1.CC
    newCCerror$ptype1.S  = CCerrorList$ptype1.S
    newCCerror$ptype2.CC = CCerrorList$ptype2.CC
    newCCerror$ptype2.S  = CCerrorList$ptype2.S
    newCCerror$ptype4.CC = CCerrorList$ptype4.CC
    newCCerror$ptype4.S  = CCerrorList$ptype4.S
    newCCerror$ptype5.CC = CCerrorList$ptype5.CC
    newCCerror$ptype5.S  = CCerrorList$ptype5.S
    if(length(CCerrorList$ptype1.CC)==1)newCCerror$ptype1.CC<-rep(CCerrorList$ptype1.CC,nAges)
    if(length(CCerrorList$ptype1.S)==1)newCCerror$ptype1.S<-rep(CCerrorList$ptype1.S,nAges)
    if(length(CCerrorList$ptype2.CC)==1)newCCerror$ptype2.CC<-rep(CCerrorList$ptype2.CC,nAges)
    if(length(CCerrorList$ptype2.S)==1)newCCerror$ptype2.S<-rep(CCerrorList$ptype2.S,nAges)
    if(length(CCerrorList$ptype4.CC)==1)newCCerror$ptype4.CC<-rep(CCerrorList$ptype4.CC,nAges)
    if(length(CCerrorList$ptype4.S)==1)newCCerror$ptype4.S<-rep(CCerrorList$ptype4.S,nAges)
    if(length(CCerrorList$ptype5.CC)==1)newCCerror$ptype5.CC<-rep(CCerrorList$ptype5.CC,nAges)
    if(length(CCerrorList$ptype5.S)==1)newCCerror$ptype5.S<-rep(CCerrorList$ptype5.S,nAges)
  }
  newCCerror
}
#######################################################################
#' @export   
make.A2A<-function(A2A,ageMin,ageMax){
  if(!is.null(A2A)){
    mat.ages<-as.integer(rownames(A2A))
    if(is.null(mat.ages))mat.ages<-0:(nrow(A2A)-1)
    if(ageMax>max(mat.ages)){
      nextra<-ageMax-max(mat.ages)
      A2A<-cbind(A2A,matrix(0,nrow=nrow(A2A),ncol=nextra))
      A2A<-rbind(A2A,cbind(matrix(0,ncol=nrow(A2A),nrow=nextra),diag(nextra)))
      mat.ages<-c(min(mat.ages):ageMax)
    }
    if(ageMin<min(mat.ages)){
      nextra<-min(mat.ages)-ageMin
      A2A<-cbind(matrix(0,nrow=nrow(A2A),ncol=nextra),A2A)
      A2A<-rbind(cbind(diag(nextra),matrix(0,ncol=nrow(A2A),nrow=nextra)),A2A)
      mat.ages<-ageMin:max(mat.ages)
    }
    rownames(A2A)<-colnames(A2A)<-mat.ages
    if(ageMin>min(mat.ages)){
      A2A<-rbind(colSums(matrix(A2A[mat.ages<=ageMin,],ncol=ncol(A2A))),
                 matrix(A2A[mat.ages>ageMin,],ncol=ncol(A2A)))
      A2A<-A2A[,mat.ages>=ageMin]
      mat.ages<-ageMin:max(mat.ages)}
    if(ageMax<max(mat.ages)){
      A2A<-rbind(matrix(A2A[mat.ages<ageMax,],ncol=ncol(A2A)),
                 colSums(matrix(A2A[mat.ages>=ageMax,],ncol=ncol(A2A))))
      A2A<-A2A[,mat.ages<=ageMax]}
    rownames(A2A)<-colnames(A2A)<-ageMin:ageMax
    
  }
  A2A
}
#######################################################################
#' @export   
make.hsz.data<-function(stoxdata){
  stoxdata<-stoxdata
  ihsz<-colnames(stoxdata$AgeLength$CovariateMatrix)=="haulsize"
  DataMatrix<-stoxdata$AgeLength$CovariateMatrix[,ihsz]
  HSZ<-stoxdata$AgeLength
  HSZ$CovariateMatrix<-HSZ$CovariateMatrix[,!ihsz]
  HSZ$random<-HSZ$random[!ihsz]
  HSZ$in.catch<-HSZ$in.catch[!ihsz]
  HSZ$spatial<-HSZ$spatial[!ihsz]
  HSZ$continuous<-HSZ$continuous[!ihsz]
  HSZ$nlev<-HSZ$nlev[!ihsz]
  HSZ$in.slopeModel<-HSZ$in.slopeModel[!ihsz]
  HSZ$interaction<-HSZ$interaction[!ihsz]
  stoxdata$HSZ<-HSZ
  stoxdata$DataMatrix<-DataMatrix
  stoxdata
}
#########################################################################################
#' @export   
make.fixed<-function(nsamples,lganew=NULL,lganew.cc=NULL,wglnew=NULL,wglnew.cc=NULL){
  lga<-list(fixed.int.cc=NULL,fixed.slp.cc=NULL,fixed.tau.cc=NULL,fixed.g.a.gamma.cc=NULL,
            fixed.int=NULL,fixed.slp=NULL,fixed.tau=NULL,fixed.g.a.gamma=NULL)
  wgl<-list(fixed.int.cc=NULL,fixed.slp.cc=NULL,fixed.tau.cc=NULL,
            fixed.int=NULL,fixed.slp=NULL,fixed.tau=NULL)
  if(!is.null(lganew)){
    lga$fixed.int=rep(lganew[,1],length=nsamples)
    lga$fixed.slp=rep(lganew[,2],length=nsamples)
    lga$fixed.tau=1/rep(lganew[,3],length=nsamples)^2
    if(ncol(lganew)==4)lga$fixed.g.a.gamma=rep(lganew[,4],length=nsamples)
    if(!is.null(lganew.cc)){
      lga$fixed.int.cc=rep(lganew.cc[,1],length=nsamples)
      lga$fixed.slp.cc=rep(lganew.cc[,2],length=nsamples)
      lga$fixed.tau.cc=1/rep(lganew.cc[,3],length=nsamples)^2
      if(ncol(lganew.cc)==4)lga$fixed.g.a.gamma.cc=rep(lganew.cc[,4],length=nsamples)
    }}
  if(!is.null(wglnew)){
    wgl$fixed.int=rep(wglnew[,1],length=nsamples)
    wgl$fixed.slp=rep(wglnew[,2],length=nsamples)
    wgl$fixed.tau=1/rep(wglnew[,3],length=nsamples)^2
    if(!is.null(wglnew.cc)){
      wgl$fixed.int.cc=rep(wglnew.cc[,1],length=common$burnin+common$thin*common$nsamples)
      wgl$fixed.slp.cc=rep(wglnew.cc[,2],length=common$burnin+common$thin*common$nsamples)
      wgl$fixed.tau.cc=1/rep(wglnew.cc[,3],length=common$burnin+common$thin*common$nsamples)^2
    }}
  list(lga=lga,wgl=wgl)
}
#######################################################################
#' @export   
make.priorList<-function(ageList,wglList,common){
  a.vec<-common$minage:common$maxage
  nAges<-length(a.vec)
  ncat<-nAges
  if(common$cc)ncat<-ncat*2
  
  age.priors<-make.prior(ageList,ncat,F)
  lga.priors<-make.prior(ageList,1,T)
  wgl.priors<-make.prior(wglList,1,T)
  list(
    age.eff.mean=age.priors$int$mean,
    age.eff.prec=age.priors$int$prec,
    age.prec.par=age.priors$int$rand.prec,
    age.ar=age.priors$int$ar,
    
    lga.eff.mean=c(lga.priors$int$mean,lga.priors$slp$mean),
    lga.eff.prec=c(lga.priors$int$prec,lga.priors$slp$prec),
    lga.prec.par=c(lga.priors$int$rand.prec,lga.priors$slp$rand.prec),
    lga.ar=c(lga.priors$int$ar,lga.priors$slp$ar),
    
    wgl.eff.mean=c(wgl.priors$int$mean,wgl.priors$slp$mean),
    wgl.eff.prec=c(wgl.priors$int$prec,wgl.priors$slp$prec),
    wgl.prec.par=c(wgl.priors$int$rand.prec,wgl.priors$slp$rand.prec),
    wgl.ar=c(wgl.priors$int$ar,wgl.priors$slp$ar)
  )
}

#######################################################################
#' @export   
make.prior<-function(data,nages,fish){
  info<-data$cell$info
  int<-slp<-vector("list")
  int$mean<-rep(0,sum(info$nlev*(1-info$random)*nages))
  int$prec<-rep(0.01,sum(1-info$random))
  int$prec[1]<-0.0000001
  int$rand.prec<-rep(0.0001,2*sum(info$random))
  int$rand.prec<-c(int$rand.prec,0.0001,0.0001)#haul effect
  int$ar<-c(1,1)
  if(fish){
    use<-info$in.slopeModel==1
    slp$mean<-rep(0,sum(info$nlev[use]*(1-info$random[use])*nages))
    slp$prec<-rep(0.01,sum(1-info$random[use]))
    slp$prec[1]<-0.0000001
    slp$rand.prec<-rep(0.0001,2*sum(info$random[use]))
    slp$ar<-c(1,1)
    int$rand.prec<-c(int$rand.prec,0.01,0.01)#fish effect
  }
  list(int=int,slp=slp)
}
#######################################################################
#' @export   
make.numpar<-function(ageList,wglList,common){
  a.vec<-common$minage:common$maxage
  nAges<-length(a.vec)
  ncat<-nAges
  if(common$cc)ncat<-ncat*2
  age.numpar<-calc.numpar(ageList,T,ncat)
  lga.numpar<-calc.numpar(ageList,F,1)
  wgl.numpar<-calc.numpar(wglList,F,1)
  if(common$lgamodel=="log-linear")lgarelnpar = 0
  if(common$lgamodel=="non-linear")lgarelnpar = 3
  if(common$cc){
    wgl.numpar = c(wgl.numpar,wgl.numpar)
    numpar1=c(age.numpar,lga.numpar,lgarelnpar,lga.numpar,lgarelnpar)
    if(common$CCerror)numpar1[4]<-numpar1[4]+2
  }
  else{
    wgl.numpar = c(wgl.numpar,0)
    numpar1 = c(age.numpar,lga.numpar,lgarelnpar,0,0)}
  list(numpar1=numpar1,wgl.numpar=wgl.numpar)
}

#######################################################################
#' @export   
calc.numpar<-function(data,age,ncat){
  #note - assumes covmat includes column of 1s as constant term (nlev=1)
  info<-data$cell$info
  numpar<-ncat*sum(info$nlev)+sum(info$random)
  numpar<-numpar+1 #loglikelihood
  numpar<-numpar+1 #haul precision
  numpar<-numpar+sum(info$CAR)
  if(!age){
    numpar<-numpar+length(data$cell$covmat[[1]]) #hauls
    numpar<-numpar+1 #obs precision
    numpar<-numpar+sum(info$in.slopeModel*info$nlev)+
      sum(info$in.slopeModel*info$random)
    numpar<-numpar-sum(info$continuous)#haul size not included in lga or wgl
  }
  numpar
}
#######################################################################
#' @export   
code.na<-function(x,code){
  type<-class(x)
  x[is.na(x)]<-code
  class(x)<-type
  x
}

#######################################################################
#' @export   
write.bin.C.list<-function(listname,dir,fName=NULL){
  if(is.null(fName))
    fileName <- paste0(dir,"/",deparse(substitute(listname)))
  else
    fileName <- paste0(dir,"/",fName)
  zz<-file(fileName,"wb")
  cat(paste("Write to binary file ",fileName,"\n",sep=""))
  nuse<-count.list(listname,0)
  writeBin(as.integer(nuse),zz)
  write.list(listname,zz)
  close(zz)
}
#######################################################################
#' @export   
count.list<-function(listname,nuse){
  nuse<-nuse
  for(i in 1:length(listname)){
    if(length(listname[[i]])>0 & is.element(typeof(listname[[i]]),c("integer","double","character","logical")))
      nuse<-nuse+1
    if(length(listname[[i]])>0 & typeof(listname[[i]])=="list"){
      nuse<-count.list(listname[[i]],nuse)}
  }
  nuse
}
#######################################################################
#' @export   
write.list<-function(listname,zz){
  names<-names(listname)
  for(i in 1:length(listname)){
    xx<-listname[[i]]
    if(is.matrix(xx))xx<-as.vector(xx)
    lx<-length(xx)
    if(lx>0){
      type<-typeof(xx)
      if(type=="list")write.list(xx,zz)   
      else{
        #print(c(names[i],type))
        writeBin(nchar(names[i]),zz)
        writeChar(names[i],zz,nchars=nchar(names[i]),eos=NULL)
        writeBin(lx,zz)
        if(type=="integer"){
          writeBin(as.integer(0),zz) #type integer
          writeBin(xx,zz)
        } else if(type=="double"){
          writeBin(as.integer(1),zz) # type double
          writeBin(xx,zz)
        } else if(type=="character"){
          writeBin(as.integer(2),zz) # type character
          for(ilength in 1:lx){
            writeBin(nchar(xx[ilength]),zz)
            writeChar(xx[ilength],zz,nchars=nchar(xx[ilength]),eos=NULL)}
        } else if(type=="logical"){
          writeBin(as.integer(0),zz) #type integer
          xx<-as.integer(xx)
          writeBin(xx,zz)
        }     }
    }}
}
#######################################################################
#######################################################################
#' @export   
write.ascii.C.list<-function(listname,dir,fName=NULL){
  if(is.null(fName))
    fileName <- paste0(dir,"/",deparse(substitute(listname)))
  else
    fileName <- paste0(dir,"/",fName)
  zz<-file(fileName,"wt")
  cat(paste("Write to ascii file ",fileName,"\n",sep=""))
  nuse<-count.list(listname,0)
  write(as.integer(nuse),zz,append=TRUE)
  write.ascii.list(listname,zz)
  close(zz)
}
#' @export   
write.ascii.list<-function(listname,zz){
  names<-names(listname)
  for(i in 1:length(listname)){
    xx<-listname[[i]]
    if(is.matrix(xx))xx<-as.vector(xx)
    lx<-length(xx)
    if(lx>0){
      type<-typeof(xx)
      if(type=="list")write.ascii.list(xx,zz)   
      else{
        #print(c(names[i],type))
        #write(nchar(names[i]),zz,append=TRUE)
        write(names[i],zz,append=TRUE)
        write(lx,zz,append=TRUE)
        if(type=="integer"){
          #write(as.integer(0),zz,append=TRUE) #type integer
          write(xx,zz,append=TRUE)
        } else if(type=="double"){
          #write(as.integer(1),zz,append=TRUE) # type double
          write(xx,zz,append=TRUE)
        } else if(type=="character"){
          #write(as.integer(2),zz,append=TRUE) # type character
          for(ilength in 1:lx){
            #write(nchar(xx[ilength]),zz,append=TRUE)
            write(xx[ilength],zz,append=TRUE)}
        } else if(type=="logical"){
          #write(as.integer(0),zz,append=TRUE) #type integer
          xx<-as.integer(xx)
          write(xx,zz,append=TRUE)
        }     }
    }}
}
