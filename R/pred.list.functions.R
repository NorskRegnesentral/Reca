###################################################################
#' @export   
make.data.catch<-function(stoxdata,pc.age,pc.wgl){
  agecell<-sum(stoxdata$AgeLength$info[,"interaction"])>1
  wglcell<-sum(stoxdata$WeightLength$info[,"interaction"])>1
  agevars<-names(stoxdata$AgeLength$CovariateMatrix)[stoxdata$AgeLength$info[,"in.landings"]==1]
  wglvars<-names(stoxdata$WeightLength$CovariateMatrix)[stoxdata$WeightLength$info[,"in.landings"]==1]
  
  agecov<-stoxdata$Landings$AgeLengthCov[,agevars,drop=FALSE]
  wglcov<-stoxdata$Landings$WeightLengthCov[,wglvars,drop=FALSE]
  fac.lga.slp<-(1:(dim(agecov)[2]))[stoxdata$AgeLength$info[,"in.slopeModel"]==1]
  fac.wgl.slp<-(1:(dim(wglcov)[2]))[stoxdata$WeightLength$info[,"in.slopeModel"]==1]
  if(agecell)agecov<-cbind(agecov,agecell=pc.age$predict.cell) # NB! agecell must be the last column
  if(wglcell)wglcov<-cbind(wglcov,wglcell=pc.wgl$predict.cell) # NB! wglcell must be the last column
  nagevar<-dim(agecov)[2]
  nwglvar<-dim(wglcov)[2]
  factors<-cbind(agecov,wglcov)
  for(i in 1:ncol(factors)){  ## 25.10.2018: Must somehow have factors as integer - constant term was not!!!
    class(factors[,i]) <- "integer"
  }
  list(nCell=as.integer(length(stoxdata$Landings$LiveWeightKG)),                  
       nfactors=as.integer(dim(factors)[2]),
       fac.age.int=as.integer(1:nagevar), 
       fac.lga.int=as.integer(1:nagevar), 
       fac.lga.slp=as.integer(fac.lga.slp), 
       fac.wgl.int=as.integer(nagevar+(1:nwglvar)),
       fac.wgl.slp=as.integer(nagevar+fac.wgl.slp), 
       fac.hsz.int=as.integer(1:nagevar),
       n.col.cov=as.integer(ncol(factors)),
       factors=factors,      
       catch=as.double(stoxdata$Landings$LiveWeightKG),
       midseason=as.double(stoxdata$Landings$AgeLengthCov$midseason))
#       season=as.integer(stoxdata$Landings$AgeLengthCov$season))
}
###################################################################
#' @export   
make.dist.cell<-function(pc.age,pc.lga,pc.wgl,coastal.cod){
  
  num.cell.u1 = pc.age$cell.u.dist$n.u
  num.cell.u2 = c(pc.lga$cell.u.dist$n.u,0,pc.wgl$cell.u.dist$n.u,0,pc.lga$cell.u.dist$n.u)
  
  num.cell.o1 = pc.age$cell.u.dist$n.o
  num.cell.o2 = c(pc.lga$cell.o.dist$n.o,0,pc.wgl$cell.u.dist$n.o,0,pc.lga$cell.u.dist$n.o)
  num.cell.u<-c(num.cell.u1,num.cell.u2)
  num.cell.o<-c(num.cell.o1,num.cell.o2)
  if(coastal.cod){
    num.cell.u = c(num.cell.u,num.cell.u2)
    num.cell.o = c(num.cell.o,num.cell.o2)
  }
  
  
  list(
    num.cell.o=as.integer(num.cell.o),
    num.cell.u=as.integer(num.cell.u),
    age.int.E=as.double(t(pc.age$cell.u.dist$E)),
    lga.int.E=as.double(t(pc.lga$cell.u.dist$E)),
    lga.slp.E=NULL,
    wgl.int.E=as.double(t(pc.wgl$cell.u.dist$E)),
    wgl.slp.E=NULL,
    hsz.int.E=as.double(t(pc.lga$cell.u.dist$E)),
    age.int.C=as.double(t(pc.age$cell.u.dist$C)),
    age.int.nC=as.integer(dim(pc.age$cell.u.dist$C)[2]),
    lga.int.C=as.double(t(pc.lga$cell.u.dist$C)),
    lga.int.nC=as.integer(dim(pc.lga$cell.u.dist$C)[2]),
    lga.slp.C=NULL,
    lga.slp.nC=NULL,
    wgl.int.C=as.double(pc.wgl$cell.u.dist$C),
    wgl.int.nC=as.integer(dim(pc.wgl$cell.u.dist$C)[2]),
    wgl.slp.C=NULL,
    wgl.slp.nC=NULL,
    hsz.int.C=as.double(t(pc.lga$cell.u.dist$C)),
    hsz.int.nC=as.integer(dim(pc.lga$cell.u.dist$C)[2]))
}
###################################################################
