#' Estimate catch-at-age model
#'
#' This function runs the C fitting program and writes results to binary files.
#' @param AgeLength A list with elements connected to the age-length model, see Details
#' @param WeightLength A list with elements connected to the weight-length model, see Details
#' @param Landings A list with elements connected to landings, not used until prediction.
#' @param GlobalParameters A list with input parameters that are given by the user, see Details
#' @details
#' AgeLength
#' \itemize{
#'  \item DataMatrix  - matrix with one row for each fish and columns \emph{age} (age in whole years), \emph{part.year} 
#'  (fraction of year, >0 and <=1, i.e. realage=age+part.year, e.g. for seasonal sampling part.year=0.25 for season 1 etc, 
#'  for daily sampling part.year=1/365 for day 1 etc), \emph{lengthCM} (length in centimeters), \emph{samplingID} (id for 
#'  sampling unit, e.g. haul number, numbered from 1 to number of sampling units), \emph{partnumber}, \emph{partcount}, 
#'  \emph{otolithtype} (classification when multiple stocks (type1 (Coastal cod) and type2 (Atlantic cod)): 
#'  1-certain type1, 2-uncertain type1, 4-uncertain type2, 5-certain type2).
#'  \item CovariateMatrix  - matrix with columns for the different covariates, e.g. \emph{constant}, \emph{season}, 
#'  \emph{gearfactor}, \emph{spatial}, \emph{haulcount}, \emph{boat}. Each row corresponds to the samplingID \emph{in DataMatrix}.
#'  \item info - matrix with one row for each covariate and with columns \emph{random} (=1 if random effect, =0 otherwise), 
#'  \emph{CAR} (=0 if spatial effect, =0 otherwise), \emph{continuous} (=1 if continuous covariate, =0 otherwise), 
#'  \emph{in.landings}, \emph{nlev} (number of levels for the different covariates), \emph{interaction} 
#'  (=1 if including a cell effect/interaction term for the specified covariates), \emph{in.slopeModel}.
#'  \item CARNeighbours - list with \emph{numNeighbours} and \emph{idNeighbours}.
#'  \item AgeErrorMatrix - matrix that defines the probability of observing each age given the true age. Must be included if 
#'  uncertainty in ageing is included in the analysis, i.e. GlobalParameters$age.error=1. The true ages are the 
#'  columns, the observed ages are the rows and the number of rows and columns must equal the number of age categories defined
#'  by minage and maxage. The element in row i, column j defines the probability of a fish of true age j being observed as 
#'  having age i. The columns must sum to 1 (though the rows need not). The dimension is nAges x nAges, where nAges is the 
#'  number of age categories.
#'  \item CCerrorList - list with 8 elements of length 1 that defines the classification error:
#'  ptype1.CC: Probability of a certain coastal cod being classified as type 1.
#'  ptype1.S: Probability of a certain Atlantic cod being classified as type 1.
#'  ptype2.CC: Probability of an uncertain coastal cod being classified as type 2.
#'  ptype2.S: Probability of an uncertain Atlantic cod being classified as type 2.
#'  ptype4.CC: Probability of an uncertain coastal cod being classified as type 4.
#'  ptype4.S: Probability of an uncertain Atlantic cod being classified as type 4.
#'  ptype5.CC: Probability of a certain coastal cod being classified as type 5.
#'  ptype5.S: Probability of a certain Atlantic cod being classified as type 5.
#'  }
#'  WeightLength
#'  \itemize{
#'  \item DataMatrix  - matrix with one row for each fish and columns \emph{weightKG} (weight in kilograms), 
#'  \emph{lengthCM} (length in centimeters), \emph{samplingID} (id for sampling unit, e.g. haul number), 
#'  \emph{partnumber}, \emph{partcount}.
#'  \item CovariateMatrix  - matrix with columns for the different covariates (see AgeLength).
#'  \item info - matrix with one row for each covariate (see AgeLength).
#'  }
#'  GlobalParameters
#'  \itemize{
#'  \item nSamples  - number of MCMC samples that are saved
#'  \item thin  - number of MCMC samples run before each sample that is saved, e.g. thin=1 saves every sample, thin=10 saves every 10th sample
#'  \item burnin  - number of MCMC samples run and discarded before any samples are saved
#'  \item lengthresCM  - gives the length interval if catch-at-age by length is to be calculated, in centimeters
#'  \item maxlength - the maximum length of fish in the data set, in centimeters
#'  \item minage - minimum age category, used to exclude unreasonable data. The model will include ages from minage to maxage and
#'  any fish with observed ages outside this range will be omitted. Note that fish with small lengths but no observed age will be
#'  included even though they may be below minage
#'  \item maxage - maximum age category
#'  \item resultdir - directory path where the results are stored. 
#'  Two directories \strong{cfiles} and \strong{resfiles} are in this folder (created if they don't exist). 
#'  In \strong{cfiles} the temporary binary input files to the c-program are stored, and in \strong{resfiles} the
#'  binary output files are stored.
#'  \item fitfile - name of the output files from estimating the model. 
#'  The files are given the extensions \strong{.mcmc1} and \strong{.mcmc2} for the age and length-given-age model 
#'  and for the weight-given-length model, respectively.
#'  \item delta.age - default 0.001 - parameter used to improve estimation when there are ages with no observations 
#'  \item age.error - TRUE if uncertainty in ageing is included in the model, FALSE otherwise
#'  \item lgamodel - default "log-linear" - The length-given-age model will be log-linear unless specified 
#'  as "non-linear" in which case the Schnute-Richards model will be used
#'  \item CC - TRUE if coastal cod is included in the model, FALSE otherwise
#'  \item CCerror - TRUE if classification error is included in coastal cod model, FALSE otherwise
#'  \item seed - random seed value that initializes the random number generator. If NULL, a random number is generated for 
#'  this value. Otherwise, if the same seed is used and the other input values are the same, the results will be the same in 
#'  different simulations.
#'  }
#' @return A list with elements
#' \item{ProportionAtAge: }{a list with \strong{Intercept} and \strong{LogLikelihood} for the age model 
#' with all the mcmc samples for the different parameters included, see below}
#' \item{LengthGivenAge: }{a list with \strong{Intercept}, \strong{Slope} and \strong{LogLikelihood} 
#' for the length-given-age model, with all the mcmc samples for the different parameters included, see below}
#' \item{WeightGivenLength: }{a list with \strong{Intercept}, \strong{Slope} and \strong{LogLikelihood}
#' for the weight-given-length model, with all the mcmc samples for the different parameters included, see below.}
#' The \strong{Intercept} and \strong{Slope} lists consist of:
#' \itemize{
#' \item \strong{cov:} a list containing (depending on the selected model):
#'   \strong{constant} - \strong{season} - \strong{gearfactor} - \strong{spatial} - \strong{cell} - \strong{catchSample}:
#'    arrays of mcmc samples of the effects with dimensions #categories x #levels x #samples
#' \item \strong{tau:} list containing (depending on the selected model): 
#'   \strong{season} - \strong{gearfactor} - \strong{spatial} - \strong{cell} - \strong{catchSample} - \strong{fish}: 
#'   vectors with mcmc samples of precision parameters for random effects
#' \item \strong{CAR:} list with if spatial model  spatial.ar: mcmc samples of the AR-coefficient
#' }
#' If a multiple stocks (coastal cod/Atlantic cod) analysis is run, the list consists of
#' \item{ProportionAtAge: }{a list for the age model as described above, except that the age category dimension 
#' is 2 times the number of age categories(ncat). 
#' The first ncat values correspond to the Atlantic cod and the last ncat values correspond to the coastal cod}
#' \item{LengthGivenAge: }{a list for the length-given-age model as described above for Atlantic cod}
#' \item{LengthGivenAgeCC: }{a list for the length-given-age model as described above for coastal cod}
#' \item{WeightGivenLength: }{a list for the weight-given-length model as described above for Atlantic cod}
#' \item{WeightGivenLengthCC: }{a list for the weight-given-length model as described above for coastal cod}
#' The \strong{LogLikelihood} is a vector with the log-likelihood value for the different models in each sample
#' @seealso \link{eca.predict} for running predictions from the catch-at-age model.
#' @keywords c-program mcmc
#' @export
## Created 27.10.2017
## Modified
## 04.01.2018  HWR  added parameters in GlobalParameters
## 12.10.2018  HWR  changed to ascii files for common parameters
eca.estimate <- function(AgeLength,WeightLength,Landings,GlobalParameters)
{
  ## Create variable for specifying specific windows commands
  sysinfo <- Sys.info()
  if(sysinfo["sysname"]=="Windows"){
    win <- 1
  } else {
    win <- 0
  }

  stoxdata<-list(AgeLength=AgeLength,WeightLength=WeightLength,
                 Landings=Landings,GlobalParameters=GlobalParameters)

  if(is.null(GlobalParameters$seed)){
    GlobalParameters$seed <- sample(1:9999,1) #Generate a positive number
  }
  if(is.null(GlobalParameters$delta.age)){
    GlobalParameters$delta.age <- 0.001
  }
  if(is.null(GlobalParameters$lgamodel)){
    GlobalParameters$lgamodel <- "log-linear"
  }
  if(GlobalParameters$CC==FALSE){
    GlobalParameters$CCerror <- FALSE
  }
  
  ## Folder for input files
  cfilesfolder<-"cfiles/"
  inputfolder<-paste0(GlobalParameters$resultdir,"/",cfilesfolder)
  if(!dir.exists(inputfolder)){
    cat(paste("Create new directory",inputfolder,"\n"))
    dir.create(inputfolder)
  }
  
  ## Folder for binary result files with mcmc runs
  fitfolder<-"resfiles/"
  resfolder<-paste0(GlobalParameters$resultdir,"/",fitfolder)
  if(!dir.exists(resfolder)){
    cat(paste("Create new directory",resfolder,"\n"))
    dir.create(resfolder)
  }
  filename.mcmc1 <- paste0(resfolder,GlobalParameters$fitfile,"_mcmc1")
  filename.mcmc2 <- paste0(resfolder,GlobalParameters$fitfile,"_mcmc2")
  filename.hsz.mcmc2 <- paste0(resfolder,GlobalParameters$fitfile,"_hsz_mcmc2")
  filename.hsz.it <- paste0(resfolder,GlobalParameters$fitfile,"_hsz_it")
  filename.hsz.hauleff <- paste0(resfolder,GlobalParameters$fitfile,"_hsz_hauleff")
  
  common<-list(nSamples=GlobalParameters$nSamples,
               thin=GlobalParameters$thin,
               burnin=GlobalParameters$burnin,
               inputfolder=inputfolder,
               filename.mcmc1=filename.mcmc1,
               filename.mcmc2=filename.mcmc2,
               filename.hsz.mcmc2=filename.hsz.mcmc2,
               filename.hsz.it=filename.hsz.it,
               filename.hsz.hauleff=filename.hsz.hauleff,
               minage=GlobalParameters$minage,
               maxage=GlobalParameters$maxage,
               cc=GlobalParameters$CC,
               CCerror=GlobalParameters$CCerror,
               A2A=stoxdata$AgeLength$AgeErrorMatrix,
               age.error=GlobalParameters$age.error,
               resolution=GlobalParameters$lengthresCM,
               seed=GlobalParameters$seed,
               lgamodel=GlobalParameters$lgamodel,
               delta.age=GlobalParameters$delta.age,
               sim.ar=T,usedebug=F,print.boat=T,
               save.age=F,save.lga=F,save.wgl=F,inc.haulsize=F,
               print.format=GlobalParameters$print.format,
               old.version=GlobalParameters$old.version)

  run.fit(stoxdata,common,win)
  stoxdata<<-stoxdata
  fit <- read.fit.bin(filename.mcmc1,filename.mcmc2,filename.hsz.mcmc2,stoxdata)
  
  return(fit)
}

#' Predict catch-at-age
#'
#' This function runs the C prediction program and writes results to binary files.
#' @param AgeLength A list with elements connected to the age-length model, see \link{eca.estimate}.
#' @param WeightLength A list with elements connected to the weight-length model, see \link{eca.estimate}.
#' @param Landings A list with elements connected to landings, see Details
#' @param GlobalParameters A list with input parameters that are given by the user, see Details
#' @details
#'  Landings
#'  \itemize{
#'  \item LiveWeightKG  - Catch data given in kilograms for each cell.
#'  \item AgeLengthCov  - Covariate matrix for age-length model specifying which level of the different covariates the 
#'  corresponding cell in \emph{LiveWeightKG} has, e.g. constant, season, gearfactor, spatial. 
#'  It also contains a column \emph{midseason} which is the fraction of the year where the catch is taken.
#'  \item WeightLengthCov  - Covariate matrix for weight-length model specifying which level of the different covariates 
#'  the corresponding cell in \emph{LiveWeightKG} has, e.g. constant, season, gearfactor, spatial.
#'  It also contains a column \emph{midseason} which is the fraction of the year where the catch is taken.
#'  }
#'  GlobalParameters
#'  \itemize{
#'  \item minage - minimum age category
#'  \item maxage - maximum age category
#'  \item CC - TRUE if coastal cod is included in the model, FALSE otherwise
#'  \item resultdir - directory path where the results are stored. 
#'  Two directories \strong{cfiles} and \strong{resfiles} are in this folder (created if they don't exist). 
#'  In \strong{cfiles} the temporary binary input files to the c-program are stored, and in \strong{resfiles} the
#'  binary output files are stored.
#'  \item fitfile - name of the output files from estimating the model. 
#'  The files are given the extensions \strong{.mcmc1} and \strong{.mcmc2} for the age and length-given-age model 
#'  and for the weight-given-length model, respectively.
#'  \item predictfile - name of the output file from prediction. The file is given the extension \strong{.pred}.
#'  \item caa.burnin - number of MCMC samples discarded before running the prediction
#'  \item seed - random seed value
#'  }
#' @return A list with elements
#' \item{TotalCount}{array of the total catch (in numbers) with dimensions #length intervals x #age categories x #samples}
#' \item{MeanLength}{matrix of the mean length in the different age categories with dimensions #age categories x #samples}
#' \item{MeanWeight}{matrix of the mean weight in the different age categories with dimensions #age categories x #samples}
#' \item{AgeCategories}{vector with age categories}
#' \item{LengthIntervalsLog}{vector with length intervals, on log scale}
#' If a multiple stocks (coastal cod/Atlantic cod) analysis is run, the function returns the same list, but where the
#' dimension of age categories in all the elements is 2 times the number of age categories in the model (ncat). The first ncat 
#' valuels correspond to Atlantic cod while the last ncat values correspond to coastal cod.
#' @seealso \link{eca.estimate} for estimating the catch-at-age model.
#' @keywords c-program mcmc
#' @export
## Created 27.10.2017
## Modified
## 12.10.2018  HWR  changed to ascii files for common parameters
##
eca.predict <- function(AgeLength,WeightLength,Landings,GlobalParameters)
{
  ## Create variable for specifying specific windows commands
  sysinfo <- Sys.info()
  if(sysinfo["sysname"]=="Windows"){
    win <- 1
  } else {
    win <- 0
  }

  stoxdata<-list(AgeLength=AgeLength,WeightLength=WeightLength,
                 Landings=Landings,GlobalParameters=GlobalParameters)
  
  if(GlobalParameters$caa.burnin>=GlobalParameters$nSamples){
    stop(paste("Error: Number of burnin samples in prediction(",GlobalParameters$caa.burnin,
               ") must be smaller than the number of samples when fitting the model(",
               GlobalParameters$nSamples,")\n",sep=""))
  }

    
  ## Folder for input files
  cfilesfolder<-"cfiles/"
  inputfolder<-paste0(GlobalParameters$resultdir,"/",cfilesfolder)
  
  ## Folder for binary result files with mcmc runs
  fitfolder<-"resfiles/"
  resfolder<-paste0(GlobalParameters$resultdir,"/",fitfolder)
  if(!dir.exists(resfolder)){
    stop(paste("Error: Directory ",resfolder," does not exist\n"))
  }
  filename.mcmc1 <- paste0(resfolder,GlobalParameters$fitfile,"_mcmc1")
  filename.mcmc2 <- paste0(resfolder,GlobalParameters$fitfile,"_mcmc2")
  filename.hsz.mcmc2 <- paste0(resfolder,GlobalParameters$fitfile,"_hsz_mcmc2")
  filename.hsz.it <- paste0(resfolder,GlobalParameters$fitfile,"_hsz_it")
  filename.predict <- paste0(resfolder,GlobalParameters$predictfile,"_predict")

  common <- list(npredictMC=1000,
                 caa.burnin=GlobalParameters$caa.burnin,
                 minage=GlobalParameters$minage,
                 maxage=GlobalParameters$maxage,
                 maxlength=GlobalParameters$maxlength,
                 lengthresCM=GlobalParameters$lengthresCM,
                 cc=GlobalParameters$CC,
                 inputfolder=inputfolder,
                 filename.mcmc1=filename.mcmc1,
                 filename.mcmc2=filename.mcmc2,
                 filename.hsz.mcmc2=filename.hsz.mcmc2,
                 filename.hsz.it=filename.hsz.it,
                 filename.predict=filename.predict)
  
  run.pred(stoxdata,common,win)
                     
  pred <- read.predict.bin(filename.predict)
  
  return(pred)
}
