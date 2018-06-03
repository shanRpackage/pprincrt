#' Bayesian general linear regression analysis of current trial data
#'
#'\code{hmodelBUGS} returns posterior samples of treatment effect
#'
#' @param hData historical trial data. It should be in a data frame.
#' @param hForm model for historical trial data. It should be a formula object.
#' @param file.dir directory to save model BUGS file and input and output files of OpenBUGS.
#' @param file.name the name of the model BUGS file.
#' @param file.rm whether to remove all files after Bayesian computations are done.
#' @param OpenBUGS.dir directory to install OpenBUGS package.
#' @param OpenBUGS.seed the random number generator state of OpenBUGS.
#' @param family outcome distribution.
#' @param nchain number of MCMC chains.
#' @param niter length of MCMC chain.
#' @param nburnin burnin period of MCMC chain.
#' @param nthin thinning rate of MCMC chain.
#' @return The posterior samples of treatment effect
#' @examples
#' \dontrun{
#' #' my historical data;
#' my.hData <- data.frame(y=rbinom(3,50,0.5),n=rep(50,3),x=c(1,2,3))
#'
#' my.hsample <- hmodelBUGS(hData=my.hData,hForm=y|n ~ x,
#' file.dir=tempdir(),file.name='hmodelfile',file.rm=T, OpenBUGS.dir=NULL,
#' OpenBUGS.seed=14, family="binomial",nchain=2,niter=5000,nburnin=2500,nthin=5)
#' }
#' @export

hmodelBUGS <- function(hData, hForm, file.dir, file.name, file.rm=F, OpenBUGS.dir,OpenBUGS.seed=1, family,nchain,niter,nburnin,nthin){

  ### extract the variable names;
  hStrForm <- as.character(hForm)
  if((length(base::strsplit(hStrForm[3], split="(\\+)|(\\|)")[[1]])==1) & (length(base::strsplit(hStrForm[2], split="\\|")[[1]])==2)){
    hOutcomeVar <- stringr::str_trim(base::strsplit(hStrForm[2],split="\\|")[[1]][1])
    hTrtVar <- hStrForm[3]
    hSumVar <- stringr::str_trim(base::strsplit(hStrForm[2],split="\\|")[[1]][2])
  } else{
    stop("There should be one fixed effect, no random effect and the number of subjects by arm in historical trial!")
  }

  ntrt <- length(unique(hData[,hTrtVar]))
  if (ntrt == 1){
    stop ("At least two arms are needed in historical trial !")
  } else {
    pprincrt::hmodelfile(ntrt,family,paste(file.dir,"\\",file.name,".txt",sep=""))

    ### data to assign;

    data.str <- list(ntrt)
    htrtdummies <- stats::model.matrix( ~ factor(hData[,hTrtVar]))[,-c(1)]

    if(ntrt==2){
      names(htrtdummies) <- NULL
      data.str[[2]] <- htrtdummies
    } else {
      colnames(htrtdummies) <- NULL
      for(i in 1:(ntrt-1)){
        data.str[[i+1]] <- htrtdummies[,i]
      }
    }

    data.str[[ntrt+1]] <- hData[,hOutcomeVar]
    data.str[[ntrt+2]] <- hData[,hSumVar]

    data.name <- c("I")
    for(i in 1:(ntrt-1)){
      data.name <- c(data.name,paste("trt",i,sep=""))
    }
    data.name <- c(data.name,"y","n")
    names(data.str) <- data.name

    ### parameters to monitor;

    par.str <- c()
    for(i in 1:(ntrt-1)){
      par.str <- c(par.str,paste("beta",i,sep=""))
    }

    bugsdata <- data.str
    bugsinits <- list()
    for(i in 1:nchain){bugsinits[[i]] <- list(beta0=stats::rnorm(1,0,sqrt(10)))}
    bugsparameters <- par.str
    Bayes.sim <- R2OpenBUGS::bugs(bugsdata,bugsinits,bugsparameters,model.file=paste(file.dir,"\\",file.name,".txt",sep=""),
                                  OpenBUGS.pgm=OpenBUGS.dir,n.chain=nchain,n.iter=niter,n.burnin=nburnin,n.thin=nthin,DIC=FALSE,
                                  clearWD=file.rm, working.directory=file.dir,bugs.seed=OpenBUGS.seed)
    BUGS.rlts <- Bayes.sim$sims.matrix

    if(file.rm){
      base::file.remove(paste(file.dir,"\\",file.name,".txt",sep=""))
    }
    return(BUGS.rlts)
  }
}
