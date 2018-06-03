#' Bayesian mixed model analysis of current trial data
#'
#'\code{cmodelBUGS} returns posterior samples of treatment effect
#'
#' @param cData current trial data. It should be in a data frame.
#' @param cForm model for current trial data. It should be in a formula object.
#' @param file.dir directory to save model BUGS file and input and output files of OpenBUGS.
#' @param file.name the name of the model BUGS file.
#' @param file.rm whether to remove all files after Bayesian computations are done.
#' @param OpenBUGS.dir directory to install OpenBUGS package.
#' @param OpenBUGS.seed random number generator state of OpenBUGS.
#' @param family outcome distribution.
#' @param nchain number of MCMC chains.
#' @param niter length of MCMC chain.
#' @param nburnin burnin period of MCMC chain.
#' @param nthin thinning rate of MCMC chain.
#' @return The posterior samples of treatment effect
#' @examples
#' \dontrun{
#'   ### my current data;
#' my.cData <- data.frame(cl=rep(1:10,each=50),y=rbinom(500,1,0.5),x=c(rep(1,250),rep(2,250)))
#'
#' my.csample <- cmodelBUGS(cData=my.cData,cForm=y ~ x + (1|cl),
#' file.dir=tempdir(),file.name='cmodelfile',file.rm=T,
#' OpenBUGS.dir="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe",
#' family="binomial",nchain=2,niter=500,nburnin=50,nthin=5)
#' }
#' @export

cmodelBUGS <- function(cData,cForm, file.dir, file.name, file.rm=F, OpenBUGS.dir,OpenBUGS.seed=1, family,nchain,niter,nburnin,nthin){

  ### extract the variable names;
  cStrForm <- as.character(cForm)
  if((length(base::strsplit(cStrForm[3], split="\\+")[[1]])==2) & (length(base::strsplit(cStrForm[3], split="\\|")[[1]])==2)){
    cOutcomeVar <- cStrForm[2]
    cTrtVar <- stringr::str_trim(base::strsplit(cStrForm[3],split="\\+")[[1]][1])
    cClusterVar <- stringr::str_replace(stringr::str_trim(base::strsplit(cStrForm[3],split="\\|")[[1]][2]),pattern="\\)",replacement="")
  } else{
    stop("There should be one fixed effect and one random effect in current trial!")
  }

  ntrt <- length(unique(cData[,cTrtVar]))
  if (ntrt ==1){
    stop("At least two arms are needed in current trial !")
  } else{
    nobs <- dim(cData)[1]
    ncluster <- length(unique(cData[,cClusterVar]))

    pprincrt::cmodelfile(ntrt,family,paste(file.dir,"\\",file.name,".txt",sep=""))

    ### data to assign;

    data.str <- list(nobs,ncluster,cData[,cClusterVar])
    ctrtdummies <- stats::model.matrix( ~ factor(cData[,cTrtVar]))[,-c(1)]

    if(ntrt==2){
      names(ctrtdummies) <- NULL
      data.str[[4]] <- ctrtdummies
    } else {
      colnames(ctrtdummies) <- NULL
      for(i in 1:(ntrt-1)){
        data.str[[i+3]] <- ctrtdummies[,i]
      }
    }

    data.str[[ntrt+3]] <- cData[,cOutcomeVar]

    data.name <- c("N","K","cluster")
    for(i in 1:(ntrt-1)){
      data.name <- c(data.name,paste("trt",i,sep=""))
    }
    data.name <- c(data.name,"y")
    names(data.str) <- data.name

    ### parameters to monitor;

    par.str <- c()
    for(i in 1:(ntrt-1)){
      par.str <- c(par.str,paste("beta",i,sep=""))
    }

    bugsdata <- data.str
    bugsinits <- list()
    for(i in 1:nchain){bugsinits[[i]] <- list(sigma.b=stats::runif(1,0,5))}
    bugsparameters <- par.str
    Bayes.sim <- R2OpenBUGS::bugs(bugsdata,bugsinits,bugsparameters,model.file=paste(file.dir,"\\",file.name,".txt",sep=""),
                                  OpenBUGS.pgm=OpenBUGS.dir,n.chain=nchain,n.iter=niter,n.burnin=nburnin,n.thin=nthin,DIC=FALSE,
                                  clearWD=file.rm, working.directory=file.dir, bugs.seed=OpenBUGS.seed)
    BUGS.rlts <- Bayes.sim$sims.matrix

    if(file.rm){
      base::file.remove(paste(file.dir,"\\",file.name,".txt",sep=""))
    }

    return(BUGS.rlts)
  }
}
