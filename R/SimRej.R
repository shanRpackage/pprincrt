#' rejection decision of null treatment effect
#'
#'\code{SimRej} returns rejection decision of null treatment effect
#'
#' @param family outcome distribution
#' @param to.do.option a value in 'sim_power', 'sim_SSD' and 'real_SSD'
#' @param cSet current trial setting. It should be in a list.
#' @param cNarmVar the variable name of arm in current trial.
#' @param cNcluster.armVar the variable name of total clusters in each arm in current trial. Its default value is 'ncluster.arm'.
#' @param cNpat.clusterVar the variable name of total patients in each cluster in current trial. Its default value is 'npat.cluster'.
#' @param cTrtVecVar the variable name of treatment in current trial. It should be in the order of first treatment, second treatment until the last one.
#' @param cBtwVar the variable name of between cluster variation in current trial.
#' @param cWthVar the variable name of within cluster variation in current trial.
#' @param hSet historical trial setting. It should be in a list. It is NULL by default.
#' @param hNarmVar the variable name of arm in historical trial.
#' @param hNpat.armVar the variable name of total patients in each arm in historical trial. Its default value is 'npat.arm'.
#' @param hTrtVecVar the variable name of treatment in historical trial. It should be in the order of first treatment, second treatment until the last one.
#' @param hWthVar the variable name of within cluster variation in historical trial.
#' @param hData Real historical trial data. It is NULL by default.
#' @param hTrtVar the variable name of treatment in historical trial data.
#' @param hOutcomeVar the variable name of outcome in historical trial data.
#' @param hSumVar the variable name of total observation in historical trial data.
#' @param cRdmSeed the random seed to generate current trial data.
#' @param hRdmSeed the random seed to generate historical trial data.
#' @param weight discounting parameter. It is 'asym' by default.
#' @param cover.level nomial coverage level of HPD interval.
#' @param nchain number of mcmc chains to run.
#' @param niter length of mcmc chains.
#' @param nburnin length of burnin of mcmc chains.
#' @param nthin thinning rate of mcmc chain.
#' @param file.dir directory to save the model BUGS file and input and output files of OpenBUGS.
#' @param cfile.name the file name of current BUGS model.
#' @param hfile.name the file name of historical BUGS model.
#' @param pprfile.name the file name of power prior BUGS model.
#' @param graphics.name the file name of mcmc convergence graphics.
#' @param file.rm whether to remove all files after Bayesian computations are done.
#' @param OpenBUGS.dir directory to install OpenBUGS package.
#' @param OpenBUGS.seed the random number generator state of OpenBUGS.
#' @return rejection decision of null treatment effect
#' @examples
#' \dontrun{
#' my.curset <- list(narm=3,ncluster.arm=rep(5,3),npat.cluster=rep(20,15),
#' 'trt1'=0.1,'trt2 v.s. trt1'=1,'trt3 v.s. trt1'=1,sigma.b=1)
#' my.histset <- list(narm=3,npat.arm=rep(20,3),
#' 'trt1'=-0.1,'trt2 v.s. trt1'=1,'trt3 v.s. trt1'=1)
#' myrej <- SimRej(family='binomial', to.do.option="sim_power",cSet=my.curset,
#' cTrtVecVar=c('trt1', paste(paste('trt', 2:3, sep=""),'v.s.', 'trt1',sep=" ")),
#' hSet=my.histset,
#' hTrtVecVar=c('trt1', paste(paste('trt', 2:3, sep=""),'v.s.', 'trt1',sep=" ")),
#' file.rm=T, nchain=2,niter=200,nburnin=50,nthin=5, OpenBUGS.seed=4)
#' }
#' @export

SimRej <- function(family, to.do.option, cSet,cNarmVar='narm', cNcluster.armVar='ncluster.arm', cNpat.clusterVar='npat.cluster',cTrtVecVar,
                   cBtwVar='sigma.b',cWthVar='sigma.e', hSet=NULL, hNarmVar='narm', hNpat.armVar='npat.arm', hTrtVecVar=NULL, hWthVar='sigma.e',
                   hData=NULL, hTrtVar='x', hOutcomeVar='y', hSumVar='n', cRdmSeed=1, hRdmSeed=1, weight='asym', cover.level=0.95,
                   nchain=1,niter=2000,nburnin=1000,nthin=1, file.dir=tempdir(), cfile.name='cmodelfile', hfile.name='hmodelfile',
                   pprfile.name='pprmodelfile',graphics.name='graphics',file.rm=F, OpenBUGS.dir=NULL, OpenBUGS.seed=1){

  ### simulated or extracted historical data;
  if(identical(to.do.option,"sim_power") || identical(to.do.option, "sim_SSD")){

    if(cSet[[cNarmVar]] != hSet[[hNarmVar]]){
      stop ('number of arms should be the same across the two trials !')
    } else if(cSet[[cNarmVar]]==1 && hSet[[hNarmVar]] ==1){
      stop ('at least two arms are needed in both trials !')
    } else{
      narm <- hSet[[hNarmVar]]
      npat.arm <- hSet[[hNpat.armVar]]

      beta.prime <- c()
      for(i in 1:length(hTrtVecVar)){
        if(length(hSet[[hTrtVecVar[i]]])==1){
          beta.prime <- c(beta.prime, hSet[[hTrtVecVar[i]]])
        } else if(length(hSet[[hTrtVecVar[i]]])==2){
          base::set.seed(hRdmSeed)
          beta.prime.tmpt <- stats::runif(1,hSet[[hTrtVecVar[i]]][1],hSet[[hTrtVecVar[i]]][2])
          beta.prime <- c(beta.prime, beta.prime.tmpt)
        } else{
          stop("treatment effect in historical trial should be fixed value or within a fixed interval !")
        }
      }

      Trt.hist <- rep(1:narm,each=1)
      X.hist <- matrix(0,ncol=(narm-1),nrow=narm)
      for(i in 1:(narm-1)){
        X.hist[(i+1),i] <- 1
      }
      pred.hist <- cbind(rep(1,narm),X.hist) %*% beta.prime

      if(identical(family,"gaussian")){

        sigma.e.hist <- c()
        if(length(hSet[[hWthVar]])==1 && hSet[[hWthVar]] > 0){
          sigma.e.hist <- c(sigma.e.hist,hSet[[hWthVar]])
        } else if(length(hSet[[hWthVar]])==2 && hSet[[hWthVar]][1] > 0 && hSet[[hWthVar]][2] >= hSet[[hWthVar]][1] ){
          base::set.seed(hRdmSeed)
          sigma.e.hist.tmpt <- stats::runif(1,hSet[[hWthVar]][1],hSet[[hWthVar]][2])
          sigma.e.hist <- c(sigma.e.hist,sigma.e.hist.tmpt)
        } else{
          stop("within-cluster variation in historical trial should be (positive) fixed value or within a fixed interval above zero ! ")
        }

        y.hist <- c()
        for(i in 1:narm){
          set.seed(hRdmSeed+i)
          y.hist.tmpt <- stats::rnorm(1,npat.arm[i]*pred.hist[i],sigma.e.hist*sqrt(npat.arm[i]))
          y.hist <- c(y.hist,y.hist.tmpt)
        }
      } else if (identical(family, "poisson")){
        lambda.hist <- exp(pred.hist)
        y.hist <- c()
        for(i in 1:narm){
          set.seed(hRdmSeed+i)
          y.hist.tmpt <- stats::rpois(1,npat.arm[i]*lambda.hist[i])
          y.hist <- c(y.hist,y.hist.tmpt)
        }
      } else if(identical(family,"binomial")){
        p.hist <- 1/(1+exp(-pred.hist))
        y.hist <- c()
        for(i in 1:narm){
          set.seed(hRdmSeed+i)
          y.hist.tmpt <- stats::rbinom(1,npat.arm[i],p.hist[i])
          y.hist <- c(y.hist,y.hist.tmpt)
        }
      } else {
        stop ("family should be a character 'gaussian', 'poisson' or 'binomial' ! ")
      }

      data.hist <- data.frame(x=Trt.hist,y=y.hist,n=npat.arm)
      hTrtVar <- 'x'
      hOutcomeVar <- 'y'
      hSumVar <- 'n'
    }
  } else if (identical(to.do.option, "real_SSD")){
    if(!is.null(hData)){
      if(cSet[[cNarmVar]] != length(unique(hData[,hTrtVar]))){
        stop('number of arms should be the same across the two trials !')
      } else if(cSet[[cNarmVar]]==1 && length(unique(hData[,hTrtVar]))==1){
        stop ('at least two arms are needed in both trials !')
      } else{
        data.hist <- hData
      }
    } else {
      stop("historical data should not be empty !")
    }
  } else {
    stop ("to.do.option should be a character 'sim_power', 'sim_SSD' or 'real_SSD' !")
  }

  ### simulate current data;
  narm <- cSet[[cNarmVar]]
  ncluster.arm <- cSet[[cNcluster.armVar]]
  npat.cluster <- cSet[[cNpat.clusterVar]]

  beta <- c()
  for(i in 1:length(cTrtVecVar)){
    if(length(cSet[[cTrtVecVar[i]]])==1){
      beta <- c(beta, cSet[[cTrtVecVar[i]]])
    } else if (length(cSet[[cTrtVecVar[i]]])==2){
      base::set.seed(cRdmSeed)
      beta.tmpt <- stats::runif(1,cSet[[cTrtVecVar[i]]][1], cSet[[cTrtVecVar[i]]][2])
      beta <- c(beta,beta.tmpt)
    } else {
      stop("treatment effect in current trial should be fixed value or  within a fixed interval !")
    }
  }

  sigma.b <- c()
  if(length(cSet[[cBtwVar]])==1 && cSet[[cBtwVar]] > 0){
    sigma.b <- c(sigma.b, cSet[[cBtwVar]])
  } else if(length(cSet[[cBtwVar]])==2 && cSet[[cBtwVar]][1] > 0 && cSet[[cBtwVar]][2] >= cSet[[cBtwVar]][1]){
    base::set.seed(cRdmSeed)
    sigma.b.tmpt <- stats::runif(1,cSet[[cBtwVar]][1],cSet[[cBtwVar]][2])
    sigma.b <- c(sigma.b,sigma.b.tmpt)
  } else{
    stop("between-cluster variation in current trial should be (positive) fixed value or within a fixed interval above zero !")
  }

  clusterID <- rep(1:sum(ncluster.arm),npat.cluster)
  Trt.cur <- rep(1,sum(npat.cluster))
  X.cur <- matrix(0,ncol=(narm-1),nrow=sum(npat.cluster))
  for(i in 1:(narm-1)){
    lwr.index <- sum(npat.cluster[1:cumsum(ncluster.arm)[i]])+1
    upr.index <- sum(npat.cluster[1:cumsum(ncluster.arm)[i+1]])

    Trt.cur[lwr.index:upr.index] <- rep(i+1,upr.index-lwr.index+1)
    X.cur[lwr.index:upr.index,i] <- rep(1,(upr.index-lwr.index+1))
  }
  base::set.seed(cRdmSeed)
  b.tmpt <- stats::rnorm(sum(ncluster.arm),0,sigma.b)
  b <- rep(b.tmpt,npat.cluster)
  pred.cur <- cbind(rep(1,sum(npat.cluster)),X.cur)%*%beta + b

  if (identical(family, "gaussian")){

    sigma.e.cur <- c()
    if(length(cSet[[cWthVar]])==1 && cSet[[cWthVar]] > 0){
      sigma.e.cur <- c(sigma.e.cur, cSet[[cWthVar]])
    } else if(length(cSet[[cWthVar]])==2 && cSet[[cWthVar]][1] > 0 && cSet[[cWthVar]][2] >= cSet[[cWthVar]][1]){
      base::set.seed(cRdmSeed)
      sigma.e.cur.tmpt <- stats::runif(1, cSet[[cWthVar]][1], cSet[[cWthVar]][2])
      sigma.e.cur <- c(sigma.e.cur,sigma.e.cur.tmpt)
    } else{
      stop ("within-cluster variation in current trial should be (positive) fixed value or within a fixed interval above zero !")
    }

    y.cur <- c()
    for(i in 1:sum(npat.cluster)){
      base::set.seed(cRdmSeed+i)
      y.cur.tmpt <- stats::rnorm(1,pred.cur[i],sigma.e.cur)
      y.cur <- c(y.cur,y.cur.tmpt)
    }
  } else if (identical(family,"poisson")){
    lambda.cur <- exp(pred.cur)
    y.cur <- c()
    for(i in 1:sum(npat.cluster)){
      base::set.seed(cRdmSeed+i)
      y.cur.tmpt <- stats::rpois(1,lambda.cur[i])
      y.cur <- c(y.cur,y.cur.tmpt)
    }
  } else if(identical(family, "binomial")){
    p.cur <- 1/(1+exp(-pred.cur))
    y.cur <- c()
    for(i in 1:sum(npat.cluster)){
      base::set.seed(cRdmSeed+i)
      y.cur.tmpt <- stats::rbinom(1,1,p.cur[i])
      y.cur <- c(y.cur,y.cur.tmpt)
    }
  } else {
    stop("family should be a character 'gaussian', 'poisson' or 'binomial' !")
  }

  data.cur <- data.frame(cl=clusterID,x=Trt.cur,y=y.cur)
  cClusterVar <- 'cl'
  cTrtVar <- 'x'
  cOutcomeVar <- 'y'

  #### make rejections;

  pprobject <- pprincrt::pprmodelBUGS(data.cur,stats::as.formula(paste(cOutcomeVar,paste(cTrtVar,paste(paste("(1",cClusterVar,sep="|"),")",sep=""),sep="+"),sep="~")),
                            data.hist,stats::as.formula(paste(paste(hOutcomeVar,hSumVar,sep="|"),hTrtVar,sep="~")),
                            weight,file.dir, cfile.name,hfile.name, pprfile.name, graphics.name, file.rm, OpenBUGS.dir,
                            OpenBUGS.seed, family,nchain,niter,nburnin,nthin,cover.level)
  cur.rej <- ifelse((pprobject$mcmc.summary[,3] >= 0 | pprobject$mcmc.summary[,4] <= 0),1,0)
  if(length(cur.rej)==1) {names(cur.rej) <- "trt2 v.s. trt1"}

  return(cur.rej)
}
