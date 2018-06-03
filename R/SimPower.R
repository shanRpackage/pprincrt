#' power in designing a new cluster randomized trial with the power prior method
#'
#'\code{SimPower} returns power in designing a new cluster randomized trial with the power prior method
#'
#' @param family outcome distribution.
#' @param to.do.option a value in 'sim_power', 'sim_SSD' and 'real_SSD'.
#' @param cSet current trial setting. It should be a list.
#' @param cNarmVar the variable name of arm in current trial.
#' @param cNcluster.armVar the variable name of total clusters in each arm in current trial. Its default value is 'ncluster.arm'.
#' @param cNpat.clusterVar the variable name of total patients in each cluster in current trial. Its default value is 'npat.cluster'.
#' @param cBtwVar the variable name of between cluster variation in current trial.
#' @param cWthVar the variable name of within cluster variation in current trial.
#' @param cTrtVecVar the variable name of treatment in current trial. It should be in the order of first treatment, second treatment until the last one.
#' @param hSet historical trial setting. It should be in a list. It is NULL by default.
#' @param hNarmVar the variable name of arm in historical trial.
#' @param hNpat.armVar the variable name of total patients in each arm in historical trial. Its default value is 'npat.arm'.
#' @param hWthVar the variable name of within cluster variation in historical trial.
#' @param hTrtVecVar the variable name of treatment in historical trial. It should be in the order of first treatment, second treatment until the last one.
#' @param hData Real historical trial data. It is NULL by default.
#' @param hTrtVar the variable name of treatment in historical trial data.
#' @param hOutcomeVar the variable name of outcome in historical trial data.
#' @param hSumVar the variable name of total observation in historical trial data.
#' @param cRdmSeed.init the initial value of random seed to generate current trial data.
#' @param hRdmSeed.init the initial value of random seed to generate historical trial data.
#' @param weight discounting parameter. It is 'asym' by default.
#' @param cover.level nomial coverage level of HPD interval.
#' @param nchain number of mcmc chains to run.
#' @param niter length of mcmc chains.
#' @param nburnin length of burnin of mcmc chains.
#' @param nthin thinning rate of mcmc chain.
#' @param file.dir.init the start of the directory to save the model BUGS file and input and output files of OpenBUGS.
#' @param cfile.name the file name of current BUGS model.
#' @param hfile.name the file name of historical BUGS model.
#' @param pprfile.name the file name of power prior BUGS model.
#' @param graphics.name the file name of mcmc convergence graphics.
#' @param file.rm whether to remove all files after Bayesian computations are done.
#' @param OpenBUGS.dir directory to install OpenBUGS package.
#' @param OpenBUGS.seed the random number generator state of OpenBUGS.
#' @param Rep number of replications to determine power.
#' @param num.core number of cores for parallel computation.
#' @return power in designing a new cluster randomized trial with the power prior method
#' @examples
#' \dontrun{
#' #### SSD with simulated historical data;
#'my.curset <- list(narm=2,ncluster.arm=rep(10,2),npat.cluster=rep(20,20),
#' 'trt1'=-0.1,'trt2 v.s. trt1'=log(4),sigma.b=1)
#'my.histset <- list(narm=2,npat.arm=rep(200,2),'trt1'=0.1,'trt2 v.s. trt1'=log(4))
#'mypower <- SimPower(family='binomial', to.do.option='sim_SSD', cSet=my.curset,
#'cTrtVecVar=c('trt1','trt2 v.s. trt1'), hSet=my.histset,
#'hTrtVecVar=c('trt1','trt2 v.s. trt1'),
#'weight="asym",nchain=2,niter=500,nburnin=100,nthin=5,Rep=1,num.core=3,file.rm=T,
#'OpenBUGS.dir="C:\\Program Files (x86)\\OpenBUGS\\OpenBUGS323\\OpenBUGS.exe",
#'OpenBUGS.seed=12)
#'  #### SSD with real historical data;
#'my.curset <- list(narm=2,ncluster.arm=rep(10,2),npat.cluster=rep(20,20),
#''trt1'=-0.1,'trt2 v.s. trt1'=log(4),sigma.b=1)
#'my.histdata <- data.frame(x=c(1,2),y=c(10,15),n=c(20,20))
#'mypower <- SimPower(family='binomial', to.do.option='real_SSD',cSet=my.curset,
#'cTrtVecVar=c('trt1','trt2 v.s. trt1'),hData=my.histdata,weight='sym',
#'nchain=2,niter=500,nburnin=100,nthin=5,file.rm=T, Rep=4,num.core=2,
#'OpenBUGS.dir="C:\\Program Files (x86)\\OpenBUGS\\OpenBUGS323\\OpenBUGS.exe")
#' }
#' @export
#' @importFrom foreach %dopar%
#' @name SimPower

SimPower <- function(family, to.do.option, cSet, cNarmVar='narm', cNcluster.armVar='ncluster.arm', cNpat.clusterVar='npat.cluster',cBtwVar='sigma.b',cWthVar='sigma.e',
                     cTrtVecVar, hSet=NULL, hNarmVar='narm',hNpat.armVar='npat.arm', hWthVar='sigma.e', hTrtVecVar=NULL,hData=NULL, hTrtVar='x', hOutcomeVar='y', hSumVar='n',
                     cRdmSeed.init=0, hRdmSeed.init=0, weight='asym', cover.level=0.95,nchain=1,niter=2000,nburnin=1000,nthin=1, file.dir.init=tempdir(),
                     cfile.name='cmodelfile', hfile.name='hmodelfile', pprfile.name='pprmodelfile', graphics.name='graphics',file.rm=T,OpenBUGS.dir=NULL,
                     OpenBUGS.seed=1, Rep=1, num.core=parallel::detectCores(logical=F)-1){

  #### create folders to store input & output of OpenBUGS;
  file.dir <- paste(file.dir.init, "\\",1:Rep, sep="")
  for(i in 1:length(file.dir)){ base::dir.create(file.dir[i])}

  #### start parallel computing;
  mycl <- snow::makeCluster(num.core,type='SOCK')
  doSNOW::registerDoSNOW(mycl)
  sim.power <- c()

  if(identical(to.do.option, 'sim_power')){
    Rej.rlts <- foreach::foreach(i=1:Rep,.combine='rbind',.packages=c("pprincrt")) %dopar% {
      pprincrt::SimRej(family, to.do.option,cSet,cNarmVar, cNcluster.armVar, cNpat.clusterVar,cTrtVecVar,
             cBtwVar,cWthVar, hSet, hNarmVar, hNpat.armVar, hTrtVecVar, hWthVar, hData, hTrtVar, hOutcomeVar, hSumVar,
             cRdmSeed.init+i, hRdmSeed.init+i, weight, cover.level, nchain,niter,nburnin,nthin,
             file.dir[i], cfile.name, hfile.name,pprfile.name, graphics.name,file.rm,OpenBUGS.dir, OpenBUGS.seed)
    }
  } else if (identical(to.do.option, 'sim_SSD')){
    Rej.rlts <- foreach::foreach(i=1:Rep,.combine='rbind',.packages=c("pprincrt")) %dopar% {
      pprincrt::SimRej(family, to.do.option,cSet,cNarmVar, cNcluster.armVar, cNpat.clusterVar,cTrtVecVar,
             cBtwVar,cWthVar, hSet, hNarmVar, hNpat.armVar, hTrtVecVar, hWthVar, hData, hTrtVar, hOutcomeVar, hSumVar,
             cRdmSeed.init+i, hRdmSeed.init, weight, cover.level, nchain,niter,nburnin,nthin,
             file.dir[i],cfile.name, hfile.name,pprfile.name, graphics.name, file.rm,OpenBUGS.dir, OpenBUGS.seed)
    }
  } else if (identical(to.do.option, 'real_SSD')){
    Rej.rlts <- foreach::foreach(i=1:Rep,.combine='rbind',.packages=c("pprincrt")) %dopar% {
      pprincrt::SimRej(family, to.do.option,cSet,cNarmVar, cNcluster.armVar, cNpat.clusterVar,cTrtVecVar,
             cBtwVar,cWthVar, hSet, hNarmVar, hNpat.armVar, hTrtVecVar, hWthVar, hData, hTrtVar, hOutcomeVar, hSumVar,
             cRdmSeed.init+i, hRdmSeed.init, weight, cover.level, nchain,niter,nburnin,nthin,
             file.dir[i],cfile.name, hfile.name,pprfile.name, graphics.name, file.rm,OpenBUGS.dir, OpenBUGS.seed)
    }
  } else{
    stop ("to.do.option should be a character 'sim_power', 'sim_SSD' or 'real_SSD' !")
  }
  snow::stopCluster(mycl)

  #### delete all folders;
  for(i in 1:length(file.dir)){base::unlink(file.dir[i],recursive=file.rm)}

  if(Rep==1){
    sim.power <- Rej.rlts
  } else {
    sim.power <- colMeans(Rej.rlts)
  }

  return(sim.power)
}
