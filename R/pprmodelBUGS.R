#' power prior method analysis
#'
#'\code{pprmodelBUGS} returns a pprMod object, including the posterior samples of the
#' treatment effect, the corresponding summary results and convergence diagnostic results
#'
#' @param cData current trial data. It should be in data frame.
#' @param cForm model for current trial data. It should be a formula object.
#' @param hData historical trial data. It should be in data frame.
#' @param hForm model for historical trial data. It should be a formula object.
#' @param weight discounting parameter. It should be a value between 0 and 1, or a character 'asym' or 'sym'.
#' @param file.dir directory to save model BUGS file and input and output files of OpenBUGS.
#' @param cfile.name the file name of curent trial BUGS model.
#' @param hfile.name the file name of historical trial BUGS model.
#' @param pprfile.name the file name of power prior BUGS model.
#' @param graphics.name the file name of the MCMC diagnostic graphics file.
#' @param file.rm whether to remove all files after all computations are done.
#' @param OpenBUGS.dir directory to install OpenBUGS package.
#' @param OpenBUGS.seed the random number generator state of OpenBUGS.
#' @param family outcome distribution.
#' @param nchain number of MCMC chains.
#' @param niter length of MCMC chain.
#' @param nburnin burnin period of MCMC chain.
#' @param nthin thinning rate of MCMC chain.
#' @param cover.level nominal coverage level of the HPD interval
#' @return The posterior samples of treatment effect
#' @examples
#' \dontrun{
#' ### gaussian data;
#' my.cData <- data.frame(cl=rep(1:20,each=50),
#' x=c(rep(1,250),rep(2,250),rep(3,250),rep(4,250)), y=rnorm(1000,0,1))
#' my.hData <- data.frame(x=c(1,2,3,4),y=rnorm(4,0,1),n=rep(10,4))
#' my.pprobject <- pprmodelBUGS(cData=my.cData,cForm=y~x+(1|cl),
#' hData=my.hData,hForm=y|n ~ x,weight='sym',file.rm=F,OpenBUGS.seed=2,
#' family="gaussian",niter=5000,nburnin=2500,nthin=5,nchain=2,cover.level=0.95)
#'  ### binomial data;
#' my.cData <- data.frame(cl=rep(1:15,each=50),x=c(rep(1,250),
#' rep(2,250),rep(3,250)),y=rbinom(750,1,0.5))
#' my.hData <- data.frame(x=c(1,2,3),y=rbinom(3,50,0.5),n=rep(50,3))
#' my.pprobject <- pprmodelBUGS(cData=my.cData,cForm=y~x+(1|cl),
#' hData=my.hData,hForm=y|n ~ x,weight="sym",family="binomial",
#' niter=200,nburnin=50,nthin=1,nchain=1)
#' }
#' @export

pprmodelBUGS <- function(cData, cForm, hData, hForm, weight,
                         file.dir=tempdir(),cfile.name='cmodelfile',hfile.name='hmodelfile',pprfile.name='pprmodelfile',
                         graphics.name='graphics', file.rm=F, OpenBUGS.dir=NULL,OpenBUGS.seed=1, family,nchain=1,niter=2000,nburnin=base::floor(niter/2),nthin=1,cover.level=0.95){

  ### extract the variable names;
  cStrForm <- as.character(cForm)
  if((length(base::strsplit(cStrForm[3], split="\\+")[[1]])==2) & (length(base::strsplit(cStrForm[3], split="\\|")[[1]])==2)){
    cOutcomeVar <- cStrForm[2]
    cTrtVar <- stringr::str_trim(base::strsplit(cStrForm[3],split="\\+")[[1]][1])
    cClusterVar <- stringr::str_replace(stringr::str_trim(base::strsplit(cStrForm[3],split="\\|")[[1]][2]),pattern="\\)",replacement="")
  } else{
    stop("There should be one fixed effect and one random effect in current trial!")
  }

  hStrForm <- as.character(hForm)
  if((length(base::strsplit(hStrForm[3], split="(\\+)|(\\|)")[[1]])==1) & (length(base::strsplit(hStrForm[2], split="\\|")[[1]])==2)){
    hOutcomeVar <- stringr::str_trim(base::strsplit(hStrForm[2],split="\\|")[[1]][1])
    hTrtVar <- hStrForm[3]
    hSumVar <- stringr::str_trim(base::strsplit(hStrForm[2],split="\\|")[[1]][2])
  } else{
    stop("There should be one fixed effect, no random effect and the number of subjects by arm in historical trial!")
  }

  if(identical(family,"gaussian") | identical(family,"poisson") | identical(family,"binomial")){

    cntrt <-  length(unique(cData[,cTrtVar]))
    hntrt <- length(unique(hData[,hTrtVar]))
    if(cntrt != hntrt){
      stop("Number of arms should be the same across the two trials !")
    } else if ((cntrt ==1) && (hntrt ==1)){
      stop("At least two arms are needed in the two trials !")
    } else{
      if(identical(weight,"sym")){
        k <- base::ceiling(0.05*nchain*(niter-nburnin))
        ctrtsample <- pprincrt::cmodelBUGS(cData,cForm, file.dir,cfile.name,file.rm, OpenBUGS.dir, OpenBUGS.seed, family,nchain,niter,nburnin,nthin)
        htrtsample <- pprincrt::hmodelBUGS(hData,hForm, file.dir,hfile.name,file.rm, OpenBUGS.dir, OpenBUGS.seed,family,nchain,niter,nburnin,nthin)
        distance.par <- FNN::KL.divergence(ctrtsample,htrtsample,k)[k] + FNN::KL.divergence(htrtsample,ctrtsample,k)[k]
        DisPar <- exp(- ifelse(distance.par >= 0,distance.par,0))
      } else if (identical(weight,"asym")){
        k <- base::ceiling(0.05*nchain*(niter-nburnin))
        ctrtsample <- pprincrt::cmodelBUGS(cData,cForm,file.dir,cfile.name,file.rm, OpenBUGS.dir, OpenBUGS.seed, family,nchain,niter,nburnin,nthin)
        htrtsample <- pprincrt::hmodelBUGS(hData,hForm,file.dir,hfile.name,file.rm, OpenBUGS.dir, OpenBUGS.seed, family,nchain,niter,nburnin,nthin)
        distance.par <- FNN::KL.divergence(ctrtsample,htrtsample,k)[k]
        DisPar <- exp(- ifelse(distance.par >= 0,distance.par,0))
      } else if (weight>=0 && weight <=1){
        DisPar <- weight
      } else{
        stop("Weight should be a value from 0 to 1 or a character 'sym' or 'asym'!")
      }

      ### run power prior analysis under appropriate dist. & discounting parameter;
      nobs <- dim(cData)[1]
      ncluster <- length(unique(cData[,cClusterVar]))
      ntrt <-  length(unique(cData[,cTrtVar]))

      pprincrt::pprmodelfile(ntrt,family,paste(file.dir,"\\",pprfile.name,".txt",sep=""))

      ### data to assign;

      data.str <- list(DisPar,nobs,ncluster,ntrt,cData[,cClusterVar])
      ctrtdummies <- stats::model.matrix( ~ factor(cData[,cTrtVar]))[,-c(1)]
      htrtdummies <- stats::model.matrix( ~ factor(hData[,hTrtVar]))[,-c(1)]

      if(ntrt == 2){
        names(ctrtdummies) <- names(htrtdummies) <- NULL
        data.str[[6]] <- c(ctrtdummies, htrtdummies)
      } else{
        colnames(ctrtdummies) <- colnames(htrtdummies) <- NULL
        for(i in 1:(ntrt-1)){
          data.str[[i+5]] <- c(ctrtdummies[,i],htrtdummies[,i])
        }
      }
      data.str[[ntrt+5]] <- c(cData[,cOutcomeVar],hData[,hOutcomeVar])
      data.str[[ntrt+6]] <- hData[,hSumVar]

      data.name <- c("a","N","K","I","cluster")
      for(i in 1:(ntrt-1)){
        data.name <- c(data.name,paste("trt",i,sep=""))
      }
      data.name <- c(data.name,"y","n")
      names(data.str) <- data.name

      ### parameters to monitor;

      par.str <- c()
      for(i in 1:(ntrt-1)){
        par.str <- c(par.str,paste("beta",i,"[1]",sep=""))
      }

      if(ntrt > 2){
        for(i in 1:(ntrt-2)){
          for(j in (i+1):(ntrt-1)){
            par.str <- c(par.str,paste("beta",i,j,sep=""))
          }
        }
      }

      bugsdata <- data.str
      bugsinits <- list()
      for(i in 1:nchain){bugsinits[[i]] <- list(sigma.b=stats::runif(1,0,5))}
      bugsparameters <- par.str
      Bayes.sim <- R2OpenBUGS::bugs(bugsdata,bugsinits,bugsparameters,model.file=paste(file.dir,"\\",pprfile.name,".txt",sep=""),
                                    OpenBUGS.pgm=OpenBUGS.dir,n.chain=nchain,n.iter=niter,n.burnin=nburnin,n.thin=nthin,DIC=FALSE,
                                    clearWD=file.rm, working.directory=file.dir, bugs.seed=OpenBUGS.seed)
      BUGS.rlts <- Bayes.sim$sims.array

      if(file.rm){
        base::file.remove(paste(file.dir,"\\",pprfile.name,".txt",sep=""))
      }

      BUGS.rlts.name <- c()
      for(i in 1:(ntrt-1)){
        BUGS.rlts.name <- c(BUGS.rlts.name,paste(paste("trt",i+1,sep=""),"v.s.",paste("trt",1,sep=""),sep=" "))
      }

      if(ntrt > 2){
        for(i in 1:(ntrt-2)){
          for(j in (i+1):(ntrt-1)){
            BUGS.rlts.name <- c(BUGS.rlts.name,paste(paste("trt",i+1,sep=""),"v.s.",paste("trt",j+1,sep=""),sep=" "))
          }
        }
      }
      dimnames(BUGS.rlts)[[3]] <- BUGS.rlts.name

      ### mcmc diagnostic results;

      npar <- dim(BUGS.rlts)[3]
      color.option <- c(2,4,3,5:(nchain+5))

      if(nchain==1){
        if(npar==1){
          mcmc.sample.obj.tmpt <- as.matrix(BUGS.rlts[,1,],ncol=1)
          colnames(mcmc.sample.obj.tmpt) <- paste("trt2", "v.s.", "trt1",sep=" ")
          mcmc.sample.obj <- coda::as.mcmc(mcmc.sample.obj.tmpt)
        } else{
          mcmc.sample.obj <- coda::as.mcmc(BUGS.rlts[,1,])
        }

        rtn.gelman <- "number of chain should be more than 1!"

        grDevices::pdf(file=paste(file.dir,"\\",graphics.name,".pdf",sep=""))
        graphics::plot(mcmc.sample.obj,trace=T,density=F,smooth=F,col=color.option[1],auto.layout=F)
        grDevices::dev.off()
        rtn.plot <- paste(file.dir,"\\",graphics.name,".pdf",sep="")

      } else{
        mcmc.sample.list1 <- list()
        if(npar==1){
          for(i in 1:nchain){
            mcmc.sample.tmpt <- as.matrix(BUGS.rlts[,i,],ncol=1)
            colnames(mcmc.sample.tmpt) <- paste("trt2", "v.s.", "trt1",sep=" ")
            mcmc.sample.list1[[i]] <- coda::as.mcmc(mcmc.sample.tmpt)
          }
        } else{
          for(i in 1:nchain){
            mcmc.sample.list1[[i]] <- coda::as.mcmc(BUGS.rlts[,i,])
          }
        }

        mcmc.sample.list2 <- coda::as.mcmc.list(mcmc.sample.list1)
        rtn.gelman <- coda::gelman.diag(mcmc.sample.list2)

        grDevices::pdf(file=paste(file.dir,"\\",graphics.name,".pdf",sep=""))
        coda::gelman.plot(mcmc.sample.list2,auto.layout=F)
        graphics::plot(mcmc.sample.list2,trace=T,density=F,smooth=F,col=color.option[1:nchain],auto.layout=F)
        grDevices::dev.off()
        rtn.plot <- paste(file.dir,"\\",graphics.name,".pdf",sep="")
      }

      if(file.rm){
        base::file.remove(paste(file.dir,"\\",graphics.name,'.pdf',sep=""))
        diagnostic.rlts <- list(rtn.gelman)
        names(diagnostic.rlts) <- c("diagnostic statistics")
      } else{
        diagnostic.rlts <- list(rtn.gelman, rtn.plot)
        names(diagnostic.rlts) <- c("diagnostic statistics", "diagnostic plots")
      }

      ### mcmc summary results;
      mcmc.sample.pool_1 <- c()
      if(npar==1){
        for(i in 1:nchain){
          mcmc.sample.tmpt <- as.matrix(BUGS.rlts[,i,],ncol=1)
          colnames(mcmc.sample.tmpt) <- paste("trt2", "v.s.", "trt1",sep=" ")

          mcmc.sample.pool_1 <- rbind(mcmc.sample.pool_1,mcmc.sample.tmpt)
        }
      } else{
        for(i in 1:nchain){
          mcmc.sample.pool_1 <- rbind(mcmc.sample.pool_1,BUGS.rlts[,i,])
        }
      }
      mcmc.sample.pool_2 <- coda::as.mcmc(mcmc.sample.pool_1)

      if(npar==1){
        mcmc.sample.median <- base::summary(mcmc.sample.pool_2)$quantiles[3]
        mcmc.sample.sd <- base::summary(mcmc.sample.pool_2)$statistics[2]
      } else{
        mcmc.sample.median <- base::summary(mcmc.sample.pool_2)$quantiles[,3]
        mcmc.sample.sd <- base::summary(mcmc.sample.pool_2)$statistics[,2]
      }
      mcmc.sample.HPD <- coda::HPDinterval(mcmc.sample.pool_2, prob=cover.level)

      summary.rlts <- cbind(mcmc.sample.median,mcmc.sample.sd,mcmc.sample.HPD)
      colnames(summary.rlts) <- c("median","sd",paste(cover.level*100,"% HPD Lower",sep=""), paste(cover.level*100,"% HPD Upper",sep=""))
      row.names(summary.rlts) <- row.names(mcmc.sample.HPD)

      rtn.rlts <- list(DisPar,BUGS.rlts,diagnostic.rlts,summary.rlts)
      names(rtn.rlts) <- c("DisPar","mcmc.sample","mcmc.diagnostic", "mcmc.summary")
      class(rtn.rlts) <- "pprMod"
      return(rtn.rlts)
    }
  } else{
    stop("Distribution should be a character 'gaussian' or 'poisson' or 'binomial'!")
  }
}
