#' Generate power prior method BUGS file
#'
#'\code{pprmodelfile} returns power prior method BUGS file
#'
#' @param ntrt number of arms.
#' @param family outcome distribution.
#' @param mypprfile the power prior method BUGS file.
#' @return power prior method BUGS file
#' @examples
#' \dontrun{
#' pprmodelfile(4,"gaussian",paste(tempdir(),"\\","pprmodelfile.txt",sep=""))
#' pprmodelfile(3,"poisson",paste(tempdir(),"\\","pprmodelfile.txt",sep=""))
#' pprmodelfile(3,"binomial",paste(tempdir(),"\\","pprmodelfile.txt",sep=""))
#' file.show(paste(tempdir(),"\\","pprmodelfile.txt",sep=""))
#' }
#' @export

pprmodelfile <- function(ntrt,family,mypprfile){

  str_1 <- ""
  str_2 <- ""
  for (i in 1:(ntrt-1))
  {
    str_1 <- paste(str_1, paste("beta", i, "[1]*trt", i, "[i] + ", sep = ""), sep = "")
    str_2 <- paste(str_2, paste(" + beta", i, "[2]*trt", i, "[i]", sep = ""), sep = "")
  }
  str_1 <- paste("beta0 + ", str_1 , "b[cluster[i]]", sep = "")
  str_2 <- paste("beta0.star", str_2, sep="")

  sink(mypprfile)

  cat("model{ \n")
  cat("   for(i in 1:N){ \n")

  if (identical(family,"gaussian")){
    cat("       y[i] ~ dnorm(mu[i],tau_1) \n")
    cat("       \r")
    cat(paste("mu[i] <- ",str_1,sep=""),fill=TRUE)
  } else if (identical(family,"binomial")){
    cat("       y[i] ~ dbin(p[i],1) \n")
    cat("       \r")
    cat(paste("logit(p[i]) <- ",str_1,sep=""),fill=TRUE)
  } else if (identical(family,"poisson")){
    cat("       y[i] ~ dpois(lambda[i]) \n")
    cat("       \r")
    cat(paste("log(lambda[i]) <- ",str_1,sep=""),fill=TRUE)
  }

  cat("       } \n")
  cat("   M <- 100000 \n")
  cat("   for(i in (N+1):(N+I)){ \n")

  if (identical(family,"gaussian")){
    cat("       \r")
    cat(paste("mu[i] <- ",str_2,sep=""),fill=TRUE)
    cat("       m[i-N] <- mu[i]*n[i-N] \n")
    cat("       t[i-N] <- tau_2*pow(n[i-N],-1) \n")
    cat("       l[i-N] <- a*0.5*(log(t[i-N])-t[i-N]*pow(y[i]-m[i-N],2)) \n")
  } else if (identical(family,"binomial")){
    cat("       \r")
    cat(paste("logit(p[i]) <- ", str_2, sep=""),fill=TRUE)
    cat("       l[i-N] <- a*(y[i]*log(p[i])+(n[i-N]-y[i])*log(1-p[i])) \n")
  } else if (identical(family,"poisson")){
    cat("       \r")
    cat(paste("log(lambda[i]) <- ", str_2, sep=""),fill=TRUE)
    cat("       lda[i-N] <- lambda[i]*n[i-N] \n")
    cat("       l[i-N] <- a*(y[i]*log(lda[i-N])-lda[i-N]-logfact(y[i])) \n")
  }

  cat("       zeros[i-N] <- 0 \n")
  cat("       phi[i-N] <- -l[i-N]+M \n")
  cat("       zeros[i-N] ~ dpois(phi[i-N]) \n")
  cat("       } \n")

  if(ntrt > 2){

    for(i in 1:(ntrt-2)){
      for(j in (i+1):(ntrt-1)){
        cat("   \r")
        cat(paste("beta",i,j," <- beta",i,"[1] - beta",j,"[1]",sep=""),fill=TRUE)
      }
    }

  }

  cat("   for(k in 1:K){ \n")
  cat("       b[k] ~ dnorm(0,tau.b) \n")
  cat("       }\n")
  cat("   beta0 ~ dnorm(0,0.1) \n")
  cat("   beta0.star ~ dnorm(0,0.1) \n")

  cat("   tau.b <- pow(sigma.b,-2) \n")
  cat("   sigma.b ~ dunif(0,5) \n")

  if (identical(family,"gaussian")){
    cat("   tau_1 <- pow(sigma_1,-2) \n")
    cat("   sigma_1 ~ dunif(0,5) \n")
    cat("   tau_2 <- pow(sigma_2,-2) \n")
    cat("   sigma_2 ~ dunif(0,5) \n")
  }

  for(i in 1:(ntrt-1)){
    cat("\n")
    cat("   \r")
    cat(paste("beta",i,"[1:2]"," ~ dmnorm(", paste("M.h",i,"[]",sep=""),",",paste("Tau.h",i,"[,]",sep=""),")",sep=""),fill=TRUE)
    cat("   \r")
    cat(paste("M.h",i,"[1]"," <- mu.h",i,sep=""),fill=TRUE)
    cat("   \r")
    cat(paste("M.h",i,"[2]"," <- mu.h",i,sep=""),fill=TRUE)
    cat("   \r")
    cat(paste("Tau.h",i,"[1,1]"," <- tau.h",i,"/(1-pow(rho.h",i,",2))",sep=""),fill=TRUE)
    cat("   \r")
    cat(paste("Tau.h",i,"[1,2]"," <- -tau.h",i,"*rho.h",i,"/(1-pow(rho.h",i,",2))",sep=""),fill=TRUE)
    cat("   \r")
    cat(paste("Tau.h",i,"[2,1]"," <- -tau.h",i,"*rho.h",i,"/(1-pow(rho.h",i,",2))",sep=""),fill=TRUE)
    cat("   \r")
    cat(paste("Tau.h",i,"[2,2]"," <- tau.h",i,"/(1-pow(rho.h",i,",2))",sep=""),fill=TRUE)

    cat("\n")
    cat("   \r")
    cat(paste("mu.h",i," ~ dnorm(0,0.01)",sep=""),fill=TRUE)
    cat("   \r")
    cat(paste("tau.h",i," <- pow(sigma.h",i,",-2)",sep=""),fill=TRUE)
    cat("   \r")

    if (identical(family,"gaussian")){
      cat(paste("sigma.h",i," ~ dunif(0,10)",sep=""),fill=TRUE)
      cat("   \r")
    } else if (identical(family,"binomial")){
      cat(paste("sigma.h",i," ~ dunif(0,10)",sep=""),fill=TRUE)
      cat("   \r")
    } else if (identical(family,"poisson")){
      cat(paste("sigma.h",i," ~ dunif(0,5)",sep=""),fill=TRUE)
      cat("   \r")
    }

    cat(paste("rho.h",i," ~ dunif(0,1)",sep=""),fill=TRUE)
  }

  cat("} \n")

  sink()
}

