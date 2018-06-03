#' Generate current model BUGS file
#'
#'\code{cmodelfile} returns current model BUGS file
#'
#' @param ntrt number of arms.
#' @param family outcome distribution.
#' @param mycfile the current model BUGS file.
#' @return current model BUGS file
#' @examples
#' \dontrun{
#' cmodelfile(4,"binomial",paste(tempdir(),"\\","cmodelfile.txt",sep=""))
#' cmodelfile(4,"gaussian",paste(tempdir(),"\\","cmodelfile.txt",sep=""))
#' cmodelfile(4,"poisson",paste(tempdir(),"\\","cmodelfile.txt",sep=""))
#' file.show(paste(tempdir(),"\\","cmodelfile.txt",sep=""))
#' }
#' @export

cmodelfile <- function(ntrt,family,mycfile){

  str <- ""
  for (i in 1:(ntrt-1)){
    str <- paste(str, paste("beta", i, "*trt", i, "[i] + ", sep = ""), sep = "")
  }
  str <- paste("beta0 + ", str , "b[cluster[i]]", sep = "")

  sink(mycfile)

  cat("model{ \n")
  cat("   for(i in 1:N){ \n")

  if (identical(family,"gaussian")){
    cat("      y[i] ~ dnorm(mu[i],tau) \n")
    cat("      \r")
    cat(paste("mu[i] <- ",str,sep=""),fill=TRUE)
  } else if (identical(family,"binomial")){
    cat("      y[i] ~ dbin(p[i],1) \n")
    cat("      \r")
    cat(paste("logit(p[i]) <- ",str,sep=""),fill=TRUE)
  } else if (identical(family,"poisson")){
    cat("      y[i] ~ dpois(lambda[i]) \n")
    cat("      \r")
    cat(paste("log(lambda[i]) <- ",str,sep=""),fill=TRUE)
  }

  cat("   } \n")

  cat("   for(k in 1:K){ \n")
  cat("      b[k] ~ dnorm(0,tau.b) \n")
  cat("   } \n")

  cat("   beta0 ~ dnorm(0,0.1) \n")
  for(i in 1:(ntrt-1)){
    cat("   \r")
    cat(paste("beta",i," ~ dnorm(0,0.1)",sep=""),fill=TRUE)
  }

  cat("   tau.b <- pow(sigma.b,-2) \n")
  cat("   sigma.b ~ dunif(0,5) \n")

  if (identical(family,"gaussian")){
    cat("   tau <- pow(sigma,-2) \n")
    cat("   sigma ~ dunif(0,5) \n")
  }

  cat("} \n")

  sink()
}

