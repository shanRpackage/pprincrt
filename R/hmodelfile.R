#' Generate historical model BUGS file
#'
#'\code{hmodelfile} returns historical model BUGS file
#'
#' @param ntrt number of arms.
#' @param family outcome distribution.
#' @param myhfile the historical model BUGS file.
#' @return historical model BUGS file
#' @examples
#' \dontrun{
#' hmodelfile(2,"binomial",paste(tempdir(),"\\","hmodelfile.txt",sep=""))
#' hmodelfile(3,"gaussian",paste(tempdir(),"\\","hmodelfile.txt",sep=""))
#' hmodelfile(4,"poisson",paste(tempdir(),"\\","hmodelfile.txt",sep=""))
#' file.show(paste(tempdir(),"\\","hmodelfile.txt",sep=""))
#' }
#' @export

hmodelfile <- function(ntrt,family,myhfile){

  str <- ""
  for (i in 1:(ntrt-1)){
    str <- paste(str, paste(" + beta", i, "*trt", i, "[i]", sep = ""), sep = "")
  }
  str <- paste("beta0", str, sep = "")

  sink(myhfile)

  cat("model{ \n")
  cat("   for(i in 1:I){ \n")

  if (identical(family,"gaussian"))
  {
    cat("        y[i] ~ dnorm(m[i],t[i]) \n")
    cat("        m[i] <- mu[i]*n[i] \n")
    cat("        t[i] <- tau*pow(n[i],-1) \n")
    cat("        \r")
    cat(paste("mu[i] <- ",str,sep=""),fill=TRUE)
  }
  else if (identical(family,"binomial"))
  {
    cat("       y[i] ~ dbin(p[i],n[i]) \n")
    cat("       \r")
    cat(paste("logit(p[i]) <- ",str,sep=""),fill=TRUE)
  }
  else if (identical(family,"poisson"))
  {
    cat("       y[i] ~ dpois(lda[i]) \n")
    cat("       lda[i] <- lambda[i]*n[i] \n")
    cat("       \r")
    cat(paste("log(lambda[i]) <- ",str,sep=""),fill=TRUE)
  }

  cat("     } \n")

  cat("     beta0 ~ dnorm(0,0.1) \n")

  for(i in 1:(ntrt-1)){
    cat("     \r")
    cat(paste("beta",i," ~ dnorm(0,0.1)",sep=""),fill=TRUE)
  }

  if (identical(family,"gaussian"))
  {
    cat("     tau <- pow(sigma,-2) \n")
    cat("     sigma ~ dunif(0,5) \n")
  }

  cat("} \n")

  sink()
}
