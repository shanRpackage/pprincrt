#' print pprMod object, returned by pprmodelBUGS function
#'
#'\code{pprmodelBUGS} returns summary and convergence diagnositc results of MCMC chains
#' @param x pprMod object. A S3 object defined by this package developer.
#' @param open.plot An indicator of opening the convergence diagnostic plots. It is TRUE by default.
#' @param \dots other arguments
#' @return summary and convergence diagnositc results of MCMC chains
#' @method print pprMod
#' @export

print.pprMod <- function(x, open.plot=TRUE, ...){
  base::cat("Discounting parameter value: \n")
  base::print(x$DisPar)
  base::cat("Summay of mcmc chain: \n")
  base::print(x$mcmc.summary)
  base::cat("Convergence diagnostic statistics for mcmc chain: \n")
  base::print(x$mcmc.diagnostic[[1]])
  if(length(x$mcmc.diagnostic)==2){
    base::cat("Convergence diagnostic plots for mcmc chain: \n")
    base::print(x$mcmc.diagnostic[[2]])
    if (open.plot) {base::system2('open',args=x$mcmc.diagnostic[[2]],wait=F)}
  }
}
