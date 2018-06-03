#' Animation Plot
#'
#'\code{AniPlot} returns animation plots
#'
#' @param x data to plot. It should be in data frame.
#' @param GroupVar the group variable name.
#' @param FrameVar the frame variable name.
#' @param XVar the x-axis variable name.
#' @param YVar the y-axis variable name.
#' @param play.int time interval to play each frame
#' @param AniPlot.dir directory to save the animation plot in PDF file
#' @param plotcaption animation plot caption
#' @param plotwidth animation plot width
#' @param plotheight animation plot height
#' @param imagename animation image file name
#' @param plotmain animation plot title
#' @return Animation plot
#' @examples
#' \dontrun{
#' AniPlot(anidata)
#' }
#' @export

AniPlot <- function(x,GroupVar='group', FrameVar='frame', XVar='x', YVar='y', play.int=0.5, AniPlot.dir=tempdir(),
                    plotcaption="My Animation Plot",plotwidth=7,plotheight=7,imagename="My_Animation_Image",plotmain=NULL){

  x.order <- x[order(x[,GroupVar],x[,FrameVar],x[,XVar]),]
  unigroupvar <- unique(x.order[,GroupVar])

  x.short <- x.order[which(x.order[,GroupVar]==unigroupvar[1]),]
  if(length(unigroupvar) > 1){
    for(k in 2:length(unigroupvar)){
      x.short <- cbind(x.short,x.order[which(x.order[,GroupVar]==unigroupvar[k]),][,YVar])
    }
  }
  x.short <- data.frame(x.short[,-c(1)])
  colnames(x.short) <- c(FrameVar,XVar,paste(YVar,1:length(unigroupvar),sep=""))

  nframe <- length(unique(x.short[,FrameVar]))
  ncolor <- dim(x.short)[2]-2
  color.option <- c(2,4,3,5:(ncolor+5))

  nstep <- length(unique(x.short[,XVar]))
  start.index <- seq(from=1,to=dim(x.short)[1],by=nstep)
  end.index <- seq(from=nstep,to=dim(x.short)[1],by=nstep)

  uniframevar <- unique(x.short[,FrameVar])

  base::setwd(AniPlot.dir)
  animation::saveLatex({
    oopt <- animation::ani.options(interval = play.int, nmax = nframe)
    for (i in 1:animation::ani.options("nmax")) {
      grDevices::dev.hold()

      if(length(unigroupvar)==1){
        graphics::plot(x.short[start.index[i]:end.index[i],XVar],x.short[start.index[i]:end.index[i],paste(YVar,1,sep="")],type='b',main=plotmain,xlab=XVar,ylab=YVar,
                       ylim=c(min(x.short[,paste(YVar,1:length(unigroupvar),sep="")]),max(x.short[,paste(YVar,1:length(unigroupvar),sep="")])),col=color.option[1])
        graphics::legend("topright",legend=paste(GroupVar,'=',unigroupvar),col=color.option[1:length(unigroupvar)],lty=rep(1,length(unigroupvar)),cex=0.8)
        graphics::legend("topleft",legend = substitute(paste(f,'=',b),list(f=FrameVar,b=sprintf("%.03f",uniframevar[i]))), bty = "n")
      } else{
        graphics::plot(x.short[start.index[i]:end.index[i],XVar],x.short[start.index[i]:end.index[i],paste(YVar,1,sep="")],type='b',main=plotmain,xlab=XVar,ylab=YVar,
                       ylim=c(min(x.short[,paste(YVar,1:length(unigroupvar),sep="")]),max(x.short[,paste(YVar,1:length(unigroupvar),sep="")])),col=color.option[1])
        for(j in 2:ncolor){
          graphics::lines(x.short[start.index[i]:end.index[i],XVar],x.short[start.index[i]:end.index[i],paste(YVar,j,sep="")],col=color.option[j])
          graphics::points(x.short[start.index[i]:end.index[i],XVar],x.short[start.index[i]:end.index[i],paste(YVar,j,sep="")],col=color.option[j])
        }

        graphics::legend("topright",legend=paste(GroupVar,'=',unigroupvar),col=color.option[1:length(unigroupvar)],lty=rep(1,length(unigroupvar)),cex=0.8)
        graphics::legend("topleft",legend = substitute(paste(f,'=',b),list(f=FrameVar,b=sprintf("%.03f",uniframevar[i]))), bty = "n")
      }
      animation::ani.pause()
    }},caption=plotcaption,ani.dev = "pdf", ani.type = "pdf", ani.width = plotwidth, ani.height = plotheight, img.name=imagename)

  animation::ani.options(oopt)
}

