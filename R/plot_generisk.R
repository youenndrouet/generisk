#' plot penetrance curves estimated with the generisk() function
#'
#' @param x an object returned by the generisk function
#' @param type absolute or relative
#' @param conf show bootstrapped 95% confidence intervals
#' @param add add to existing graph
#' @param multiple multiple graphs in case of multiple phenotypes
#' @param ymax maximum value for y axis
#' @param xmax maximum value for x axis
#' @param allcurves show all bootstrapped penetrance curves
#' @param legendloc location of the legend
#' @param legendtext text of the legend
#' @param cols colors
#' @param main main title of the graph
#' @param small_labels logical, small or large labels for penetrance curves
#'
#' @import graphics
#' @importFrom stats qnorm
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom utils tail
#'
#' @return to be completed
#'
#' @examples
#' #to be completed
#'
#' @export
#'
plot_generisk <- function(x,
                          type = "absolute",
                          conf = 0.95,
                          add = FALSE,
                          multiple = FALSE,
                          ymax = 1,
                          xmax = 100,
                          allcurves = FALSE,
                          legendloc = "right",
                          legendtext = NULL,
                          cols = rep(c("#87CEFA","#FF0000","#00CD00",
                                     "#00CED1","#E9967A","#006400",
                                     "#FF8C00","#00008B","#8B008B",
                                     "#008B8B"),3),
                          main = "",
                          small_labels = TRUE){

  is.bootstrap <- ('boot' %in% names(x))
  qn <- qnorm(1-(1-conf)/2)
  mask  <- x$fit$paramsmask

  ncurves <- 1

  if(legendloc == "right"){
    marg <- c(5,6,5,20)
  }else{
    marg <- c(5,4,4,1)
    TEXTleg <- NULL
  }

  if(!add & !multiple){

    empty.plot(type, ymax=ymax, xmax=xmax, marg=marg, main=main)
    LKL <- round(x$fit$obj,1)
    suppl.par.np <- grep("120",rownames(x$fit$paramsmask),fixed = TRUE)

    if(any(suppl.par.np)){
      pars <-  x$fit$paramsmask[-suppl.par.np,]
    }else{
      pars <- x$fit$paramsmask
    }

    npar <- length(unique(pars[pars != 0]))
    message <- x$fit$message
    AIC <- 2*LKL + 2*npar
    lab <- paste("-2logL= ", 2*LKL," with ", npar," parameters (AIC= ",AIC,")\n ", message, sep="")

  }

  ll  <- nrow(mask)
  dis <- length(x$par$FIT.pars)

  ng  <- ncol(mask)/2
  sexe <- rep(c(1,2),each=ng)
  geno <- rep(1:ng, times = 2)

  while(ll >= 1){ #loop over diseases

    paramll  <- unique(mask[ll,])[-1]
    penetmod <- x$par$FIT.pars[[dis]]$penet.model
    disname  <- names(x$par$FIT.pars)[dis]
    ageminmax <- x$par$AgeDef[[dis]]

    seqx <- ageminmax[1]:ageminmax[2]
    seqx.all <- 0:120

    if(penetmod == "np"){

      agenodes <- x$par$FIT.pars[[dis]]$agenodes
      np <- length(agenodes) + 1
      penetmod <- paste0("NP nodes at: ", paste(agenodes, collapse = "/"))

    }else{
      if(penetmod == "Weibull"){
        np <- 3
      }else{
        if(penetmod == "Cox"){
          np <- 1
        }else{
          warning("Unknown penetrance model")
        }
      }
    }

    for (cc in match(paramll, mask[ll,])){

      labcurve  <- paste(disname, paste(colnames(mask)[which(mask[ll,cc] == mask[ll,])], collapse="/"),sep=": ")
      colk      <- ncurves
      Ft        <- x$fit$ft[,sexe[cc],dis,geno[cc]][seqx + 1]
      Ftpop     <- x$fit$ft[,sexe[cc],dis,1][seqx + 1]
      Ft.all    <- x$fit$ft[,sexe[cc],dis,geno[cc]][seqx.all + 1]
      Ftpop.all <- x$fit$ft[,sexe[cc],dis,1][seqx.all + 1]
      names(Ft) <- names(Ftpop) <- as.character(seqx)
      names(Ft.all) <- names(Ftpop.all) <- as.character(seqx.all)

      if(!add & multiple){
        empty.plot(type, ymax=ymax, xmax=xmax, marg=marg, main=main)
      }

      if(is.bootstrap){    # draw bootstrap IC
        Ftb    <- sapply(x$boot$ft.boot,FUN=function(x) x[,sexe[cc],dis,geno[cc]])
        tra.all    <- sapply(seqx.all,FUN = function(x){
          meanx  <- exp(mean(log(Ftb[x+1,])))
          medix  <- median(Ftb[x+1,])
          # sdx <- sd(Ftb[x+1,])
          # qlo.param    <- meanx - qn*sdx
          # qup.param    <- meanx + qn*sdx
          qlo    <- quantile(Ftb[x+1,], (1-conf)/2)
          qup    <- quantile(Ftb[x+1,], 1- ((1-conf)/2))
          return(c("mean" = meanx,"median" = medix, "qlo" = as.double(qlo), "qup" = as.double(qup)))
        })

        colnames(tra.all) <- seqx.all
        tra <- tra.all[, as.character(seqx)]

        CI.U <- tra['qup',]
        CI.L <- tra['qlo',]
        CI.U.all <- tra.all['qup',]
        CI.L.all <- tra.all['qlo',]
        X.Vec <- c(seqx, tail(seqx, 1), rev(seqx), seqx[1])
        X.Vec.all <- c(seqx.all, tail(seqx.all, 1), rev(seqx.all), seqx.all[1])
        if(type=="absolute"){
          lines(seqx.all, tra.all["median",], col=cols[colk], lwd=1, lty=2)
          lines(seqx, tra["median",], col=cols[colk], lwd=2)
          Y.Vec <- c(CI.L, tail(CI.U, 1), rev(CI.U), CI.L[1])
          Y.Vec.all <- c(CI.L.all, tail(CI.U.all, 1), rev(CI.U.all), CI.L.all[1])
        }else{
          Y.Vec <- c(CI.L/Ftpop, tail(CI.U/Ftpop, 1), rev(CI.U/Ftpop), (CI.L/Ftpop)[1])
          Y.Vec.all <- c(CI.L.all/Ftpop.all, tail(CI.U.all/Ftpop.all, 1), rev(CI.U.all/Ftpop.all), (CI.L.all/Ftpop.all)[1])
        }
        if(!allcurves){
          polygon(X.Vec.all, Y.Vec.all, col=colpct(cols[colk],30), border = NA)
        }
      }


      if(type=="absolute"){
        if(allcurves){

          apply(Ftb,2,FUN=function(x)lines(seqx.all, x, col=colpct("#000000",20)))
          lines(seqx.all, Ft.all, col=cols[colk], lwd=1, lty=2)
          lines(seqx, Ft, col="black", lwd=2)
          cols[colk] <- "#000000"

        }else{

          if(!is.bootstrap){

              lines(seqx.all, Ft.all, col=cols[colk], lwd=1, lty=2)
              lines(seqx, Ft, col=cols[colk], lwd=2)

          }

        }

        if(is.null(legendtext)){
          if(small_labels){
            textleg <- labcurve
          }else{
            textleg <- paste(labcurve,paste('(',penetmod,')',sep=""), sep=' ')
          }

        }else{
          textleg <- legendtext[colk]
        }

        if(legendloc=="right"){
          mtext(side = 4, line=0.5,
                at = ifelse(Ft.all[xmax + 1]>ymax, ymax, Ft.all[xmax + 1]),
                text= textleg,
                las=1,
                col=cols[colk],
                cex=0.8)
        }else{
          TEXTleg <- c(TEXTleg, textleg)
        }

      }else{ #type != "absolute"

        if(allcurves){

          apply(Ftb,2,FUN=function(x)lines(seqx.all, x/Ftpop, col=colpct("#000000",20)))
          lines(seqx.all, Ft.all/Ftpop.all, col=cols[colk], lwd=1, lty=2)
          lines(seqx, Ft/Ftpop, col="black", lwd=2)
          cols[colk] <- "#000000"

        }else{

          if(is.bootstrap){
            lines(seqx.all, tra.all["median",]/Ftpop.all, col=cols[colk], lwd=1, lty=2)
            lines(seqx, tra["median",]/Ftpop, col=cols[colk], lwd=2)
          }else{
            lines(seqx.all, Ft.all/Ftpop.all, col=cols[colk], lwd=1, lty=2)
            lines(seqx, Ft/Ftpop, col=cols[colk], lwd=2)
          }

        }

        if(is.null(legendtext)){
          if(small_labels){
            textleg <- labcurve
          }else{
            textleg <- paste(labcurve,paste('(',penetmod,')',sep=""), sep=' ')
          }
        }else{
          textleg <- legendtext[colk]
        }
        if(legendloc=="right"){
          mtext(side = 4,
                line=0.5,
                at = (Ft.all/Ftpop.all)[xmax + 1],
                text= textleg,
                las=1,
                col=cols[colk],
                cex=0.8)
        }else{
          TEXTleg <- c(TEXTleg, textleg)
        }

        abline(h=1,lty=2)
      }
      ncurves <- ncurves + 1
    }
    ll <- ll - np
    dis <- dis - 1
  }

  if(legendloc!="right"){
    legend(x="topleft",legend=TEXTleg, lty=1, col=cols[1:colk],lwd=2)
  }

}


empty.plot <- function(type="absolute", ymax, xmax, marg, main){

  par(mar=marg,las=1)
  if(type=="absolute"){
    plot(NA,NA, type="n", xlim=c(0,xmax), ylim=c(0,ymax), ann=FALSE, bty="n",axes=FALSE)
    #rect(0,0,100,1, col=gray(.96))
    abline(h = seq(0,1,0.02), v = seq(10,90,10),col=gray(.95),lty=3 )
    abline(h = seq(0,1,0.1) , v = seq(0,90,20), col=gray(.85),lty=3)
    axis(2,lwd=1,col = gray(.7), at = seq(0,1,0.1), labels=seq(0,1,0.1)*100)
    title(ylab="",main="Cumulative risk (%)", cex.lab = 1.3)
  }else{
    plot(NA,NA, type="n", xlim=c(0,xmax), ylim=c(0.1,100000), log="y", xlab="", ylab ="", bty="n",axes=FALSE)
    abline(h = c(seq(0.2,0.9,0.1),
                 seq(2,9,1),
                 seq(20,90,10),
                 seq(200,900,100),
                 seq(2000,9000,1000),
                 seq(20000,90000,10000)),
           v = seq(10,90,10),
           col=gray(.95),lty=3 )
    abline(h = c(.01,.1,10,100,1000,10000), v = seq(0,90,20),col=gray(.85),lty=3)
    axis(2, at=c(0.01,0.1,1,10,100,1000,1e4,1e5 ), labels=c(".01", ".1","1","10","100","1000", 1e4, 1e5),lwd=1,col = gray(.7))
    title(ylab="",main="Risk ratio to pop. incidence \n (log-scale)", cex.lab = 1.3)
  }
  title(xlab="Age (y)", cex.lab = 1.3)
  axis(1,lwd=1,col = gray(.7))

}

colpct <- function(col,pct) return(paste(col,pct,sep=""))

