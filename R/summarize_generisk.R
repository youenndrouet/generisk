#' summarize object estimated with the generisk() function
#'
#' @param x an object returned by the generisk function
#' @param conf show bootstrapped 95% confidence intervals
#' @param ages at which ages to summarize ?
#'
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
summarize_generisk <- function(x,
                          conf = 0.95,
                          ages = seq(10,70,10)){

  is.bootstrap <- ('boot' %in% names(x))
  qn <- qnorm(1-(1-conf)/2)

  TAB.absolute <- tab.absolute <- NULL
  TAB.relative <- tab.relative <- NULL

  mask  <- x$fit$paramsmask

  LKL <- round(x$fit$obj,1)
  npar <- length(unique(mask[mask != 0]))
  message <- x$fit$message
  AIC <- 2*LKL + 2*npar
  lab <- paste("-2logL= ", 2*LKL," with ", npar," parameters (AIC= ",AIC,")\n ", message, sep="")

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

      Ft        <- x$fit$ft[,sexe[cc],dis,geno[cc]][seqx + 1]
      Ftpop     <- x$fit$ft[,sexe[cc],dis,1][seqx + 1]
      Ft.all    <- x$fit$ft[,sexe[cc],dis,geno[cc]][seqx.all + 1]
      Ftpop.all <- x$fit$ft[,sexe[cc],dis,1][seqx.all + 1]
      names(Ft) <- names(Ftpop) <- as.character(seqx)
      names(Ft.all) <- names(Ftpop.all) <- as.character(seqx.all)

      if(is.bootstrap){    # compute bootstrap IC

         Ftb    <- sapply(x$boot$ft.boot, FUN=function(x) x[,sexe[cc],dis,geno[cc]])
         tra.all  <- sapply(seqx.all, FUN = function(x){
                              meanx  <- exp(mean(log(Ftb[x+1,])))
                              medix  <- median(Ftb[x+1,])
                              qlo    <- quantile(Ftb[x+1,], (1-conf)/2)
                              qup    <- quantile(Ftb[x+1,], 1- ((1-conf)/2))
                              return(c("mean" = meanx,
                                       "median" = medix,
                                       "qlo" = as.double(qlo),
                                       "qup" = as.double(qup)))
                            })

         colnames(tra.all) <- seqx.all
      }

      # absolute

      if(!is.bootstrap){

        tab.absolute <- rbind(tab.absolute, c(paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '), round(Ft.all[as.character(ages)]*100,1)))

        TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),
                                                       "age" = ages,
                                                       "abs_risk" = Ft.all[as.character(ages)]))

      }else{

        icaa <- paste(round(tra.all['qlo',as.character(ages)]*100,1), round(tra.all['qup',as.character(ages)]*100,1), sep='-')
        stat <- paste(round(tra.all["median",as.character(ages)]*100,1)," (",icaa,")",sep="")
        tab.absolute <- rbind(tab.absolute,c(paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),stat ))

        TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),
                                                       "age" = ages,
                                                       "median_abs_risk" =  tra.all["median",as.character(ages)],
                                                       "lower_abs_risk" = tra.all['qlo',as.character(ages)],
                                                       "upper_abs_risk" = tra.all['qup',as.character(ages)]))

      }


      # relative

      if(is.bootstrap){
          icaa <- paste(round((tra.all['qlo',]/Ftpop.all)[as.character(ages)],1), round((tra.all['qup',]/Ftpop.all)[as.character(ages)],1), sep='-')
          stat <- paste(round((tra.all["median",]/Ftpop.all)[as.character(ages)],1)," (",icaa,")",sep="")
          tab.relative <- rbind(tab.relative,c(paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),stat ))
          TAB.relative <- rbind(TAB.relative, data.frame("outcome" = paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),
                                                         "age" = ages,
                                                         "median_relative_risk" = (tra.all["median",]/Ftpop.all)[as.character(ages)],
                                                         "lower_relative_risk" = (tra.all['qlo',]/Ftpop.all)[as.character(ages)],
                                                         "upper_relative_risk" = (tra.all['qup',]/Ftpop.all)[as.character(ages)]))


      }else{
          tab.relative <- rbind(tab.relative,c(paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '), round((Ft.all/Ftpop.all)[as.character(ages)],1)))
          TAB.relative <- rbind(TAB.relative, data.frame("outcome" = paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),
                                                         "age" = ages,
                                                         "relative_risk" = (Ft.all/Ftpop.all)[as.character(ages)]))

      }

    }
    ll <- ll - np
    dis <- dis - 1
  }

  colnames(tab.absolute) <- colnames(tab.relative) <- c("strata", paste0("Age=",ages))

  cat('\n')
  cat('Model overall statistics: \n', lab, '\n')

  cat('\n')
  cat('Estimated risks: \n')

  return(list("ABSOLUTE_CUM_RISKS" = tab.absolute,
              "RELATIVE_CUM_RISKS" = tab.relative,
              "TAB_data_abs" = TAB.absolute,
              "TAB_data_relative" = TAB.relative))

}
