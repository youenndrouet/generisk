#' summarize object estimated with the generisk() function
#'
#' @param x an object returned by the generisk function
#' @param conf show bootstrapped 95% confidence intervals
#' @param ages at which ages to summarize ?
#' @param digit_abs digit for absolute risks
#' @param digit_rr digit for relative risks
#' @param digit_abs_agec digit for absolute risks by agec
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
                          ages = seq(10,70,10),
                          digit_abs = 1,
                          digit_rr = 1,
                          digit_abs_agec = 1){

  is.bootstrap <- ('boot' %in% names(x))
  qn <- qnorm(1-(1-conf)/2)

  TAB.absolute <- tab.absolute <- TAB.absolute.agec <- tab.absolute.agec <- NULL
  TAB.relative <- tab.relative <- TAB.relative.agec <- tab.relative.agec <- NULL

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

  age_class <- paste(c("0", as.character(ages)), c(as.character(ages), "120"), sep = "-")
  seqx.all <- 0:120

  while(ll >= 1){ #loop over diseases

    paramll  <- unique(mask[ll,])[-1]
    penetmod <- x$par$FIT.pars[[dis]]$penet.model
    disname  <- names(x$par$FIT.pars)[dis]
    ageminmax <- x$par$AgeDef[[dis]]
    seqx <- ageminmax[1]:ageminmax[2]

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

      mask_matched <- which(mask[ll,cc] == mask[ll,])

      labcurve  <- paste(disname, paste(colnames(mask)[mask_matched], collapse="/"),sep=": ")

      Ft        <- x$fit$ft[,sexe[cc],dis,geno[cc]][seqx + 1]

      if(length(unique(sexe[mask_matched]))==2){
        # 1 parameter for both sex
        Ftpop     <- rowMeans(x$fit$ft[,,dis,1])[seqx + 1]
        Ftpop.all <- rowMeans(x$fit$ft[,,dis,1])[seqx.all + 1]

      }else{
        Ftpop     <- x$fit$ft[,sexe[cc],dis,1][seqx + 1]
        Ftpop.all <- x$fit$ft[,sexe[cc],dis,1][seqx.all + 1]
      }

      Ft.all    <- x$fit$ft[,sexe[cc],dis,geno[cc]][seqx.all + 1]
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

         Ftb.agec <- NULL
         aa <- c(0, ages, 120)

         for (ii in (2:length(aa))){

           ftb.agec <- Ftb[aa[ii]+1,] - Ftb[aa[ii-1]+1,]
           Ftb.agec <- rbind(Ftb.agec, ftb.agec)
         }

         tra.agec  <- sapply(1:nrow(Ftb.agec), FUN = function(x){
           medix  <- median(Ftb.agec[x,])
           qlo    <- quantile(Ftb.agec[x,], (1-conf)/2)
           qup    <- quantile(Ftb.agec[x,], 1- ((1-conf)/2))
           return(c("median" = medix,
                    "qlo" = as.double(qlo),
                    "qup" = as.double(qup)))
         })

         colnames(tra.agec) <- age_class

      }

      labc <- paste(labcurve,paste('(',penetmod,')',sep=""), sep=' ')

      # absolute

      if(!is.bootstrap){

        # cumulative risks
        cumabs <- Ft.all[as.character(ages)]

        tab.absolute <- rbind(tab.absolute, c(labc, round(cumabs*100, digit_abs)))

        TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = labc,
                                                       "age" = ages,
                                                       "abs_risk" = cumabs))
        # risks by age class
        agecabs <- Ft.all[c(as.character(ages), "120")] - Ft.all[c("0", as.character(ages))]

        tab.absolute.agec <- rbind(tab.absolute.agec, c(labc, round(agecabs*100, digit_abs_agec)))

        TAB.absolute.agec <- rbind(TAB.absolute.agec, data.frame("outcome" = labc,
                                                                 "age_class" = age_class,
                                                                 "abs_risk" = agecabs))

      }else{

        rmed <- tra.all["median",as.character(ages)]
        binf <- tra.all['qlo',as.character(ages)]
        bsup <- tra.all['qup',as.character(ages)]

        icaa <- paste(round(binf*100, digit_abs), round(bsup*100, digit_abs), sep='-')
        stat <- paste(round(rmed*100, digit_abs)," (",icaa,")",sep="")

        tab.absolute <- rbind(tab.absolute, c(labc, stat))

        TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = labc,
                                                       "age" = ages,
                                                       "median_abs_risk" = rmed,
                                                       "lower_abs_risk" = binf,
                                                       "upper_abs_risk" = bsup))

        # risks by age class
        rmed <- tra.agec["median",] #tra.all["median", c(as.character(ages), "120")] - tra.all["median", c("0", as.character(ages))]
        binf <- tra.agec["qlo",]#tra.all["qlo", c(as.character(ages), "120")] - tra.all["qlo", c("0", as.character(ages))]
        bsup <- tra.agec["qup",]#tra.all["qup", c(as.character(ages), "120")] - tra.all["qup", c("0", as.character(ages))]

        icaa <- paste(round(binf*100, digit_abs_agec), round(bsup*100, digit_abs_agec), sep='-')
        stat <- paste(round(rmed*100, digit_abs_agec)," (",icaa,")",sep="")

        tab.absolute.agec <- rbind(tab.absolute.agec, c(labc, stat))

        TAB.absolute.agec <- rbind(TAB.absolute.agec, data.frame("outcome" = labc,
                                                                 "age_class" = age_class,
                                                                 "median_abs_risk" = rmed,
                                                                 "lower_abs_risk" = binf,
                                                                 "upper_abs_risk" = bsup))


      }


      # relative

      if(is.bootstrap){

          rrmed <- (tra.all["median",]/Ftpop.all)[as.character(ages)]
          binf <- (tra.all['qlo',]/Ftpop.all)[as.character(ages)]
          bsup <- (tra.all['qup',]/Ftpop.all)[as.character(ages)]

          icaa <- paste(round(binf, digit_rr), round(bsup, digit_rr), sep='-')
          stat <- paste(round(rrmed, digit_rr)," (",icaa,")",sep="")

          tab.relative <- rbind(tab.relative,c(paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),stat ))
          TAB.relative <- rbind(TAB.relative, data.frame("outcome" = paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),
                                                         "age" = ages,
                                                         "median_relative_risk" = rrmed,
                                                         "lower_relative_risk" = binf,
                                                         "upper_relative_risk" = bsup))

          # relative risks by age class
          rmed <- tra.agec["median",]
          binf <- tra.agec["qlo",]
          bsup <- tra.agec["qup",]

          agecabs.pop <- Ftpop.all[c(as.character(ages), "120")] - Ftpop.all[c("0", as.character(ages))]

          rrmed <- rmed/agecabs.pop
          rrinf <- binf/agecabs.pop
          rrsup <- bsup/agecabs.pop

          icaa <- paste(round(rrinf, digit_rr), round(rrsup, digit_rr), sep='-')
          stat <- paste(round(rrmed, digit_rr)," (",icaa,")",sep="")

          tab.relative.agec <- rbind(tab.relative.agec, c(labc,stat ))
          TAB.relative.agec <- rbind(TAB.relative.agec, data.frame("outcome" = labc,
                                                         "age_class" = age_class,
                                                         "median_relative_risk" = rrmed,
                                                         "lower_relative_risk" = rrinf,
                                                         "upper_relative_risk" = rrsup))


      }else{

          rr <- (Ft.all/Ftpop.all)[as.character(ages)]

          tab.relative <- rbind(tab.relative, c(labc, round(rr, digit_rr)))
          TAB.relative <- rbind(TAB.relative, data.frame("outcome" = labc,
                                                         "age" = ages,
                                                         "relative_risk" = rr))
          # relative risks by age class
          agecabs <- Ft.all[c(as.character(ages), "120")] - Ft.all[c("0", as.character(ages))]
          agecabs.pop <- Ftpop.all[c(as.character(ages), "120")] - Ftpop.all[c("0", as.character(ages))]

          rr <- agecabs/agecabs.pop

          tab.relative.agec <- rbind(tab.relative.agec, c(labc, round(rr, digit_rr)))
          TAB.relative.agec <- rbind(TAB.relative.agec, data.frame("outcome" = labc,
                                                         "age_class" = age_class,
                                                         "relative_risk" = rr))


      }

    }
    ll <- ll - np
    dis <- dis - 1
  }

  colnames(tab.absolute) <- colnames(tab.relative) <- c("strata", paste0("Age=",ages))
  colnames(tab.absolute.agec) <- colnames(tab.relative.agec) <- c("strata", paste0("Age=",age_class))

  cat('\n')
  cat('Model overall statistics: \n', lab, '\n')

  cat('\n')
  cat('Estimated risks: \n')

  return(list("ABSOLUTE_CUM_RISKS" = tab.absolute,
              "ABSOLUTE_RISKS_AGEC" = tab.absolute.agec,
              "RELATIVE_CUM_RISKS" = tab.relative,
              "RELATIVE_CUM_RISKS_AGEC" = tab.relative.agec,
              "data_ABSOLUTE_CUM_RISKS" = TAB.absolute,
              "data_ABSOLUTE_RISKS_AGEC" = TAB.absolute.agec,
              "data_RELATIVE_CUM_RISKS" = TAB.relative,
              "data_RELATIVE_RISKS_AGEC" = TAB.relative.agec))

}
