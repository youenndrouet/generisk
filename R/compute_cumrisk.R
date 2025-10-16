#' compute cumulative risk over multiple traits from a bootstrapped generisk call object
#'
#' @param x an object returned by the generisk function
#' @param conf show bootstrapped 95% confidence intervals
#' @param ages at which ages ?
#' @param digit_abs digit for absolute risks
#' @param digit_rr digit for relative risks
#' @param digit_abs_agec digit for absolute risks by agec
#' @param traits over which traits ?
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
compute_cumrisk <- function(x,
                          conf = 0.95,
                          ages = seq(30,70,10),
                          traits = names(x$params$FIT.pars),
                          digit_abs = 1,
                          digit_rr = 1,
                          digit_abs_agec = 1){

  is.bootstrap <- ('boot' %in% names(x))
  qn <- qnorm(1-(1-conf)/2)

  TAB.absolute <- tab.absolute <- TAB.absolute.agec <- tab.absolute.agec <- NULL
  TAB.relative <- tab.relative <- TAB.relative.agec <- tab.relative.agec <- NULL

  mask  <- x$fit$paramsmask

  ll  <- nrow(mask)
  dis <- length(x$params$FIT.pars)
  ng  <- ncol(mask)/2
  sexe <- rep(c(1,2),each=ng)
  geno <- rep(1:ng, times = 2)

  age_class <- paste(c("0", as.character(ages)), c(as.character(ages), "120"), sep = "-")
  seqx.all <- 0:120

  if(!is.bootstrap){
    stop("x must be a bootstrapped generisk object")

  }else{

    # compute cumulative risk over selected traits

    if(length(traits)<=1){

      stop("traits length must be > 1")

    }else{

      one_minus_Ftb_prod_males <- one_minus_Ftb_prod_females <- 1
      one_minus_Ftpop_prod_males <- one_minus_Ftpop_prod_females  <- 1

      traits_male <- FALSE
      traits_female <- FALSE

      while(ll >= 1){ #loop over diseases (starting from end)

        paramll  <- unique(mask[ll,])[-1]
        penetmod <- x$params$FIT.pars[[dis]]$penet.model
        disname  <- names(x$params$FIT.pars)[dis]

        if(penetmod == "np"){

          agenodes <- x$params$FIT.pars[[dis]]$agenodes
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
         labc <- paste(labcurve,paste('(',penetmod,')',sep=""), sep=' ')

         Ftb  <- sapply(x$boot$ft.boot, FUN=function(boot) boot[,sexe[cc],dis,geno[cc]])

         if(disname %in% traits){

           # compute Pr(at least one trait) = prod_over_traits(1-pr(risk_trait_j))

           if(length(unique(sexe[mask_matched]))==2){

             # 1 parameter for both sex
             # Ftb is identical for males and females for this trait

             traits_male <- TRUE
             traits_female <- TRUE
             Ftpop_males <- x$fit$ft[,1,dis,1][seqx.all + 1]
             Ftpop_females <- x$fit$ft[,2,dis,1][seqx.all + 1]

             one_minus_Ftb_prod_males <- one_minus_Ftb_prod_males*(1-Ftb)
             one_minus_Ftpop_prod_males <- one_minus_Ftpop_prod_males * (1-Ftpop_males)
             one_minus_Ftb_prod_females <- one_minus_Ftb_prod_females*(1-Ftb)
             one_minus_Ftpop_prod_females <- one_minus_Ftpop_prod_females * (1-Ftpop_females)

           }else{
             # sex specific risks..
             if(sexe[cc] == 1){
               #males
               traits_male <- TRUE
               Ftpop_males <- x$fit$ft[,1,dis,1][seqx.all + 1]
               one_minus_Ftb_prod_males <- one_minus_Ftb_prod_males*(1-Ftb)
               one_minus_Ftpop_prod_males <- one_minus_Ftpop_prod_males * (1-Ftpop_males)
             }else{
               #females
               traits_female <- TRUE
               Ftpop_females <- x$fit$ft[,2,dis,1][seqx.all + 1]
               one_minus_Ftb_prod_females <- one_minus_Ftb_prod_females*(1-Ftb)
               one_minus_Ftpop_prod_females <- one_minus_Ftpop_prod_females * (1-Ftpop_females)
             }
           }

         }else{

           # summarize_generisk

           tra.all  <- sapply(seqx.all, FUN = function(x){
                              medix  <- median(Ftb[x+1,])
                              qlo    <- quantile(Ftb[x+1,], (1-conf)/2)
                              qup    <- quantile(Ftb[x+1,], 1- ((1-conf)/2))
                              return(c("median" = medix,
                                       "qlo" = as.double(qlo),
                                       "qup" = as.double(qup)))
                            })

           colnames(tra.all) <- seqx.all

           ## agec

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


           ## fin agec

           rmed <- tra.all["median",as.character(ages)]
           binf <- tra.all['qlo',as.character(ages)]
           bsup <- tra.all['qup',as.character(ages)]

           icaa <- paste(round(binf*100, digit_abs), round(bsup*100, digit_abs), sep='-')
           stat <- paste(round(rmed*100, digit_abs)," (",icaa,")",sep="")

           tab.absolute <- rbind(tab.absolute, c(labc, stat))

          if(length(unique(sexe[mask_matched]))==2){
            # 1 parameter for both sex
            Ftpop.all <- rowMeans(x$fit$ft[,,dis,1])
          }else{
            Ftpop.all <- x$fit$ft[,sexe[cc],dis,1][seqx.all + 1]
          }

          names(Ftpop.all) <- seqx.all

          TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = labc,
                                                         "age" = ages,
                                                         "median_abs_risk" = rmed,
                                                         "lower_abs_risk" = binf,
                                                         "upper_abs_risk" = bsup,
                                                         "pop_risk" = Ftpop.all[ages + 1]))
          # risks by age class
          rmed <- tra.agec["median",]
          binf <- tra.agec["qlo",]
          bsup <- tra.agec["qup",]

          icaa <- paste(round(binf*100, digit_abs_agec), round(bsup*100, digit_abs_agec), sep='-')
          stat <- paste(round(rmed*100, digit_abs_agec)," (",icaa,")",sep="")

          tab.absolute.agec <- rbind(tab.absolute.agec, c(labc, stat))

          TAB.absolute.agec <- rbind(TAB.absolute.agec, data.frame("outcome" = labc,
                                                                   "age_class" = age_class,
                                                                   "median_abs_risk" = rmed,
                                                                   "lower_abs_risk" = binf,
                                                                   "upper_abs_risk" = bsup))




      # relative
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



      }

    }
    ll <- ll - np
    dis <- dis - 1
  }

  #now we compute the prod..


      if(traits_male){

        ## risk of at leat one trait is equal to 1-prod(Pr(no trait))
        Ftb_prod_males <- (1-one_minus_Ftb_prod_males)
        Ftpop_prod_males <- (1-one_minus_Ftpop_prod_males)

        names(Ftpop_prod_males) <- seqx.all

        tra.all_males  <- sapply(seqx.all, FUN = function(x){
          medix  <- median(Ftb_prod_males[x+1,])
          qlo    <- quantile(Ftb_prod_males[x+1,], (1-conf)/2)
          qup    <- quantile(Ftb_prod_males[x+1,], 1- ((1-conf)/2))
          return(c("median" = medix,
                   "qlo" = as.double(qlo),
                   "qup" = as.double(qup)))
        })

        colnames(tra.all_males) <- seqx.all

        ##agec
        Ftb.agec <- NULL
        aa <- c(0, ages, 120)

        for (ii in (2:length(aa))){

          ftb.agec <- Ftb_prod_males[aa[ii]+1,] - Ftb_prod_males[aa[ii-1]+1,]
          Ftb.agec <- rbind(Ftb.agec, ftb.agec)
        }

        tra.agec_males  <- sapply(1:nrow(Ftb.agec), FUN = function(x){
          medix  <- median(Ftb.agec[x,])
          qlo    <- quantile(Ftb.agec[x,], (1-conf)/2)
          qup    <- quantile(Ftb.agec[x,], 1- ((1-conf)/2))
          return(c("median" = medix,
                   "qlo" = as.double(qlo),
                   "qup" = as.double(qup)))
        })

        colnames(tra.agec_males) <- age_class

        ## fin agec

        label_cum_males <- paste(paste(traits, collapse = "|"), "males", sep = ' ')

        rmed <- tra.all_males["median",as.character(ages)]
        binf <- tra.all_males['qlo',as.character(ages)]
        bsup <- tra.all_males['qup',as.character(ages)]

        icaa <- paste(round(binf*100, digit_abs), round(bsup*100, digit_abs), sep='-')
        stat <- paste(round(rmed*100, digit_abs)," (",icaa,")",sep="")

        tab.absolute <- rbind(tab.absolute,c(label_cum_males, stat ))

        TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = label_cum_males,
                                                       "age" = ages,
                                                       "median_abs_risk" =  tra.all_males["median",as.character(ages)],
                                                       "lower_abs_risk" = tra.all_males['qlo',as.character(ages)],
                                                       "upper_abs_risk" = tra.all_males['qup',as.character(ages)],
                                                       "pop_risk" = Ftpop_prod_males[ages + 1]))

        # risks by age class
        rmed <- tra.agec_males["median",]
        binf <- tra.agec_males["qlo",]
        bsup <- tra.agec_males["qup",]

        icaa <- paste(round(binf*100, digit_abs_agec), round(bsup*100, digit_abs_agec), sep='-')
        stat <- paste(round(rmed*100, digit_abs_agec)," (",icaa,")",sep="")

        tab.absolute.agec <- rbind(tab.absolute.agec, c(label_cum_males, stat))

        TAB.absolute.agec <- rbind(TAB.absolute.agec, data.frame("outcome" = label_cum_males,
                                                                 "age_class" = age_class,
                                                                 "median_abs_risk" = rmed,
                                                                 "lower_abs_risk" = binf,
                                                                 "upper_abs_risk" = bsup))

        # relative

        rrmed <- (tra.all_males["median",]/Ftpop_prod_males)[as.character(ages)]
        binf <- (tra.all_males['qlo',]/Ftpop_prod_males)[as.character(ages)]
        bsup <- (tra.all_males['qup',]/Ftpop_prod_males)[as.character(ages)]

        icaa <- paste(round(binf, digit_rr), round(bsup, digit_rr), sep='-')
        stat <- paste(round(rrmed, digit_rr)," (",icaa,")",sep="")

        tab.relative <- rbind(tab.relative,c(label_cum_males, stat ))
        TAB.relative <- rbind(TAB.relative, data.frame("outcome" = label_cum_males,
                                                       "age" = ages,
                                                       "median_relative_risk" = (tra.all_males["median",]/Ftpop_prod_males)[as.character(ages)],
                                                       "lower_relative_risk" = (tra.all_males['qlo',]/Ftpop_prod_males)[as.character(ages)],
                                                       "upper_relative_risk" = (tra.all_males['qup',]/Ftpop_prod_males)[as.character(ages)]))


        # relative risks by age class

        rmed <- tra.agec_males["median",]
        binf <- tra.agec_males["qlo",]
        bsup <- tra.agec_males["qup",]

        agecabs.pop <- Ftpop_prod_males[c(as.character(ages), "120")] - Ftpop_prod_males[c("0", as.character(ages))]

        rrmed <- rmed/agecabs.pop
        rrinf <- binf/agecabs.pop
        rrsup <- bsup/agecabs.pop

        icaa <- paste(round(rrinf, digit_rr), round(rrsup, digit_rr), sep='-')
        stat <- paste(round(rrmed, digit_rr)," (",icaa,")",sep="")

        tab.relative.agec <- rbind(tab.relative.agec, c(label_cum_males,stat ))
        TAB.relative.agec <- rbind(TAB.relative.agec, data.frame("outcome" = label_cum_males,
                                                                 "age_class" = age_class,
                                                                 "median_relative_risk" = rrmed,
                                                                 "lower_relative_risk" = rrinf,
                                                                 "upper_relative_risk" = rrsup))



      }


      if(traits_female){


        ## risk of at leat one trait is equal to 1-prod(Pr(no trait))
        Ftb_prod_females <- (1-one_minus_Ftb_prod_females)
        Ftpop_prod_females <- (1-one_minus_Ftpop_prod_females)

        names(Ftpop_prod_females) <- seqx.all

        tra.all_females  <- sapply(seqx.all, FUN = function(x){
          medix  <- median(Ftb_prod_females[x+1,])
          qlo    <- quantile(Ftb_prod_females[x+1,], (1-conf)/2)
          qup    <- quantile(Ftb_prod_females[x+1,], 1- ((1-conf)/2))
          return(c("median" = medix,
                   "qlo" = as.double(qlo),
                   "qup" = as.double(qup)))
        })

        colnames(tra.all_females) <- seqx.all

        ##agec
        Ftb.agec <- NULL
        aa <- c(0, ages, 120)

        for (ii in (2:length(aa))){

          ftb.agec <- Ftb_prod_females[aa[ii]+1,] - Ftb_prod_females[aa[ii-1]+1,]
          Ftb.agec <- rbind(Ftb.agec, ftb.agec)
        }

        tra.agec_females  <- sapply(1:nrow(Ftb.agec), FUN = function(x){
          medix  <- median(Ftb.agec[x,])
          qlo    <- quantile(Ftb.agec[x,], (1-conf)/2)
          qup    <- quantile(Ftb.agec[x,], 1- ((1-conf)/2))
          return(c("median" = medix,
                   "qlo" = as.double(qlo),
                   "qup" = as.double(qup)))
        })

        colnames(tra.agec_females) <- age_class

        ## fin agec


        label_cum_females <- paste(paste(traits, collapse = "|"), "females", sep = ' ')

        rmed <- tra.all_females["median",as.character(ages)]
        binf <- tra.all_females['qlo',as.character(ages)]
        bsup <- tra.all_females['qup',as.character(ages)]

        icaa <- paste(round(binf*100, digit_abs), round(bsup*100, digit_abs), sep='-')
        stat <- paste(round(rmed*100, digit_abs)," (",icaa,")",sep="")

        tab.absolute <- rbind(tab.absolute,c(label_cum_females, stat ))

        TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = label_cum_females,
                                                       "age" = ages,
                                                       "median_abs_risk" =  tra.all_females["median",as.character(ages)],
                                                       "lower_abs_risk" = tra.all_females['qlo',as.character(ages)],
                                                       "upper_abs_risk" = tra.all_females['qup',as.character(ages)],
                                                       "pop_risk" = Ftpop_prod_females[ages + 1]))

        # risks by age class
        rmed <- tra.agec_females["median",]
        binf <- tra.agec_females["qlo",]
        bsup <- tra.agec_females["qup",]

        icaa <- paste(round(binf*100, digit_abs_agec), round(bsup*100, digit_abs_agec), sep='-')
        stat <- paste(round(rmed*100, digit_abs_agec)," (",icaa,")",sep="")

        tab.absolute.agec <- rbind(tab.absolute.agec, c(label_cum_females, stat))

        TAB.absolute.agec <- rbind(TAB.absolute.agec, data.frame("outcome" = label_cum_females,
                                                                 "age_class" = age_class,
                                                                 "median_abs_risk" = rmed,
                                                                 "lower_abs_risk" = binf,
                                                                 "upper_abs_risk" = bsup))

        # relative

        rrmed <- (tra.all_females["median",]/Ftpop_prod_females)[as.character(ages)]
        binf <- (tra.all_females['qlo',]/Ftpop_prod_females)[as.character(ages)]
        bsup <- (tra.all_females['qup',]/Ftpop_prod_females)[as.character(ages)]

        icaa <- paste(round(binf, digit_rr), round(bsup, digit_rr), sep='-')
        stat <- paste(round(rrmed, digit_rr)," (",icaa,")",sep="")

        tab.relative <- rbind(tab.relative,c(label_cum_females, stat ))
        TAB.relative <- rbind(TAB.relative, data.frame("outcome" = label_cum_females,
                                                       "age" = ages,
                                                       "median_relative_risk" = (tra.all_females["median",]/Ftpop_prod_females)[as.character(ages)],
                                                       "lower_relative_risk" = (tra.all_females['qlo',]/Ftpop_prod_females)[as.character(ages)],
                                                       "upper_relative_risk" = (tra.all_females['qup',]/Ftpop_prod_females)[as.character(ages)]))

        # relative risks by age class
        rmed <- tra.agec_females["median",]
        binf <- tra.agec_females["qlo",]
        bsup <- tra.agec_females["qup",]

        agecabs.pop <- Ftpop_prod_females[c(as.character(ages), "120")] - Ftpop_prod_females[c("0", as.character(ages))]

        rrmed <- rmed/agecabs.pop
        rrinf <- binf/agecabs.pop
        rrsup <- bsup/agecabs.pop

        icaa <- paste(round(rrinf, digit_rr), round(rrsup, digit_rr), sep='-')
        stat <- paste(round(rrmed, digit_rr)," (",icaa,")",sep="")

        tab.relative.agec <- rbind(tab.relative.agec, c(label_cum_females,stat ))
        TAB.relative.agec <- rbind(TAB.relative.agec, data.frame("outcome" = label_cum_females,
                                                                 "age_class" = age_class,
                                                                 "median_relative_risk" = rrmed,
                                                                 "lower_relative_risk" = rrinf,
                                                                 "upper_relative_risk" = rrsup))



      }

  colnames(tab.absolute) <- colnames(tab.relative) <- c("strata", paste0("Age=",ages))
  colnames(tab.absolute.agec) <- colnames(tab.relative.agec) <- c("strata", paste0("Age=",age_class))


  return(list("ABSOLUTE_CUM_RISKS" = tab.absolute,
              "ABSOLUTE_RISKS_AGEC" = tab.absolute.agec,
              "RELATIVE_CUM_RISKS" = tab.relative,
              "RELATIVE_RISKS_AGEC" = tab.relative.agec,
              "data_ABSOLUTE_CUM_RISKS" = TAB.absolute,
              "data_ABSOLUTE_RISKS_AGEC" = TAB.absolute.agec,
              "data_RELATIVE_CUM_RISKS" = TAB.relative,
              "data_RELATIVE_RISKS_AGEC" = TAB.relative.agec))


    }
  }


}


