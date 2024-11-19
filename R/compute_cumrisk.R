#' compute cumulative risk over multiple traits from a bootstrapped generisk call object
#'
#' @param x an object returned by the generisk function
#' @param conf show bootstrapped 95% confidence intervals
#' @param ages at which ages ?
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
                          traits = names(x$params$FIT.pars)){

  is.bootstrap <- ('boot' %in% names(x))
  qn <- qnorm(1-(1-conf)/2)

  TAB.absolute <- tab.absolute <- NULL
  TAB.relative <- tab.relative <- NULL

  mask  <- x$fit$paramsmask

  ll  <- nrow(mask)
  dis <- length(x$params$FIT.pars)
  ng  <- ncol(mask)/2
  sexe <- rep(c(1,2),each=ng)
  geno <- rep(1:ng, times = 2)


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
         Ftb  <- sapply(x$boot$ft.boot, FUN=function(boot) boot[,sexe[cc],dis,geno[cc]])

         if(disname %in% traits){

           # compute Pr(at least one trait) = prod_over_traits(1-pr(risk_trait_j))

           if(length(unique(sexe[mask_matched]))==2){

             # 1 parameter for both sex
             # Ftb is identical for males and females for this trait

             traits_male <- TRUE
             traits_female <- TRUE
             Ftpop_males <- x$fit$ft[,1,dis,1][1:121]
             Ftpop_females <- x$fit$ft[,2,dis,1][1:121]

             one_minus_Ftb_prod_males <- one_minus_Ftb_prod_males*(1-Ftb)
             one_minus_Ftpop_prod_males <- one_minus_Ftpop_prod_males * (1-Ftpop_males)
             one_minus_Ftb_prod_females <- one_minus_Ftb_prod_females*(1-Ftb)
             one_minus_Ftpop_prod_females <- one_minus_Ftpop_prod_females * (1-Ftpop_females)

           }else{
             # sex specific risks..
             if(sexe[cc] == 1){
               #males
               traits_male <- TRUE
               Ftpop_males <- x$fit$ft[,1,dis,1][1:121]
               one_minus_Ftb_prod_males <- one_minus_Ftb_prod_males*(1-Ftb)
               one_minus_Ftpop_prod_males <- one_minus_Ftpop_prod_males * (1-Ftpop_males)
             }else{
               #females
               traits_female <- TRUE
               Ftpop_females <- x$fit$ft[,2,dis,1][1:121]
               one_minus_Ftb_prod_females <- one_minus_Ftb_prod_females*(1-Ftb)
               one_minus_Ftpop_prod_females <- one_minus_Ftpop_prod_females * (1-Ftpop_females)
             }
           }


         }else{

           tra.all  <- sapply(0:120, FUN = function(x){
                              medix  <- median(Ftb[x+1,])
                              qlo    <- quantile(Ftb[x+1,], (1-conf)/2)
                              qup    <- quantile(Ftb[x+1,], 1- ((1-conf)/2))
                              return(c("median" = medix,
                                       "qlo" = as.double(qlo),
                                       "qup" = as.double(qup)))
                            })

           colnames(tra.all) <- 0:120


          icaa <- paste(round(tra.all['qlo',as.character(ages)]*100,1), round(tra.all['qup',as.character(ages)]*100,1), sep='-')
          stat <- paste(round(tra.all["median",as.character(ages)]*100,1)," (",icaa,")",sep="")
          tab.absolute <- rbind(tab.absolute,c(paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),stat ))


          if(length(unique(sexe[mask_matched]))==2){
            # 1 parameter for both sex
            Ftpop.all <- rowMeans(x$fit$ft[,,dis,1])
          }else{
            Ftpop.all <- x$fit$ft[,sexe[cc],dis,1][1:121]
          }

          TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),
                                                         "age" = ages,
                                                         "median_abs_risk" =  tra.all["median",as.character(ages)],
                                                         "lower_abs_risk" = tra.all['qlo',as.character(ages)],
                                                         "upper_abs_risk" = tra.all['qup',as.character(ages)],
                                                         "pop_risk" = Ftpop.all[ages + 1]))


      # relative

          icaa <- paste(round((tra.all['qlo',]/Ftpop.all)[as.character(ages)],1), round((tra.all['qup',]/Ftpop.all)[as.character(ages)],1), sep='-')
          stat <- paste(round((tra.all["median",]/Ftpop.all)[as.character(ages)],1)," (",icaa,")",sep="")
          tab.relative <- rbind(tab.relative,c(paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),stat ))
          TAB.relative <- rbind(TAB.relative, data.frame("outcome" = paste(labcurve,paste('(',penetmod,')',sep=""), sep=' '),
                                                         "age" = ages,
                                                         "median_relative_risk" = (tra.all["median",]/Ftpop.all)[as.character(ages)],
                                                         "lower_relative_risk" = (tra.all['qlo',]/Ftpop.all)[as.character(ages)],
                                                         "upper_relative_risk" = (tra.all['qup',]/Ftpop.all)[as.character(ages)]))

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

        tra.all_males  <- sapply(0:120, FUN = function(x){
          medix  <- median(Ftb_prod_males[x+1,])
          qlo    <- quantile(Ftb_prod_males[x+1,], (1-conf)/2)
          qup    <- quantile(Ftb_prod_males[x+1,], 1- ((1-conf)/2))
          return(c("median" = medix,
                   "qlo" = as.double(qlo),
                   "qup" = as.double(qup)))
        })

        colnames(tra.all_males) <- 0:120
        label_cum_males <- paste(paste(traits, collapse = "|"), "males", sep = ' ')

        icaa <- paste(round(tra.all_males['qlo',as.character(ages)]*100,1), round(tra.all_males['qup',as.character(ages)]*100,1), sep='-')
        stat <- paste(round(tra.all_males["median",as.character(ages)]*100,1)," (",icaa,")",sep="")

        tab.absolute <- rbind(tab.absolute,c(label_cum_males, stat ))

        TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = label_cum_males,
                                                       "age" = ages,
                                                       "median_abs_risk" =  tra.all_males["median",as.character(ages)],
                                                       "lower_abs_risk" = tra.all_males['qlo',as.character(ages)],
                                                       "upper_abs_risk" = tra.all_males['qup',as.character(ages)],
                                                       "pop_risk" = Ftpop_prod_males[ages + 1]))

        # relative

        icaa <- paste(round((tra.all_males['qlo',]/Ftpop_prod_males)[as.character(ages)],1), round((tra.all_males['qup',]/Ftpop_prod_males)[as.character(ages)],1), sep='-')
        stat <- paste(round((tra.all_males["median",]/Ftpop_prod_males)[as.character(ages)],1)," (",icaa,")",sep="")
        tab.relative <- rbind(tab.relative,c(label_cum_males, stat ))
        TAB.relative <- rbind(TAB.relative, data.frame("outcome" = label_cum_males,
                                                       "age" = ages,
                                                       "median_relative_risk" = (tra.all_males["median",]/Ftpop_prod_males)[as.character(ages)],
                                                       "lower_relative_risk" = (tra.all_males['qlo',]/Ftpop_prod_males)[as.character(ages)],
                                                       "upper_relative_risk" = (tra.all_males['qup',]/Ftpop_prod_males)[as.character(ages)]))


      }


      if(traits_female){

        Ftb_prod_females <- (1-one_minus_Ftb_prod_females)
        Ftpop_prod_females <- (1-one_minus_Ftpop_prod_females)

        tra.all_females  <- sapply(0:120, FUN = function(x){
          medix  <- median(Ftb_prod_females[x+1,])
          qlo    <- quantile(Ftb_prod_females[x+1,], (1-conf)/2)
          qup    <- quantile(Ftb_prod_females[x+1,], 1- ((1-conf)/2))
          return(c("median" = medix,
                   "qlo" = as.double(qlo),
                   "qup" = as.double(qup)))
        })

        colnames(tra.all_females) <- 0:120

        label_cum_females <- paste(paste(traits, collapse = "|"), "females", sep = ' ')

        icaa <- paste(round(tra.all_females['qlo',as.character(ages)]*100,1), round(tra.all_females['qup',as.character(ages)]*100,1), sep='-')
        stat <- paste(round(tra.all_females["median",as.character(ages)]*100,1)," (",icaa,")",sep="")
        tab.absolute <- rbind(tab.absolute,c(label_cum_females, stat ))

        TAB.absolute <- rbind(TAB.absolute, data.frame("outcome" = label_cum_females,
                                                       "age" = ages,
                                                       "median_abs_risk" =  tra.all_females["median",as.character(ages)],
                                                       "lower_abs_risk" = tra.all_females['qlo',as.character(ages)],
                                                       "upper_abs_risk" = tra.all_females['qup',as.character(ages)],
                                                       "pop_risk" = Ftpop_prod_females[ages+1]))

        #relative

        icaa <- paste(round((tra.all_females['qlo',]/Ftpop_prod_females)[as.character(ages)],1), round((tra.all_females['qup',]/Ftpop_prod_females)[as.character(ages)],1), sep='-')
        stat <- paste(round((tra.all_females["median",]/Ftpop_prod_females)[as.character(ages)],1)," (",icaa,")",sep="")
        tab.relative <- rbind(tab.relative,c(label_cum_females, stat ))
        TAB.relative <- rbind(TAB.relative, data.frame("outcome" = label_cum_females,
                                                       "age" = ages,
                                                       "median_relative_risk" = (tra.all_females["median",]/Ftpop_prod_females)[as.character(ages)],
                                                       "lower_relative_risk" = (tra.all_females['qlo',]/Ftpop_prod_females)[as.character(ages)],
                                                       "upper_relative_risk" = (tra.all_females['qup',]/Ftpop_prod_females)[as.character(ages)]))


      }

  colnames(tab.absolute) <- colnames(tab.relative) <- c("strata", paste0("Age=",ages))

  return(list("ABSOLUTE_CUM_RISKS" = tab.absolute,
              "RELATIVE_CUM_RISKS" = tab.relative,
              "TAB_data_abs" = TAB.absolute,
              "TAB_data_relative" = TAB.relative))

    }
  }


}


