#'
#' @importFrom stats splinefun
#' @importFrom stats approx
#'
ftvproc <- function(x, args){

    #reconstruct params from x, with help of PARAMS.mask
    paramsmask <- args$PARAMS.mask
    params <- matrix(NA, ncol=ncol(paramsmask), nrow=nrow(paramsmask), dimnames =list(rownames(paramsmask),colnames(paramsmask)))
    for (k in 1:length(x)) params[paramsmask == k] <- x[k]

    ndis  <- length(args$Ft.pop)
    ng    <- args$G$ngen

    # compute Ft with help of PARAMS
    Ft <- array(0,c(121,2,ndis,ng))
    ll <- 1
    for (dis in 1:ndis){
       cc <- 1
       if(args$FIT.pars[[dis]]$penet.model == "Weibull"){
            for (i in 1:2){
                for(j in 1:ng){
                   if(j == 1){
                     if(sum(paramsmask[ll:(ll+2), cc:(cc+ng-1)]) == 0){
                        #risk is null (for all genotypes) for that gender
                        Ft[,i,dis,] <- 0.
                     }else{
                        #homozygote at all loci, risk is deduced from Pop incidence
                        if(args$approxFt.aa){
                           #assumes risk for aa genotype equal pop incidence
                           Ft[,i,dis,1] <- args$Ft.pop[[dis]][,i]
                        }else{
                            #not yet implemented
                        }
                     }
                   }else{
                       if(sum(paramsmask[ll:(ll+2), cc]) > 0){
                         #risk function is defined by weibull params
                         y <- params[ll:(ll+2), cc]
                         names(y) <- c("d", "H1", "H2")
                         Ft[,i,dis,j] <- Weibull.modified(dH1H2toakl(as.list(y)), t = 0:120)
                         #which.infpop <- which(Ft[,i,dis,j] < args$Ft.pop[[dis]][,i])
                         #if(any(which.infpop)) Ft[which.infpop,i,dis,j]  <- args$Ft.pop[[dis]][which.infpop,i]
                         #Ft[2:6,i,dis,j]  <- args$Ft.pop[[dis]][2:6,i]
                       }else{
                         #not estimable genotypic configuration !
                         Ft[,i,dis,j] = 0.
                       }
                   }
                   cc <- cc + 1
                }
            }
            ll <- ll + 3
        }else{
            if(args$FIT.pars[[dis]]$penet.model == "Cox"){
              for (i in 1:2){
                for(j in 1:ng){
                    if(j == 1){
                     if(sum(paramsmask[ll, cc:(cc+ng-1)]) == 0){
                        #risk is null (for all genotypes) for that gender
                        Ft[,i,dis,] <- 0.
                     }else{
                        #homozygote at all loci, risk is deduced from Pop incidence
                        if(args$approxFt.aa)  Ft[,i,dis,1] <- args$Ft.pop[[dis]][,i]
                        else{
                        #not yet implemented
                        }
                     }
                    }else{
                      if(paramsmask[ll, cc] > 0){
                         #risk function is defined by cox HR
                         y <- params[ll, cc]
                         Ft[,i,dis,j] <- Constrain01(y * args$Ft.pop[[dis]][,i])
                      }else{
                         #not estimable genotypic configuration !
                         Ft[,i,dis,j] = 0.
                      }
                    }
                    cc <- cc + 1
                }
              }
              ll <- ll + 1
            }else{
               #non parametric n-points
               for (i in 1:2){
                for(j in 1:ng){
                    if(j == 1){
                     if(sum(paramsmask[ll, cc:(cc+ng-1)]) == 0){
                        #risk is null (for all genotypes) for that gender
                        Ft[,i,dis,] <- 0.
                     }else{
                        #homozygote at all loci, risk is deduced from Pop incidence
                        if(args$approxFt.aa)  Ft[,i,dis,1] <- args$Ft.pop[[dis]][,i]
                        else{
                        #not yet implemented
                        }
                     }
                    }else{
                      if(paramsmask[ll, cc] > 0){
                         #risk function is defined by non parametric n-points

                         xx <- args$FIT.pars[[dis]]$agenodes
                         n.points <- length(xx)
                         xx <- c(xx, 120)

                         y <- Constrain01(cumsum(exp(params[ll:(ll + n.points +1 -1), cc])))
                         #y <- Constrain01(cumsum(params[ll:(ll + args$n.points +1 -1), cc]))

                         if(any(is.na(y))) stop("ERROR in the ftprocv function !")

                         if(args$approx.np == "spline"){
                             yy <- splinefun(x=c(0,xx),y=c(0,y), method="hyman")
                             #yyy <- c(approx(x=c(0,20), y=c(0,yy(20)), xout = 0:20)$y, yy(21:90), rep(yy(90),30))
                             yyy <- yy(0:120)
                             #splinecoef <- get("z", envir = environment(yy))
                             #plot(c(0,20,50,90,120),c(0,y,y[args$n.points]), ylim=c(0,1))
                             #lines(yyy)
                         }else{
                               if(args$approx.np == "linear"){
                                    yyy <- approx(x=c(0,xx), y=c(0,y), xout = 0:120)$y
                               }else{}
                         }
                         Ft[,i,dis,j] <- Constrain01(yyy)
                         #which.infpop <- which(Ft[,i,dis,j] < args$Ft.pop[[dis]][,i])
                         #if(any(which.infpop)) Ft[which.infpop,i,dis,j]  <- args$Ft.pop[[dis]][which.infpop,i]
                         #Ft[2:6,i,dis,j]  <- args$Ft.pop[[dis]][2:6,i]

                      }else{
                         #not estimable genotypic configuration !
                         Ft[,i,dis,j] = 0
                      }
                    }
                    cc <- cc + 1
                }
              }
              ll <- ll + n.points + 1
            }
       }
    }

  # at age =0 risk is always 0
  Ft[1,,,] = 0
  ftv  <- as.vector(as.matrix(Ft))
  return(ftv)

}
