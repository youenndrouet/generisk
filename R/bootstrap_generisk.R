#' bootstrap a generisk model
#'
#' @param B number of bootstrap
#' @param ncores number of workers multi-core computing
#' @param generisk_obj object from generisk
#' @param start init values for parameters, if "fit" starts from the initial fit, otherwise starts from population incidence
#'
#' @import parallel
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats nlminb
#'
#' @return to be completed
#'
#' @export
#'
#' @examples
#' #to be completed
bootstrap_generisk <- function(generisk_obj, B, ncores, start = "fit"){

  if(start == "fit"){
    pars.init    <- generisk_obj$fit$par
  }else{
    pars.init    <- generisk_obj$params$pars.init
  }

  Ft.pop       <- generisk_obj$params$Ft.pop
  LIK.method   <- generisk_obj$params$LIK.method
  FIT.pars     <- generisk_obj$params$FIT.pars
  approxFt.aa  <- generisk_obj$params$approxFt.aa
  approx.np    <- generisk_obj$params$approx.np

  penetmodels  <- generisk_obj$params$penetmodels
  pars.lower   <- generisk_obj$params$pars.lower
  pars.upper   <- generisk_obj$params$pars.upper
  rel.tol.boot <- generisk_obj$params$rel.tol.boot
  fA           <- generisk_obj$params$fA

  PARAMS.mask  <- generisk_obj$fit$paramsmask
  X     <- generisk_obj$X
  G     <- generisk_obj$G

  ndis  <- length(Ft.pop)
  nloci <- length(fA)
  ng    <- 3^nloci

  cat("\n  -> Bootstrap analysis  \n")
  pb <- txtProgressBar(char = "*", style =3)
  Sys.sleep(0.1)
  RESboot <- NULL
  ftaboot <- NULL

  if(ncores > 1){

    available_cores <- detectCores()

    if(available_cores >= ncores){
      cat("Parallel computing using", ncores, "/", available_cores, "available cores. \n")

      # set up each worker
      cl <- makeCluster(ncores)
      clusterEvalQ(cl, {
        library(generisk)
        NULL
      })

    }else{
      cl <- NULL
      cat("Available cores < ", ncores, ", please consider decreasing ncores to perform parallel computing. \n")
    }

  }else{
    cl <- NULL
  }


  for (b in seq(B)){

    # directly use the pre-computed X object
    X.boot <- bootfam(X)

    fit.boot <- nlminb(start = pars.init,
                       objective = function(x) nloglik(x,
                                                       return_lklbyfam = FALSE,
                                                       cl = cl,
                                                       'G' = G,
                                                       'X' = X.boot,
                                                       'Ft.pop' = Ft.pop,
                                                       'LIK.method' = LIK.method,
                                                       'FIT.pars' = FIT.pars,
                                                       'approxFt.aa' = approxFt.aa,
                                                       'approx.np' = approx.np,
                                                       'PARAMS.mask' = PARAMS.mask,
                                                       'penetmodels' = penetmodels),
                       lower = pars.lower,
                       upper = pars.upper,
                       control = list(trace=10, rel.tol = rel.tol.boot))

    ftab <- array(ftvproc(x = fit.boot$par,
                          args = list('Ft.pop' = Ft.pop,
                                      'FIT.pars' = FIT.pars,
                                      'approxFt.aa' = approxFt.aa,
                                      'approx.np' = approx.np,
                                      'PARAMS.mask' = PARAMS.mask,
                                      'G' =G)),
                  dim = c(121, 2, ndis, ng))

    ftaboot <- c(ftaboot, list(ftab))
    RESboot <- rbind(RESboot, c(convergence = fit.boot$convergence,
                                #message = fit.boot$message,
                                #iterations = fit.boot$iterations,
                                fit.boot$par,
                                nloglik = fit.boot$objective))
    setTxtProgressBar(pb, b/B)
  }
  close(pb)
  cat("\n  job done !\n")

  if(!is.null(cl)){stopCluster(cl)}

  is.bootstrap <- ('boot' %in% names(generisk_obj))
  out <- generisk_obj

  if(!is.bootstrap){
    #first bootstrap
    out$boot <- list("fit.boot" = RESboot,
                     "ft.boot" = ftaboot)
  }else{
    previous_fit.boot <- generisk_obj$boot$fit.boot
    new_fit.boot <- rbind(previous_fit.boot, RESboot)

    previous_ft.boot <- generisk_obj$boot$ft.boot
    new_ft.boot <- c(previous_ft.boot, ftaboot)

    out$boot <- list("fit.boot" = new_fit.boot,
                     "ft.boot" = new_ft.boot)

  }

  return(out)

}
