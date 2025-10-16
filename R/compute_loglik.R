#' compute the log-likelihood from a bootstrapped model
#'
#' @param generisk_obj an object returned by the generisk function
#'
#' @return to be completed
#' @export
#'
#' @examples
#' #to be completed
compute_loglik_boot <- function(generisk_obj){

  X     <- generisk_obj$X
  G     <- generisk_obj$G
  LIK.method   <- generisk_obj$params$LIK.method
  ndis <- length(generisk_obj$par$FIT.pars)
  ng  <- G$ngen

  ## calcul de ftv

  seqx.all <- 0:120

  Ft <- array(0,c(121,2,ndis,ng))

  for (dis in 1:ndis){
    for (i in 1:2){
      for(j in 1:ng){
          Ftb    <- sapply(generisk_obj$boot$ft.boot, FUN=function(x) x[,i,dis,j])
          med.all  <- sapply(seqx.all, FUN = function(x){return(median(Ftb[x+1,]))})
          Ft[,i,dis,j] <- med.all
      }
    }
  }
  # at age =0 risk is always 0
  Ft[1,,,] = 0
  ftv  <- as.vector(as.matrix(Ft))

  out <- nloglik(x = NULL,
          return_lklbyfam = TRUE,
          'ftv' = ftv,
          'G' = G,
          'X' = X,
          'LIK.method' = LIK.method)

  return(out)

}
