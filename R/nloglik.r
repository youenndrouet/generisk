
#' Title
#'
#' @param x parameters of the penetrance curve
#' @param return_lklbyfam return a table with LKL by family
#' @param ... internal args
#'
#' @return to be completed
#'
#' @examples
#' #to be completed
#'
#' @importFrom utils write.table

nloglik <- function(x, return_lklbyfam = FALSE, ...){

  args  <- list(...)

  # STEP 3 : calculation of risk function(s) from disease model(s) parameters
  if('ftv' %in% names(args)){
    ftv <- args$ftv
  }else{
    ftv <- ftvproc(x, args)
  }


  # STEP 4 : likelihood calculation using the Elston-Stewart algo

  loglikbyfam <- data.frame(matrix(nrow = length(args$X), ncol = 3))
  rownames(loglikbyfam) <- names(args$X)
  colnames(loglikbyfam) <- c("mLKL_numerator", "mLKL_denominator", "mLKL_optim")

  loglikbyfam <- t(sapply(X = args$X,
                          FUN = loglikf,
                          G = args$G,
                          ftv = ftv,
                          LIK.method = args$LIK.method)
  )

  if(return_lklbyfam){
    return(loglikbyfam)
  }else{
    return(sum(loglikbyfam[,"mLKL_optim"]))
  }


}



loglikf <- function(f, G, ftv, LIK.method){

  ## STEP 4.1 : Init of likelihood matrices
  L <- likprocFor(X = f, ftv = ftv)

  ## STEP 4.2 : likelihood by Elston Stewart algorithm
  cid <- f$ni #last individual, to begin likelihood at the bottom of the tree
  logliknum   <- peelingFor(X = f, G = G, LIKv=L$liknumv,   counselee.id = cid)$loglik

  if(LIK.method == "GRL"){
     loglikdenom <- peelingFor(X = f, G = G, LIKv=L$likdenomv, counselee.id = cid)$loglik
  }else{
     if(LIK.method %in% c("PEL", "PL")){
       loglikdenom <- 0
     }
  }


  if(any(!is.finite(c(logliknum, loglikdenom)))){
    loglikoptim <- 0
    loglik_all <-  c(0,0,0)
    warning("family excluded for unknown reason \n")
  }else{
    loglikoptim   <- logliknum - loglikdenom
    loglik_all  <- c("mLKL_numerator" = -logliknum,
                     "mLKL_denominator" = -loglikdenom,
                     "mLKL_optim" = -loglikoptim)
                     #"LKL_method" = LIK.method)
  }

  return(loglik_all)

}
