
#' Title
#'
#' @param x parameters of the penetrance curve
#' @param write_lklbyfam write a csv file with LKL by family
#' @param cl cluster object for parallel computing
#' @param ... internal args
#'
#' @return to be completed
#'
#' @examples
#' #to be completed
#'
#' @importFrom utils write.table

nloglik <- function(x, write_lklbyfam = FALSE, cl = NULL, ...){

  args  <- list(...)

  # ETAPE 3 : calcul de(s) fonction(s) de risque(s) a partir des parametres de(s) modele(s) maladie(s)
  if('ftv' %in% names(args)){
    ftv <- args$ftv
  }else{
    ftv <- ftvproc(x, args)
  }


  # ETAPE 4 : calcul de la vraisemblance par l'algo d'Elston-Stewart

  loglikbyfam <- matrix(0, nrow=length(args$X), ncol=3)
  rownames(loglikbyfam) <- names(args$X)
  colnames(loglikbyfam) <- c("mLKL_Phe_Gen","mLKL_Phe_Gen_Index","mLKL_GRL")

  if(args$LIK.method == "GRL"){

    if(!is.null(cl)){

      loglikbyfam <- t(parSapply(cl = cl,
                              X = args$X,
                              FUN = loglikf,
                              G = args$G,
                              ftv = ftv)
                       )

    }else{

      loglikbyfam <- t(sapply(X = args$X,
                              FUN = loglikf,
                              G = args$G,
                              ftv = ftv)
                       )

    }

    if(write_lklbyfam){
      write.table(loglikbyfam, "loglikbyfam.csv", sep = ";")
      cat('\n',"The GRL log-likelihood per family (-LKL) is written in 'loglikbyfam.csv' file at your current workspace.\n")
    }


  }

  return(-sum(loglikbyfam[,"mLKL_GRL"]))

}



loglikf <- function(f, G, ftv){

  ## STEP 4.1 : Init of likelihood matrices
  L <- likprocFor(X = f, ftv = ftv)

  ## STEP 4.2 : likelihood by Elston Stewart algorithm
  cid <- f$ni #last individual, to begin likelihood at the bottom of the tree
  logliknum   <- peelingFor(X = f, G = G, LIKv=L$liknumv,   counselee.id = cid)$loglik
  loglikdenom <- peelingFor(X = f, G = G, LIKv=L$likdenomv, counselee.id = cid)$loglik

  if(any(!is.finite(c(logliknum,loglikdenom)))){
    loglikgrl <- 0
    loglik_all <-  c(0,0,0)
    warning("family excluded for unknown reason \n")
  }else{
    loglikgrl   <- logliknum - loglikdenom
    loglik_all  <- c("mLKL_Phe_Gen" = logliknum,
                     "mLKL_Phe_Gen_Index" = loglikdenom,
                     "mLKL_GRL" = loglikgrl)
  }

  return(loglik_all)

}
