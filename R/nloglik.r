
#' Title
#'
#' @param x parameters of the penetrance curve
#' @param write_lklbyfam write a csv file with LKL by family
#' @param ... internal args
#'
#' @return to be completed
#'
#' @examples
#' #to be completed
#'
#' @importFrom utils write.table

nloglik <- function(x, write_lklbyfam = FALSE, ...){

  args  <- list(...)

  # ETAPE 3 : calcul de(s) fonction(s) de risque(s) a partir des parametres de(s) modele(s) maladie(s)
  if('ftv' %in% names(args)){
    ftv <- args$ftv
  }else{
    ftv <- ftvproc(x, args)
  }


  # ETAPE 4 : calcul de la vraisemblance par l'algo d'Elston-Stewart

  loglikbyfam <- matrix(0,nrow=length(args$X),ncol=3)
  rownames(loglikbyfam) <- names(args$X)
  colnames(loglikbyfam) <- c("mLKL_Phe_Gen","mLKL_Phe_Gen_Index","mLKL_GRL")

  if(args$LIK.method == "GRL"){

    # chargement des objets et fonctions dans les coeurs de calcul du cluster
    #clusterExport(cl = cl, varlist = c("args"))

    for (f in names(args$X)){

      # # ETAPE 4.1 : Init of likelihood matrices
       L <- likprocFor(X = args$X[[f]], ftv = ftv)
      #
      # # ETAPE 4.2 : likelihood by Elston Stewart algorithm
       cid <- args$X[[f]]$ni #last individual, to begin likelihood at the bottom of the tree
      #
       logliknum   <- peelingFor(X = args$X[[f]], G = args$G, LIKv=L$liknumv,   counselee.id = cid)$loglik
       loglikdenom <- peelingFor(X = args$X[[f]], G = args$G, LIKv=L$likdenomv, counselee.id = cid)$loglik
       loglikgrl   <- logliknum - loglikdenom
       loglikbyfam[f,1:3] <- c(logliknum,loglikdenom,loglikgrl)
       if(any(!is.finite(c(logliknum,loglikdenom)))){
         loglikgrl <- 0
         loglikbyfam[f,1:3] <-  0
         warning("family :",f,"excluded for unknown reason \n")
       }

    }

    # loglikbyfam <- clusterApply(cl = cl,
    #                             x = args$X,
    #                             fun = function(x){
    #                               #  STEP 4.1 : Init of likelihood matrices
    #                               L <- likprocFor(X = x, ftv = ftv)
    #
    #                               #  STEP 4.2 : likelihood computing by Elston Stewart algorithm
    #                               cid <- x$ni #last individual, to begin likelihood at the bottom of the tree
    #
    #                               logliknum   <- peelingFor(X = x, G = args$G, LIKv=L$liknumv,   counselee.id = cid)$loglik
    #                               loglikdenom <- peelingFor(X = x, G = args$G, LIKv=L$likdenomv, counselee.id = cid)$loglik
    #                               loglikgrl   <- logliknum - loglikdenom
    #                               loglikbyfam[f,1:3] <- c(logliknum,loglikdenom,loglikgrl)
    #                               if(any(!is.finite(c(logliknum,loglikdenom)))){
    #                                 loglikgrl <- 0
    #                                 #loglikbyfam[f,1:3] <-  0
    #                                 warning("family :",f,"excluded for unknown reason \n")
    #                               }
    #                             })

    if(write_lklbyfam){
      loglikbyfam$fid <- rownames(loglikbyfam)
      write.table(loglikbyfam[,c(4,1,2,3)], "loglikbyfam.csv", sep = ";", row.names = FALSE)
      cat('\n',"The GRL log-likelihood per family (-LKL) is written in 'loglikbyfam.csv' file at your current workspace.\n")
    }


  }

  return(-sum(loglikbyfam[,3]))

}
