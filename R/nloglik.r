
nloglik <- function(x, ...){

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

  if(args$LIK.method == "GRL"){

    # chargement des objets et fonctions dans les coeurs de calcul du cluster
    # clusterExport(cl = cl, varlist = c("args"))

    for (f in names(args$X)){

      # ETAPE 4.1 : Init of likelihood matrices
      L <- likprocFor(X = args$X[[f]], ftv = ftv)
      #array(L$liknumv, dim=c(L$ni,L$ngen))

      # ETAPE 4.2 : likelihood by Elston Stewart algorithm
      cid <- args$X[[f]]$ni #last individual, to begin likelihood at the bottom of the tree

      logliknum   <- peelingFor(X = args$X[[f]], G = args$G, LIKv=L$liknumv,   counselee.id = cid)$loglik
      loglikdenom <- peelingFor(X = args$X[[f]], G = args$G, LIKv=L$likdenomv, counselee.id = cid)$loglik
      loglikgrl   <- logliknum - loglikdenom
      loglikbyfam[f,1:3] <- c(logliknum,loglikdenom,loglikgrl)
      if(any(!is.finite(c(logliknum,loglikdenom)))){
        loglikgrl <- 0
        #loglikbyfam[f,1:3] <-  0
        warning("family :",f,"excluded for unknown reason \n")
      }

     }

    #write.csv2(loglikbyfam, "loglikbyfam.csv")

  }

  return(-sum(loglikbyfam[,3]))

}
