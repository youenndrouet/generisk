#' estimate age-risk (penetrance) functions from family data
#'
#' generisk() is the main function of the generisk package.
#' It performs penetrance estimation using maximum likelihood with the GRL method.
#' It relies on the Elston-Stewart algorithm to factorizing family data.
#' Different options are implemented for parametrization.
#'
#' @param Ft.pop a list of 121*2 matrices of diseases incidences (length = ndis)
#' @param LIK.method if GRL : first element of fA must be the allele frequency for the gene under selection
#' @param DATA a N*(6 + nloci + 2*ndis) matrix
#' @param FIT.pars a list of diseases parameters (length = ndis) (see details)
#' @param fA a vector of pop. allele frequencies of tested loci (order must corresponds to that of DATA columns)
#' @param B option to calculate 95% Bootstrap estimates # number of bootstrap resample (B=0 means no bootstrap)
#' @param rel.tol precision for the likelihood
#' @param approxFt.aa option to calculate Ft.aa
#' @param init.weibull initial values for Weibull model
#' @param approx.np linear or spline
#' @param multi.pheno strategy in case of multiple diseases (can be "first" or "all")
#' @param ncores number of cores used for parallel computing
#' @param imput_missing_age_last_news logical: imputation of missing ages at last news ?
#'
#' @import parallel
#' @importFrom stats runif
#' @importFrom stats nlminb
#' @importFrom stats na.omit
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom utils packageVersion
#' @importFrom dplyr bind_rows
#'
#' @return to be completed
#'
#' @examples
#' #to be completed
#'
#' @export
#'
generisk <- function(

  Ft.pop,
  LIK.method = "GRL",
  DATA,
  FIT.pars = list(                  # list of diseases (phenotypes)
    "cancer" = list(                       # each phenotype is described by following parameters
      penet.model   = "np",                # can be "Weibull" and "cox"
      agenodes = c(30,50,70),
      implic.loci   = c(MMR = TRUE),       # named TRUE/FALSE vector of implicated loci (gene with deletere mutation)
      inheritance   = c(MMR = "dominant"), # named vector of inheritance at each implicated locus can be dominant/recessive/free
      gender.effect = c(MMR = TRUE)        # named TRUE/FALSE vector of gender effect at each locus
    )
  ),
  fA = c(MMR = 0.001),
  B = 0,
  rel.tol = 1e-6,
  approxFt.aa = TRUE,
  approx.np = "spline",
  multi.pheno = "first",
  init.weibull = NULL,
  ncores = 1,                   # number of cores (default = 1)
  imput_missing_age_last_news = FALSE
){

  minrisk = 1e-16   #default minimum absolute risk is almost 0 (prevent singularities)
  maxrisk = 0.99999 #default maximum absolute risk is almost 1 (prevent singularities)

  ndis  <- length(Ft.pop)
  nloci <- length(fA)
  ng <- 3^nloci

  # check and recode DATA in DATA_generisk
  DATA_generisk <- DATA

  cat(' ________________________________________\n')
  cat('                                         \n')
  cat('  The generisk R program \n')
  cat('  version:', paste(packageVersion("generisk"), collapse = "."), '\n')
  cat('                                        \n')

  cat('  Number of families: ', length(unique(DATA_generisk[,1])), '\n')

  ## we should check here that diseases are coded 0/1
  ## and accept no NA

  cc_event <- seq((6+nloci+2)-1, (6+nloci+2*ndis)-1, 2) # columns coding for diseases events 0/1
  cc_ages  <- seq((6+nloci+2), (6+nloci+2*ndis), 2)  # columns coding for diseases ages at diagnosis or last news

  ## individuals are either:
  # unaffected (wi_unaff) of affected in at least one disease (wi_aff)
  if(ndis == 1){
    wi_unaff <- which(DATA_generisk[,cc_event] == 0)
    wi_aff <- which(DATA_generisk[,cc_event] != 0)
    minages_unaff <- DATA[wi_unaff, cc_ages]

  }else{
    wi_unaff <- which(rowSums(DATA_generisk[,cc_event]) == 0)
    wi_aff <- which(rowSums(DATA_generisk[,cc_event]) != 0)
    minages_unaff <- apply(DATA[wi_unaff, cc_ages], 1, pmin_na)

  }

  ## unaffected -> set min age at last news at all diseases
  ## note : individual with all NA at least news stay NA (see next imputation)

  for (cc in cc_ages){DATA_generisk[wi_unaff, cc] <- minages_unaff}

  ## affected (one or multiple diseases)

  if(ndis == 1){

      ## one disease
      ages <- DATA[wi_aff, cc_ages]
      ages[is.na(ages)] <- 0 ## missing age at diagnosis
      DATA_generisk[wi_aff, cc_ages] <- ages

  }else{

      # multiple diseases

      if(multi.pheno == "first"){

        ## since we restrict analyses to the "competitive risk" framework,
        ## only the first disease is kept for individuals with multiple diseases

        temp <- apply(as.matrix(DATA[wi_aff, min(cc_event):max(cc_ages)]), 1,
                      FUN = function(x){
                        status <- x[seq(1, length(x)-1, 2)]
                        ages <- x[seq(2, length(x), 2)]
                        ages[is.na(ages)] <- 0
                        ages[status == 0 | ages == 0] <- 999
                        wmin <- which.min(ages)
                        minage <- min(ages)
                        if(minage == 999){ ## missing ages at diagnosis
                          ages.out <- rep(0, length(ages))
                          status.out <- rep(0, length(status))
                        }else{ # at least one disease with age at diagnosis
                          ages.out <- rep(minage, length(ages))
                          status.out <- rep(0, length(status))
                          status.out[wmin] <- 1
                        }
                        out <-  as.vector(matrix(c(status.out, ages.out),
                                                 ncol = length(ages.out),
                                                 byrow = TRUE))
                        return(out)

                      }
        )

        DATA_generisk[wi_aff, min(cc_event):max(cc_ages)] <- t(temp)

    }else{

      if(multi.pheno == "all"){

        # disease with missing age at diagnosis are discarded
        DATA_generisk[wi_aff, cc_ages][is.na(DATA_generisk[wi_aff, cc_ages])] <- 0
        DATA_generisk[wi_aff, cc_event][DATA_generisk[wi_aff, cc_ages] == 0] <- 0

      }else{

        stop("unknown multi.pheno option")
      }

    }
  }

    for (cc in cc_ages){
      nb_aff <- sum(DATA[,cc-1])
      nb_aff_analyzed <- sum((DATA_generisk[,cc-1] == 1) & (DATA_generisk[,cc] != 0 ))
      cat('   -',names(DATA)[cc-1],": ", nb_aff_analyzed, "/",nb_aff, " affected individuals will be analyzed","\n")
    }

  ##save the DATA before imputation (for bootstrap)
  DATA_generisk_before_imputed_step <- DATA_generisk

  wi_to_impute <- which(is.na(DATA_generisk[,cc_ages[1]]))

  if(imput_missing_age_last_news){
    ## replace with median age at last news for unaff
    DATA_generisk[wi_to_impute, cc_ages] <- median(minages_unaff, na.rm = TRUE)
    cat('   -', "Unaffected individuals: missing age at last news imputed for", length(wi_to_impute),"/",length(wi_unaff), "individuals", "\n")
  }else{
    ## replace NA by O
    DATA_generisk[wi_to_impute, cc_ages] <- 0
    cat('   -', "Unaffected individuals:", length(wi_to_impute),"/",length(wi_unaff), " individuals with missing age at last news removed from analysis (age --> 0)", "\n")
  }

  cat('\n')
  # if(ncores > 1){
  #   # set up each worker
  #   cl <- makeCluster(ncores)
  #   clusterEvalQ(cl, {
  #     library(generisk)
  #     NULL
  #   })
  # }


  #ETAPE 0 : vérification de la cohérence des données familiales

  # checkit <- checkpedigrees(famid = DATA[,1],
  #                           id = paste(DATA[,1], DATA[,2], sep="/"),
  #                           father.id = paste(DATA[,1], DATA[,4], sep="/"),
  #                           mother.id = paste(DATA[,1], DATA[,5], sep="/"))

  #table(checkit$split) #should be all 1's
  #table(checkit$unrelated) #should be all 0's
  #table(checkit$join) #should be all 0's
  # a faire : traiter les probl?mes de boucle + ind. unrelated ...

  # ETAPE 1 : modèle génétique à n loci

  cat("  -> Program initialization. \n")
  G <- genetFor(allef=fA)

  #ETAPE 2.1 : matrice des paramètres

  #intervalles de definition des pénétrances
  AgeDef <- NULL
  for (cc in seq((6+nloci +2),(6+nloci + 2*ndis),2)){
    AgeDef <- c(AgeDef, list(c(min(DATA_generisk[(DATA_generisk[,cc-1] == 1) & (DATA_generisk[,cc] > 0),cc]),
                               max(DATA_generisk[(DATA_generisk[,cc-1] == 1) & (DATA_generisk[,cc] > 0),cc]))))
  }
  names(AgeDef) <- names(FIT.pars)

  PARAMS <- NULL
  parcolnames <- c(0,1,2)
  parrownames <- NULL
  penetmodels <- integer(ndis)

  for (dis in 1:ndis){
    disname <- names(FIT.pars)[dis]
    if(FIT.pars[[dis]]$penet.model == "Weibull"){
      penetmodels[dis] <- 1
    }else{
      if(FIT.pars[[dis]]$penet.model == "Cox"){
        penetmodels[dis] <- 2
      }else{
        penetmodels[dis] <- 3
      } #non parametric
    }
    for (loci in 1:nloci){
      if(FIT.pars[[dis]]$implic.loci[loci]){
        inher <- FIT.pars[[dis]]$inheritance[loci]
        gender <- FIT.pars[[dis]]$gender.effect[loci]
        if(gender){ #gender effect
          if(inher == "dominant"){
            param.m <- c("0m","1m","1m")
            param.f <- c("0f","1f","1f")
          }else{
            if(inher == "recessive"){
              param.m <- c("0m","0m","1m")
              param.f <- c("0f","0f","1f")
            }else{
              param.m <- c("0m","1m","2m")
              param.f <- c("0f","1f","2f")
            }
          }
        }else{ #no gender effect
          if(inher == "dominant"){
            param.m <- param.f <- c("0f/m","1f/m","1f/m")
          }else{
            if(inher == "recessive"){
              param.m <- param.f <- c("0f/m","0f/m","1f/m")
            }else{
              param.m <- param.f <- c("0f/m","1f/m","2f/m")
            }
          }
        }
      }else{ #loci not implicated
        param.m <- param.f <- c("0f/m","0f/m","0f/m")
      }
      if(loci==1){
        PARAM.f <- param.f
        PARAM.m <- param.m
        PARnames <- parcolnames
      }else{
        PARAM.f <- paste(rep(PARAM.f, each = 3), rep(param.f, times = 3^(loci-1)), sep ="")
        PARAM.m <- paste(rep(PARAM.m, each = 3), rep(param.m, times = 3^(loci-1)), sep ="")
        PARnames <- paste(rep(PARnames, each = 3), rep(parcolnames, times = 3^(loci-1)), sep ="")
      }
    }
    PARAM <- c(PARAM.m,PARAM.f)
    if(penetmodels[dis] == 1){
      PARAMS <- rbind(PARAMS,paste(dis,"d",PARAM,sep=":"),paste(dis,"H1",PARAM,sep=":"),paste(dis,"H2",PARAM,sep=":"))
      parrownames <- c(parrownames, paste(disname,c("d","H1","H2"), sep=":"))
    }else{
      if(penetmodels[dis] == 2){ #cox
        PARAMS <- rbind(PARAMS,paste(dis,"HR",PARAM,sep=":"))
        parrownames <- c(parrownames, paste(disname,"HR", sep=":"))
      }else{ #np

        agenodes <- FIT.pars[[dis]]$agenodes
        aa <- c(agenodes,120)
        PARAMS <- rbind(PARAMS, matrix(paste(dis,rep(aa, each = length(PARAM)),PARAM,sep=":"), ncol=length(PARAM), byrow=TRUE))
        parrownames <- c(parrownames, paste(disname, aa, sep=":"))
      }
    }
  }

  colnames(PARAMS) <- c(paste("m",PARnames, sep=":"), paste("f",PARnames, sep=":"))
  rownames(PARAMS) <- parrownames

  #remove parameters corresponding to disease with risk=0 for one sex
  ll <- 1
  for (dis in 1:ndis){
    if(penetmodels[dis] == 1){
      nl <- 3 #weibull
    }else{
      if(penetmodels[dis] == 2){
        nl <- 1 #cox
      }else{
        nl <- length(FIT.pars[[dis]]$agenodes) + 1
      }
    }
    if(sum(Ft.pop[[dis]][,1]) == 0){ #no risk for men
      PARAMS[ll:(ll+nl-1), 1:ng] <- ""
    }else{
      if(sum(Ft.pop[[dis]][,2]) == 0){#no risk for women
        PARAMS[ll:(ll+nl-1), (ng+1):(2*ng)] <- ""
      }
    }
    if(penetmodels[dis] == 1){
      ll <- ll + 3 #weibull
    }else{
      if(penetmodels[dis] == 2){
        ll <- ll + 1 #cox
      }else{
        ll <- ll + nl #np
      }
    }

  }

  #remove one genotype (wild type homozygote at all loci), because the associated risk is deduced (pop incidence)
  PARAMS[,c(1,ng+1)] <- ""

  #remove the genotypes never observed in data
  POSS.genotypes <- names(table(rowPaste(DATA_generisk[,7:(7+length(fA)-1)])))
  POSS.genotypes <- c(paste("m",POSS.genotypes, sep=":"), paste("f",POSS.genotypes, sep=":"))
  for (k in 1:ncol(PARAMS)) if(!(colnames(PARAMS)[k] %in% POSS.genotypes)) PARAMS[,k] <- ""

  PARAMS.mask <- matrix(0,
                        ncol = ncol(PARAMS),
                        nrow = nrow(PARAMS),
                        dimnames = list(rownames(PARAMS),colnames(PARAMS)))

  parameters <- unique(as.vector(PARAMS[PARAMS != ""]))

  for (k in 1:length(parameters)) PARAMS.mask[PARAMS == parameters[k]] <- k

  # ETAPE 2.2 : calculs preliminaires sur les donn?es familiales (pour ne pas refaire ces calculs dans la boucle d'optimisation et de bootstrap)

  cat("  -> Pre-calculations to speed-up likelihood algorithm. \n")
  fids  <- unique(DATA_generisk[,1])
  X <- lapply(fids, function(x) preprocFor(ped = DATA_generisk[DATA_generisk[,1]==x,-1],allef=fA, ndis=ndis))
  names(X) <- fids

  #init values for parameters
  pars.init <- NULL

  ll <- 1
  for (dis in 1:ndis){

    if(penetmodels[dis] == 1){

      nl <- 3 #weibull

      for (g in 2:ng){

        if(PARAMS.mask[ll,g] == 0 & PARAMS.mask[ll,g+ng] == 0){
        }else{

          if(PARAMS.mask[ll,g] == PARAMS.mask[ll,g+ng]){ #equal for male and female

            if(is.null(init.weibull)){
              fitmf <- nlminb(start = c(d=3/5, H1=2/3, H2=1/4),
                              objective = function(x) SCEweibull(x, y=rowMeans(Ft.pop[[dis]])),
                              lower = list(d=0.001, H1=0.001, H2=0.001),
                              upper = list(d=0.999, H1=0.999, H2=0.999))
              pars.init <- c(pars.init, fitmf$par)
            }else{

            }

          }else{

            if((PARAMS.mask[ll,g] != 0) & (PARAMS.mask[ll,g+ng] == 0)){ #only men
              fitm <- nlminb(start = c(d=3/5, H1=2/3, H2=1/4),
                             objective = function(x) SCEweibull(x, y=Ft.pop[[dis]][,"m"]),
                             lower = list(d=0.001, H1=0.001, H2=0.001),
                             upper = list(d=0.999, H1=0.999, H2=0.999))
              pars.init <- c(pars.init, fitm$par)
            }else{
              if((PARAMS.mask[ll,g] == 0) & (PARAMS.mask[ll,g+ng] != 0)){#only women
                fitf <- nlminb(start = c(d=3/5, H1=2/3, H2=1/4),
                               objective = function(x) SCEweibull(x, y=Ft.pop[[dis]][,"f"]),
                               lower = list(d=0.001, H1=0.001, H2=0.001),
                               upper = list(d=0.999, H1=0.999, H2=0.999))
                pars.init <- c(pars.init, fitf$par)
              }else{ #diff for male and female
                fitm <- nlminb(start = c(d=3/5, H1=2/3, H2=1/4),
                               objective = function(x) SCEweibull(x, y=Ft.pop[[dis]][,"m"]),
                               lower = list(d=0.001, H1=0.001, H2=0.001),
                               upper = list(d=0.999, H1=0.999, H2=0.999))
                fitf <- nlminb(start = c(d=3/5, H1=2/3, H2=1/4),
                               objective = function(x) SCEweibull(x, y=Ft.pop[[dis]][,"f"]),
                               lower = list(d=0.001, H1=0.001, H2=0.001),
                               upper = list(d=0.999, H1=0.999, H2=0.999))
                pars.init <- c(pars.init, fitm$par, fitf$par)
              }
            }
          }
        }
      }

    }else{
      if(penetmodels[dis] == 2){
        nl <- 1 #cox
        for (g in 2:ng){
          if(PARAMS.mask[ll,g] == 0 & PARAMS.mask[ll,g+ng] == 0){
          }else{
            if(PARAMS.mask[ll,g] == PARAMS.mask[ll,g+ng]){ #equal for male and female
              pars.init <- c(pars.init, 1)
            }else{
              if((PARAMS.mask[ll,g] != 0) & (PARAMS.mask[ll,g+ng] == 0)){ #only men
                pars.init <- c(pars.init, 1)
              }else{
                if((PARAMS.mask[ll,g] == 0) & (PARAMS.mask[ll,g+ng] != 0)){#only women
                  pars.init <- c(pars.init, 1)
                }else{ #diff for male and female
                  pars.init <- c(pars.init, 1, 1)
                }
              }
            }
          }
        }

      }else{  #np

        n.points <- length(FIT.pars[[dis]]$agenodes)
        nl <- n.points + 1
        aa <- c(FIT.pars[[dis]]$agenodes,120)

        for (g in 2:ng){
          if(PARAMS.mask[ll,g] == 0 & PARAMS.mask[ll,g+ng] == 0){
          }else{
            if(PARAMS.mask[ll,g] == PARAMS.mask[ll,g+ng]){ #equal for male and female

              a <- rowMeans(Ft.pop[[dis]])[aa[1:(n.points+1)] + 1]
              b <- rowMeans(Ft.pop[[dis]])[c(1,aa[1:(n.points)] + 1)]
              ## log-transformation to speed-up and improve convergence with nlminb
              ##1e-6 to deal with 0 values
              lamb <- log(a-b + 1e-6)
              pars.init <- c(pars.init, lamb)
            }else{
              if((PARAMS.mask[ll,g] != 0) & (PARAMS.mask[ll,g+ng] == 0)){ #only men
                a <- Ft.pop[[dis]][aa[1:(n.points+1)] + 1, "m"]
                b <- Ft.pop[[dis]][c(1,aa[1:(n.points)] + 1), "m"]
                lamb <- log(a-b + 1e-6)
                pars.init <- c(pars.init, lamb)
              }else{
                if((PARAMS.mask[ll,g] == 0) & (PARAMS.mask[ll,g+ng] != 0)){#only women
                  a <- Ft.pop[[dis]][aa[1:(n.points+1)] + 1, "f"]
                  b <- Ft.pop[[dis]][c(1,aa[1:(n.points)] + 1), "f"]
                  lamb <- log(a-b + 1e-6)
                  pars.init <- c(pars.init, lamb)
                }else{ #diff for male and female
                  a <- Ft.pop[[dis]][aa[1:(n.points+1)] + 1, "m"]
                  b <- Ft.pop[[dis]][c(1,aa[1:(n.points)] + 1), "m"]
                  lamb <- log(a-b + 1e-6)
                  pars.init <- c(pars.init, lamb)
                  a <- Ft.pop[[dis]][aa[1:(n.points+1)] + 1, "f"]
                  b <- Ft.pop[[dis]][c(1,aa[1:(n.points)] + 1), "f"]
                  lamb <- log(a-b + 1e-6)
                  pars.init <- c(pars.init, lamb)
                }
              }
            }
          }
        }
      }
    }
    ll <- ll + nl
  }

  #lower.weibull <- c(d=0.001, H1=0.001, H2=0.001)
  upper.weibull <- c(d=0.999, H1=0.999, H2=0.999)
  #lower.cox     <- 1/100
  upper.cox     <- 100
  #lower.np <- log(1e-6)
  upper.np <- log(1)

  pars.lower <- pars.upper <- numeric(length(parameters))

  ## by default, set lower = population incidence
  pars.lower <- pars.init

   for (a in as.character(1:100)){
  #    pars.lower[grep(a,parameters,fixed = TRUE)] <- lower.np
      pars.upper[grep(a,parameters,fixed = TRUE)] <- upper.np
   }

  for (a in c("d","H1","H2")){
    #pars.lower[grep(a,parameters,fixed = TRUE)] <- lower.weibull[a]
    pars.upper[grep(a,parameters,fixed = TRUE)] <- upper.weibull[a]

    if(!is.null(init.weibull)){
      pars.init[grep(a,parameters,fixed = TRUE)] <- init.weibull[a]
    }

  }

  #pars.lower[grep("HR",parameters,fixed = TRUE)] <- lower.cox
  pars.upper[grep("HR",parameters,fixed = TRUE)] <- upper.cox

#  nloglik(pars.init, 'G' =G, 'X' = X, 'Ft.pop' = Ft.pop, 'LIK.method' = LIK.method, 'FIT.pars' = FIT.pars, 'approxFt.aa' = approxFt.aa, 'approx.np' = approx.np, 'PARAMS.mask' = PARAMS.mask, 'penetmodels' = penetmodels)

  cat("  -> ML optimization by nlminb. \n")
  fit <- nlminb(start = pars.init,  # ETAPES 3 et 4
                objective = function(x) nloglik(x, 'G' =G, 'X' = X, 'Ft.pop' = Ft.pop, 'LIK.method' = LIK.method, 'FIT.pars' = FIT.pars, 'approxFt.aa' = approxFt.aa, 'approx.np' = approx.np, 'PARAMS.mask' = PARAMS.mask, 'penetmodels' = penetmodels),
                lower = pars.lower,
                upper = pars.upper,
                control = list(trace=2,rel.tol=rel.tol))

  fta <- array(ftvproc(x = fit$par, args = list('Ft.pop' = Ft.pop, 'FIT.pars' = FIT.pars, 'approxFt.aa' = approxFt.aa, 'approx.np' = approx.np, 'PARAMS.mask' = PARAMS.mask, 'G' =G)), dim=c(121,2,ndis,ng))

  nloglik(fit$par, write_lklbyfam = TRUE, 'G' =G, 'X' = X, 'Ft.pop' = Ft.pop, 'LIK.method' = LIK.method, 'FIT.pars' = FIT.pars, 'approxFt.aa' = approxFt.aa, 'approx.np' = approx.np, 'PARAMS.mask' = PARAMS.mask, 'penetmodels' = penetmodels)

  if(B > 0){

    cat("\n  -> Bootstrap analysis for Confidence Intervals \n")
    pb <- txtProgressBar(char = "*", style =3)
    Sys.sleep(0.1)
    RESboot <- NULL
    ftaboot <- NULL

    for (b in seq(B)){

      if(!imput_missing_age_last_news){

        # directly use the pre-computed X object
        X.boot <- bootfam(X)

      }else{

        # have to re-compute from DATA_generisk_before_imputed_step
        fids.boot <- sample(fids, replace = TRUE)
        DATA_generisk.boot <- lapply(fids.boot, FUN=function(f) DATA_generisk_before_imputed_step[DATA_generisk_before_imputed_step[,1] == f,][,-1])

        # numfam is changed to fix the case of families selected multiple times..
        DATA_generisk.boot <- bind_rows(DATA_generisk.boot, .id = "fid")

        # imputation of missing ages
        wi_to_impute.boot <- which(is.na(DATA_generisk.boot[,cc_ages[1]]))
        n_to_impute <- length(wi_to_impute.boot)
        ages_imputed <- sample(x = na.omit(minages_unaff), size = n_to_impute, replace = TRUE)

        for (cc in cc_ages){DATA_generisk.boot[wi_to_impute.boot, cc] <- ages_imputed}

        ## pré-calculs pour FORTRAN
        fids.boot  <- unique(DATA_generisk.boot$fid)
        X.boot <- lapply(fids.boot, function(f) preprocFor(ped = DATA_generisk.boot[DATA_generisk.boot$fid==f,-1],allef=fA, ndis=ndis))
        names(X.boot) <- fids.boot

      }


      fit.boot <- nlminb(start = fit$par, #init ? partir du max LIK pour les bootstrap
                         objective = function(x) nloglik(x,
                                                         'G' =G,
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
                         control = list(trace=0,rel.tol=rel.tol))
      ftab <- array(ftvproc(x = fit.boot$par, args = list('Ft.pop' = Ft.pop, 'FIT.pars' = FIT.pars, 'approxFt.aa' = approxFt.aa,'approx.np' = approx.np, 'PARAMS.mask' = PARAMS.mask, 'G' =G)), dim=c(121,2,ndis,ng))
      ftaboot <- c(ftaboot,list(ftab))
      RESboot <- rbind(RESboot,c(convergence = fit.boot$convergence, fit.boot$par, nloglik = fit.boot$objective))
      setTxtProgressBar(pb, b/B)
    }
    close(pb)
    cat("\n  job done !\n")
    return(list(fit=c(fit,paramsmask = list(PARAMS.mask), ft=list(fta)), boot = list(fit.boot = RESboot, ft.boot = ftaboot), params=list(Ft.pop = Ft.pop,'approx.np' = approx.np ,FIT.pars = FIT.pars, LIK.method=LIK.method, fA=fA, rel.tol=rel.tol, AgeDef = AgeDef), X=X, G=G))
  }else{
    cat("\n  job done !\n")
    return(list(fit=c(fit,paramsmask = list(PARAMS.mask), ft=list(fta)),params=list(Ft.pop = Ft.pop,'approx.np' = approx.np, FIT.pars = FIT.pars, LIK.method=LIK.method, fA=fA, rel.tol=rel.tol, AgeDef = AgeDef), X=X, G=G))
  }

  #  stopCluster(cl)
}
