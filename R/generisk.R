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
#' @param n.points number of points used for the non parametric penetrance function
#' @param approx.np linear or spline
#' @param ncores number of cores used for parallel computing
#' @param imput_ageNA logical: imputation of missing ages at last news ?
#'
#' @importFrom stats runif
#' @importFrom stats nlminb
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
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
  n.points = 5,
  approx.np = "spline",
  ncores = 4,
  imput_ageNA = FALSE
){

  minrisk = 1e-16   #default minimum absolute risk is almost 0 (prevent singularities)
  maxrisk = 0.99999 #default maximum absolute risk is almost 1 (prevent singularities)

  ndis  <- length(Ft.pop)
  nloci <- length(fA)
  ng <- 3^nloci

  if(imput_ageNA){
    ## imputation des âges manquants pour age_t1 (et age_t1b)
    n_toimpute <- sum(DATA$toimpute)
    DATA$age_T1[DATA$toimpute == 1] <- floor(runif(n = n_toimpute, min = 0, max = 100))
    DATA$age_T1b <- DATA$age_T1
  }

  #cl <- makeCluster(ncores)
  #clusterExport(cl = cl, varlist = c("likprocFor","peelingFor"))
  # chargement des librairies et programmes FORTRAN utilisées dans les coeurs
  # clusterEvalQ(cl, {
  #   source("./Rfunctions/GENERISK_compileFORTRAN.R")
  # })


  #ETAPE 0 : vérification de la cohérence des données familiales

  #checkit <- checkpedigrees(famid = DATA[,1], id = paste(DATA[,1], DATA[,2], sep="/"), father.id = paste(DATA[,1], DATA[,4], sep="/"), mother.id = paste(DATA[,1], DATA[,5], sep="/"))

  #table(checkit$split) #should be all 1's
  #table(checkit$unrelated) #should be all 0's
  #table(checkit$join) #should be all 0's
  # a faire : traiter les probl?mes de boucle + ind. unrelated ...

  # ETAPE 1 : modèle génétique à n loci

  cat("program initialization. \n")
  G <- genetFor(allef=fA)

  # if(n.points == 3){
  #   aa <- c(20,35,50,120)
  # }else{
  #   if(n.points == 5){
  #     aa <- c(20,35,50,70,90,120)
  #   }else{
  #     if(n.points == 8){
  #       #aa <- c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,120)
  #       aa <- c(5,10,15,20,30,40,50,60,70,80,120)
  #       # aa <- c(30,40,50,60,70,80,120)
  #
  #       n.points <- 10
  #
  #     }else{stop("The number of non parametric points must be specified with the n.points parameter")}
  #   }
  # }

  #ETAPE 2.1 : matrice des paramètres

  #intervalles de definition des pénétrances
  AgeDef <- NULL
  for (cc in seq((6+nloci +2),(6+nloci + 2*ndis),2)){
    AgeDef <- c(AgeDef, list(c(min(DATA[(DATA[,cc-1] == 1) & (DATA[,cc] > 0),cc]), max(DATA[(DATA[,cc-1] == 1) & (DATA[,cc] > 0),cc]))))
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
      }else{ penetmodels[dis] <- 3} #non parametric
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

        aa <- c(FIT.pars[[dis]]$agenodes,120)

        PARAMS <- rbind(PARAMS, matrix(paste(dis,rep(aa, each = length(PARAM)),PARAM,sep=":"), ncol=length(PARAM), byrow=TRUE))
        parrownames <- c(parrownames, paste(disname,aa, sep=":"))
      }
    }
  }

  colnames(PARAMS) <- c(paste("m",PARnames, sep=":"), paste("f",PARnames, sep=":"))
  rownames(PARAMS) <- parrownames

  #remove parameters corresponding to disease with risk=0 for one sex
  ll <- 1
  for (dis in 1:ndis){
    if(penetmodels[dis] == 1){ nl <- 3 #weibull
    }else{
      if(penetmodels[dis] == 2) nl <- 1 #cox
      else{
        n.points <- length(FIT.pars[[dis]]$agenodes)
        nl <- n.points + 1 #np
      }
    }
    if(sum(Ft.pop[[dis]][,1]) == 0){ #no risk for men
      PARAMS[ll:(ll+nl-1),1:ng] <- ""
    }else{
      if(sum(Ft.pop[[dis]][,2]) == 0){#no risk for women
        PARAMS[ll:(ll+nl-1),(ng+1):(2*ng)] <- ""
      }
    }
    if(penetmodels[dis] == 1){ ll <- ll + 3 #weibull
    }else{
      if(penetmodels[dis] == 2) ll <- ll + 1 #cox
      else ll <- ll + n.points + 1 #np
    }

  }

  #remove one genotype (wild type homozygote at all loci), because the associated risk is deduced (pop incidence)
  PARAMS[,c(1,ng+1)] <- ""

  #remove the genotypes never observed in data
  POSS.genotypes <- names(table(rowPaste(DATA[,7:(7+length(fA)-1)])))
  POSS.genotypes <- c(paste("m",POSS.genotypes, sep=":"), paste("f",POSS.genotypes, sep=":"))
  for (k in 1:ncol(PARAMS)) if(!(colnames(PARAMS)[k] %in% POSS.genotypes)) PARAMS[,k] <- ""

  PARAMS.mask <- matrix(0,ncol=ncol(PARAMS), nrow=nrow(PARAMS), dimnames =list(rownames(PARAMS),colnames(PARAMS)))
  parameters <- unique(as.vector(PARAMS[PARAMS != ""]))
  for (k in 1:length(parameters)) PARAMS.mask[PARAMS == parameters[k]] <- k

  # ETAPE 2.2 : calculs preliminaires sur les donn?es familiales (pour ne pas refaire ces calculs dans la boucle d'optimisation et de bootstrap)

  cat("pre-calculations to speed-up likelihood algorithm. \n")
  fids  <- unique(DATA[,1])
  X <- lapply(fids, function(x) preprocFor(ped = DATA[DATA[,1]==x,-1],allef=fA, ndis=ndis))
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
            fitmf <- nlminb(start = c(d=3/5, H1=2/3, H2=1/4),
                            objective = function(x) SCEweibull(x, y=rowMeans(Ft.pop[[dis]])),
                            lower = list(d=0.001, H1=0.001, H2=0.001),
                            upper = list(d=0.999, H1=0.999, H2=0.999))
            pars.init <- c(pars.init, fitmf$par)
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
        for (g in 2:ng){
          if(PARAMS.mask[ll,g] == 0 & PARAMS.mask[ll,g+ng] == 0){
          }else{
            if(PARAMS.mask[ll,g] == PARAMS.mask[ll,g+ng]){ #equal for male and female
              a <- rowMeans(Ft.pop[[dis]])[aa[1:(n.points+1)] + 1]
              b <- rowMeans(Ft.pop[[dis]])[c(1,aa[1:(n.points)] + 1)]
              pars.init <- c(pars.init, a-b)
            }else{
              if((PARAMS.mask[ll,g] != 0) & (PARAMS.mask[ll,g+ng] == 0)){ #only men
                a <- Ft.pop[[dis]][aa[1:(n.points+1)] + 1, "m"]
                b <- Ft.pop[[dis]][c(1,aa[1:(n.points)] + 1), "m"]
                pars.init <- c(pars.init, a-b)
              }else{
                if((PARAMS.mask[ll,g] == 0) & (PARAMS.mask[ll,g+ng] != 0)){#only women
                  a <- Ft.pop[[dis]][aa[1:(n.points+1)] + 1, "f"]
                  b <- Ft.pop[[dis]][c(1,aa[1:(n.points)] + 1), "f"]
                  pars.init <- c(pars.init, a-b)
                }else{ #diff for male and female
                  a <- Ft.pop[[dis]][aa[1:(n.points+1)] + 1, "m"]
                  b <- Ft.pop[[dis]][c(1,aa[1:(n.points)] + 1), "m"]
                  pars.init <- c(pars.init, a-b)
                  a <- Ft.pop[[dis]][aa[1:(n.points+1)] + 1, "f"]
                  b <- Ft.pop[[dis]][c(1,aa[1:(n.points)] + 1), "f"]
                  pars.init <- c(pars.init, a-b)
                }
              }
            }
          }
        }

        ## log-transformation to speed-up and improve convergence with nlminb
        pars.init <- log(pars.init  + 1e-6) ##1e-6 to deal with 0 values

      }
    }
    ll <- ll + nl
  }

  lower.weibull <- c(d=0.001, H1=0.001, H2=0.001)
  upper.weibull <- c(d=0.999, H1=0.999, H2=0.999)
  lower.cox     <- 1/100
  upper.cox     <- 100
  lower.np <- log(1e-6)
  upper.np <- log(1)

  pars.lower <- pars.upper <- numeric(length(parameters))

  for (a in as.character(1:100)){
    pars.lower[grep(a,parameters,fixed = TRUE)] <- lower.np
    pars.upper[grep(a,parameters,fixed = TRUE)] <- upper.np
  }

  for (a in c("d","H1","H2")){
    pars.lower[grep(a,parameters,fixed = TRUE)] <- lower.weibull[a]
    pars.upper[grep(a,parameters,fixed = TRUE)] <- upper.weibull[a]
  }
  pars.lower[grep("HR",parameters,fixed = TRUE)] <- lower.cox
  pars.upper[grep("HR",parameters,fixed = TRUE)] <- upper.cox


  cat("ML optimization by nlminb. \n")
  fit <- nlminb(start = pars.init,  # ETAPES 3 et 4  c(d=3/5, H1=2/3, H2=1/4)
                objective = function(x) nloglik(x, 'G' =G, 'X' = X, 'Ft.pop' = Ft.pop, 'LIK.method' = LIK.method, 'FIT.pars' = FIT.pars, 'approxFt.aa' = approxFt.aa,'approx.np' = approx.np, 'PARAMS.mask' = PARAMS.mask, 'penetmodels' = penetmodels),
                lower = pars.lower,
                upper = pars.upper,
                control = list(trace=2,rel.tol=rel.tol))

  fta <- array(ftvproc(x = fit$par, args = list('Ft.pop' = Ft.pop, 'FIT.pars' = FIT.pars, 'approxFt.aa' = approxFt.aa, 'approx.np' = approx.np, 'PARAMS.mask' = PARAMS.mask, 'G' =G)), dim=c(121,2,ndis,ng))

  if(B > 0){
    cat("Bootstrap in progress... \n")
    pb <- txtProgressBar(char = "*", style =3)
    Sys.sleep(0.1)
    RESboot <- NULL
    ftaboot <- NULL
    for (b in seq(B)){

      # ech bootstrap

      if(imput_ageNA){

        ## imputation des âges NA

        fids.boot <- sample(fids, replace = TRUE)
        DATA.boot <- lapply(fids.boot, FUN=function(f) DATA[DATA$numFam == f,])

        # on change le numfam (cas des familles sélectionnées plusieurs fois...)
        DATA.boot <- bind_rows(DATA.boot, .id = "fid")
        DATA.boot$numFam <- NULL

        ## imputation des âges manquants pour age_t1 et age_t1b
        n_toimpute <- sum(DATA.boot$toimpute)
        DATA.boot$age_T1[DATA.boot$toimpute == 1] <- floor(runif(n = n_toimpute, min = 0, max = 100))
        DATA.boot$age_T1b <- DATA.boot$age_T1

        ## pré-calculs pour FORTRAN
        fids.boot  <- unique(DATA.boot$fid)
        X.boot <- lapply(fids.boot, function(f) preprocFor(ped = DATA.boot[DATA.boot$fid==f,-1],allef=fA, ndis=ndis))
        names(X.boot) <- fids.boot


      }else{

        X.boot <- bootfam(X)

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
    return(list(fit=c(fit,paramsmask = list(PARAMS.mask), ft=list(fta)), boot = list(fit.boot = RESboot, ft.boot = ftaboot), params=list(Ft.pop = Ft.pop,'approx.np' = approx.np ,FIT.pars = FIT.pars, LIK.method=LIK.method, fA=fA, rel.tol=rel.tol, AgeDef = AgeDef), X=X, G=G))
  }else{
    return(list(fit=c(fit,paramsmask = list(PARAMS.mask), ft=list(fta)),params=list(Ft.pop = Ft.pop,'approx.np' = approx.np, FIT.pars = FIT.pars, LIK.method=LIK.method, fA=fA, rel.tol=rel.tol, AgeDef = AgeDef), X=X, G=G))
  }

  #  stopCluster(cl)
}
