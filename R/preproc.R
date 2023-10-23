
#' @useDynLib generisk preproc_f
#'
preproc_Fortran <- function(pedbrutv , pedv, idsv, disv, agev, tesv, ascv, ngen, ni, ndi, nl, km){
  .Fortran(preproc_f,
           pedbrutv = as.integer(pedbrutv),
           pedv = as.integer(pedv),
           idsv = as.integer(idsv),
           disv = as.integer(disv),
           agev = as.integer(agev),
           tesv = as.double(tesv),
           ascv = as.integer(ascv),
           ngen = as.integer(ngen),
           ni = as.integer(ni),
           ndi = as.integer(ndi),
           nl = as.integer(nl),
           km = as.integer(km))
}

#' @export
#'
preprocFor <- function(ped, allef, ndis){

    kmax  <- 50
    n     <- nrow(ped)
    nloci <- length(allef)
    ng    <- 3^nloci

    res <- preproc_Fortran(
        pedbrutv = as.numeric(as.matrix(ped)),
        pedv = integer(n*4),
        idsv = integer(n*5*kmax),
        disv = integer(n*ndis),
        agev = integer(n*ndis),
        tesv = integer(n*ng),
        ascv = integer(n*ng),
        ngen = ng,
        nl   = nloci,
        ni   = n,
        ndi  = ndis,
        km   = kmax)

    return(res)
}




