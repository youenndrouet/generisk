
#' @useDynLib generisk genet_f

genet_Fortran <- function(hwpr, gprv, ngen, af, nl) {
  .Fortran(genet_f,
           hwpr = as.double(hwpr),
           gprv = as.double(gprv),
           ngen = as.integer(ngen),
           af = as.double(af),
           nl = as.integer(nl))
}

#' @export
#'
genetFor <- function(allef){ # R wrapper

  nloci <- length(allef)
  ng    <- 3^nloci
  res <- genet_Fortran(
    hwpr = numeric(ng),
    gprv = numeric(ng*ng*ng),
    ngen = ng,
    af   = allef,
    nl   = nloci)

  return(res)

}


