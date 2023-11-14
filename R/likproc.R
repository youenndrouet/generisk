
#' @useDynLib generisk likproc_f

likproc_Fortran <- function(pedbrutv, ftv, disv, tesv, ascv, agev, liknumv,
                            likdenomv, ngen, ni, ndi, nl){
  .Fortran(likproc_f,
          pedbrutv = as.integer(pedbrutv),
          ftv = as.double(ftv),
          disv = as.integer(disv),
          tesv = as.integer(tesv),
          ascv = as.integer(ascv),
          agev = as.integer(agev),
          liknumv = as.double(liknumv),
          likdenomv = as.double(likdenomv),
          ngen = as.integer(ngen),
          ni = as.integer(ni),
          ndi = as.integer(ndi),
          nl = as.integer(nl))

}


likprocFor <- function(X, ftv){
    res <- likproc_Fortran(
        pedbrutv = X$pedbrutv,
        ftv  = ftv,
        disv = X$disv,
        agev = X$agev,
        tesv = X$tesv,
        ascv = X$ascv,
        liknumv    = double(X$ni*X$ngen),
        likdenomv  = double(X$ni*X$ngen),
        ngen = X$ngen,
        ni   = X$ni,
        ndi  = X$ndi,
        nl   = X$nl)

    return(res)
}



