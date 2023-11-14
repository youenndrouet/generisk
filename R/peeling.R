
#' @useDynLib generisk peeling_f

peeling_Fortran <- function(counsid, pedv, likv, counspr, hwpr,
                            gprv, idsv, w, ll, ngen, ni, km){
  .Fortran(peeling_f,
           counsid = as.integer(counsid),
           pedv = as.integer(pedv),
           likv = as.double(likv),
           counspr = as.double(counspr),
           hwpr = as.double(hwpr),
           gprv = as.double(gprv),
           idsv = as.integer(idsv),
           w = as.double(w),
           ll = as.double(ll),
           ngen = as.integer(ngen),
           ni = as.integer(ni),
           km = as.integer(km))

}

peelingFor <- function(X, G, LIKv, counselee.id){

      res <- peeling_Fortran(
        counsid = counselee.id,
        pedv    = X$pedv,
        idsv    = X$idsv,
        km      = X$km,
        counspr = numeric(X$ngen),
        hwpr    = G$hwpr,
        gprv    = G$gprv,
        w       = numeric(1),
        ll      = numeric(1),
        likv    = LIKv,
        ni      = X$ni,
        ngen    = X$ngen)

    return(list("counselee.Pr" = res$counspr, "loglik" = res$ll))
}
