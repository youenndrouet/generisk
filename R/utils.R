#' @import graphics
#' @importFrom stats pchisq
#' @importFrom grDevices gray

Constrain01 <- function(f, maxrisk = 0.99999, minrisk = 1e-16) return(ifelse(f>=maxrisk,maxrisk,ifelse(f<=minrisk,minrisk,f)))

Weibull.modified <- function(par,t=0:120){
    out <- Constrain01((1 - par$k)*(1 - exp(-(par$l*t)^par$a)))
    out[1] <- 0.
    names(out) <- t
    return(out)
}

pmin_na <- function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = T))

#find weibull param
SCEweibull <- function(x,y) return(sum((y-Weibull.modified(dH1H2toakl(list("d"=x[1],"H1"=x[2],"H2"=x[3]))))^2))

akltodH1H2 <- function(par, t1 = 40, t2 = 70){
    d <- 1 - par$k
    #P1 <- Weibull.modified(a,k,l,t1)
    #P2 <- Weibull.modified(a,k,l,t2)
    #H1 = (d - P1)/d
    #H2 = (d - P2)/(d - P1)
    return(list(d = d, H1 = exp(-(par$l*t1)^par$a), H2 = exp(-(par$l*t2)^par$a)/exp(-(par$l*t1)^par$a)))
    #d is the lifetime penetrance
    #H1 is the proportion of unaffected at t1 years among at risk individuals
    #H2 is the proportion of unaffected at t2 years among unaffected at t1 but at risk at t2.
    #These 3 parameters are independent between them and are defined between 0 and 1.
    #Practical experience shows that it is always much easier to find solutions with d, H1 and H2 than with a, k and l.
}

dH1H2toakl <- function(par, t1 = 40, t2 = 70){
    d  <- par$d
    H1 <- par$H1
    H2 <- par$H2
    k  <- 1 - d
    a  <- (log(-log(H1*H2)) - log(-log(H1)))/(log(t2) - log(t1))
    l  <- exp(log(-log(H1))/a - log(t1))
    return(list(a = a, k = k, l = l))
}

rowPaste <- function(x){
  if(is.vector(x)) return(x)
  else{
    res <- x[,1]
    if(ncol(x)>=2) for (k in 2:ncol(x)) res <- paste(res,x[,k], sep="")
    return(res)
  }
}


bootfam <- function(x, n=NULL){

  if(is.null(n)) size = length(x)
  else size = n

  pos.f.resample <- sample(1:size, size, replace = TRUE, prob = NULL)
  x.resample <- x[pos.f.resample]

  names(x.resample) <- 1:size

  return(x.resample)

}


constrain1 <- function(x) ifelse(x>1,1,x)




