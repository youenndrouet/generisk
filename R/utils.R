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

compareModels <- function(x1,x2){
  nlkl1 <- x1$fit$objective
  nlkl2 <- x2$fit$objective
  npar1 <- length(x1$fit$par)
  npar2 <- length(x2$fit$par)
  lr <- 2*abs(nlkl1 - nlkl2)
  pvalue <- 1-pchisq(lr,df=abs(npar2-npar1))
  return(list('lr'=lr,'pvalue' = pvalue))
}

colpct <- function(col,pct) return(paste(col,pct,sep=""))

empty.plot <- function(type="absolute", ymax, xmax, marg, main){

  par(mar=marg,las=1)
  if(type=="absolute"){
    plot(NA,NA, type="n", xlim=c(0,xmax), ylim=c(0,ymax), ann=FALSE, bty="n",axes=FALSE)
    #rect(0,0,100,1, col=gray(.96))
    abline(h = seq(0,1,0.02), v = seq(10,90,10),col=gray(.95),lty=3 )
    abline(h = seq(0,1,0.1) , v = seq(0,90,20), col=gray(.85),lty=3)
    axis(2,lwd=1,col = gray(.7), at = seq(0,1,0.1), labels=seq(0,1,0.1)*100)
    title(ylab="",main="Cumulative risk (%)", cex.lab = 1.3)
  }else{
    plot(NA,NA, type="n", xlim=c(0,xmax), ylim=c(0.1,100000), log="y", xlab="", ylab ="", bty="n",axes=FALSE)
    abline(h = c(seq(0.2,0.9,0.1),
                 seq(2,9,1),
                 seq(20,90,10),
                 seq(200,900,100),
                 seq(2000,9000,1000),
                 seq(20000,90000,10000)),
           v = seq(10,90,10),
           col=gray(.95),lty=3 )
    abline(h = c(.01,.1,10,100,1000,10000), v = seq(0,90,20),col=gray(.85),lty=3)
    axis(2, at=c(0.01,0.1,1,10,100,1000,1e4,1e5 ), labels=c(".01", ".1","1","10","100","1000", 1e4, 1e5),lwd=1,col = gray(.7))
    title(ylab="",main="Risk ratio to pop. incidence \n (log-scale)", cex.lab = 1.3)
  }
  title(xlab="Age (y)", cex.lab = 1.3)
  axis(1,lwd=1,col = gray(.7))

}

