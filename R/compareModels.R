#' function to compare two models
#'
#' @param model1 a model object returned by the generisk() function
#' @param model2 a model object returned by the generisk() function
#'
#' @return a list of two: likelihood ratio and p-value
#' @export
#'
#' @examples
#' #to be completed
#'

compareModels <- function(model1,model2){
  nlkl1 <- model1$fit$objective
  nlkl2 <- model2$fit$objective
  npar1 <- length(model1$fit$par)
  npar2 <- length(model2$fit$par)
  lr <- 2*abs(nlkl1 - nlkl2)
  pvalue <- 1-pchisq(lr,df=abs(npar2-npar1))
  return(list('lr'=lr,'pvalue' = pvalue))
}
