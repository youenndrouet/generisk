#' function to compute the cumulative risk function Ft from incidence rates
#'
#' @param rates Populational incidence rates by 1-year from 0y to 120y
#' @param rates_denom Number "at risk" for the incidence rates (eg 100,000)
#' @param smo Logical, smooth the Ft function or not ? Without smoothing the Ft function can reflect the age classes.
#' @param smo_df degree of freedom for the smooth function
#'
#' @importFrom stats smooth.spline
#'
#' @return the cumulative risk function Ft from incidence rate data
#' @export
#'
#' @examples
#' #to be completed
Ft_compute <- function(rates, rates_denom = 100000, smo = TRUE, smo_df = 15){

  if(length(rates) != 121){
    stop("rates must be a vector of length=121, going from age=0 to age=120.")
  }

  if(smo){
    rates_smo <- smooth.spline(rates, df=smo_df)$y
    wiout <- which(rates_smo <= 0)
    rates_smo[wiout] <- rates[wiout]
  }else{
    rates_smo <- rates
  }

  Ft <- numeric(121)
  for(t in (2:121)){
    Ft[t] = Ft[t-1] + (1 - Ft[t-1])*(rates_smo[t]/rates_denom)
  }

  names(Ft) <- as.character(0:120)
  Ft["0"] <- 0

  return(Ft)

}
