#' function checking the consistency of family data with help of the pedtools::ped() function
#'
#' @param fams a dataframe with columns famid, id, fatherid, motherid, sex
#'
#' @importFrom pedtools ped
#'
#' @return to be completed
#'
#' @examples
#' #to be completed
#'
#' @export
#'
checkpedigrees <- function(fams){

  fids <- unique(fams$famid)
  peds <- vector(mode = "list", length = length(fids))
  names(peds) <- as.character(fids)

  for (fi in fids){
    fami <- subset(fams, fams$famid == fi)
    pedi <- pedtools::ped(id = fami$id,
                fid = fami$fatherid,
                mid = fami$motherid,
                sex = fami$sex,
                famid = fi)

    peds[[as.character(fi)]] <- pedi
  }

  return(peds)

}

