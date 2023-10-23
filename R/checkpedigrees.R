#' function largely inspired from the "familycheck" fct of the Kinship package
#'
#' @param famid family identification number (vector)
#' @param id individual identification number (vector)
#' @param father.id father id
#' @param mother.id mother id
#'
#' @importFrom kinship2 makefamid
#'
#' @return to be completed
#'
#' @examples
#' #to be completed
#'
#' @export
#'
checkpedigrees <- function (famid, id, father.id, mother.id){
  nfam <- length(unique(famid))
  newfam <- makefamid(id, father.id, mother.id)
  xtab <- table(famid, newfam)
  if (any(newfam == 0)) {
    print(paste("family",famid[newfam == 0],"unrelated individual", id[newfam == 0]))
    unrelated <- xtab[, 1]
    xtab <- xtab[, -1, drop = FALSE]
  }
  else unrelated <- rep(0, nfam)
  splits <- apply(xtab > 0, 1, sum)
  joins <- apply(xtab > 0, 2, sum)
  temp <- apply((xtab > 0) * outer(rep(1, nfam), joins - 1),
                1, sum)
  out <- data.frame(famid = dimnames(xtab)[[1]], n = as.vector(table(famid)),
                    unrelated = as.vector(unrelated), split = as.vector(splits),
                    join = temp, row.names = 1:nfam)
  if (any(joins > 1)) {
    tab1 <- xtab[temp > 0, ]
    tab1 <- tab1[, apply(tab1 > 0, 2, sum) > 0]
    dimnames(tab1) <- list(dimnames(tab1)[[1]], NULL)
    attr(out, "join") <- tab1
  }
  out
}

