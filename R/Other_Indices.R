#' PUC
#'
#' \code{PUC} computes the proportion of uncontaminated correlations for a bifactor mode
#'
#' \code{PUC} is called by \code{\link{bifactorIndices}} and \code{\link{bifactorIndicesMPlus}}, which are the only functions in this package intended for casual users
#'
#' @param Lambda is a matrix of factor loadings
#'
#' @return \code{numeric}
#'
#' @seealso \code{\link{bifactorIndices}}
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' ARPB(MplusAutomation::readModels(file.choose), MplusAutomation::readModels(file.choose))
#' }
#'
PUC <- function(Lambda) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  if (!isBifactor(Lambda)) return(NULL)
  ### Count how many items are on each factor
  numItemsOnFactor <- colSums(Lambda != 0)
  ## contaminated correlations: add up n*(n-1)/2 for all factors, then subtract off the n(n-1)/2 for general factor
  specificCorrelationCount <- sum(sapply(numItemsOnFactor, function (x) {x*(x-1)/2})) - (nrow(Lambda)*(nrow(Lambda)-1)/2)
  ## PUC = 1 - PCC
  1 - specificCorrelationCount/(nrow(Lambda)*(nrow(Lambda)-1)/2)
}



#' ARPB
#'
#' \code{ARPB} computes absolute relative bias in factor loadings between the general factor of a bifactor model and a unidimensional model.
#'
#'\code{ARPB} is called by \code{\link{bifactorIndices}} and \code{\link{bifactorIndicesMPlus}}, which are the only functions in this package intended for casual users
#'
#' @param Lambda is a matrix of factor loadings
#' @param UniLambda is a matrix of factor loadings
#'
#' @return a list where the first element is a vector of absolute relative bias by item, and the second
#' item is
#'
#' @seealso \code{\link{bifactorIndices}}
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' ARPB(MplusAutomation::readModels(file.choose), MplusAutomation::readModels(file.choose))
#' }
#'
#'
## put another example in there from HS data using lavaan (do a 1+3 bifactor, compare to unidimensional)
ARPB <- function(Lambda, UniLambda = NULL) {
  if (is.null(UniLambda)) return(NULL)
  if (is.null(getGen(Lambda))) return(NULL)
  genFac <- getGen(Lambda)
  genLambda <- Lambda[,genFac]
  relBias <- abs((UniLambda - genLambda)/genLambda)
  rownames(relBias) <- rownames(Lambda)
  colnames(relBias) <- c("Abs Rel Bias")
  list(ARPB = mean(relBias), relBias = relBias)
}
