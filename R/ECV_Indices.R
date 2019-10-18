#' ECV_SS_C
#'
#' Computes ECV index for a single factor. Here, ECV is compute with respect to items which load on the factor. In the Excel version of the bifactor indices calculator, this index is referred to as ECV (NEW)
#'
#' @param Fac is the factor for which ECV is being requested
#' @param Lambda is a matrix of factor loadings or an object the function can convert to a matrix of factor loadings
#' @param standardized if the factor loading matrix needs to be extracted, \code{standardized} tells the function whether to extract standardized or unstandardized loadings
#'
#' @return A \code{numeric}, the ECV of the factor with respect to the items in that factor
#'
#' @seealso \code{\link{ECV_SS_All}}
#'
#'
### Compute ECV within-domain for the specific factor (see S&E p. 201)   [ECV_S.S]
  ECV_SS_C <- function(Fac, Lambda, standardized = TRUE) {
    Lambda <- getLambda(Lambda, standardized = standardized)
    ## Isolate the items which load on the chosen factor FAC
    inFactor <- Lambda[,Fac] != 0
    ## Square the loadings
    L2 = Lambda^2
    ## Compute the appropriate ratio of sums
    sum(L2[,Fac]*inFactor)/sum(L2*inFactor)
  }

#' ECV_SS_All
#'
#' Computes ECV index for a all factors. Here, ECV is computed with respect to items which load on the factor. In the Excel version of the bifactor indices calculator, this index is referred to as ECV (NEW)
#'
#' @param Lambda is a matrix of factor loadings or an object the function can convert to a matrix of factor loadings
#' @param standardized if the factor loading matrix needs to be extracted, \code{standardized} tells the function whether to extract standardized or unstandardized loadings
#'
#' @return A vector of ECVs for all factors
#'
#' @export
#'
#' @seealso \code{\link{ECV_SS}}, \code{\link{bifactorIndices}}, \code{\link{bifactorIndicesMPlus}}
#'
ECV_SS_All <- function(Lambda, standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  sapply(colnames(Lambda), ECV_SS_C, Lambda, standardized = standardized)
}

### Compute ECV for each factor (S&E Interpretation, p. 199)   [ECV_S.G]
ECV_SG_C <- function(Fac, Lambda, standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  if (is.null(getGen(Lambda))) return(NULL)
  ## Make a matrix of logical vectors for non-zero elements of Lambda.
  inFactor <- Lambda[,Fac] != 0
  ## Square the loadings
  L2 <- Lambda^2
  ## Compute the appropriate ratio of sums
  sum(L2[,Fac]*inFactor)/sum(L2)
}

ECV_SG_All <- function(Lambda, standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  if (is.null(getGen(Lambda))) return(NULL)
  sapply(colnames(Lambda), ECV_SG_C, Lambda, standardized = standardized)
}

## Compute ECV within-domain version (see S&E p. 201)   [ECV_G.S]
ECV_GS_C <- function(Fac, Lambda, standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  if (is.null(getGen(Lambda))) return(NULL)
  genFac <- getGen(Lambda)
  ## Make a matrix of logical vectors for non-zero elements of Lambda.
  inFactor <- Lambda[,Fac] != 0
  ## Square the loadings
  L2 <- Lambda^2
  ## Compute the appropriate ratio of sums
  sum(L2[,genFac]*inFactor)/sum(L2*inFactor)
}


ECV_GS_All <- function(Lambda, standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  if (is.null(getGen(Lambda))) return(NULL)
  sapply(colnames(Lambda), ECV_GS_C, Lambda, standardized = standardized)
}

IECV_C <- function(Item, Lambda, standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  if (is.null(getGen(Lambda))) return(NULL)
  genFac <- getGen(Lambda)
  ## Square the loadings
  L2 <- Lambda^2
  IECV <- L2[Item, genFac]/sum(L2[Item,])
  names(IECV) <- rownames(Lambda)[Item]
  IECV
}

IECV_All <- function(Lambda, standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  if (is.null(getGen(Lambda))) return(NULL)
  genFac <- getGen(Lambda)
  sapply(1:nrow(Lambda), IECV_C, Lambda, standardized = standardized)
}
