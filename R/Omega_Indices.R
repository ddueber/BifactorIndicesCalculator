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

### Compute Omega, OmegaS
Omega_S_C <- function(Fac, Lambda, Theta = getTheta(Lambda, standardized = standardized), standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  ## Make a matrix of logical vectors for non-zero elements of Lambda.
  inFactor <- Lambda[,Fac] != 0
  ## Compute the appropriate ratio of sums
  sum(colSums(Lambda*inFactor)^2)/(sum(colSums(Lambda*inFactor)^2) + sum(Theta*inFactor))
}


Omega_S_All <- function(Lambda, Theta = getTheta(Lambda, standardized = TRUE), fac_var = rep(1, nrow(lambda)) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  sapply(colnames(Lambda), Omega_S_C, Lambda = Lambda, Theta = Theta, standardized = standardized)
}

  ### Compute Omega_H, Omega_HS
  Omega_H_C <- function(Fac, Lambda, Theta = getTheta(Lambda, standardized = standardized), standardized = TRUE) {
    Lambda <- getLambda(Lambda, standardized = standardized)
    ## Make a matrix of logical vectors for non-zero elements of Lambda. Let's replace NA with zero at the start!!
    inFactor <- Lambda[,Fac] != 0
    ## Compute the appropriate ratio of sums
    sum(Lambda[,Fac])^2/(sum(colSums(Lambda*inFactor)^2) + sum(Theta*inFactor))
  }

  Omega_H_All <- function(Lambda, Theta = getTheta(Lambda), standardized = TRUE) {
    Lambda <- getLambda(Lambda, standardized = standardized)
    sapply(colnames(Lambda), Omega_H_C, Lambda = Lambda, Theta = Theta, standardized = standardized)
  }

}
