#' Omega_S
#'
#' Computes omega reliability estimate for all factors as described in Rodriguez, Reise, and Haviland (2016).
#'
#' \code{Omega_S} is called by \code{\link{bifactorIndices}} and \code{\link{bifactorIndicesMPlus}}, which are the only functions in this package intended for casual users
#'
#' @param Lambda is a matrix of factor loadings
#' @param Theta is a vector of indicator error variances
#'
#' @return A \code{numeric}, the omega reliability estimate for all factors.
#'
#' @examples
#'
#' @section References:
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models: calculating and interpreting statistical indices. Psychological methods, 21(2), 137 <doi:10.1037/met0000045>.
#'
#' @export
#'
#' @seealso \code{\link{Omega_H}}, \code{\link{bifactorIndices}}, \code{\link{bifactorIndicesMPlus}}
#'
#'

Omega_S <- function(Lambda, Theta) {
  Omega_S_C <- function(Fac, Lambda, Theta) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Compute the appropriate ratio of sums
    sum(colSums(Lambda*inFactor)^2)/(sum(colSums(Lambda*inFactor)^2) + sum(Theta*inFactor))
  }
  sapply(colnames(Lambda), Omega_S_C, Lambda = Lambda, Theta = Theta, standardized = standardized)
}


#' Omega_H
#'
#' Computes hierarchical omega reliability estimate for all factors as described in Rodriguez, Reise, and Haviland (2016).
#'
#' \code{Omega_H} is called by \code{\link{bifactorIndices}} and \code{\link{bifactorIndicesMPlus}}, which are the only functions in this package intended for casual users
#'
#' @param Lambda is a matrix of factor loadings
#' @param Theta is a vector of indicator error variances
#'
#' @return A \code{numeric}, the omega reliability estimate for all factors.
#'
#' @examples
#'
#' @section References:
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models: calculating and interpreting statistical indices. Psychological methods, 21(2), 137 <doi:10.1037/met0000045>.
#'
#' @export
#'
#' @seealso \code{\link{Omega_S}}, \code{\link{bifactorIndices}}, \code{\link{bifactorIndicesMPlus}}
#'
#'

Omega_H <- function(Lambda, Theta) {
  Omega_H_C <- function(Fac, Lambda, Theta) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda. Let's replace NA with zero at the start!!
    inFactor <- Lambda[,Fac] != 0
    ## Compute the appropriate ratio of sums
    sum(Lambda[,Fac])^2/(sum(colSums(Lambda*inFactor)^2) + sum(Theta*inFactor))
  }
  sapply(colnames(Lambda), Omega_H_C, Lambda = Lambda, Theta = Theta)
}

