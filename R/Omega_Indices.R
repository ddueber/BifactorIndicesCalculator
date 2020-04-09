#' Omega_S
#'
#' Computes an omega reliability estimate for all factors as described in Rodriguez, Reise, and
#' Haviland (2016).
#'
#' \code{Omega_S} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of factor loadings
#' @param Theta is a vector of indicator error variances
#'
#' @return A \code{numeric}, the omega reliability estimate for all factors.
#'
#' @examples
#' Lambda <- matrix(c(.82, .10,   0,   0,
#'                    .77, .35,   0,   0,
#'                    .79, .32,   0,   0,
#'                    .66, .39,   0,   0,
#'                    .51,   0, .71,   0,
#'                    .56,   0, .43,   0,
#'                    .68,   0, .13,   0,
#'                    .60,   0, .50,   0,
#'                    .83,   0,   0, .47,
#'                    .60,   0,   0, .27,
#'                    .78,   0,   0, .28,
#'                    .55,   0,   0, .75),
#'                    ncol = 4, byrow = TRUE)
#' colnames(Lambda) <- c("General", "SF1", "SF2", "SF3")
#' Theta <- rep(1, nrow(Lambda)) - rowSums(Lambda^2)
#' Omega_S(Lambda, Theta)
#'
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models:
#' calculating and interpreting statistical indices. \emph{Psychological Methods, 21}(2),
#' 137 \doi{10.1037/met0000045}.
#'
#' @export
#'
#' @seealso \code{\link{Omega_H}}, \code{\link{bifactorIndices}}
#'
#'

Omega_S <- function(Lambda, Theta) {
  Omega_S_C <- function(Fac, Lambda, Theta) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Compute the appropriate ratio of sums
    sum(colSums(Lambda*inFactor)^2)/(sum(colSums(Lambda*inFactor)^2) + sum(Theta*inFactor))
  }
  if (is.null(Theta)) return(NULL)
  omega_results <- sapply(1:ncol(Lambda), Omega_S_C, Lambda = Lambda, Theta = Theta)
  names(omega_results) <- colnames(Lambda)
  omega_results
}


#' OmegaH
#'
#' Computes hierarchical omega reliability estimate for all factors as described in
#' Rodriguez, Reise, and Haviland (2016).
#'
#' \code{Omega_H} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of factor loadings
#' @param Theta is a vector of indicator error variances
#'
#' @return A \code{numeric}, the omega reliability estimate for all factors.
#'
#' @examples
#' Lambda <- matrix(c(.82, .10,   0,   0,
#'                    .77, .35,   0,   0,
#'                    .79, .32,   0,   0,
#'                    .66, .39,   0,   0,
#'                    .51,   0, .71,   0,
#'                    .56,   0, .43,   0,
#'                    .68,   0, .13,   0,
#'                    .60,   0, .50,   0,
#'                    .83,   0,   0, .47,
#'                    .60,   0,   0, .27,
#'                    .78,   0,   0, .28,
#'                    .55,   0,   0, .75),
#'                    ncol = 4, byrow = TRUE)
#' colnames(Lambda) <- c("General", "SF1", "SF2", "SF3")
#' Theta <- rep(1, nrow(Lambda)) - rowSums(Lambda^2)
#' Omega_H(Lambda, Theta)
#'
#' @section References:
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models:
#' Calculating and interpreting statistical indices. Psychological methods, 21(2),
#' 137 \doi{10.1037/met0000045}.
#'
#' @export
#'
#' @seealso \code{\link{Omega_S}}, \code{\link{bifactorIndices}}
#'
#'

Omega_H <- function(Lambda, Theta) {
  Omega_H_C <- function(Fac, Lambda, Theta) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Compute the appropriate ratio of sums
    sum(Lambda[,Fac])^2/(sum(colSums(Lambda*inFactor)^2) + sum(Theta*inFactor))
  }
  if (is.null(Theta)) return(NULL)
  omega_results <- sapply(1:ncol(Lambda), Omega_H_C, Lambda = Lambda, Theta = Theta)
  names(omega_results) <- colnames(Lambda)
  omega_results
}
