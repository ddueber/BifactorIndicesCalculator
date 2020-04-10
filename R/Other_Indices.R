#' PUC
#'
#' \code{PUC} computes the proportion of uncontaminated correlations for a bifactor mode
#'
#' \code{PUC} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
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
#' PUC(Lambda)
#'
#'
PUC <- function(Lambda) {
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
#' \code{ARPB} computes absolute relative bias in factor loadings between the general factor of a
#' bifactor model and a unidimensional model.
#'
#'\code{ARPB} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of factor loadings
#' @param UniLambda is a matrix of factor loadings
#'
#' @return a list where the first element is the average absolute relative paramter bias, and the second
#' element is a vector of absolute relative bias by item
#'
#' @seealso \code{\link{bifactorIndices}}
#'
#' @export
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
#' UniLambda <- c(.78, .84, .82, .77, .69, .62, .69, .66, .82, .56, .74, .65)
#' ARPB(Lambda, UniLambda)
#'
#'
ARPB <- function(Lambda, UniLambda) {
  if (is.null(UniLambda)) return(NULL)
  if (is.null(getGen(Lambda))) return(NULL)
  genFac <- getGen(Lambda)
  genLambda <- Lambda[,genFac]
  relBias <- abs((UniLambda - genLambda)/genLambda)
  names(relBias) <- rownames(Lambda)
  list(ARPB = mean(relBias), AbsRelBias = relBias)
}



#' Factor Determinacy
#'
#' \code{FD} computes factor determinacies for all factors provided
#' standardized factor loadings and an interfactor correlation matrix.
#'
#' \code{FD} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of standardized factor loadings
#' @param Phi is the matrix of factor intercorrelations. For bifactor models
#' \code{Phi} is diagonal with ones on the diagonal.
#'
#' @return a vector of factor determinacies.
#'
#' @seealso \code{\link{bifactorIndices}}
#'
#' @export
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
#' Phi <- matrix(c(1, 0, 0, 0,
#'                 0, 1, 0, 0,
#'                 0, 0, 1, 0,
#'                 0, 0, 0, 1), ncol = 4)
#' colnames(Lambda) <- c("General", "SF1", "SF2", "SF3")
#' FD(Lambda, Phi)
#'

FD <- function(Lambda, Phi) {
  Psi   <- getTheta(Lambda)
  Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(Psi)
  FacDet <- sqrt(diag(Phi %*% t(Lambda) %*% solve(Sigma) %*% Lambda %*% Phi))
  names(FacDet) <- colnames(Lambda)
  FacDet
}


#' Construct Replicability
#'
#' \code{H} computes construct replicability for all factors given
#' standardized factor loadings.
#'
#' \code{H} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of standardized factor loadings
#'
#' @return a vector of construct reliabilities.
#'
#' @seealso \code{\link{bifactorIndices}}
#'
#' @export
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
#' H(Lambda)
#'

H <- function(Lambda) {
  1/(1+1/(colSums(Lambda^2/(1-Lambda^2))))
}

