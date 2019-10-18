#' bifactorIndices
#'
#' Computes all available bifactor indices for the input given.
#'
#' @param Lambda is a matrix of factor loadings or an object that can be converted to a matrix of factor loadings by \code{\link{getLambda}}. Currently fitted lavaan objects and fitted mirt objects are supported in addition to raw factor loading matrix input. For Mplus output files, use \code{\link{bifactorIndicesMplus}}.
#' @param Theta is a vector of residual variances. If omitted, Theta will be computed from input for Lambda
#' @param UniLambda is a matrix of factor loadings or an object that can be converted to a matrix of factor loadings by \code{\link{getLambda}}
#' @param standardized lets the function know whether to look for standardized or unstandardized results from Mplus or lavaan. If lambda is not a lavaan object, then \code{standardized} will be ignored.
#'
#' @return a list of bifactor indices, including three different ECV indices, IECV, PUC, Omega, OmegaH, and ARPB. Please note that many of these indices are interpretable even when the model being used is not a bifactor model; some indices may be useful for two-tier, trifactor, correlated traits, and even unidimensional models.
#'
#' @details Currently, factor loading matrices, fitted lavaan objects, and fitted mirt objects are supported. For Mplus output, see \code{\link{bifactorIndicesMplus}}. IRT parameters from mirt are converted to standardized factor loadings via the correspondence described in Kamata & Bauer (2008). If you wish to use standardized coefficients, item error variance will be computed directly from standardized factor loadings. ARPB will only be computed if the factor loadings from a unidimensional model are included, while PUC, ECV_GS, and ECV_SG will only be computed if the the model is a true bifactor model. Note that if a correlated traits model is provided, the omega indices will simply be the regular omega values for those factors. Interpretations for individual indices as well as details about their computation can be found in the man page for the individual indices.
#'
#' @export
#'
#' @seealso \code{\link{bifactorIndicesMplus}},
#'          \code{\link{ECV_SS}},
#'          \code{\link{ECV_SG}},
#'          \code{\link{ECV_GS}},
#'          \code{\link{IECV}},
#'          \code{\link{PUC}},
#'          \code{\link{Omega_S}},
#'          \code{\link{Omega_H}},
#'          \code{\link{ARPB}}
#'
#' @examples
#'
#' GENERATE MATRICES, SEND IT TO BINDICES
#'
#' USE SOME LAVAAN DATA, CREATE BIFACTOR MODEL, FEED IT TO BINDICES
#'
#' USE BINDICES TO COMPUTE OMEGA
#'


bifactorIndices <- function(Lambda, Theta = getTheta(Lambda), UniLambda = NULL, standardized = TRUE) {
  Lambda <- getLambda(Lambda, standardized = standardized)
  UniLambda <- getLambda(UniLambda, standardized = standardized)
  indicesList <- list(ECV_SS  = ECV_SS(Lambda),
                      ECV_SG  = ECV_SG(Lambda),
                      ECV_GS  = ECV_GS(Lambda),
                      IECV    = IECV(Lambda),
                      PUC     = PUC(Lambda),
                      Omega   = Omega_S(Lambda),
                      Omega_H = Omega_H(Lambda),
                      ARPB    = ARPB(Lambda, UniLambda)
  )
  indicesList[which(!sapply(indicesList, is.null))]
}

#' bifactorIndicesMPlus
#'
#' Computes all available bifactor indices given an Mplus .out file for a bifactor model
#'
#' @param Lambda is an Mplus .out file. Defaults to an open file dialog box
#' @param UniLambda is an object that the function can convert to a matrix of factor loadings. The expected behavior is to store an Mplus output file as a variable and pass that variable as UniLambda. Defaults to \code{NULL}, as UniLambda is only required if you wish to compute \code{\link{ARPB}}
#' @param standardized lets the function know whether it should be looking in
#'   the unstandardized results or the STDYX results from the Mplus output.
#'
#' @return a list of bifactor indices, including three different ECV indices, IECV, PUC, Omega, OmegaH, and ARPB. Please note that many of these indices are interpretable even when the model being used is not a bifactor model; some indices may be useful for two-tier, trifactor, correlated traits, and even unidimensional models.
#'
#' @details ARPB will only be compute if the factor loadings from a unidimensional model (as a vector or as the result of using MplusAutomation::readModels() on an Mplus .out file) are included. Note that if a correlated traits model is provided, the omega indices will simply be the regular omega values for those factors. Interpretations for individual indices as well as details about their computation can be found in the man page for the individual indices.
#'
#' @seealso \code{\link{bifactorIndices}},
#'          \code{\link{ECV_SS}},
#'          \code{\link{ECV_SG}},
#'          \code{\link{ECV_GS}},
#'          \code{\link{IECV}},
#'          \code{\link{PUC}},
#'          \code{\link{Omega_S}},
#'          \code{\link{Omega_H}},
#'          \code{\link{ARPB}}
#'
#'
#' @export
#'
#' @examples
#' FIGURE OUT A WAY TO DO AN EXAMPLE OR THREE
#'
bifactorIndicesMPlus <- function(Lambda = file.choose(), UniLambda = NULL, standardized = TRUE) {
  Lambda <- MplusAutomation::readModels(Lambda)
  Theta <- getTheta(Lambda, standardized = standardized)
  Lambda <- getLambda(Lambda, standardized = standardized)
  bifactorIndices(Lambda, Theta, UniLambda, standardized)
}
