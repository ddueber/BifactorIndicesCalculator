#' bifactorIndices
#'
#' Computes all available bifactor indices for the input given.
#'
#' @param Lambda is a matrix of factor loadings or an object that can be converted to a
#' matrix of factor loadings by \code{\link{getLambda}}. Currently fitted \pkg{lavaan}
#' objects and fitted \pkg{mirt} objects are supported in addition to raw factor loading
#' matrix input. For \code{Mplus} output files, use \code{\link{bifactorIndicesMplus}}.
#' @param Theta is a vector of residual variances. If omitted, \code{Theta} will be computed from
#' input for \code{Lambda}.
#' @param UniLambda is a matrix of factor loadings or an object that can be converted to
#' a matrix of factor loadings such as a fitted \pkg{lavaan} objects or fitted \pkg{mirt}
#' object. Defaults to \code{NULL}, as \code{UniLambda} is only required if you wish to
#' compute \code{\link{ARPB}}
#' @param standardized lets the function know whether to look for standardized or
#' unstandardized results from \pkg{lavaan}. If \code{Lambda} is not a \pkg{lavaan} object,
#' then \code{standardized} will be ignored.
#'
#' @return A list of bifactor indices, including three different ECV indices, IECV, PUC,
#' Omega, OmegaH, and ARPB. Please note that many of these indices are interpretable even
#' when the model being used is not a bifactor model; some indices may be useful for
#' two-tier, trifactor, correlated traits, and even unidimensional models.
#'
#' @details Currently, factor loading matrices, fitted \pkg{lavaan} objects, and fitted \pkg{mirt}
#' objects are supported. For \code{Mplus} output, see \code{\link{bifactorIndicesMplus}}.
#' IRT parameters from \pkg{mirt} are converted to standardized factor loadings via the
#' correspondence described in Kamata & Bauer (2008). If you wish to use standardized
#' coefficients, item error variance will be computed directly from standardized factor
#' loadings. \code{\link{ARPB}} will only be computed if the factor loadings from a unidimensional model
#' are included, while \code{\link{ECV_GS}} and \code{\link{ECV_SG}} will only be computed for
#' models with a general factor, and \code{\link{PUC}} will only be conputed for a true bifactor
#' model. Note that if a correlated traits model is provided, the omega indices
#' will simply be the regular omega values for those factors. Interpretations for individual
#' indices as well as details about their computation can be found in the man page for the
#' individual indices.
#'
#' @references
#' Kamata, A., & Bauer, D. J. (2008). A note on the relation between factor analytic and item
#' response theory models. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 15}
#' (1), 136-153.
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
#' # Computing bifactor indices from fitted lavaan object
#' # (using mirt object is similar)
#' HS_model_bifactor <- "visual  =~ x1 + x2 + x3
#'                       textual =~ x4 + x5 + x6
#'                       speed   =~ x7 + x8 + x9
#'                       general =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
#' # lavaan cannot find a good solution, but that's ok since this is just for illustration
#' fit <- suppressWarnings(lavaan::cfa(HS_model_bifactor,
#'                    data = lavaan::HolzingerSwineford1939,
#'                    orthogonal = TRUE))
#' bifactorIndices(fit)
#'
#'
#' # Computing bifactor indices from standardized factor loading matrices
#' Lambda <-  matrix(c(.82, .10,   0,   0,
#'                     .77, .35,   0,   0,
#'                     .79, .32,   0,   0,
#'                     .66, .39,   0,   0,
#'                     .51,   0, .71,   0,
#'                     .56,   0, .43,   0,
#'                     .68,   0, .13,   0,
#'                     .60,   0, .50,   0,
#'                     .83,   0,   0, .47,
#'                     .60,   0,   0, .27,
#'                     .78,   0,   0, .28,
#'                     .55,   0,   0, .75),
#'                     ncol = 4, byrow = TRUE)
#' colnames(Lambda) <- c("General", "SF1", "SF2", "SF3")
#' UniLambda <- c(.78, .84, .82, .77, .69, .62, .69, .66, .82, .56, .74, .65)
#' bifactorIndices(Lambda, UniLambda = UniLambda)
#'
bifactorIndices <- function(Lambda, Theta = NULL, UniLambda = NULL, standardized = TRUE) {
  ## if fitted mirt object, then throw warning about Omegas probably being meaningless
  if ("SingleGroupClass" %in% class(Lambda)) {
    warning("Interpreting omega indices for IRT models is not recommended at this time")
  }

  ## Make Lambda, Theta, and UniLambda matrices. Do Theta first because getTheta needs the
  ## Lambda object, not the Lambda matrix
  if (is.null(Theta)) {Theta = getTheta(Lambda, standardized = standardized)}

  Lambda <- getLambda(Lambda, standardized = standardized)

  if (!is.null(UniLambda)) {UniLambda <- getLambda(UniLambda, standardized = standardized)}

  ## Build up the lists of indices. FactorLevelIndices first
  FactorLevelIndices = list(ECV_SS  = ECV_SS(Lambda),
                            ECV_SG  = ECV_SG(Lambda),
                            ECV_GS  = ECV_GS(Lambda),
                            Omega   = Omega_S(Lambda, Theta),
                            Omega_H = Omega_H(Lambda, Theta))
  ## Remove any NULL values and convert to dataframe
  FactorLevelIndices <- FactorLevelIndices[which(!sapply(FactorLevelIndices, is.null))]
  FactorLevelIndices <- as.data.frame(FactorLevelIndices)


  ## Item level indices next
  ARPB_indices <- ARPB(Lambda, UniLambda)
  ItemLevelIndices <- list(IECV             = IECV(Lambda),
                           AbsParameterBias = ARPB_indices[[1]],
                           RelParameterBias = ARPB_indices[[2]])

  ## Remove any NULL values and convert to dataframe
  ItemLevelIndices <- ItemLevelIndices[which(!sapply(ItemLevelIndices, is.null))]
  ItemLevelIndices <- as.data.frame(ItemLevelIndices)
  if (isTRUE(all.equal(dim(ItemLevelIndices), c(0, 0)))) {ItemLevelIndices <- NULL}

  ## Model level indices next
  if (is.null(getGen(Lambda))) { # No general factor
    ECV <- NULL
    Omega <- NULL
    OmegaH <- NULL
  } else {
    Gen <- getGen(Lambda)
    ECV <- ECV_SG(Lambda)[Gen]
    Omega <- Omega_S(Lambda, Theta)[Gen]
    OmegaH <- Omega_H(Lambda, Theta)[Gen]
  }

  ModelLevelIndices <- c(ECV = ECV, PUC = PUC(Lambda), Omega = Omega, OmegaH  = OmegaH, ARPB = ARPB_indices[[1]])

  ## Now put them all together
  indicesList <- list(FactorLevelIndices = FactorLevelIndices,
                      ItemLevelIndices   = ItemLevelIndices,
                      ModelLevelIndices  = ModelLevelIndices)
  ## if any index type is emtirely missing, remove that index type entirely (e.g., no model or item level indices if not bifactor)
  indicesList[which(!sapply(indicesList, is.null))]

}

#' bifactorIndicesMplus
#'
#' Computes all available bifactor indices given an \code{Mplus} .out file for a bifactor model
#'
#' @param Lambda is an Mplus .out file. Defaults to an open file dialog box
#' @param UniLambda is an object that the function can convert to a matrix of factor loadings.
#' The expected behavior is to store an Mplus output file as a variable and pass that variable
#' as \code{UniLambda}. Defaults to \code{NULL}, as \code{UniLambda} is only required if you wish to
#' compute \code{\link{ARPB}}.
#' @param standardized lets the function know whether it should be looking in
#'   the unstandardized results or the STDYX results from the Mplus output.
#'
#' @return A list of bifactor indices, including three different ECV indices, IECV, PUC, Omega,
#' OmegaH, and ARPB. Please note that many of these indices are interpretable even when the
#' model being used is not a bifactor model; some indices may be useful for two-tier, trifactor,
#' correlated traits, and even unidimensional models.
#'
#' @details ARPB will only be compute if the factor loadings from a unidimensional model
#' (as a vector or as the result of using \code{\link[MplusAutomation]{readModels}} on an
#' \code{Mplus} .out file) are included. Note that if a correlated traits model is provided,
#' the omega indices will simply be the regular omega values for those factors. Interpretations
#' for individual indices as well as details about their computation can be found in the
#' man page for the individual indices.
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


bifactorIndicesMplus <- function(Lambda = file.choose(), UniLambda = NULL, standardized = TRUE) {
  ## Expectation is that UniLambda is either a .out file, which will have a class of "character"
  if ("character" %in% class(UniLambda)) {UniLambda <- MplusAutomation::readModels(UniLambda)}

  if (!("mplus.model" %in% class(Lambda))) {Lambda <- MplusAutomation::readModels(Lambda)}
  ## if categorical, then error if standardized = FALSE and manually compute Theta if standardized - TRUE
  categorical <- !is.null(Lambda$input$variable$categorical)
  if (categorical) {
    if (standardized) {
      Lambda <- getLambda(Lambda, standardized = standardized)
      Theta <- getTheta(Lambda, standardized = standardized)
    } else {
      stop("Bifactor indices based on unstandardized coefficients with categorical variables is not available")
    }
  } else {
    Theta <- getTheta(Lambda, standardized = standardized)
    Lambda <- getLambda(Lambda, standardized = standardized)
  }

  bifactorIndices(Lambda, Theta, UniLambda, standardized)
}

