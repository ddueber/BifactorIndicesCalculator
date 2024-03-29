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
#' compute \code{\link{ARPB}}.
#' @param standardized lets the function know whether to look for standardized or
#' unstandardized results from \pkg{lavaan} and defaults to \code{TRUE}. If \code{Lambda} is not a
#'  \pkg{lavaan} object, then \code{standardized} will be ignored.
#' @param Phi is the correlation matrix of factors and defaults to \code{NULL}. User should generally ignore this
#' parameter. If not provided, \code{bifactorIndices} will try to determine \code{Phi} from \code{Lambda} when \code{Lambda}
#' is a fitted lavaan model or will assume it is the identity matrix otherwise.
#' @param Thresh is a list of vectors of item thresholds, used only when items are categorical.\code{bifactorIndices}
#'  will try to determine \code{Thresh} from \code{Lambda} when \code{Lambda}
#' is a fitted lavaan model and the indicators are categorical.
#' \code{Thresh} defaults to null, which indicates items are continuous.
#'
#' @return A list of bifactor indices, including three different ECV indices, IECV, PUC,
#' Omega, OmegaH, Factor Determinacy (FD), Construct Replicability (H) and ARPB.
#' Please note that many of these indices are interpretable even
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
#' Formulas for all indices can be found in Rodriguez et al. (2016). When indicators are categorical,
#' the methodology of Green and Yang (2009) is used for computing Omega and OmegaH.
#'
#' @references
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
#' structural equation modeling: An alternative to coefficient alpha.
#' \emph{Psychometrika, 74}(1), 155-167 \doi{10.1007/s11336-008-9099-3}.
#'
#' Kamata, A., & Bauer, D. J. (2008). A note on the relation between factor analytic and item
#' response theory models. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 15}
#' (1), 136-153.
#'
#' #' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models:
#' calculating and interpreting statistical indices. \emph{Psychological Methods, 21}(2),
#' 137 \doi{10.1037/met0000045}.
#'

#' @export
#'
#' @seealso
#'          \code{\link{bifactorIndicesMplus}},
#'          \code{\link{bifactorIndices_expl}},
#'          \code{\link{bifactorIndicesMplus_expl}},
#'          \code{\link{bifactorIndicesMplus_ESEM}},
#'          \code{\link{ECV_SS}},
#'          \code{\link{ECV_SG}},
#'          \code{\link{ECV_GS}},
#'          \code{\link{IECV}},
#'          \code{\link{PUC}},
#'          \code{\link{Omega_S}},
#'          \code{\link{Omega_H}},
#'          \code{\link{cat_Omega_S}},
#'          \code{\link{cat_Omega_H}},
#'          \code{\link{H}},
#'          \code{\link{FD}},
#'          \code{\link{ARPB}}
#'
#' @examples
#'
#' # Computing bifactor indices from fitted lavaan object
#' # (using mirt object is similar). Use of the unidimensional
#' # model is optional; it is only used to compute ARPB.
#'
#'\donttest{
#'SRS_UnidimensionalModel <-
#'   "SRS =~ SRS_1  + SRS_2  + SRS_3  + SRS_4  + SRS_5  +
#'           SRS_6  + SRS_7  + SRS_8  + SRS_9  + SRS_10 +
#'           SRS_11 + SRS_12 + SRS_13 + SRS_14 + SRS_15 +
#'           SRS_16 + SRS_17 + SRS_18 + SRS_19 + SRS_20"
#'
#'SRS_Unidimensional <- lavaan::cfa(SRS_UnidimensionalModel,
#'                                  SRS_data,
#'                                  ordered = paste0("SRS_", 1:20),
#'                                  orthogonal = TRUE)
#'
#'
#' SRS_BifactorModel <-
#' "SRS =~ SRS_1  + SRS_2  + SRS_3  + SRS_4  + SRS_5  +
#'         SRS_6  + SRS_7  + SRS_8  + SRS_9  + SRS_10 +
#'         SRS_11 + SRS_12 + SRS_13 + SRS_14 + SRS_15 +
#'         SRS_16 + SRS_17 + SRS_18 + SRS_19 + SRS_20
#'  Function     =~ SRS_5  + SRS_9  + SRS_12 + SRS_15 + SRS_18
#'  Pain         =~ SRS_1  + SRS_2  + SRS_8  + SRS_11 + SRS_17
#'  SelfImage    =~ SRS_4  + SRS_6  + SRS_10 + SRS_14 + SRS_19
#'  MentalHealth =~ SRS_3  + SRS_7  + SRS_13 + SRS_16 + SRS_20"
#'
#' SRS_bifactor <- lavaan::cfa(SRS_BifactorModel,
#'                             SRS_data,
#'                             ordered = paste0("SRS_", 1:20),
#'                             orthogonal = TRUE)
#'
#' bifactorIndices(SRS_bifactor, UniLambda = SRS_Unidimensional)
#'}
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
#' bifactorIndices(Lambda)
#'
#'
#' # bifactorIndices can also be used on two-tier models
#' MTMM_model <- "
#' Trait1  =~ T1M1_1 + T1M1_2 + T1M1_3 +
#'            T1M2_1 + T1M2_2 + T1M2_3 +
#'            T1M3_1 + T1M3_2 + T1M3_3
#' Trait2  =~ T2M1_1 + T2M1_2 + T2M1_3 +
#'            T2M2_1 + T2M2_2 + T2M2_3 +
#'            T2M3_1 + T2M3_2 + T2M3_3
#' Trait3  =~ T3M1_1 + T3M1_2 + T3M1_3 +
#'            T3M2_1 + T3M2_2 + T3M2_3 +
#'            T3M3_1 + T3M3_2 + T3M3_3
#'
#' Method1  =~ T1M1_1 + T1M1_2 + T1M1_3 +
#'             T2M1_1 + T2M1_2 + T2M1_3 +
#'             T3M1_1 + T3M1_2 + T3M1_3
#' Method2  =~ T1M2_1 + T1M2_2 + T1M2_3 +
#'             T2M2_1 + T2M2_2 + T2M2_3 +
#'             T3M2_1 + T3M2_2 + T3M2_3
#' Method3  =~ T1M3_1 + T1M3_2 + T1M3_3 +
#'             T2M3_1 + T2M3_2 + T2M3_3 +
#'             T3M3_1 + T3M3_2 + T3M3_3
#'
#' Trait1 ~~ 0*Method1
#' Trait1 ~~ 0*Method2
#' Trait1 ~~ 0*Method3
#' Trait2 ~~ 0*Method1
#' Trait2 ~~ 0*Method2
#' Trait2 ~~ 0*Method3
#' Trait3 ~~ 0*Method1
#' Trait3 ~~ 0*Method2
#' Trait3 ~~ 0*Method3
#'
#' Method1 ~~ 0*Method2
#' Method1 ~~ 0*Method3
#' Method2 ~~ 0*Method3"
#'
#' MTMM_fit <- lavaan::cfa(MTMM_model, MTMM_data)
#' bifactorIndices(MTMM_fit)
#'

bifactorIndices <- function(Lambda, Theta = NULL, UniLambda = NULL, standardized = TRUE, Phi = NULL, Thresh = NULL) {
  ## If categorical and lavaan, then force standardized = TRUE and message about it
  if (!standardized & "lavaan" %in% class(Lambda) ) {
    if (length(lavaan::lavInspect(Lambda, "ordered")) > 0) {
      standardized <- TRUE
      message("Only bifactor indices based on standardized coefficients make sense for categorical indicators.")
    }
  }


  ## if fitted mirt object, then throw warning about Omegas probably being meaningless
  if ("SingleGroupClass" %in% class(Lambda)) {
    message("Interpreting omega indices for IRT models is not recommended at this time")
  }

  ## Make Lambda, Theta, Phi, and UniLambda matrices. Do Theta and Phi first because
  ## they need the Lambda object, not the Lambda matrix
  if (is.null(Theta)) {Theta <- getTheta(Lambda, standardized = standardized)}

  # Can do Phi for SingleGroupClass and for lavaan
  if (is.null(Phi)) {
    if ("SingleGroupClass" %in% class(Lambda)) {
      Phi <- mirt::summary(Lambda)$fcor
    } else if ("lavaan" %in% class(Lambda)) {
      Phi <- lavaan::lavInspect(Lambda, "std.lv")$psi
      # I hate that dumb symmetric matrix print method
      class(Phi) <- "matrix"
    } else {
      message("Latent variable covariance matrix is assumed to be the identity. This influences
            FD and, in the case of categorical models, the omega indices.")
      Phi <- diag(nrow = ncol(Lambda))  ## we can only get here if Lambda was provided as a matrix
    }
  }

  # Can do Thresh for lavaan. Maybe we should do it for mirt as well.
  if (is.null(Thresh) & ("lavaan" %in% class(Lambda))) {
    # Check to see if items are ordered; if so, rip out the thresholds
    if (length(lavaan::lavInspect(Lambda, "ordered")) > 0) {
      thresh_long <- lavaan::lavInspect(Lambda, "std")$tau  ## if categorical, then standardized.
      rownames(thresh_long) <- sapply(strsplit(rownames(thresh_long), "[|]"), "[[", 1)
      items <- unique(rownames(thresh_long))
      Thresh <- lapply(items, function (i) {
        thresh_long[rownames(thresh_long) == i]
      })
    }
  }


  Lambda <- getLambda(Lambda, standardized = standardized)


  if (!is.null(UniLambda)) {UniLambda <- getLambda(UniLambda, standardized = standardized)}

  ## Build up the lists of indices. FactorLevelIndices first

  FactorLevelIndices = list(ECV_SS  = ECV_SS(Lambda),
                            ECV_SG  = ECV_SG(Lambda),
                            ECV_GS  = ECV_GS(Lambda))
  if (is.null(Thresh)) {
    FactorLevelIndices[["Omega"]]  <- Omega_S(Lambda, Theta)
    FactorLevelIndices[["OmegaH"]] <- Omega_H(Lambda, Theta)
  } else {
    FactorLevelIndices[["Omega"]]  <- cat_Omega_S(Lambda, Thresh)
    FactorLevelIndices[["OmegaH"]] <- cat_Omega_H(Lambda, Thresh)
  }

  if (standardized) {
    FactorLevelIndices[["H"]] <- H(Lambda)
    FactorLevelIndices[["FD"]] <- FD(Lambda, Phi)
  } else {
    message("H and FD are currently only available when standardized = TRUE")
  }

  ## Remove any NULL values and convert to dataframe
  FactorLevelIndices <- FactorLevelIndices[which(!sapply(FactorLevelIndices, is.null))]
  FactorLevelIndices <- as.data.frame(FactorLevelIndices)


  ## Item level indices next. Figure out label on ARPB later
  ARPB_indices <- ARPB(Lambda, UniLambda)
  ItemLevelIndices <- list(IECV             = IECV(Lambda),
                           RelParameterBias = ARPB_indices[[2]])

  ## Remove any NULL values and convert to dataframe
  ItemLevelIndices <- ItemLevelIndices[which(!sapply(ItemLevelIndices, is.null))]
  ItemLevelIndices <- as.data.frame(ItemLevelIndices)
  if (isTRUE(all.equal(dim(ItemLevelIndices), c(0, 0)))) {ItemLevelIndices <- NULL}
  if (!is.null(ItemLevelIndices) && ncol(ItemLevelIndices) == 2) {colnames(ItemLevelIndices)[2] <- "RelParBias"}

  ## Model level indices next
  if (is.null(getGen(Lambda))) { # No general factor
    ECV <- NULL
    Omega <- NULL
    OmegaH <- NULL
  } else {
    Gen <- getGen(Lambda)
    ECV <- ECV_SG(Lambda)[Gen]
    if (is.null(Thresh)) {
      Omega <- Omega_S(Lambda, Theta)[Gen]
      OmegaH <- Omega_H(Lambda, Theta)[Gen]
    } else {
      Omega <- cat_Omega_S(Lambda, Thresh)[Gen]
      OmegaH <- cat_Omega_H(Lambda, Thresh)[Gen]
    }
  }

  ModelLevelIndices <- c(ECV = ECV, PUC = PUC(Lambda), Omega = Omega, OmegaH  = OmegaH, ARPB = ARPB_indices[[1]])

  ## Now put them all together
  indicesList <- list(ModelLevelIndices  = ModelLevelIndices,
                      FactorLevelIndices = FactorLevelIndices,
                      ItemLevelIndices   = ItemLevelIndices
                      )
  ## if any index type is entirely missing, remove that index type entirely (e.g., no model or item level indices if not bifactor)
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
#' @details To use this function, simply call it without any arguments and a dialog box
#' will pop up for you to select a .out file of a confirmatory bifactor model.
#'
#' ARPB will only be computed if the factor loadings from a unidimensional model
#' (as a vector or as the result of using \code{\link[MplusAutomation]{readModels}} on an
#' \code{Mplus} .out file) are included. Note that if a correlated traits model is provided,
#' the omega indices will simply be the regular omega values for those factors. Interpretations
#' for individual indices as well as details about their computation can be found in the
#' man page for the individual indices.
#'
#' @seealso \code{\link{bifactorIndices}},
#'          \code{\link{bifactorIndices_expl}},
#'          \code{\link{bifactorIndicesMplus_expl}},
#'          \code{\link{bifactorIndicesMplus_ESEM}},
#'          \code{\link{ECV_SS}},
#'          \code{\link{ECV_SG}},
#'          \code{\link{ECV_GS}},
#'          \code{\link{IECV}},
#'          \code{\link{PUC}},
#'          \code{\link{Omega_S}},
#'          \code{\link{Omega_H}},
#'          \code{\link{H}},
#'          \code{\link{FD}},
#'          \code{\link{ARPB}}
#'
#' @export
#'
bifactorIndicesMplus <- function(Lambda = file.choose(), UniLambda = NULL, standardized = TRUE) {
  ## Expectation is that UniLambda is either a .out file, which will have a class of "character"
  if ("character" %in% class(UniLambda)) {UniLambda <- MplusAutomation::readModels(UniLambda)}

  if (!("mplus.model" %in% class(Lambda))) {Lambda <- MplusAutomation::readModels(Lambda)}

  ## Check if categorical indicators.
  categorical <- !is.null(Lambda$input$variable$categorical)
  ## categorical -> standardized
  if (!standardized & categorical) {
    standardized <- TRUE
    message("Only bifactor indices based on standardized coefficients make sense for categorical indicators.")
  }

  # Fetch Phi matrix. Why isn't there a getPhi function?
  if (standardized) {
    # need to have standardized parameters if standardized!
    if (is.null(Lambda$parameters$stdyx.standardized)) {
      stop("You must request STDYX output for computing bifactor indices based on standardized coefficients.")
    }
    params <- Lambda$parameters$stdyx.standardized
  } else {
    params <- Lambda$parameters$unstandardized
  }

  ## We need factor names to be in the same order as factor loading matrix
  facNames <- params[grep(".BY", params$paramHeader), "paramHeader"]
  facNames <- gsub(".BY", "", facNames, fixed = TRUE)
  facNames <- unique(facNames)
  facVar <- sapply(facNames, function (fac) {
    params[params$paramHeader == "Variances" & substring(params$param, 1, 8) == fac,"est"]
  })

  # grab factor correlations, make them more easily parsed, then grab them
  factorCorrs <- params[grep(".WITH", params$paramHeader), ]
  factorCorrs$paramHeader <- gsub(".WITH", "", factorCorrs$paramHeader, fixed = TRUE)
  Phi <- lapply(1:length(facNames), function (x) {
    fac1 <- facNames[x]
    sapply(1:length(facNames), function (y) {
      fac2 <- facNames[y]
      ## Factor variances are different
      if (x == y) {
        facVar[x]
      } else {
        ## Look for (x,y) and if that's not there look for (y,x)
        if (length(factorCorrs[substring(factorCorrs$paramHeader, 1, 8) == fac1 & substring(factorCorrs$param, 1, 8) == fac2, "est"]) == 1) {
          factorCorrs[substring(factorCorrs$paramHeader, 1, 8) == fac1 & substring(factorCorrs$param, 1, 8) == fac2, "est"]
        } else {
          factorCorrs[substring(factorCorrs$paramHeader, 1, 8) == fac2 & substring(factorCorrs$param, 1, 8) == fac1, "est"]
        }
      }
    })
  })
  Phi <- matrix(unlist(Phi), byrow=TRUE, nrow=length(Phi) )

  if (categorical) {
    Lambda <- getLambda(Lambda, standardized = standardized)
    Theta <- getTheta(Lambda, standardized = standardized)
    ## now get thresholds
    items <- rownames(Lambda)
    thresh_long <- params[params$paramHeader == "Thresholds",]
    thresh_long$itemName <- sapply(strsplit(thresh_long$param, "[$]"), "[[", 1)
    Thresh <- lapply(items, function (i) {
      thresh_long[thresh_long$itemName == i, "est"]
    })
  } else {
    Theta <- getTheta(Lambda, standardized = standardized)
    Lambda <- getLambda(Lambda, standardized = standardized)
    Thresh <- NULL
  }

  # if diag(Phi) is not all ones, then we are in trouble. Let's deal with that
  D <- diag(x = sqrt(diag(Phi)))
  Phi <- solve(D) %*% Phi %*% t(solve(D))
  Lambda <- Lambda %*% solve(D)
  colnames(Lambda) <- facNames



  bifactorIndices(Lambda, Theta, UniLambda, standardized, Phi, Thresh)
}

