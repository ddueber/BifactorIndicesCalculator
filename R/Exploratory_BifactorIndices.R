#' bifactorIndices_expl
#'
#' Computes all available bifactor indices for the input given.
#'
#' @param Lambda is a factor loading matrix from EFA or an object which can be converted to such.
#' Currently only \code{psych::fa()} objects are supported.
#'
#' @param ItemsBySF is a list, indexed by factor, of vectors of item names belonging to each
#' factor. You must include the general factor in this list, and the list must have names which
#' match the factor names in \code{Lambda}. It is recommended you look at the EFA solution first
#' to see which factor is which. Defaults to \code{NULL}, in which case composition of specific
#' factors is automated by comparing loadings to \code{LoadMin}
#'
#' @param LoadMin is the minimum loading size so that an item is considered to "belong" to a factor.
#' If \code{ItemsBySF} is not provided, then items are assigned to factors based on whether their
#' loading on that factor is greater than \code{LoadMin}. If \code{ItemsBySF} is provided, then
#' warnings are issued whenever items load above \code{LoadMin} on factors to which they do not belong,
#' or do not load above \code{LoadMin} on factors to which they do belong, \code{LoadMin} defaults to 0.2.
#'
#' @return A list of bifactor indices, including three different ECV indices, Omega, and
#' OmegaH.
#'
#' @details Only standardized models are considered for exploratory models. PUC and ARPB are not
#' supported for exploratory models currently, although that may change.
#'
#' @seealso \code{\link{bifactorIndices}},
#'          \code{\link{bifactorIndicesMplus}},
#'          \code{\link{bifactorIndicesMplus_expl}},
#'          \code{\link{bifactorIndicesMplus_ESEM}},
#'          \code{\link{ECV_SS}},
#'          \code{\link{ECV_SG}},
#'          \code{\link{ECV_GS}},
#'          \code{\link{IECV}},
#'          \code{\link{Omega_S}},
#'          \code{\link{Omega_H}},
#'          \code{\link{H}},
#'          \code{\link{FD}}
#'
#'
#' @export
#'
#' @examples
#'
#'# psych::fa() can not access the rotations We have to load the library.
#'library(psych)
#'SRS_BEFA <- fa(SRS_data, nfactors = 5, rotate = "bifactor")
#'
#'# inspect the solution to see which exploratory factors belong to which subdomain
#'SRS_BEFA$loadings
#'ItemsBySF = list(MR4 = paste0("SRS_", c(5, 9, 12, 15, 18)),
#'                 MR2 = paste0("SRS_", c(1, 2, 8, 11, 17)),
#'                 MR3 = paste0("SRS_", c(4, 6, 10, 14, 19)),
#'                 MR5 = paste0("SRS_", c(3, 7, 13, 16, 20)))
#'
#'bifactorIndices_expl(SRS_BEFA, ItemsBySF = ItemsBySF)

bifactorIndices_expl <- function(Lambda, ItemsBySF = NULL, LoadMin = 0.2) {
  ## I'll make this into S3 methods once MplusAutomation supports EFA
  ## This is the method for pscyh::fa
  ## Actually, since I use a separate function for Mplus, that's not a problem
  ## Leaving this as is until I have a reason not to.
  getLambdaExploratory <- function (Lambda) {
    Lambda <- Lambda$loadings
    class(Lambda) <- "matrix"
    Lambda
  }
  if ("psych" %in% class(Lambda)) Lambda  <- getLambdaExploratory(Lambda)

  Items   <- rownames(Lambda)
  names(Items) <- Items
  Factors <- colnames(Lambda)
  names(Factors) <- Factors

  if (is.null(ItemsBySF)) {
    ItemsBySF <- lapply(Factors, function (Fac) {
      Items[Lambda[,Fac] > LoadMin]
    })
    names(ItemsBySF) <- Factors
    SmallLambda <- round(Lambda, 3)
    SmallLambda[SmallLambda < LoadMin] <- 0
    cat("This matrix describes assignemnt of items to factors \n")
    print(ifelse(SmallLambda == 0, "", SmallLambda), quote = FALSE)
    cat("\n \n")
  } else { # issue a warning for each loading above LoadMin on the wrong factor or loading below LoadMin on the right factor
    for (I in Items) {
      for (Fac in Factors) {
        if (!(I %in% ItemsBySF[[Fac]]) & (Lambda[I,Fac] > LoadMin)) {
          warning(paste0("Item ", I, " loads on factor ", Fac, "above ", LoadMin))
        }
        if ((I %in% ItemsBySF[[Fac]]) & (Lambda[I,Fac] < LoadMin)) {
          warning(paste0("Item ", I, " loads on factor ", Fac, "below ", LoadMin))
        }
      }
    }
  }

  # Is there single factor that pervades all items
  FactorLengths <- sapply(ItemsBySF, length)

  # Issue a warning if no true gneral factor
  if (max(FactorLengths) != nrow(Lambda)) warning("The exploratory model has no general factor")

  ## Some of the indices we want involve all items
  GlobalIndices <- bifactorIndices(Lambda)

  ## For specific factor indices, we only use the items on the specific factor
  SpecificIndicesList <- lapply(Factors, function (Fac) {
    bifactorIndices(Lambda[ItemsBySF[[Fac]],])
  })

  SpecificIndices <- as.data.frame(t(sapply(Factors, function (Fac) {
    SpecificIndicesList[[Fac]]$FactorLevelIndices[Fac,]
  })))

  if (max(FactorLengths) == nrow(Lambda)) {
    ModelIndices <- GlobalIndices[["FactorLevelIndices"]][1,]
    names(ModelIndices) <- c("ECV", "Omega", "OmegaH")

    # ECV_SG taken from version with all items
    SpecificIndices$ECV_SG <- GlobalIndices$FactorLevelIndices$ECV_SS
    # ECV_GS is the general factor's ECV_SS when only items on the specific are included
    SpecificIndices$ECV_GS <- sapply(Factors, function (Fac) {
      SpecificIndicesList[[Fac]]$FactorLevelIndices[1,"ECV_SS"]
    })
    # Reorder the columns
    SpecificIndices <- SpecificIndices[,c("ECV_SS", "ECV_SG", "ECV_GS", "Omega", "Omega_H")]

    # If only one factor is general, then we can do I-ECV
    if (sum(FactorLengths == nrow(Lambda)) == 1) {
      # The I-ECV function cannot be used because there is no "true" general factor
      GenFac <- which(FactorLengths == nrow(Lambda))
      L2 <- Lambda^2
      IECV <- L2[,GenFac] / rowSums(L2)
    }

    return(list(ModelLevelIndices  = ModelIndices,
                FactorLevelIndices = SpecificIndices,
                ItemLevelIndices   = IECV)
                )
  } else {
    return(SpecificIndices)
  }

}


#' bifactorIndicesMplus_expl
#'
#' Computes all available bifactor indices given an \code{Mplus} .out file for a bifactor EFA
#'
#' @param Lambda is an Mplus .out file. Defaults to an open file dialog box
#'
#' @param ItemsBySF is a list, indexed by factor, of vectors of item names belonging to each
#' factor. You must include the general factor in this list, and the list must have names which
#' match the factor names in Mplus. Defaults to \code{NULL}, in which case composition of specific
#' factors in automated by comparing loadings to \code{LoadMin}
#'
#' @param LoadMin is the minimum loading size so that an item is considered to "belong" to a factor.
#' If \code{ItemsBySF} is not provided, then items are assigned to factors based on whether their
#' loading on that factor is greater than \code{LoadMin}. If \code{ItemsBySF} is provided, then
#' warnings are issued whenever items load above \code{LoadMin} on factors to which they do not belong,
#' or do not load above \code{LoadMin} on factors to which they do belong,
#'
#' @return A list of bifactor indices, including three different ECV indices, Omega, and
#' OmegaH.
#'
#' @details To use this function, simply call it without any arguments and a dialog box
#' will pop up for you to select a .out file of an exploratory bifactor model.
#'
#' EFA models are not currently (3/3/2020) supported by \code{MplsuAutomation::ReadModels()},
#' but they will be in the very near future, at which time this function will be completed.
#'
#' @seealso \code{\link{bifactorIndices}},
#'          \code{\link{bifactorIndicesMplus}},
#'          \code{\link{bifactorIndicesMplus_ESEM}},
#'          \code{\link{bifactorIndices_expl}}
#'
#'
#' @export
#'
bifactorIndicesMplus_expl <- function(Lambda = file.choose(), ItemsBySF = NULL, LoadMin = 0.2) {
  ## If Lambda hasn't been put through MplusAutomation::readModels, then we need to do that
  if (!("mplus.model" %in% class(Lambda))) {Lambda <- MplusAutomation::readModels(Lambda)}

  stop("MplusAutomation does not support EFA output yet, but should soon!")
}


#' bifactorIndicesMplus_ESEM
#'
#' Computes all available bifactor indices given an \code{Mplus} .out file for a bifactor ESEM
#'
#' @param Lambda is an Mplus .out file. Defaults to an open file dialog box
#'
#' @param ItemsBySF is a list, indexed by factor, of vectors of item names belonging to each
#' factor. You must include the general factor in this list, and the list must have names which
#' match the factor names in Mplus. Defaults to \code{NULL}, in which case composition of specific
#' factors in automated by comparing loadings to \code{LoadMin}
#'
#' @param LoadMin is the minimum loading size so that an item is considered to "belong" to a factor.
#' If \code{ItemsBySF} is not provided, then items are assigned to factors based on whether their
#' loading on that factor is greater than \code{LoadMin}. If \code{ItemsBySF} is provided, then
#' warnings are issued whenever items load above \code{LoadMin} on factors to which they do not belong,
#' or do not load above \code{LoadMin} on factors to which they do belong,
#'
#' @return A list of bifactor indices, including three different ECV indices, Omega, and
#' OmegaH.
#'
#' @details To use this function, simply call it without any arguments and a dialog box
#' will pop up for you to select a .out file for an ESEM model.
#'
#' Only standardized models are considered for exploratory models. PUC and ARPB are not
#' supported for exploratory models currently, although that may change.
#'
#' @seealso \code{\link{bifactorIndices}},
#'          \code{\link{bifactorIndicesMplus}},
#'          \code{\link{bifactorIndicesMplus_expl}},
#'          \code{\link{bifactorIndices_expl}}
#'
#'
#' @export
#'
bifactorIndicesMplus_ESEM <- function(Lambda = file.choose(),
                                      ItemsBySF = NULL,
                                      LoadMin = 0.2) {

  ## If Lambda hasn't been put through MplusAutomation::readModels, then we need to do that
  if (!("mplus.model" %in% class(Lambda))) {Lambda <- MplusAutomation::readModels(Lambda)}

  ## Now we need to fish out the factor loading matrix
  Lambda <- getLambda(Lambda)

  bifactorIndices_expl(Lambda, ItemsBySF = NULL, LoadMin = 0.2)
}
