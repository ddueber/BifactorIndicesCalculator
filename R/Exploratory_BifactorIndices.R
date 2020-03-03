#' bifactorIndices_EFA
#'
#' Computes all available bifactor indices for the input given.
#'
#' @param Lambda is a factor loading matrix from EFA or an object which can be converted to such.
#' Currently only \code{psych} objects are supported.
#'
#' @param ItemsBySF is a list, indexed by factor, of vectors of item names belonging to each
#' factor. You must include the general factor in this list, and the list must have names which
#' match the factor names in \code{Lambda}. It is recommended you look at the EFA solution first
#' to see which factor is which. Defaults to \code{NULL}, in which case composition of specific
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
#' @details Only standardized models are considered for exploratory models. PUC and ARPB are not
#' supported for exploratory models currently, although that may change.
#'
#' @seealso \code{\link{bifactorIndices}},
#'          \code{\link{bifactorIndicesMplus}},
#'          \code{\link{bifactorIndicesMplus_EFA}},
#'          \code{\link{bifactorIndices_EFA}},
#'          \code{\link{ECV_SS}},
#'          \code{\link{ECV_SG}},
#'          \code{\link{ECV_GS}},
#'          \code{\link{Omega_S}},
#'          \code{\link{Omega_H}},
#'
#'
#' @export
#'
bifactorIndices_EFA <- function(Lambda, ItemsBySF = NULL, LoadMin = 0.2) {
  ## I'll make this into S3 methods once MplusAutomation supports EFA
  ## This is the method for pscyh::fa
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
    print("This matrix describes assignemnt of items to factors")
    print(ifelse(SmallLambda == 0, "", SmallLambda), quote = FALSE)
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

  # Is the first factor a true general factor?
  hasGen <- length(ItemsBySF[[1]]) == length(Items)

  # Issue a warning if no true gneral factor
  if (!hasGen) warning("The exploratory model has no general factor")

  ## Some of the indices we want involve all items
  GlobalIndices <- bifactorIndices(Lambda)

  ## For specific factor indices, we only use the items on the specific factor
  SpecificIndicesList <- lapply(Factors, function (Fac) {
    bifactorIndices(Lambda[ItemsBySF[[Fac]],])
  })

  SpecificIndices <- as.data.frame(t(sapply(Factors, function (Fac) {
    SpecificIndicesList[[Fac]]$FactorLevelIndices[Fac,]
  })))

  if (hasGen) {
    # ECV_SG taken from version with all items
    SpecificIndices$ECV_SG <- GlobalIndices$FactorLevelIndices$ECV_SS
    # ECV_GS is the general factor's ECV_SS when only items on the specific are included
    SpecificIndices$ECV_GS <- sapply(Factors, function (Fac) {
      SpecificIndicesList[[Fac]]$FactorLevelIndices[1,"ECV_SS"]
    })
    # Reorder the columns
    SpecificIndices <- SpecificIndices[,c("ECV_SS", "ECV_SG", "ECV_GS", "Omega", "Omega_H")]
    return(list(FactorLevelIndices = SpecificIndices,
                ModelLevelIndices = GlobalIndices[["FactorLevelIndices"]][1,]))
  } else {
    return(SpecificIndices)
  }

}


#' bifactorIndicesMplus_EFA
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
#' @details EFA models are nt currently (3/3/2020) supported by \code{MplsuAutomation::ReadModels()},
#' but they will be in the very near future, at which time this function will be completed.
#'
#' @seealso \code{\link{bifactorIndices}},
#'          \code{\link{bifactorIndicesMplus}},
#'          \code{\link{bifactorIndicesMplus_EFA}},
#'          \code{\link{bifactorIndices_EFA}},
#'          \code{\link{ECV_SS}},
#'          \code{\link{ECV_SG}},
#'          \code{\link{ECV_GS}},
#'          \code{\link{Omega_S}},
#'          \code{\link{Omega_H}},
#'
#'
#' @export
#'
bifactorIndicesMplus_EFA <- function(Lambda, ItemsBySF = NULL, LoadMin = 0.2) {
  stop("MplusAutomation does not support EFA output yet, but should within a month or so.")
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
#' @details Only standardized models are considered for exploratory models. PUC and ARPB are not
#' supported for exploratory models currently, although that may change.
#'
#' @seealso \code{\link{bifactorIndices}},
#'          \code{\link{bifactorIndicesMplus}},
#'          \code{\link{bifactorIndicesMplus_EFA}},
#'          \code{\link{bifactorIndices_EFA}},
#'          \code{\link{ECV_SS}},
#'          \code{\link{ECV_SG}},
#'          \code{\link{ECV_GS}},
#'          \code{\link{Omega_S}},
#'          \code{\link{Omega_H}},
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

  bifactorIndices_EFA(Lambda, ItemsBySF = NULL, LoadMin = 0.2)
}
