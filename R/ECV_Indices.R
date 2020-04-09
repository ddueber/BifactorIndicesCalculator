#' ECV_SS
#'
#' Computes an ECV index for all factors which can be interpreted as the proportion of common
#' variance of the items in each factor which is due to that factor;
#' \code{ECV_SS} should be read 'ECV of a specific factor with respect to itself.' Here, ECV is computed
#' only with respect to items which load on the factor. Note that \code{ECV_SS} of the general factor
#' is simply the ECV. Stucky and Edelen (2015, p. 201) do not refer to this form of ECV. In the Excel
#' version of the bifactor indices calculator (Dueber, 2017), this index is referred to as
#' 'ECV (NEW).' \code{ECV_SS} is useful in that it can be computed when there is no general factor, such
#' as in a two-tier model, and interpreted in the same way as ECV for general factors.
#'
#' \code{ECV_SS} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.``
#'
#' @param Lambda is a matrix of factor loadings. Be sure that all factors have the same variance
#'  before calling this function.
#'
#' @return A vector of ECVs for all factors
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
#' ECV_SS(Lambda)
#'
#'
#' @references
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to calculate various indices relevant to bifactor CFA models. \doi{10.13023/edp.tool.01}
#'
#' Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create unidimensional measures from multidimensional data. In S. P. Reise & D. A. Revicki (Eds.), \emph{Handbook of item response theory modeling: Applications to typical performance assessment} (pp.183-206). New York: Routledge.
#'
#' @export
#'
#' @seealso \code{\link{ECV_SG}}, \code{\link{ECV_GS}}, \code{\link{bifactorIndices}}
#'
ECV_SS <- function(Lambda) {
  ECV_SS_C <- function(Fac, Lambda) {
    ## Isolate the items which load on the chosen factor FAC
    inFactor <- Lambda[,Fac] != 0
    ## Square the loadings
    L2 = Lambda^2
    ## Compute the appropriate ratio of sums
    sum(L2[,Fac]*inFactor)/sum(L2*inFactor)
  }
  ECV_results <- sapply(1:ncol(Lambda), ECV_SS_C, Lambda)
  names(ECV_results) <- colnames(Lambda)
  ECV_results
}


#' ECV_SG
#'
#' Computes an ECV index for all factors which can be interpreted as the proportion of
#' common variance of all items which is due to the specific factor;
#' \code{ECV_SG} should be read 'ECV of a specific factor with respect to the general
#' factor.' Here,
#' ECV is computed with respect to the items of the general factor using the specific factor loadings in
#' the numerator; Stucky and Edelen (2015, p. 199)
#' refer to this index simply as 'specific-dimension ECV.' Note that \code{ECV_SG} of the general factor
#' is simply the ECV. In the Excel version of the Bifactor
#' Indices Calculator (Dueber, 2017), this form of ECV is referred to as 'ECV (S&E).'
#
#' \code{ECV_SG} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of factor loadings. Be sure that all factors have the same
#' variance before calling this function.
#'
#' @return A vector of ECVs for all factors
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
#' ECV_GS(Lambda)
#'
#' @references
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to
#' calculate various indices relevant to bifactor CFA models. \doi{10.13023/edp.tool.01}
#'
#' Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create
#' unidimensional measures from multidimensional data. In S. P. Reise & D. A. Revicki
#' (Eds.), \emph{Handbook of item response theory modeling: Applications to typical
#' performance assessment} (pp.183-206). New York: Routledge.
#'
#' @export
#'
#' @seealso \code{\link{ECV_SS}}, \code{\link{ECV_GS}}, \code{\link{bifactorIndices}}
#'

ECV_SG <- function(Lambda) {
  if (is.null(getGen(Lambda))) return(NULL)
  ECV_SG_C <- function(Fac, Lambda) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Square the loadings
    L2 <- Lambda^2
    ## Compute the appropriate ratio of sums
    sum(L2[,Fac]*inFactor)/sum(L2)
  }
  ECV_results <- sapply(1:ncol(Lambda), ECV_SG_C, Lambda)
  names(ECV_results) <- colnames(Lambda)
  ECV_results
}



#' ECV_GS
#'
#' Computes an ECV index for all factors which can be interpreted as the proportion of common
#' variance of the items in each specific factor which is due to the general factor;
#' \code{ECV_GS} should be read 'ECV of the general factor with respect to a specific
#' factor.' Here, ECV is
#' computed only with respect to the items of a specific factor using the general factor
#' loadings in the numerator;
#' Stucky and Edelen (2015, p. 201) refer to this index as the 'within-domain ECV' for the
#' specific factor. In the
#' Excel version of the bifactor indices calculator (Dueber, 2017), this index is not computed.
#'
#' \code{ECV_GS} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of factor loadings. Be sure that all factors have the same variance
#' before calling this function.
#'
#' @return A vector of ECVs for all factors
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
#' ECV_GS(Lambda)
#'
#' @references
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to calculate
#' various indices relevant to bifactor CFA models. \doi{10.13023/edp.tool.01}
#'
#' Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create unidimensional
#' measures from multidimensional data. In S. P. Reise & D. A. Revicki (Eds.), \emph{Handbook of item
#' response theory modeling: Applications to typical performance assessment} (pp.183-206).
#' New York: Routledge.
#'
#' @export
#'
#' @seealso \code{\link{ECV_SS}}, \code{\link{ECV_SG}}, \code{\link{bifactorIndices}}
#'

ECV_GS <- function(Lambda) {

  ECV_GS_C <- function(Fac, Lambda, genFac) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Square the loadings
    L2 <- Lambda^2
    ## Compute the appropriate ratio of sums
    sum(L2[,genFac]*inFactor)/sum(L2*inFactor)
  }
  ## ECV_GS only makes sense when there is a general factor
  if (is.null(getGen(Lambda))) return(NULL)
  genFac <- getGen(Lambda)

  ECV_results <- sapply(1:ncol(Lambda), ECV_GS_C, Lambda, genFac)
  names(ECV_results) <- colnames(Lambda)
  ECV_results
}



#' IECV
#'
#' Computes an ECV index for each item which can be interpreted as the proportion of common
#' variance of that item due to the general factor. Stucky and Edelen (2015, p. 201) define
#' I-ECV, which is also computed in the Excel version of the bifactor indices calculator
#' (Dueber, 2017).
#'
#' \code{IECV} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of factor loadings. Be sure that all factors have the same variance
#' before calling this function.
#'
#' @return A vector of item ECVs
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
#' IECV(Lambda)
#'
#' @references
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to calculate
#' various indices relevant to bifactor CFA models. \doi{10.13023/edp.tool.01}
#'
#' Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create unidimensional
#' measures from multidimensional data. In S. P. Reise & D. A. Revicki (Eds.), \emph{Handbook of item
#' response theory modeling: Applications to typical performance assessment} (pp.183-206).
#' New York: Routledge.
#'
#' @export
#'
#' @seealso \code{\link{ECV_SS}}, \code{\link{ECV_SG}}, \code{\link{ECV_GS}}, \code{\link{bifactorIndices}}
#'

IECV <- function(Lambda) {

  ## Compute IECV for single item
  IECV_C <- function(Item, Lambda, genFac) {
    ## Square the loadings
    L2 <- Lambda^2
    IECV <- L2[Item, genFac]/sum(L2[Item,])
    names(IECV) <- rownames(Lambda)[Item]
    IECV
  }

  ## I-ECV only makes sense when there is a general factor
  if (is.null(getGen(Lambda))) return(NULL)
  genFac <- getGen(Lambda)
  IECV_results <- sapply(1:nrow(Lambda), IECV_C, Lambda, genFac)
  names(IECV_results) <- rownames(Lambda)
  IECV_results
}
