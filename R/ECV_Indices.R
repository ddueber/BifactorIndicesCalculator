#' ECV_SS
#'
#' Computes an ECV index for all factors which can be interpreted as the proportion of common variance of the items in each factor which is due to that factor. Here, ECV is computed only with respect to items which load on the factor; Stucky and Edelen (2015, p. 201) refer to this index as "within-domain ECV". In the Excel version of the bifactor indices calculator (Dueber, 2017), this index is referred to as ECV (NEW). ECV_SS is useful in that it can be computed when there is no specific factor, such as in a two-tier model.
#'
#' \code{ECV_SS} is called by \code{\link{bifactorIndices}} and \code{\link{bifactorIndicesMPlus}}, which are the only functions in this package intended for casual users
#'
#' @param Lambda is a matrix of factor loadings. Be sure that all factors have the same variance before calling this function.
#'
#' @return A vector of ECVs for all factors
#'
#' @examples
#' Hey I need some examples here!!
#'
#' @section References:
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to calculate various indices relevant to bifactor CFA models. <doi:10.13023/edp.tool.01>
#' Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create unidimensional measures from multidimensional data. In S. P. Reise & D. A. Revicki (Eds.), Handbook of item response theory modeling: Applications to typical performance assessment (pp.183-206). New York: Routledge.
#'
#' @export
#'
#' @seealso \code{\link{ECV_SG}}, \code{\link{ECV_GS}}, \code{\link{bifactorIndices}}, \code{\link{bifactorIndicesMPlus}}
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
  sapply(colnames(Lambda), ECV_SS_C, Lambda, standardized = standardized)
}


#' ECV_SG
#'
#' Computes an ECV index for all factors which can be interpreted as the proportion of common variance of the items in each factor which is due to the general factor. Here, ECV is computed only with respect to all items using the specific factor loadings in the numerator. Neither Stucky and Edelen (2015) or the Excel version of the Bifactor Indices Calculator (Dueber, 2017) use this form of ECV.
#'
#' \code{ECV_SG} is called by \code{\link{bifactorIndices}} and \code{\link{bifactorIndicesMPlus}}, which are the only functions in this package intended for casual users
#'
#' @param Lambda is a matrix of factor loadings. Be sure that all factors have the same variance before calling this function.
#'
#' @return A vector of ECVs for all factors
#'
#' @examples
#' Hey I need some examples here!!
#'
#' @section References:
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to calculate various indices relevant to bifactor CFA models. <doi:10.13023/edp.tool.01>
#' Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create unidimensional measures from multidimensional data. In S. P. Reise & D. A. Revicki (Eds.), Handbook of item response theory modeling: Applications to typical performance assessment (pp.183-206). New York: Routledge.
#'
#' @export
#'
#' @seealso \code{\link{ECV_SS}}, \code{\link{ECV_GS}}, \code{\link{bifactorIndices}}, \code{\link{bifactorIndicesMPlus}}
#'

ECV_SG <- function(Lambda) {
  if (is.null(getGen(Lambda))) return(NULL)
  ECV_SG_C <- function(Fac, Lambda) {
    Lambda <- getLambda(Lambda, standardized = standardized)
    if (is.null(getGen(Lambda))) return(NULL)
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Square the loadings
    L2 <- Lambda^2
    ## Compute the appropriate ratio of sums
    sum(L2[,Fac]*inFactor)/sum(L2)
  }
  sapply(colnames(Lambda), ECV_SG_C, Lambda, standardized = standardized)
}

#' ECV_GS
#'
#' Computes an ECV index for all factors which can be interpreted as the proportion of common variance of the items in each factor which is due to the general factor. Here, ECV is computed only with respect to all items using the general factor loadings in the numerator; Stucky and Edelen (2015, p. 199) refer to this index as ECV for the specific factor. In the Excel version of the bifactor indices calculator (Dueber, 2017), this index is referred to as ECV (S&E).
#'
#' \code{ECV_GS} is called by \code{\link{bifactorIndices}} and \code{\link{bifactorIndicesMPlus}}, which are the only functions in this package intended for casual users
#'
#' @param Lambda is a matrix of factor loadings. Be sure that all factors have the same variance before calling this function.
#'
#' @return A vector of ECVs for all factors
#'
#' @examples
#' Hey I need some examples here!!
#'
#' @section References:
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to calculate various indices relevant to bifactor CFA models. <doi:10.13023/edp.tool.01>
#' Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create unidimensional measures from multidimensional data. In S. P. Reise & D. A. Revicki (Eds.), Handbook of item response theory modeling: Applications to typical performance assessment (pp.183-206). New York: Routledge.
#'
#' @export
#'
#' @seealso \code{\link{ECV_SS}}, \code{\link{ECV_SG}}, \code{\link{bifactorIndices}}, \code{\link{bifactorIndicesMPlus}}
#'

ECV_GS <- function(Lambda) {
  if (is.null(getGen(Lambda))) return(NULL)
  ECV_GS_C <- function(Fac, Lambda) {
    if (is.null(getGen(Lambda))) return(NULL)
    genFac <- getGen(Lambda)
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Square the loadings
    L2 <- Lambda^2
    ## Compute the appropriate ratio of sums
    sum(L2[,genFac]*inFactor)/sum(L2*inFactor)
  }
  sapply(colnames(Lambda), ECV_GS_C, Lambda)
}



#' IECV
#'
#' Computes an ECV index for each item which can be interpreted as the proportion of common variance of that item due to the general factor. Stucky and Edelen (2015, p. 201) define I-ECV, which is also computed in the Excel version of the bifactor indices calculator (Dueber, 2017).
#'
#' \code{IECV} is called by \code{\link{bifactorIndices}} and \code{\link{bifactorIndicesMPlus}}, which are the only functions in this package intended for casual users
#'
#' @param Lambda is a matrix of factor loadings. Be sure that all factors have the same variance before calling this function.
#'
#' @return A vector of ECVs for all factors
#'
#' @examples
#' Hey I need some examples here!!
#'
#' @section References:
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to calculate various indices relevant to bifactor CFA models. <doi:10.13023/edp.tool.01>
#' Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create unidimensional measures from multidimensional data. In S. P. Reise & D. A. Revicki (Eds.), Handbook of item response theory modeling: Applications to typical performance assessment (pp.183-206). New York: Routledge.
#'
#' @export
#'
#' @seealso \code{\link{ECV_SS}}, \code{\link{ECV_SG}}, \code{\link{bifactorIndices}}, \code{\link{bifactorIndicesMPlus}}
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
  sapply(1:nrow(Lambda), IECV_C, Lambda, genFac)
}
