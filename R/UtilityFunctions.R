#' getLambda
#'
#' getLambda computes or extracts a matrix of factor loadings given some input. Methods exist to support an input of
#' a dataframe, an mplus.model from MplusAutomation, a SingleGroupClass object from mirt, and a
#' lavaan object from lavaan. Please do not use tibbles, as they do not support row names, and it is
#' best if your items are given names.
#'
#' @param x an object to be converted into a factor loading matrix, or an object containing a fitted model from which
#' a factor loading matrix will be extracted. Supported classes are data.frame, matrix, mplus.model, lavaan, and SingleGroupClass
#'
#' @param ... can be used to specify whether a standardized or unstandardized factor loading matrix
#' should be returned. Only relevant for lavaan and mplus.model input. The standardized matrix for mplus.model
#' is taken from stdyx results.
#'
#' @return a matrix of factor loadings
#'
getLambda <- function(x, ...) {
  UseMethod("getLambda")
}

    getLambda.default <- function(x, ...) {
      x[is.na(x)] <- 0
      as.matrix(x)
    }

    getLambda.lavaan <- function(x, standardized = TRUE, ...) {
      if (standardized) {
        x <- lavaan::lavInspect(x, "std")$lambda
        x[is.na(x)] <- 0
        as.matrix(x)
      } else {
        x <- lavaan::lavInspect(x, "est")$lambda
        x[is.na(x)] <- 0
        as.matrix(x)
      }
    }

    getLambda.SingleGroupClass <- function(x, ...) {
      FitSum <- mirt::summary(x)
      x <- FitSum$rotF
      x[is.na(x)] <- 0
      as.matrix(x)
    }

    getLambda.mplus.model <- function(x, standardized = TRUE, ...) {
      if (standardized) {
        getLambda(x$parameters$stdyx.standardized)
      } else {
        getLambda(x$parameters$unstandardized)
      }
    }

    getLambda.mplus.params <- function(x, ...) {
      ## I am not proud of this function, but it works...
      ## This line throws warnings because not every row has a period. But, all the rows we care about *do* have a period. So, I am suppressing the warnings
      x <- suppressWarnings(tidyr::separate(x, col = paramHeader, into = c("Fac", "op"), sep = "\\."))
      loadings <- na.omit(x[x$op == "BY",])
      Facs <- unique(loadings$Fac)
      Items <- unique(loadings$param)
      Lambda <- matrix(ncol = length(Facs), nrow = length(Items))
      for (i in 1:length(Facs)) {
        for (j in 1:length(Items)) {
          if (length(loadings[loadings$Fac == Facs[i] & loadings$param == Items[j], "est"]) == 0) {
            Lambda[j,i] <- 0
          } else {
            Lambda[j,i] <- loadings[loadings$Fac == Facs[i] & loadings$param == Items[j], "est"]
          }
        }
      }
      rownames(Lambda) <- Items
      colnames(Lambda) <- Facs
      Lambda
    }



#' getTheta
#'
#' getTheta extracts or computes a vector of residual variance for items. If a
#' factor loading matrix is provided, then the vector of residual variances is
#' computed from that matrix if \code{standardized} is \code{TRUE} or an error
#' is thrown for unstandardized models.
#'
#' @param x an object that can be converted into a factor loading matrix, or an
#' object containing a fitted model from which a vector of residual variances
#' can be extracted. Supported classes are data.frame, matrix, mplus.model,
#' lavaan, and SingleGroupClass
#' @param ... can be used to specify whether a standardized or unstandardized factor loading matrix
#' should be returned. Only relevant for lavaan and mplus.model input. The standardized matrix for mplus.model
#' is taken from stdyx results.
#'
#' @return a vector of residual variances for items. If x is a fitted model, then
#' the residual variances are extracted from the fitted model. Lavaan, mirt
#' (SingleGroupClass), and Mplus (mplus.model) models are supported.
#' If Mplus does not report residual variances for categorical variables, then
#' factor loadings are used to compute the residual variance for standardized models
#' and an error is thrown for unstandardized models. In both cases, the user is
#' alerted that residual variances could not be found in the input and perhaps the
#' model should be rerun.
#'
#' @seealso \code{\link{getLambda}}
#'
#' @examples
#'
#' MAKE SOME EXAMPLES
#'
getTheta <- function(x, standardized = TRUE, ...) {
  UseMethod("getTheta")
}

getTheta.default <- function(x, standardized = TRUE) {
  if(!standardized) {
    stop("Not enough information is provided to compute indicator residual variances")
  } else {
    ## This is excessive. There's not way to get here unless x is a data frame or matrix. But, better safe than sorry
    Lambda <- getLambda(x)
    Ones <- rep(1, nrow(Lambda))
    Ones - rowSums(Lambda^2)    }
}

getTheta.SingleGroupClass <- function(x, ...) {
  FitSum <- summary(x)
  Theta <- 1 - FitSum$h2
  as.vector(Theta)
}

    getTheta.lavaan <- function(x, standardized = TRUE, ...) {
      if (standardized) {
        diag(lavInspect(x, "std")$theta)
      } else {
        diag(lavInspect(x, "est")$theta)
      }
    }

    getTheta.mplus.model <- function(x, ...) {
      pars <- x$parameters$unstandardized
      ## This line throws warning because not every row has a period. But, all the rows we care about *do* have a period. So, I am suppressing the warnings
      loadings <- suppressWarnings(tidyr::separate(pars, col = paramHeader, into = c("Fac", "op"), sep = "\\."))
      loadings_2 <- na.omit(loadings[loadings$op == "BY",])
      items <- unique(loadings_2$param)
      ## Item names are not preserved below.
      if (length(x$input$variable$categorical) == 0) {
        Theta <- c()
        thetaOutput <- pars[pars$paramHeader == "Residual.Variances",]
        for (i in 1:length(items)) {
          Theta <- c(Theta, thetaOutput[thetaOutput$param == items[i], "est"])
        }} else {
          Theta <- c()
          for (i in 1:length(items)) {
            Theta <- c(Theta, x$parameters$r2[x$parameters$r2$param == items[i], "resid_var"])
          }}
      names(Theta) <- items
      Theta
    }


#' getGen
#'
#' getGen detects whether or not a single factor loads on all items, and returns the column index of that factor if it exists.
#'
#' @param Lambda is a factor loading matrix or an object which can be converted to one using getLambda
#'
#' @return The index of the general factor, or \code{NULL} if there is no general factor
#'
#' @seealso isBifactor, getLambda
#'
  getGen <- function(Lambda) {
    if (!("matrix" %in% class(Lambda))) {Lambda <- getLambda(Lambda, standardized = standardized)}
    ## Make a matrix of logical vectors for non-zero elements of Lambda. Let's replace NA with zero at the start!!
    inFactorMat <- Lambda != 0
    ## Now, compute the column sums. [[ If colSum is nrow, then we have a general factor. We cannot have more than one ]]
    itemsOnFactor <- colSums(inFactorMat == TRUE)
    if (sum(itemsOnFactor == nrow(Lambda)) == 1 ) {
      which(sum(itemsOnFactor == nrow(Lambda)) == 1)
    } else {
      NULL
    }
  }


#' isBifactor
#'
#' Determines whether a model has bifactor structure
#'
#' @param Lambda Matrix of factor loadings
#'
#' @return Logical. If each item loads on a general factor and at most one specific factor, returns TRUE. Otherwise FALSE.
#' @export
#'
  isBifactor <- function(Lambda) {
    if (!("matrix" %in% class(Lambda))) {Lambda <- getLambda(Lambda, standardized = standardized)}
    if (is.null(getGen(Lambda))) return(FALSE)
    inFactorMat <- Lambda != 0
    return(sum(rowSums(inFactorMat) > 2) == 0)
  }
