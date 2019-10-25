#' getLambda
#'
#' getLambda computes or extracts a matrix of factor loadings given some input. Methods exist to
#' support an input of
#' a \code{dataframe}, an \code{mplus.model} from \pkg{MplusAutomation}, a \code{SingleGroupClass} object from \pkg{mirt}, and a
#' \code{lavaan} object from \pkg{lavaan}. Please do not use a \code{tibble}, as they do not support
#' row names, and it is best if your items are given names.
#'
#' @param x an object to be converted into a factor loading matrix, or an object containing a fitted
#'  model from which a factor loading matrix will be extracted. Supported classes are
#'  \code{data.frame}, \code{matrix}, \code{mplus.model}, \code{lavaan}, and \code{SingleGroupClass}.
#'
#' @param standardized can be used to specify whether a standardized or unstandardized factor
#' loading matrix should be returned. Only relevant for \code{lavaan} and \code{mplus.model} input. The
#' standardized matrix for \code{mplus.model} is taken from stdyx results.
#'
#' @return A matrix of factor loadings
#'
getLambda <- function(x, standardized = TRUE) {
  UseMethod("getLambda")
}

getLambda.default <- function(x, standardized = TRUE) {
  x[is.na(x)] <- 0
  as.matrix(x)
}

getLambda.lavaan <- function(x, standardized = TRUE) {
  if (standardized) {
    x <- lavaan::lavInspect(x, "std")$lambda
    x[is.na(x)] <- 0
    as.matrix(x)
  } else {
    ## Make sure all factors have a variance of one.
    x <- lavaan::lavInspect(x, "std.lv")$lambda
    x[is.na(x)] <- 0
    as.matrix(x)
  }
}

getLambda.SingleGroupClass <- function(x, standardized = TRUE) {
  ## the summary method for mirt likes to print to screen, so this next line very awkwardly suppresses that printing
  temp <- capture.output(FitSum <- mirt::summary(x))
  x <- FitSum$rotF
  x[is.na(x)] <- 0
  as.matrix(x)
}

getLambda.mplus.model <- function(x, standardized = TRUE) {
  if (standardized) {
    ## check to make sure standardized output was requested
    if (is.null(x$parameters$stdyx.standardized)) stop("You must request standardized output from Mplus when standardized = TRUE")

    getLambda(x$parameters$stdyx.standardized)
  } else {
    ## check to make sure factor variances are all one.
    if (!all(x$parameters$unstandardized[x$parameters$unstandardized$paramHeader == "Variances", "est"] == 1)) stop("All factor variances must be one when standardized = FALSE. Please respecify your model.")

    getLambda(x$parameters$unstandardized)
  }
}

getLambda.mplus.params <- function(x) {
  ## I am not proud of this function, but it works...
  ## This line throws warnings because not every row has a period. But, all the rows we care about *do* have a period. So, I am suppressing the warnings
  x <- suppressWarnings(tidyr::separate(x, col = "paramHeader", into = c("Fac", "op"), sep = "\\."))
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
#' \code{getTheta} extracts or computes a vector of residual variance for items. If a
#' factor loading matrix is provided, then the vector of residual variances is
#' computed from that matrix if \code{standardized} is \code{TRUE}.
#'
#' @param x an object that can be converted into a factor loading matrix, or an
#' object containing a fitted model from which a vector of residual variances
#' can be extracted. Supported classes are \code{data.frame}, \code{matrix}, \code{mplus.model},
#' \code{lavaan}, and \code{SingleGroupClass}
#' @param standardized can be used to specify whether a standardized or unstandardized factor
#' loading matrix should be returned. Only relevant for \code{lavaan} and \code{mplus.model}
#' input. The standardized matrix for \code{mplus.model} is taken from stdyx results.
#'
#' @return a vector of residual variances for items. If x is a fitted model, then
#' the residual variances are extracted from the fitted model. \pkg{lavaan}, \pkg{mirt}
#' (\code{SingleGroupClass}), and \code{Mplus} (\code{mplus.model}) models are supported.
#' If \code{Mplus} does not report residual variances for categorical variables, then
#' factor loadings are used to compute the residual variance for standardized models
#' and an error is thrown for unstandardized models. In both cases, the user is
#' alerted that residual variances could not be found in the input and perhaps the
#' model should be rerun.
#'
#' @seealso \code{\link{getLambda}}
#'
#'
getTheta <- function(x, standardized = TRUE) {
  UseMethod("getTheta")
}

getTheta.default <- function(x, standardized = TRUE) {
  if(!standardized) {
    stop("Not enough information is provided to compute indicator residual variances. Either provide indicator residual variances or use a standardized solution.")
  } else {
    ## This is excessive. There's no way to get here unless x is a data frame or matrix. But, better safe than sorry
    Lambda <- getLambda(x)
    Ones <- rep(1, nrow(Lambda))
    Ones - rowSums(Lambda^2)    }
}

getTheta.SingleGroupClass <- function(x, standardized = TRUE) {
  ## the summary method for mirt likes to print to screen, so this next line very awkwardly suppresses that printing
  temp <- capture.output(FitSum <- mirt::summary(x))
  Theta <- 1 - FitSum$h2
  as.vector(Theta)
}

getTheta.lavaan <- function(x, standardized = TRUE) {
  if (standardized) {
    diag(lavaan::lavInspect(x, "std")$theta)
  } else {
    diag(lavaan::lavInspect(x, "est")$theta)
  }
}

getTheta.mplus.model <- function(x, standardized = TRUE) {
  if (standardized) {
    if (is.null(x$parameters$stdyx.standardized)) stop("You must request standardized output from Mplus when standardized = TRUE")
    pars <- x$parameters$stdyx.standardized
  } else {
    pars <- x$parameters$unstandardized
  }
  ## This line throws warning because not every row has a period. But, all the rows we care about *do* have a period. So, I am suppressing the warnings
  loadings <- suppressWarnings(tidyr::separate(pars, col = "paramHeader", into = c("Fac", "op"), sep = "\\."))
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
#' \code{getGen} detects whether or not a single factor loads on all items, and returns the column
#' index of the general factor if it exists.
#'
#' @param Lambda is a factor loading matrix
#'
#' @return The index of the general factor, or \code{NULL} if there is no general factor
#'
#'
getGen <- function(Lambda) {
  ## Make a matrix of logical vectors for non-zero elements of Lambda. Let's replace NA with zero at the start!!
  inFactorMat <- Lambda != 0
  ## Now, compute the column sums. [[ If colSum is nrow, then we have a general factor. We cannot have more than one ]]
  itemsOnFactor <- colSums(inFactorMat == TRUE)
  if (length(which(itemsOnFactor == nrow(Lambda))) == 1 ) {
    which(itemsOnFactor == nrow(Lambda))
  } else {
    NULL
  }
}


#' isBifactor
#'
#' Determines whether a model has bifactor structure.
#'
#' @param Lambda Matrix of factor loadings
#'
#' @return Logical. If each item loads on a general factor and at most one specific factor, returns TRUE. Otherwise FALSE.
#'
  isBifactor <- function(Lambda) {
    if (is.null(getGen(Lambda))) return(FALSE)
    inFactorMat <- Lambda != 0
    return(sum(rowSums(inFactorMat) > 2) == 0)
  }
