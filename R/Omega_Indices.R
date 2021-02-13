#' Omega_S
#'
#' Computes an omega reliability estimate for all factors as described in Rodriguez, Reise, and
#' Haviland (2016).
#'
#' \code{Omega_S} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of factor loadings
#' @param Theta is a vector of indicator error variances
#'
#' @return A \code{numeric}, the omega reliability estimate for all factors.
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
#' Theta <- rep(1, nrow(Lambda)) - rowSums(Lambda^2)
#' Omega_S(Lambda, Theta)
#'
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models:
#' calculating and interpreting statistical indices. \emph{Psychological Methods, 21}(2),
#' 137 \doi{10.1037/met0000045}.
#'
#' @export
#'
#' @seealso \code{\link{Omega_H}}, \code{\link{bifactorIndices}}
#'

Omega_S <- function(Lambda, Theta) {
  Omega_S_C <- function(Fac, Lambda, Theta) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Compute the appropriate ratio of sums
    sum(colSums(Lambda*inFactor)^2)/(sum(colSums(Lambda*inFactor)^2) + sum(Theta*inFactor))
  }
  if (is.null(Theta)) return(NULL)
  omega_results <- sapply(1:ncol(Lambda), Omega_S_C, Lambda = Lambda, Theta = Theta)
  names(omega_results) <- colnames(Lambda)
  omega_results
}


#' cat_Omega_S
#'
#' Computes an omega reliability estimate for all factors as described in Rodriguez, Reise, and
#' Haviland (2016).
#'
#' \code{cat_Omega_S} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of standardized factor loadings
#' @param Thresh is a list (indexed by items) of vectors of item thresholds
#'
#' @return A \code{numeric}, the omega reliability estimate for all factors using the technique of
#' Green and Yang (2009).
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
#' Thresh # Gotta define Thresh
#' catOmega_S(Lambda, Thresh)
#'
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models:
#' calculating and interpreting statistical indices. \emph{Psychological Methods, 21}(2),
#' 137 \doi{10.1037/met0000045}.
#'
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
#' structural equation modeling: An alternative to coefficient alpha.
#' \emph{Psychometrika, 74}(1), 155-167 \doi{10.1007/s11336-008-9099-3}.
#'
#' @export
#'
#' @seealso \code{\link{Omega_H}}, \code{\link{bifactorIndices}}
#'
#'

cat_Omega_S <- function(Lambda, Thresh, Phi = NULL, denom = NULL) {
  cat_Omega_S_C <- function(Fac, Lambda, Thresh, Phi = NULL, denom = NULL) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## subset Lambda to only include appropriate items
    subLambda <- Lambda[inFactor,]

    ## create Phi matrix
    if (is.null(Phi)) {Phi <- diag(ncol(subLambda))}

    num_items <- nrow(subLambda)

    ## latent item covariances (correlations, because they're standardized!!)
    lat_cov <- subLambda %*% Phi %*% t(subLambda)
    if (is.null(denom)) {
      poly_cor <- lat_cov
      diag(poly_cor) <- rep(1, num_items)
    }

    ## Create covariance matrix of parallel items
    par_cov_mat <- sapply(1:num_items, function (j) {
      t_j <- Thresh[[j]]
      sapply(1:num_items, function(jp) {
        t_jp <- Thresh[[jp]]
        ## cov(x_j, x_jp)

        ## first, the left half of the expression in equation 19 in Green and Yang (2009)
        left <- sum(sapply(1:length(t_j), function(c) {
          sapply(1:length(t_jp), function(cp) {
            r <- lat_cov[j, jp]
            mnormt::pmnorm(c(t_j[c], t_jp[cp]), c(0, 0), matrix(c(1, r, r, 1), 2))
          }) ## end cp
        })) ## end c

        ## now the two expression in the right half of the expression in equation 19 in G&Y
        right_j  <- sum(pnorm(t_j))
        right_jp <- sum(pnorm(t_jp))

        ## put it together and what do you get? Bibbidi-Bobbidi-Boo
        left - right_j * right_jp

      }) ## end jp
    }) ## end j

    ## create covariance matrix of items... copy/paste of par_cov_mat,
    ## but switch from lat_cov_mat to poly_cor
    item_cov_mat <- sapply(1:num_items, function (j) {
      t_j <- Thresh[[j]]
      sapply(1:num_items, function(jp) {
        t_jp <- Thresh[[jp]]
        ## cov(x_j, x_jp)

        ## first, the left half of the expression in equation 19 in Green and Yang (2009)
        left <- sum(sapply(1:length(t_j), function(c) {
          sapply(1:length(t_jp), function(cp) {
            r <- poly_cor[j, jp]
            mnormt::pmnorm(c(t_j[c], t_jp[cp]), c(0, 0), matrix(c(1, r, r, 1), 2))
          }) ## end cp
        })) ## end c

        ## now the two expression in the right half of the expression in equation 19 in G&Y
        right_j  <- sum(pnorm(t_j))
        right_jp <- sum(pnorm(t_jp))

        ## put it together and what do you get? Bibbidi-Bobbidi-Boo
        left - right_j * right_jp

      }) ## end jp
    }) ## end j

    numer <- sum(par_cov_mat)
    denom <- sum(item_cov_mat)

    numer/denom
  }


  omega_results <- sapply(1:ncol(Lambda), cat_Omega_S_C, Lambda = Lambda, Theta = Theta)
  names(omega_results) <- colnames(Lambda)
  omega_results
}





#' OmegaH
#'
#' Computes hierarchical omega reliability estimate for all factors as described in
#' Rodriguez, Reise, and Haviland (2016).
#'
#' \code{Omega_H} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of factor loadings
#' @param Theta is a vector of indicator error variances
#'
#' @return A \code{numeric}, the omega reliability estimate for all factors.
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
#' Theta <- rep(1, nrow(Lambda)) - rowSums(Lambda^2)
#' Omega_H(Lambda, Theta)
#'
#' @section References:
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models:
#' Calculating and interpreting statistical indices. Psychological Methods, 21(2),
#' 137 \doi{10.1037/met0000045}.
#'
#' @export
#'
#' @seealso \code{\link{Omega_S}}, \code{\link{bifactorIndices}}
#'
#'

Omega_H <- function(Lambda, Theta) {
  Omega_H_C <- function(Fac, Lambda, Theta) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## Compute the appropriate ratio of sums
    sum(Lambda[,Fac])^2/(sum(colSums(Lambda*inFactor)^2) + sum(Theta*inFactor))
  }
  if (is.null(Theta)) return(NULL)
  omega_results <- sapply(1:ncol(Lambda), Omega_H_C, Lambda = Lambda, Theta = Theta)
  names(omega_results) <- colnames(Lambda)
  omega_results
}



#' cat_Omega_S_H
#'
#' Computes hierarchical omega reliability estimate for all factors as described in Rodriguez, Reise, and
#' Haviland (2016).
#'
#' \code{cat_Omega_S_H} is called by \code{\link{bifactorIndices}} and the various convenience functions
#' for exploratory models and/or Mplus output,
#' which are the only functions in this package intended for casual users.
#'
#' @param Lambda is a matrix of standardized factor loadings
#' @param Thresh is a list (indexed by items) of vectors of item thresholds
#' @param Phi is the latent variable covariance matrix. Defaults to \code{NULL}, and
#' the identity matrix will be used. No other options are currently available.
#' @param denom specifies how the variance of the total score will be computed. Defaults
#' to \code{NULL}, and the model implied total score variance will be used. No other options
#' are currently available.
#'
#' @return A \code{numeric}, the hierarchical omega reliability estimate for all factors using
#' the technique of Green and Yang (2009).
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
#' Thresh # Gotta define Thresh
#' catOmega_S_H(Lambda, Thresh)
#'
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models:
#' calculating and interpreting statistical indices. \emph{Psychological Methods, 21}(2),
#' 137 \doi{10.1037/met0000045}.
#'
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
#' structural equation modeling: An alternative to coefficient alpha.
#' \emph{Psychometrika, 74}(1), 155-167 \doi{10.1007/s11336-008-9099-3}.
#'
#' @export
#'
#' @seealso \code{\link{Omega_H}}, \code{\link{bifactorIndices}}
#'
#'

cat_Omega_S <- function(Lambda, Thresh, Phi = NULL, denom = NULL) {
  cat_Omega_S_C <- function(Fac, Lambda, Thresh, Phi = NULL, denom = NULL) {
    ## Make a matrix of logical vectors for non-zero elements of Lambda.
    inFactor <- Lambda[,Fac] != 0
    ## subset Lambda to only include appropriate items
    subLambda <- Lambda[inFactor,]

    ## create Phi matrix
    if (is.null(Phi)) {Phi <- diag(ncol(subLambda))}

    num_items <- nrow(subLambda)

    ## latent item covariances (correlations, because they're standardized!!)
    lat_cov <- subLambda %*% Phi %*% t(subLambda)
    if (is.null(denom)) {
      poly_cor <- lat_cov
      diag(poly_cor) <- rep(1, num_items)
    }

    # Now that we have the poly_cor matrix, we need to restrict to just the factor of interest
    #for computing the numerator

    subLambda <- subLambda[, Fac]
    Phi <- Phi[Fac, Fac]
    lat_cov <- subLambda %*% t(subLambda) * Phi

    ## Create covariance matrix of parallel items
    par_cov_mat <- sapply(1:num_items, function (j) {
      t_j <- Thresh[[j]]
      sapply(1:num_items, function(jp) {
        t_jp <- Thresh[[jp]]
        ## cov(x_j, x_jp)

        ## first, the left half of the expression in equation 19 in Green and Yang (2009)
        left <- sum(sapply(1:length(t_j), function(c) {
          sapply(1:length(t_jp), function(cp) {
            r <- lat_cov[j, jp]
            mnormt::pmnorm(c(t_j[c], t_jp[cp]), c(0, 0), matrix(c(1, r, r, 1), 2))
          }) ## end cp
        })) ## end c

        ## now the two expression in the right half of the expression in equation 19 in G&Y
        right_j  <- sum(pnorm(t_j))
        right_jp <- sum(pnorm(t_jp))

        ## put it together and what do you get? Bibbidi-Bobbidi-Boo
        left - right_j * right_jp

      }) ## end jp
    }) ## end j

    ## create covariance matrix of items... copy/paste of par_cov_mat,
    ## but switch from lat_cov_mat to poly_cor
    item_cov_mat <- sapply(1:num_items, function (j) {
      t_j <- Thresh[[j]]
      sapply(1:num_items, function(jp) {
        t_jp <- Thresh[[jp]]
        ## cov(x_j, x_jp)

        ## first, the left half of the expression in equation 19 in Green and Yang (2009)
        left <- sum(sapply(1:length(t_j), function(c) {
          sapply(1:length(t_jp), function(cp) {
            r <- poly_cor[j, jp]
            mnormt::pmnorm(c(t_j[c], t_jp[cp]), c(0, 0), matrix(c(1, r, r, 1), 2))
          }) ## end cp
        })) ## end c

        ## now the two expression in the right half of the expression in equation 19 in G&Y
        right_j  <- sum(pnorm(t_j))
        right_jp <- sum(pnorm(t_jp))

        ## put it together and what do you get? Bibbidi-Bobbidi-Boo
        left - right_j * right_jp

      }) ## end jp
    }) ## end j

    numer <- sum(par_cov_mat)
    denom <- sum(item_cov_mat)

    numer/denom
  }


  omega_results <- sapply(1:ncol(Lambda), cat_Omega_S_C, Lambda = Lambda, Theta = Theta)
  names(omega_results) <- colnames(Lambda)
  omega_results
}
