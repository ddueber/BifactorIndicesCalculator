% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECV_Indices.R
\name{ECV_SS}
\alias{ECV_SS}
\title{ECV_SS}
\usage{
ECV_SS(Lambda)
}
\arguments{
\item{Lambda}{is a matrix of factor loadings. Be sure that all factors have the same variance
before calling this function.}
}
\value{
A vector of ECVs for all factors
}
\description{
Computes an ECV index for all factors which can be interpreted as the proportion of common
variance of the items in each factor which is due to that factor;
\code{ECV_SS} should be read 'ECV of a specific factor with respect to itself.' Here, ECV is computed
only with respect to items which load on the factor. Note that \code{ECV_SS} of the general factor
is simply the ECV. Stucky and Edelen (2015, p. 201) do not refer to this form of ECV. In the Excel
version of the bifactor indices calculator (Dueber, 2017), this index is referred to as
'ECV (NEW).' \code{ECV_SS} is useful in that it can be computed when there is no general factor, such
as in a two-tier model, and interpreted in the same way as ECV for general factors.
}
\details{
\code{ECV_SS} is called by \code{\link{bifactorIndices}} and the various convenience functions
for exploratory models and/or Mplus output,
which are the only functions in this package intended for casual users.``
}
\examples{
Lambda <- matrix(c(.82, .10,   0,   0,
                   .77, .35,   0,   0,
                   .79, .32,   0,   0,
                   .66, .39,   0,   0,
                   .51,   0, .71,   0,
                   .56,   0, .43,   0,
                   .68,   0, .13,   0,
                   .60,   0, .50,   0,
                   .83,   0,   0, .47,
                   .60,   0,   0, .27,
                   .78,   0,   0, .28,
                   .55,   0,   0, .75),
                   ncol = 4, byrow = TRUE)
colnames(Lambda) <- c("General", "SF1", "SF2", "SF3")
ECV_SS(Lambda)


}
\references{
Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based tool to calculate various indices relevant to bifactor CFA models. \doi{10.13023/edp.tool.01}

Stucky, B. D., & Edelen, M. O. (2015). Using hierarchical IRT models to create unidimensional measures from multidimensional data. In S. P. Reise & D. A. Revicki (Eds.), \emph{Handbook of item response theory modeling: Applications to typical performance assessment} (pp.183-206). New York: Routledge.
}
\seealso{
\code{\link{ECV_SG}}, \code{\link{ECV_GS}}, \code{\link{bifactorIndices}}
}
