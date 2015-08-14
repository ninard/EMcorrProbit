\name{summary.emcorrprobit}
\alias{summary}
\alias{summary.emcorrprobit}
\alias{print.summary.emcorrprobit}
\title{Summarizing an EMcorrProbit fit}
\usage{

\method{summary}{emcorrprobit}(object, ...)
\method{print}{summary.emcorrprobit}(x, ...)

}
\arguments{
\item{object}{an object of class \code{'emcorrprobit'}}

\item{x}{an object of class \code{'summary.emcorrprobit'}}

\item{\dots}{arguments passed to or from other methods.\code{bootstrap.samples}: the number of samples used in bootsrap method, between 50 and 100 is recommended, 50 by default; \code{doParallel}: logical, if TRUE \code{\link[foreach:foreach]{foreach}} function from \code{\link[doParallel:doParallel-package]{doParallel}} package is used to speed up the calculation of standard errors, FALSE by default; \code{cores}: the number of cores used in parallel computation, if NULL (by default), the number of cores is set by the \code{\link[doParallel:doParallel-package]{doParallel}} package; \code{epsilon}: a value for the stopping criterion. If not specified, it is defined as 10 times \code{epsilon} used for \code{emcorrprobit} fit of \code{object}.}

}
\value{
The function \code{summary.emcorrprobit} computes and returns the covariance matrix
of the parameters' estimates via bootstrap method. Standard errors, z-scores and p-values
of the \code{emcorrprobit} estimates are presented via \code{print} function.

\code{vcov} the covariance matrix of the parameters' estimates
}
\description{
\code{summary} of an emcorrprobit object
}
\details{
\code{print.summary.emcorrprobit} uses smart formatting of the coefficients,
standard errors, etc.

Standard errors are obtained via bootstrap method and a summary table with
respective z-scores and p-values is printed.
}
\examples{
### using doParallel package
### for standard errors and respective P-values
#ex1.se=summary(example1, doParallel=TRUE, bootstrap.samples=50)
### print
#ex1.se
### variance-covariance matrix of the estimates
#ex1.se$vcov

### without parallel computations - very slow
#ex2.se=summary(example1, bootstrap.samples=50)
#ex2.se$vcov
}

