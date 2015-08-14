\name{print.emcorrprobit}
\alias{print}
\alias{print.emcorrprobit}
\title{Printing an EMcorrProbit fit}
\usage{

\method{print}{emcorrprobit}(x, ...)

}
\arguments{
\item{x}{an object of class \code{'emcorrprobit'}}

\item{\dots}{not used.}

}
\value{
The \code{'emcorrprobit'} is a list.

}
\description{
\code{print} of an emcorrprobit object
}
\details{
Sigma.rand.effects=sigma.rand.new,
       regression.coefficients=betanew,
       differences.in.thresholds=deltanew,
       thresholds=c(0,cumsum(deltanew)),
       random.effects=firstmomentb,
       loglikelihood=loglikelihood,
       AIC=AIC,
       BIC=BIC,
       number.iterations=number.it,
       data.ordinal=data.ordinal,
       predictors.fixed=predictors.fixed,
       predictors.random=predictors.random,
       exact=exact,
       montecarlo=montecarlo,
       epsilon=epsilon}
\examples{
}

