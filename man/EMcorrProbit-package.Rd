\name{EMcorrProbit}
\alias{EMcorrProbit}
\docType{package}
\title{
Maximum likelihood estimation of correlated probit models via EM algorithm
}
\description{
The package is designed for maximum likelihood estimation of correlated probit models, especially for ordinal longitudinal or clustered data. The fitting approach is EM algorithm. Two main functions are implemented. 

Function \code{emcorrprobit} provides estimates of the parameters of the model, random effects, log-likelihood, AIC and BIC.

The \code{summary} method provides standard error estimation via bootstrap method. 
}
\details{
\tabular{ll}{
Package: \tab EMcorrProbit\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-08-21\cr
License: \tab GPL-2\cr
}
First main parameter in the \code{emcorrprobit} function is the \code{model} parameter. For the present moment only correlated probit model for one longitudinal ordinal variable is implemented (\code{model="oneord"}). Joint model for two longitudinal ordinal variables and joint model for one longitudinal ordinal and continuous variables will be added. For more details see \code{\link{emcorrprobit}}.

The \code{summary} function allows for parallel computation which is faster. The \code{doParallel} parameter should be set to \code{TRUE}. The recommended number of bootstrap samples is between 50 and 100. For more details see \code{\link{summary.emcorrprobit}}.

}
\author{
Denitsa Grigorova, Nina Daskalova

Maintainer: Denitsa Grigorova <dpgrigorova@abv.bg>

}
\references{
R. V. Gueorguieva. Correlated probit model. In Encyclopedia of Biopharmaceutical
Statistics, chapter 59, pages 355-362. 2006. doi: 10.3109/9781439822463.057. URL
http://informahealthcare.com/doi/abs/10.3109/9781439822463.057

D. Grigorova and R. Gueorguieva. Implementation of the EM algorithm for maximum
likelihood estimation of a random effects model for one longitudinal ordinal outcome.
Pliska Stud. Math. Bulgar., 22:41-56, 2013
}

\keyword{ 
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
package }
\seealso{
\code{\link[doParallel:doParallel-package]{doParallel}}
\code{\link[tmvtnorm:rtmvnorm]{rmvtnorm}}
}
\examples{
### data simulation
############################################################
### Random intercept model for 3-level ordinal variable ####
### 750 individuals with 2 observations per subject ########
### Predcitors - intercept and time ########################
############################################################
rm(list=ls())
random.int=rnorm(750,0,0.1)
l=length(random.int)
int=-0.5
mult.obs=2
y1=sapply(1:mult.obs,function(i) random.int+int+i+rnorm(l), simplify="array")
data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.5,2,3))
table(data.ordinal)
head(data.ordinal)
time=sapply(1:mult.obs, function(i) rep(i,l), simplify="array")
predictors.fixed=sapply(1:mult.obs, function(i) cbind(1,time[,i]), simplify="array")
predictors.random=sapply(1:mult.obs, function(i) matrix(rep(1,l),ncol=1), simplify="array")

sigma.rand=matrix(.01)
beta=c(-0.55,0.95)
delta=c(1.5)
mc=500
e=TRUE

### estimation
###should work
example1=emcorrprobit(model = "oneord", y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta=beta,start.values.delta=delta,
                    start.values.sigma.rand=sigma.rand,
                    exact=e,montecarlo=mc,epsilon=.0002)

###doesn't work
#example2=emcorrprobit(model = "1ord", y=data.ordinal,xfixed=predictors.fixed,
#                    xrand=predictors.random,
#                    start.values.beta=beta,start.values.delta=delta,
#                    start.values.sigma.rand=sigma.rand,
#                    exact=e,montecarlo=mc,epsilon=.0002)


###to see the estimates
example1
## the same as 
print(example1)


###for standard errors and respective P-values
ex1.se=summary(example1, doParallel=TRUE, bootstrap.samples=50)
###
ex1.se
### variance-covariance matrix of the estimates
ex1.se$vcov

### Monte Carlo approach for estimation of moments of truncated normal distribution - slower in this 
### case
example2=emcorrprobit(model = "oneord", y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta=beta,start.values.delta=delta,
                    start.values.sigma.rand=sigma.rand,
                    exact=FALSE,montecarlo=mc,epsilon=.0002)




### example with missing data 
data.ordinal[1,2]=NA
head(data.ordinal)

example3=emcorrprobit(model = "oneord", y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta=beta,start.values.delta=delta,
                    start.values.sigma.rand=sigma.rand,
                    exact=TRUE,montecarlo=mc,epsilon=.0002)




}
