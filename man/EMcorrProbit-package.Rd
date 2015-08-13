\name{EMcorrProbit}
\alias{EMcorrProbit}
\docType{package}
\title{
Maximum likelihood estimation of correlated probit models via EM algorithm
}
\description{
The package is designed for maximum likelihood estimation of correlated probit models especially for ordinal longitudinal or clustered data. The fitting approach is EM algorithm. The main function \code{emcorrprobit} provides estimates of the parameters of the model, random effects, log-likelihood, AIC and BIC. Standard error estimation is done via bootstrap method.

}
\details{
\tabular{ll}{
Package: \tab EMcorrProbit\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-08-21\cr
License: \tab GPL-2\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Denitsa Grigorova, Nina Daskalova

Maintainer: Denitsa Grigorova <dpgrigorova@abv.bg>

}
\references{
R. V. Gueorguieva. Correlated probit model. In Encyclopedia of Biopharmaceutical Statistics, chapter 59, pages 355-362. 2006. doi: 10.3109/9781439822463.057. URL http://informahealthcare.com/doi/abs/10.3109/9781439822463.057. 


D. Grigorova and R. Gueorguieva. Implementation of the EM algorithm for maximum likelihood estimation of a random effects model for one longitudinal ordinal outcome. Pliska Stud. Math. Bulgar., 22:41-56, 2013


}

\keyword{ 
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
package }
\seealso{
\code{\link[tmvtnorm:tmvtnorm-package]{tmvtnorm}}
\code{\link[doParallel:doParallel-package]{doParallel}}
}
\examples{
### data simulation
############################################################
## RI model, 2 predictors, 3-level ordinal variable #######
### 750 individuals - 2 observations per subject ###########
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
mc=200
e=T

example.complete.cases=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
                      xrand=predictors.random,
                      start.values.beta=beta,start.values.delta=delta,
                      start.values.sigma.rand=sigma.rand,
                      exact=e,montecarlo=mc,epsilon=.0001)

###to see the estimates
print(example.complete.cases) 

###for standard errors and respective P-values
summary(example.complete.cases) 


}
