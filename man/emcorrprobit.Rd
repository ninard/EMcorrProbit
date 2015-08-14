\name{emcorrprobit}
\alias{emcorrprobit}
\title{Fitting Correlated Probit Models for Ordinal Data}
\usage{
emcorrprobit(model, y, xfixed, xrand, start.values.beta,
  start.values.delta = NULL, start.values.sigma.rand, exact, montecarlo,
  epsilon, ...)
}
\arguments{
\item{model}{specifies the response model; \code{"oneord"} is defined
for a one ordinal response variable}

\item{y}{2-way array for the the response variable with dimensions:
individuals and multiple observations. The ordinal data should be
represented by numeric values in the following way: the first level is
denoted by the number 1, second by the number 2 and so on.}

\item{xfixed}{3-way array for the predictors for the fixed effects with
dimensions: individuals, dimension of the fixed effects and multiple
observations. The intercept should be included as well.}

\item{xrand}{3-way array for the predictors for the random effects with
dimensions: individuals, dimension of the random effects and multiple
observations.}

\item{start.values.beta}{start values for the regression parameters
\eqn{\beta}}

\item{start.values.delta}{start values for the differences in the consecutive
thresholds \eqn{\delta}. Default NULL for binary data. Otherwise it should be specified.}

\item{start.values.sigma.rand}{a matrix with the start values for the
covariance matrix of the random effects \eqn{\Sigma}}

\item{exact}{logical; if TRUE analytical calculation of the moments of
truncated normal distribution is obtained (\code{\link[tmvtnorm:mtmvnorm]{mmvtnorm}} function is used), otherwise a Monte Carlo approach
for estimation is used  (\code{\link[tmvtnorm:rtmvnorm]{rmvtnorm}})}

\item{montecarlo}{numeric; the number of generated values used for the
estimation of the first two moments of truncated normal distribution. If
exact=TRUE this parameter is not needed.}

\item{epsilon}{a value for the stopping criterion.}

\item{\dots}{not used.}
}
\value{
An object of class \code{emcorrprobit}. List with following components:

\item{Sigma.rand.effects}{The estimated covariance matrix of the random effects
\eqn{\Sigma}.}
\item{regression.coefficients}{The estimated regression coefficients \eqn{\beta}.}
\item{differences.in.thresholds}{The estimated differences in the consecutive
thresholds \eqn{\delta}.}
\item{thresholds}{Estimated thresholds \eqn{\alpha}. By definition the first
threshold \eqn{\alpha_1} is zero.}
\item{random.effects}{The estimated random effects for each individual \eqn{b_i}.}
\item{loglikelihood}{Log-likelihood of the model.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{number.iterations}{The number of iterations.}
}
\description{
Maximum likelihood estimation of the parameters of a correlated probit model
via EM algorithm. The function works with wide format of the response data.
The function allows for NA values in the outcome.
}
\details{
The function fits the latent class probit model:
\deqn{y_{ij} = x'_{ij}\beta+ z'_{ij}b_i+\epsilon_{ij},}
where \eqn{y*_{ij} = k is observed, if y_{ij} < \alpha_k} and
\eqn{y*_{ij} = m}, if \eqn{y_{ij} > \alpha_{m-1},} the response variable
y*_{ij} may take a value from 1 to m. We assume \eqn{b_i ~ N(0,\Sigma)} and
\eqn{\epsilon_{ij} ~ N(0,1)}.

The model is fitted using re-parametrisation where new parameters are defined as:
\eqn{\delta_k=\alpha_k-\alpha_{k-1}, k=2,...,m-1}.

The stopping criterion of the algorithm is when the differences between the
estimates from two successive iterations of the algorithm are less than
\code{epsilon} for each parameter.

One should choose carefully the starting values for the parameters (especially
for the covariance matrix of the random effects) and the value of \code{epsilon}.
It is possible that the algorithm stops before convergence and over- or
underestimate the parameters. We recommend using different starting values for
the parameters and only after getting similar results, it can be assumed that
obtained estimates are the MLEs.

When the data consists of 2 or 3 observations per subject it is recommended using the
analytical calculation of the moments of truncated normal distribution
(\code{exact=TRUE}).
}
\examples{
### data simulation
############################################################
### Random intercept model for 3-level ordinal variable ####
### 750 individuals with 2 observations per subject ########
### Predcitors - intercept and time ########################
############################################################
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
#                      xrand=predictors.random,
#                      start.values.beta=beta,start.values.delta=delta,
#                      start.values.sigma.rand=sigma.rand,
#                      exact=e,montecarlo=mc,epsilon=.0002)
#example2=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
#                      xrand=predictors.random,
#                      start.values.beta=beta,start.values.delta=delta,
#                      start.values.sigma.rand=sigma.rand,
#                      exact=e,montecarlo=mc,epsilon=.0002)


###to see the estimates
example1
## the same as
print(example1)

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
\references{
R. V. Gueorguieva. Correlated probit model. In Encyclopedia of Biopharmaceutical
Statistics, chapter 59, pages 355-362. 2006. doi: 10.3109/9781439822463.057. URL
http://informahealthcare.com/doi/abs/10.3109/9781439822463.057

D. Grigorova and R. Gueorguieva. Implementation of the EM algorithm for maximum
likelihood estimation of a random effects model for one longitudinal ordinal outcome.
Pliska Stud. Math. Bulgar., 22:41-56, 2013
}

