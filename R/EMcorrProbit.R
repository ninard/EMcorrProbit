#' Fitting Correlated Probit Models for Ordinal Data
#' 
#' Maximum likelihood estimation of the parameters of a correlated probit model via EM algorithm. The function works with wide format of the response data. The function allows for NA values in the outcome. 
#' @param y 2-way array for the the response variable with dimensions: individuals and multiple observations. The ordinal data should be represented by numeric values in the following way: the first level is denoted by the number 1, second by the number 2 and so on. 
#' @param xfixed 3-way array for the predictors for the fixed effects with dimensions: individuals, dimension of the fixed effects and multiple observations. The intercept should be included as well.
#' @param xrand 3-way array for the predictors for the random effects with dimensions: individuals, dimension of the random effects and multiple observations
#' @param exact logical. If TRUE analytical calculation of the moments of truncated normal distribution is obtained, otherwise a Monte Carlo approach for estimation is used.
#' @param montecarlo numeric. The number of generated values used for the estimation of the first two moments of truncated normal distribution. If exact=T this parameter is not needed.
#' @param start.values.delta start values for the differences in the consecutive thresholds \eqn{\delta} 
#' @param start.values.beta start values for the regression parameters \eqn{\beta}
#' @param start.values.sigma.rand a matrix with the start values for the covariance matrix of the random effects \eqn{\Sigma}
#' @param epsilon a value for the stopping criterion
#' @details The function fits the latent class probit model:
#' \deqn{y_{ij} = x'_{ij}\beta+ z'_{ij}b_i+\epsilon_{ij},}{y_ij = x'_ij*\beta+ z'_ij*b_i+\epsilon_ij,} 
#' where we observe \eqn{y*_{ij} = k, if y_{ij} < \alpha_k} and \eqn{y*_{ij} = m}, if \eqn{y_{ij} > \alpha_{m-1},} the response variable y*_{ij} may take a value from 1 to m. We assume \eqn{b_i ~ N(0,\Sigma)} and \eqn{\epsilon_{ij} ~ N(0,1)}.
#' 
#' The model is fitted using re-parametrisation where new parameters are defined: \eqn{\delta_k=\alpha_k-\alpha_{k-1}, k=2,...,m-1}
#' 
#' The stopping criterion of the algorithm is when the differences between the estimates from two successive iterations of the algorithm are less than \code{epsilon} for each parameter.
#' 
#' One should choose carefully the starting values for the parameters (especially for the covariance matrix of the random effects) and the value of \code{epsilon}. It is possible that the algorithm stops before convergence and over- or underestimate the parameters. We recommend using different starting values for the parameters and only after getting similar results, it can be assumed that obtained estimates are the MLEs.
#' 
#' When the data consists of 2 or 3 observations per subject we recommend using the analytical calculation of the moments of truncated normal distribution (\code{exact=T}).
#' @return An object of class \code{"emcorrprobit"}. List with following components:
#'   
#' \item{Sigma.rand.effects}{The estimated covariance matrix of the random effects \eqn{\Sigma}}
#' \item{regression.coefficients}{The estimated regression coefficients \eqn{\beta}}
#' \item{differences.in.thresholds}{The estimated differences in the consecutive thresholds \eqn{\delta}}
#' \item{thresholds}{Estimated thresholds \eqn{\alpha}. By definition the first threshold \eqn{\alpha_1} is zero}
#' \item{random.effects}{The estimated random effects for each individual \eqn{b_i}}
#' \item{loglikelihood}{Log-likelihood of the model}
#' \item{AIC}{Akaike information criterion}
#' \item{BIC}{Bayesian information criterion}
#' \item{number.iterations}{The number of iterations}
#' @examples
#' ### data simulation
#' ############################################################
#' ## RI model, 2 predictors, 3-level ordinal variable #######
#' ### 750 individuals - 2 observations per subject ###########
#' ############################################################
#' rm(list=ls())
#' random.int=rnorm(750,0,0.1)
#' l=length(random.int)
#' int=-0.5
#' mult.obs=2
#' y1=sapply(1:mult.obs,function(i) random.int+int+i+rnorm(l), simplify="array")
#' data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.5,2,3))
#' table(data.ordinal)
#' head(data.ordinal)
#' time=sapply(1:mult.obs, function(i) rep(i,l), simplify="array")
#' predictors.fixed=sapply(1:mult.obs, function(i) cbind(1,time[,i]), simplify="array")
#' predictors.random=sapply(1:mult.obs, function(i) matrix(rep(1,l),ncol=1), simplify="array")
#'
#' sigma.rand=matrix(.01)
#' beta=c(-0.55,0.95)
#' delta=c(1.5)
#' mc=200
#' e=T
#'
#' example.complete.cases=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
#'                        xrand=predictors.random,
#'                        start.values.beta=beta,start.values.delta=delta,
#'                        start.values.sigma.rand=sigma.rand,
#'                        exact=e,montecarlo=mc,epsilon=.0001)
#' print(example.complete.cases)
#'
#' data.ordinal[1,2]=NA
#' head(data.ordinal)
#'
#' example.missings=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
#'                    xrand=predictors.random,
#'                    start.values.beta=beta,start.values.delta=delta,
#'                    start.values.sigma.rand=sigma.rand,
#'                    exact=e,montecarlo=mc,epsilon=.0001)
#' print(example.missings)

emcorrprobit <- function(y, ...) UseMethod("emcorrprobit")

emcorrprobit.default <- function(y, xfixed, xrand, start.values.beta, 
                                 start.values.delta=NULL,  start.values.sigma.rand, 
                                 exact, montecarlo=100, epsilon=.001, ...)
{
  #cat("Please, be patient, the function is working ... \n")
  
  #######################
  ## check
  ############################
  if(!is.matrix(y)) warning("y must be a matrix")
  if(!is.matrix(start.values.sigma.rand)) warning("start.values.sigma.rand must be a matrix")
  if(any(sort(unique(c(y)))!=1:max(y,na.rm=T))) stop("Missing levels in the response variable.")
  if(length(start.values.delta)!=(max(y,na.rm=T)-2)) stop("Incorrect dimension of start.values.delta.")
  #if(exact==F & is.na(montecarlo)) {warning("montecarlo parameter undefined, 100 taken by default")
  #montecarlo=100
  #}
  if(exact==F & is.na(montecarlo)) stop("montecarlo parameter undefined")
  
  xfixed <- as.array(xfixed)
  xrand <- as.array(xrand)
  y <- as.array(y)
  
  if (dim(xfixed)[1] != dim(xrand)[1] | dim(xrand)[1] != dim(y)[1]) stop("Incompatible dimensions!")
  if (dim(xfixed)[3] != dim(xrand)[3] | dim(xrand)[3] != dim(y)[2]) stop("Incompatible dimensions!")
  
  est <- ecm.one.ordinal(y,xfixed,xrand,start.values.beta,start.values.delta,
                         start.values.sigma.rand,
                         exact,montecarlo,epsilon, additional=T)
  
  
  est$call <-match.call()
  
 # cat("It's ready! Use print or summary to see the result.\n")
  
  class(est) <- "emcorrprobit"
  est
}

print.emcorrprobit <- function(x, ...)
{
  cat("Call: \n")
  print(x$call)
  cat("\nNumber of iterations: ",x$number.it,"\n")
  cat("\nCovariance matrix of the random effects: \n")
  print(x$Sigma.rand.effects)
  cat("\nRegression coefficients: \n")
  cat(x$regression.coefficients,"\n")
  if(length(x$differences.in.thresholds)>0) 
    {cat("\nDifferences in thresholds: \n")
     cat(x$differences.in.thresholds,"\n")
    } 
  cat("\nThresholds: \n")
  cat(x$thresholds,"\n")
  cat("\nRandom efects: \n")
  print(head(x$random.effects))
  cat("..............\n")
  cat("\nLog-likelihood: \n")
  cat(x$loglikelihood,"\n")
  cat("\nAIC: \n")
  cat(x$AIC,"\n")
  cat("\nBIC: \n")
  cat(x$BIC,"\n")
}

summary.emcorrprobit <- function(x, ...)
{ cat(" Please, be very patient ... \n")
  vcov <- standard.error.bootstrap.one.ordinal(x, ...) 
  
  se <- sqrt(diag(vcov))
  
  se.sigma <- matrix(se[1:length(x$Sigma.rand.effects)],
                     ncol=sqrt(length(x$Sigma.rand.effects)))
  
  se.sigma <- se.sigma[lower.tri(se.sigma,diag=T)]
  
  TAB.sigma <- cbind(Sigma= x$Sigma.rand.effects[lower.tri(x$Sigma.rand.effects,diag=T)], 
                    StdErr = se.sigma,
                    z.score=x$Sigma.rand.effects[lower.tri(x$Sigma.rand.effects,diag=T)]/se.sigma)
  nam <- matrix(paste("sigma ", outer(1:dim(x$Sigma.rand.effects),1:dim(x$Sigma.rand.effects),
                                  function(x,y) paste(x,y,sep="")),sep=""),ncol=dim(x$Sigma.rand.effects))
  nam <- nam[lower.tri(nam,diag=T)]
  rownames(TAB.sigma)=nam
  
  
  se.regr.coeff <- se[(length(x$Sigma.rand.effects)+1):
                        (length(x$Sigma.rand.effects)+length(x$regression.coefficients))]
  
  TAB.regr.coeff <- cbind(Regression.coeff= x$regression.coefficients, 
               StdErr = se.regr.coeff,
               z.score = x$regression.coefficients/se.regr.coeff,
               p.value = 2*pnorm(abs(x$regression.coefficients/se.regr.coeff), lower=F))
  rownames(TAB.regr.coeff)=paste("Predictor",1:length(x$regression.coefficients))
    
  if(length(x$differences.in.thresholds)>0)
    
  {se.diff.thresholds =se[(length(se)-length(x$differences.in.thresholds)+1):length(se)]
  
  TAB.diff.thresholds <- cbind(Threshold.differences= x$differences.in.thresholds, 
               StdErr = se.diff.thresholds,
               z.score = x$differences.in.thresholds/se.diff.thresholds,
               p.value = 2*pnorm(abs(x$differences.in.thresholds/se.diff.thresholds), lower=F))
  rownames(TAB.diff.thresholds)=paste("Diff",1:length(x$differences.in.thresholds))
  } else TAB.diff.thresholds=NULL
  
  res <- list(call = x$call, 
              regr.coefficients = TAB.regr.coeff,
              diff.thresholds=TAB.diff.thresholds,
              sigma=TAB.sigma,
              vcov=vcov)
  
  class(res) <- "summary.emcorrprobit"
  res
}

print.summary.emcorrprobit <- function(x, ...)
{
  cat("Call: \n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$sigma)
  cat("\n")
  printCoefmat(x$regr.coefficients,P.values=T, has.Pvalue=T)
  cat("\n")
  printCoefmat(x$diff.thresholds,P.values=T, has.Pvalue=T)
  
}

# formula.emcorrprobit <- function(formula, data=list(), ...)
# {
#   mf <- model.frame(formula=formula, data=data)
#   x <- model.matrix(attr(mf, "terms"), data = mf)
#   y <- model.response(mf)
#   
#   est <- emcorrprobit.default(x, y, ...)
#   est$call <-match.call
#   est$formula <-formula
#   est
# }

ecm.one.ordinal <- function(data.ordinal,predictors.fixed,predictors.random,start.values.beta,start.values.delta,start.values.sigma.rand,
                         exact,montecarlo,epsilon, additional) 
  {
  ########################################
  ### wide format of longitudinal data: first level denoted by 1, second - 2 and so on
  ### predictors.fixed is a 3-dimentional array: individuals x dimentions predictros fixed x time points 
  ### predictors.random is a 3-dimentional array: individuals x dimentions predictros fixed x time points  
  ### sigma.rand is the covaraince matrix of the random effects
  ##################################
  #library(MASS)
  #library(tmvtnorm)
  #library(doParallel)
  #registerDoParallel()
  
  
  ##############################################################
  
  betanew <- start.values.beta
  deltanew <- start.values.delta
  sigma.rand.new <- start.values.sigma.rand
  
  num.categories <- max(data.ordinal,na.rm=T)
  if(num.categories==2) deltanew <- NULL
  
  
  new.est <- c(start.values.sigma.rand,start.values.beta,start.values.delta)
  old.est <- new.est+1.5*epsilon
  
  number.it <- 0
  
  # definition of lower and upper boundary of the truncated normal distr. 
  # (the new variable given observed data)
  trunc.lower <- ifelse(data.ordinal==1,-Inf,0)
  trunc.upper <- ifelse(data.ordinal==1,0,ifelse(data.ordinal<num.categories,1,Inf))
  
  #number of observations in each category
  n <- table(data.ordinal)
  individuals <- length(data.ordinal[,1])
  mult.obs <- dim(predictors.fixed)[3]
  dimension.sigma <- ncol(sigma.rand.new)
  
  z <- sapply(1:individuals, function(i) matrix(t(predictors.random[i,,]), ncol=dimension.sigma,byrow=F), 
           simplify="array")
  ## mult.obs x dimension.sigma x individuals
  
  ### definition na k, knot i kk
  k <- lapply(1:individuals, function(i) which(is.na(data.ordinal[i,])))
  ### list: missings for each individual
  knot <- lapply(1:individuals, function(i) which(!is.na(data.ordinal[i,])))
  ### list: observed for each individual
  kk <- sapply(1:individuals, function(i) length(k[[i]]), simplify=T)
  ### individuals: length of missings for each individual
  
  ######################################################################
  ############             while loop           ##########################
  ##########################################################################
  while(any(abs(new.est-old.est)>epsilon)) {
    
    number.it <- number.it+1
    
    old.est <- new.est
    
    #######################################################################
    delta.exp <- c(1,deltanew,1)
    delta <- matrix(delta.exp[data.ordinal],ncol=mult.obs)
    alpha.exp <- c(0,0,cumsum(deltanew))
    alpha <- matrix(alpha.exp[data.ordinal],ncol=mult.obs)
    
    ###expectation of ynew
    if(length(start.values.beta)>1) pred <- sapply(1:mult.obs, function(i) predictors.fixed[,,i]%*%betanew, simplify=TRUE)  else pred <- betanew
    ## makes pred a matrix: observations x mult.obs
    mean.ynew <- (pred-alpha)/delta
    
    variance.ynew <- sapply(1:individuals, function(i)
      (z[,,i]%*%sigma.rand.new%*%t(z[,,i])+diag(mult.obs))/(delta[i,]%*%t(delta[i,])),
      simplify="array")
    ###
    
    #######################################################
    if(exact==F) {
      simulations <- lapply(1:individuals, function(i)  
        rtmvnorm(montecarlo, mean=mean.ynew[i,knot[[i]]], 
                 sigma=matrix(as.vector(variance.ynew[knot[[i]],knot[[i]],i]),ncol=mult.obs-kk[i]),
                 lower=trunc.lower[i,knot[[i]]],upper=trunc.upper[i,knot[[i]]], algorithm="gibbs") )  
      ## simulations is a list
      moments1<-lapply(simulations, function(i) apply(as.matrix(i),2,mean))
      ## list
      moments2<-lapply(simulations, var)
      ## list
      ## montecarlo x mult.obs x observations
    } else {
      simulations.exact <- lapply(1:individuals, function(i) 
        mtmvnorm(mean=mean.ynew[i,knot[[i]]], 
                 sigma=matrix(as.vector(variance.ynew[knot[[i]],knot[[i]],i]),ncol=mult.obs-kk[i]),
                 lower=trunc.lower[i,knot[[i]]],upper=trunc.upper[i,knot[[i]]]) )
      ##  mult.obs x observations
      moments1 <- lapply(simulations.exact, function(i) i$tmean)
      ## list
      moments2 <- lapply(simulations.exact, function(i) i$tvar)
      ## the variance of y latent
    }

    moments1A <- array(NA,dim=c(individuals,mult.obs))
    for(i in 1:individuals) moments1A[i,knot[[i]]] <- moments1[[i]]
    
    moments2A <- list(0)
    for(i in 1:individuals) {moments2A[[i]] <- array(NA, dim=c(mult.obs, mult.obs))
                             moments2A[[i]][knot[[i]],knot[[i]]] <- moments2[[i]]}
    
    
    ######### to calculate the expectation of the random effects
    
    sigmab <- lapply(1:individuals, function(i) 
      if(kk[i]!=(mult.obs-1)) t(t(sigma.rand.new%*%t(z[knot[[i]],,i]))/delta[i,knot[[i]]])%*%solve(variance.ynew[knot[[i]],knot[[i]],i]) else
        t(t(sigma.rand.new%*%z[-k[[i]],,i])/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]))
    ### list: dimension.sigma x mult.obs x individuals
    
    firstmomentb <- sapply(1:individuals, function(i) {
      if(kk[i]!=(mult.obs-1)) sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,knot[[i]]]) else 
        sigmab[[i]]%*%as.matrix(moments1[[i]]-mean.ynew[i,-k[[i]]])}, simplify="array") 
    ### array: dimension.sigma x 1 x individuals
    ### when dimension.sigma=1: vector : individuals 
    
    
    if (dimension.sigma>1)  pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
      pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
    ### individuals x mult.obs
    
    
    
    ##################################################
    #########################  beta estimate
    ##################################################
    
    
    ywave <- delta*moments1A-pred.random+alpha
    if(length(start.values.beta)==1) betanew <- lm(c(ywave)~1)$coef else {
      pred.beta <- predictors.fixed[,,1]
      for(i in 2:mult.obs) pred.beta <- rbind(pred.beta,predictors.fixed[,,i])
      betanew <- lm(c(ywave)~pred.beta-1)$coef 
      rm(pred.beta)}
    
    ############ update
    
    ###expectation of ynew
    if(length(start.values.beta)>1) pred <- sapply(1:mult.obs, function(i) predictors.fixed[,,i]%*%betanew, simplify=TRUE)  else pred <- betanew
    ## makes pred a matrix: observations x mult.obs
    mean.ynew <- (pred-alpha)/delta
    
    
    firstmomentb <- sapply(1:individuals, function(i) {
      if(kk[i]!=(mult.obs-1)) sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,knot[[i]]]) else 
        sigmab[[i]]%*%as.matrix(moments1[[i]]-mean.ynew[i,-k[[i]]])}, simplify="array") 
    ### array: dimension.sigma x 1 x individuals
    ### when dimension.sigma=1: vector : individuals 
    
    
    if (dimension.sigma>1)  pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
      pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
    ### individuals x mult.obs
    
    
    
    #######################################
    ################### delta estimates
    ##########################################
    if(num.categories>2) {
      
      ### where NA - with 0 - for first and second1
      first <- t(sapply(1:individuals, function(i) diag(moments2A[[i]])+(moments1A[i,])^2, simplify="array"))
      ### individuals x mult.obs
      
      for (delta.index in 1:(num.categories-2)) {
        
        second1 <- moments1A*(pred-alpha)
        ### individuals x mult.lbs
        
        alpha.exp.term <- c(0,0)
        for(i in 1:length(deltanew)) alpha.exp.term[i+2] <- ifelse(delta.index<=i, alpha.exp[i+2]-deltanew[delta.index], 0)
        alpha.term <- matrix(alpha.exp.term[data.ordinal],ncol=mult.obs)
        
        
        #### to check it carefully - if I have only one observation for some individual
        second2 <- lapply(1:individuals, function(i)
          sigmab[[i]]%*%(moments2[[i]]+moments1[[i]]%*%t(moments1[[i]])-mean.ynew[i,knot[[i]]]%*%t(moments1[[i]])) )
        ## for each i: q x n_i matrix 
        
        if(dimension.sigma==1)  second2 <- lapply(1:individuals, function(i)  
          z[knot[[i]],,i]*as.vector(second2[[i]]) )  else
            second2 <- sapply(1:individuals, function(i) diag(z[knot[[i]],,i]%*%second2[[i]]))
        
        second2new <- array(NA, dim=c(individuals,mult.obs))
        for(i in 1:individuals) second2new[i,knot[[i]]] <- second2[[i]]
        
        
        third <- moments1A*delta-pred-pred.random+alpha.term
        ## observations x mult.obs
        
        a <- sum(first[data.ordinal==(delta.index+1)],na.rm=T)+sum(n[(delta.index+2):num.categories])
        b <- sum((second1+second2new)[data.ordinal==delta.index+1], na.rm=T)-sum(third[data.ordinal>delta.index+1], na.rm=T)
        c <- -n[(delta.index+1)]
        
        deltanew[delta.index] <- (b+sqrt(b^2-4*a*c))/(2*a) 
        
        #################### update
        
        delta.exp <- c(1,deltanew,1)
        delta <- matrix(delta.exp[data.ordinal],ncol=mult.obs)
        alpha.exp <- c(0,0,cumsum(deltanew))
        alpha <- matrix(alpha.exp[data.ordinal],ncol=mult.obs)
        
        mean.ynew <- (pred-alpha)/delta
        
        variance.ynew <- sapply(1:individuals, function(i)
          (z[,,i]%*%sigma.rand.new%*%t(z[,,i])+diag(mult.obs))/(delta[i,]%*%t(delta[i,])),
          simplify="array")
        
        
        sigmab <- lapply(1:individuals, function(i) 
          if(kk[i]!=(mult.obs-1)) t(t(sigma.rand.new%*%t(z[knot[[i]],,i]))/delta[i,knot[[i]]])%*%solve(variance.ynew[knot[[i]],knot[[i]],i]) else
            t(t(sigma.rand.new%*%z[-k[[i]],,i])/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]))
        ### list: dimension.sigma x mult.obs x individuals
        
        firstmomentb <- sapply(1:individuals, function(i) {
          if(kk[i]!=(mult.obs-1)) sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,knot[[i]]]) else 
            sigmab[[i]]%*%as.matrix(moments1[[i]]-mean.ynew[i,-k[[i]]])}, simplify="array") 
        ### array: dimension.sigma x 1 x individuals
        ### when dimension.sigma=1: vector : individuals 
        
        
        if (dimension.sigma>1)  pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
          pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
        ### individuals x mult.obs
      } 
      
    } #if
    #############################################################################################
    ############################### sigma.rand estimate
    #############################################################################################
    
    sigma.rand <- sapply(1:individuals, function(i) {
      if(kk[[i]]!=(mult.obs-1)) sigma.rand.new-sigmab[[i]]%*%((z[knot[[i]],,i]%*%sigma.rand.new)/delta[i,knot[[i]]])+
        sigmab[[i]]%*%(moments2[[i]]+moments1[[i]]%*%t(moments1[[i]])
                       -(moments1[[i]])%*%t(mean.ynew[i,knot[[i]]])
                       -(mean.ynew[i,knot[[i]]])%*%t(moments1[[i]])
                       +(mean.ynew[i,knot[[i]]])%*%t(mean.ynew[i,knot[[i]]])
        )%*%t(sigmab[[i]]) else    sigma.rand.new-sigmab[[i]]%*%((z[-k[[i]],,i]%*%sigma.rand.new)/delta[i,-k[[i]]])+
        sigmab[[i]]%*%(moments2[[i]]+moments1[[i]]%*%t(moments1[[i]])-moments1[[i]]%*%
                         t(mean.ynew[i,-k[[i]]])-(mean.ynew[i,-k[[i]]])%*%t(moments1[[i]])+(mean.ynew[i,-k[[i]]])%*%t(mean.ynew[i,-k[[i]]])
        )%*%t(sigmab[[i]])    
    }, simplify=T)
    
    
    if(dimension.sigma==1) sigma.rand.new <- matrix(mean(sigma.rand)) else {sigma.rand.new <- matrix(apply(sigma.rand,1,mean), ncol=dimension.sigma)
                                                                            sigma.rand.new[lower.tri(sigma.rand.new)]  <-  t(sigma.rand.new)[lower.tri(sigma.rand.new)] }
    ########################
    new.est <- c(sigma.rand.new,betanew,deltanew)
    ########################
       } # while

  if(additional) {
  ###########################################################################
  ################### random effects
  ###########################################################################
  delta.exp <- c(1,deltanew,1)
  delta <- matrix(delta.exp[data.ordinal],ncol=mult.obs)
  alpha.exp <- c(0,0,cumsum(deltanew))
  alpha <- matrix(alpha.exp[data.ordinal],ncol=mult.obs)
  
  ###expectation of ynew
  if(length(start.values.beta)>1) pred <- sapply(1:mult.obs, function(i) predictors.fixed[,,i]%*%betanew, simplify=TRUE)  else pred <- betanew
  ## makes pred a matrix: observations x mult.obs
  mean.ynew <- (pred-alpha)/delta
  
  variance.ynew <- sapply(1:individuals, function(i)
    (z[,,i]%*%sigma.rand.new%*%t(z[,,i])+diag(mult.obs))/(delta[i,]%*%t(delta[i,])),
    simplify="array")
  ###
  
 ##################################################################################
  if(exact==F) {
  simulations <- lapply(1:individuals, function(i)  
    rtmvnorm(montecarlo, mean=mean.ynew[i,knot[[i]]], 
             sigma=matrix(as.vector(variance.ynew[knot[[i]],knot[[i]],i]),ncol=mult.obs-kk[i]),
             lower=trunc.lower[i,knot[[i]]],upper=trunc.upper[i,knot[[i]]], algorithm="gibbs") )  
  ## simulations is a list
  moments1<-lapply(simulations, function(i) apply(as.matrix(i),2,mean))
  ## list
  moments2<-lapply(simulations, var)
  ## list
  ## montecarlo x mult.obs x observations
  } else {
  simulations.exact <- lapply(1:individuals, function(i) 
    mtmvnorm(mean=mean.ynew[i,knot[[i]]], 
             sigma=matrix(as.vector(variance.ynew[knot[[i]],knot[[i]],i]),ncol=mult.obs-kk[i]),
             lower=trunc.lower[i,knot[[i]]],upper=trunc.upper[i,knot[[i]]]) )
  ##  mult.obs x observations
  moments1 <- lapply(simulations.exact, function(i) i$tmean)
  ## list
  moments2 <- lapply(simulations.exact, function(i) i$tvar)
  ## the variance of y latent
}

  moments1A <- array(NA,dim=c(individuals,mult.obs))
  for(i in 1:individuals) moments1A[i,knot[[i]]] <- moments1[[i]]

  moments2A <- list(0)
  for(i in 1:individuals) {moments2A[[i]] <- array(NA, dim=c(mult.obs, mult.obs))
                         moments2A[[i]][knot[[i]],knot[[i]]] <- moments2[[i]]}


  ######### to calculate the expectation of the random effects

  sigmab <- lapply(1:individuals, function(i) 
  if(kk[i]!=(mult.obs-1)) t(t(sigma.rand.new%*%t(z[knot[[i]],,i]))/delta[i,knot[[i]]])%*%solve(variance.ynew[knot[[i]],knot[[i]],i]) else
    t(t(sigma.rand.new%*%z[-k[[i]],,i])/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]))
  ### list: dimension.sigma x mult.obs x individuals

  firstmomentb <- lapply(1:individuals, function(i) {
  if(kk[i]!=(mult.obs-1)) sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,knot[[i]]]) else 
    sigmab[[i]]%*%as.matrix(moments1[[i]]-mean.ynew[i,-k[[i]]])}) 
  ### list: dimension.sigma x 1 x individuals
  ### when dimension.sigma=1: vector : individuals 
  firstmomentb=t(do.call(cbind,firstmomentb))
  
  
  #########################################################################################
  #########################################################################################
  #########################################################################################
  ################ log-likelihood #########################################################
  #########################################################################################
  #########################################################################################
  
  ordinal.part <- sapply(1:individuals, function(i) log(pmvnorm(mean=mean.ynew[i,knot[[i]]],
                                                             sigma=matrix(variance.ynew[knot[[i]],knot[[i]],i], ncol=mult.obs-kk[i]),
                                                             lower=trunc.lower[i,knot[[i]]], upper=trunc.upper[i,knot[[i]]])),
                      simplify="array")
  
  loglikelihood <- sum(ordinal.part)
  
  
  ###################################################################################
  ############## AIC and BIC ########################################################
  ###################################################################################
  AIC <- -2*loglikelihood+2*(length(betanew)+length(deltanew)+choose(dimension.sigma+1,2))
  
  BIC <- -2*loglikelihood+log(length(which(!is.na(data.ordinal))))*(length(betanew)+length(deltanew)+choose(dimension.sigma+1,2))
  
  
  } else {
  firstmomentb=NULL
  loglikelihood=NULL
  AIC=NULL
  BIC=NULL
  } ### if (additional) 

 rownames(sigma.rand.new) <- paste("sigma ", 1:dimension.sigma, ".", sep="")
 colnames(sigma.rand.new) <- paste("sigma .", 1:dimension.sigma, sep="")
  #########################################################################################
  #########################################################################################
  list(Sigma.rand.effects=sigma.rand.new,
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
       epsilon=epsilon)
  
} #function ecm.one.ordinal


standard.error.bootstrap.one.ordinal=function(x, bootstrap.samples, 
                                              doParallel, cores=NULL) {
  beta.estimate=x$regression.coefficients
  delta.estimates=x$differences.in.thresholds
  sigma.rand.estimate=x$Sigma.rand.effects
  
  start.values.beta=beta.estimate
  start.values.delta=delta.estimates
  start.values.sigma.rand=sigma.rand.estimate
  
  l=dim(x$data.ordinal)[1]
  mult.obs=ncol(x$data.ordinal)
  miss=is.na(x$data.ordinal)
  
  if(doParallel) {if (length(cores)>0) registerDoParallel(cores=cores) else registerDoParallel()
                  boot=foreach(i=1:bootstrap.samples, .packages=c("MASS","EMcorrProbit","tmvtnorm")) %dopar% {
                    
                    random.int=mvrnorm(n=l, mu=rep(0,ncol(sigma.rand.estimate)), Sigma=sigma.rand.estimate)
                    
                    pred.fixed=sapply(1:mult.obs,function(obs) as.matrix(x$predictors.fixed[,,obs])%*%beta.estimate, simplify=T)
                    if(ncol(sigma.rand.estimate)==1) pred.rand=sapply(1:mult.obs,function(obs) x$predictors.random[,,obs]*random.int,simplify=T) else
                      pred.rand=sapply(1:mult.obs,function(obs) apply(x$predictors.random[,,obs]*random.int,1,sum),simplify=T)
                    
                    y1=pred.fixed+pred.rand+matrix(rnorm(l*mult.obs),ncol=mult.obs)
                    
                    data.ordinal.new=matrix(cut(y1,c(min(y1)-1,x$thresholds,max(y1)+1), labels=F),ncol=mult.obs)
                    data.ordinal.new=ifelse(miss==T, NA, data.ordinal.new)
                    
                    ecm.one.ordinal(data.ordinal.new,x$predictors.fixed,
                                    x$predictors.random,start.values.beta,
                                    start.values.delta,start.values.sigma.rand,
                                    exact=x$exact,montecarlo=x$montecarlo,epsilon=x$epsilon*50, additional=F)
                  } #foreach
  } else {boot=list(NA)
          
          for(i in 1:bootstrap.samples) {
            random.int=mvrnorm(n=l, mu=rep(0,ncol(sigma.rand.estimate)), Sigma=sigma.rand.estimate)
            
            pred.fixed=sapply(1:mult.obs,function(obs) as.matrix(x$predictors.fixed[,,obs])%*%beta.estimate, simplify=T)
            if(ncol(sigma.rand.estimate)==1) pred.rand=sapply(1:mult.obs,function(obs) x$predictors.random[,,obs]*random.int, simplify=T) else
              pred.rand=sapply(1:mult.obs,function(obs) apply(x$predictors.random[,,obs]*random.int,1,sum),simplify="array")
            
            y1=pred.fixed+pred.rand+matrix(rnorm(l*mult.obs),ncol=mult.obs)
            
            data.ordinal.new=matrix(cut(y1,c(min(y1)-1,x$thresholds,max(y1)+1), labels=F),ncol=mult.obs)
            data.ordinal.new=ifelse(miss==T, NA, data.ordinal.new)
            
            boot[[i]]=ecm.one.ordinal(data.ordinal.new,x$predictors.fixed,
                                      x$predictors.random,start.values.beta,
                                      start.values.delta,start.values.sigma.rand,
                                      exact=x$exact,montecarlo=x$montecarlo,epsilon=x$epsilon*50, additional=F)
            
          }  #for
  } #else  
  est=t(sapply(1:bootstrap.samples, function(i) c(boot[[i]][[1]], boot[[i]][[2]],boot[[i]][[3]]), simplify=T))
  res=var(est)
  colnames(res)=c(paste("sigma ",outer(1:dim(sigma.rand.estimate),1:dim(sigma.rand.estimate),function(x,y) paste(x,y,sep="")),sep=""),
                  paste("Predictor",1:length(beta.estimate)), paste("Diff",1:length(delta.estimates)))
  rownames(res)=c(paste("sigma ",outer(1:dim(sigma.rand.estimate),1:dim(sigma.rand.estimate),function(x,y) paste(x,y,sep="")),sep=""),
                  paste("Predictor",1:length(beta.estimate)), paste("Diff",1:length(delta.estimates)))
  res
} # function standard error




