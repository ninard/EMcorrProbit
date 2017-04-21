
###model: oneord, ordcont
emcorrprobit <- function(model, y, xfixed, xrand, start.values.beta, 
                         start.values.delta=NULL,  start.values.sigma.rand, 
                         start.values.sigma22, start.values.lambda,
                         exact, montecarlo, epsilon, ...)
{
  obj <- list(model = model, y=y, xfixed=xfixed, xrand=xrand, 
              start.values.beta=start.values.beta, 
              start.values.delta=start.values.delta,  
              start.values.sigma.rand=start.values.sigma.rand, 
              start.values.sigma22 = start.values.sigma22, 
              start.values.lambda = start.values.lambda,
              exact=exact, montecarlo=montecarlo, epsilon=epsilon, ...)
  obj$call <-match.call()
  class(obj) <- c(model)
  emcorrprobitFit(obj)
}

emcorrprobitFit <- function(obj, ...) UseMethod("emcorrprobitFit", obj)

emcorrprobitFit.default <- function(obj, ...) stop("incorrect or no model specified")

emcorrprobitFit.oneord <- function(obj, ...)
{
  #cat("Please, be patient, the function is working ... \n")
  
  y <- obj$y
  xfixed <- obj$xfixed
  xrand <-obj$xrand
  exact <- obj$exact
  montecarlo <-obj$montecarlo
  start.values.beta <- obj$start.values.beta 
  start.values.delta <-obj$start.values.delta
  start.values.sigma.rand <- obj$start.values.sigma.rand
  epsilon <- obj$epsilon
  
  #######################
  ## check
  ############################
  if(!is.matrix(y)) warning("y must be a matrix")
  if(!is.matrix(start.values.sigma.rand)) warning("start.values.sigma.rand must be a matrix")
  if(any(sort(unique(c(y)))!=1:max(y,na.rm=TRUE))) stop("Missing levels in the response variable.")
  if(length(start.values.delta)!=(max(y,na.rm=TRUE)-2)) stop("Incorrect dimension of start.values.delta.")
  if(!exact & is.na(montecarlo)) stop("montecarlo parameter undefined")
  
  xfixed <- as.array(xfixed)
  xrand <- as.array(xrand)
  y <- as.array(y)
  
  if (dim(xfixed)[1] != dim(xrand)[1] | dim(xrand)[1] != dim(y)[1]) stop("Incompatible dimensions!")
  if (dim(xfixed)[3] != dim(xrand)[3] | dim(xrand)[3] != dim(y)[2]) stop("Incompatible dimensions!")
  
  est <- ecm.one.ordinal(y,xfixed,xrand,start.values.beta, start.values.delta,
                         start.values.sigma.rand,
                         exact,montecarlo,epsilon, additional)
  
  
  est$call <- obj$call
  
  # cat("It's ready! Use print or summary to see the result.\n")
  
  class(est) <- c("emcorrprobit")
  est
}

emcorrprobitFit.ordcont <- function(obj, ...)
{
  ### y should be list(data.ordinal, data.continuous) and the same for xfixed and xrand
  data.ordinal <- obj$y[[1]]
  data.continuous <- obj$y[[2]]
  predictors.fixed.ordinal <- obj$xfixed[[1]]
  predictors.fixed.continous <- obj$xfixed[[2]]
  predictors.random.ordinal <-obj$xrand[[1]]
  predictors.random.continuous <-obj$xrand[[2]]
  exact <- obj$exact
  montecarlo <-obj$montecarlo
  start.values.beta.ordinal <- obj$start.values.beta[[1]] 
  start.values.beta.continuous <- obj$start.values.beta[[2]] 
  start.values.delta <-obj$start.values.delta
  start.values.sigma.rand <- obj$start.values.sigma.rand
  start.values.sigma22 <- obj$start.values.sigma22
  start.values.lambda <- obj$start.values.lambda
  epsilon <- obj$epsilon
  additional <-obj$additional
  
  #######################
  ## check
  ############################
  if(class(data.ordinal)!="matrix") warning("data.ordinal must be a matrix")
  if(class(start.values.sigma.rand)!="matrix") warning("start.values.sigma.rand must be a matrix")
  if(any(sort(unique(c(data.ordinal)))!=1:max(data.ordinal,na.rm=T))) stop("Missing levels in the response variable")
  if(length(start.values.delta)!=max(data.ordinal,na.rm=T)-2) stop("Incorrect dimension of start.values.delta")
  #if(exact==F & !exists("montecarlo")) {stop("Montecarlo parameter isn't defined, taken by default 100")
  #montecarlo=100
  #}
  
  
  est <- ecm.ord.plus.cont(data.ordinal,data.continuous,predictors.fixed.ordinal,
                           predictors.fixed.continous,predictors.random.ordinal,
                           predictors.random.continuous, 
                           start.values.beta.ordinal,start.values.beta.continuous,
                           start.values.delta,start.values.sigma.rand,
                           start.values.sigma22, start.values.lambda,
                           exact,montecarlo, additional,
                           epsilon)
  
  
  est$call <- obj$call
  
  # cat("It's ready! Use print or summary to see the result.\n")
  
  class(est) <- c("emcorrprobit")
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
  print(x$regression.coefficients)
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

summary.emcorrprobit <- function(object, ...)
{ #cat(" Please, be very patient ... \n")
  vcov <- standard.error.bootstrap.one.ordinal(object, ...) 
  
  se <- sqrt(diag(vcov))
  
  se.sigma <- matrix(se[1:length(object$Sigma.rand.effects)],
                     ncol=sqrt(length(object$Sigma.rand.effects)))
  
  se.sigma <- se.sigma[lower.tri(se.sigma,diag=TRUE)]
  
  TAB.sigma <- cbind(Sigma= object$Sigma.rand.effects[lower.tri(object$Sigma.rand.effects,diag=TRUE)], 
                    StdErr = se.sigma,
                    z.score=object$Sigma.rand.effects[lower.tri(object$Sigma.rand.effects,diag=TRUE)]/se.sigma)
  nam <- matrix(paste("sigma ", outer(1:ncol(object$Sigma.rand.effects),1:ncol(object$Sigma.rand.effects),
                                  function(x,y) paste(x,y,sep="")),sep=""),ncol=ncol(object$Sigma.rand.effects))
  nam <- nam[lower.tri(nam,diag=TRUE)]
  rownames(TAB.sigma)=nam
  
  
  se.regr.coeff <- se[(length(object$Sigma.rand.effects)+1):
                        (length(object$Sigma.rand.effects)+length(object$regression.coefficients))]
  
  TAB.regr.coeff <- cbind(Regression.coeff= object$regression.coefficients, 
               StdErr = se.regr.coeff,
               z.score = object$regression.coefficients/se.regr.coeff,
               p.value = 2*pnorm(abs(object$regression.coefficients/se.regr.coeff), lower.tail=FALSE))
  rownames(TAB.regr.coeff)=paste("Predictor",1:length(object$regression.coefficients))
    
  if(length(object$differences.in.thresholds)>0)
    
  {se.diff.thresholds =se[(length(se)-length(object$differences.in.thresholds)+1):length(se)]
  
  TAB.diff.thresholds <- cbind(Threshold.differences= object$differences.in.thresholds, 
               StdErr = se.diff.thresholds,
               z.score = object$differences.in.thresholds/se.diff.thresholds,
               p.value = 2*pnorm(abs(object$differences.in.thresholds/se.diff.thresholds), lower.tail=FALSE))
  rownames(TAB.diff.thresholds)=paste("Diff",1:length(object$differences.in.thresholds))
  } else TAB.diff.thresholds=NULL
  
  res <- list(call = object$call, 
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
  printCoefmat(x$regr.coefficients,P.values=TRUE, has.Pvalue=TRUE)
  cat("\n")
  printCoefmat(x$diff.thresholds,P.values=TRUE, has.Pvalue=TRUE)
  
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

ecm.one.ordinal <- function(data.ordinal,predictors.fixed,predictors.random,
                            start.values.beta,start.values.delta,
                            start.values.sigma.rand,
                         exact=FALSE,montecarlo=100,epsilon=0.001, additional=FALSE) 
  {
  ########################################
  ### wide format of longitudinal data: first level denoted by 1, second - 2 and so on
  ### predictors.fixed is a 3-dimentional array: individuals x dimentions predictros fixed x time points 
  ### predictors.random is a 3-dimentional array: individuals x dimentions predictros fixed x time points  
  ### sigma.rand is the covaraince matrix of the random effects
  ##################################
  
  ##############################################################
  
  betanew <- start.values.beta
  deltanew <- start.values.delta
  sigma.rand.new <- start.values.sigma.rand
  
  num.categories <- max(data.ordinal,na.rm=TRUE)
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
  
  z <- sapply(1:individuals, function(i) matrix(t(predictors.random[i,,]), ncol=dimension.sigma,byrow=FALSE), 
           simplify="array")
  ## mult.obs x dimension.sigma x individuals
  
  ### definition na k, knot i kk
  k <- lapply(1:individuals, function(i) which(is.na(data.ordinal[i,])))
  ### list: missings for each individual
  knot <- lapply(1:individuals, function(i) which(!is.na(data.ordinal[i,])))
  ### list: observed for each individual
  kk <- sapply(1:individuals, function(i) length(k[[i]]), simplify=TRUE)
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
    if(!exact) {
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
    
    
    if (dimension.sigma>1)  pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=TRUE))       else 
      pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=TRUE))  
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
    
    
    if (dimension.sigma>1)  pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=TRUE))       else 
      pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=TRUE))  
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
        
        a <- sum(first[data.ordinal==(delta.index+1)],na.rm=TRUE)+sum(n[(delta.index+2):num.categories])
        b <- sum((second1+second2new)[data.ordinal==delta.index+1], na.rm=TRUE)-sum(third[data.ordinal>delta.index+1], na.rm=TRUE)
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
        
        
        if (dimension.sigma>1)  pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=TRUE))       else 
          pred.random <- t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=TRUE))  
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
    }, simplify=TRUE)
    
    
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
  if(!exact) {
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


standard.error.bootstrap.one.ordinal=function(x, bootstrap.samples = 50, epsilon=NULL,
                                              doParallel = FALSE, cores=NULL) {

  beta.estimate=x$regression.coefficients
  delta.estimates=x$differences.in.thresholds
  sigma.rand.estimate=x$Sigma.rand.effects
  
  start.values.beta=beta.estimate
  start.values.delta=delta.estimates
  start.values.sigma.rand=sigma.rand.estimate
  
  l=dim(x$data.ordinal)[1]
  mult.obs=ncol(x$data.ordinal)
  miss=is.na(x$data.ordinal)
  if (length(epsilon)==0) epsilon=x$epsilon*10
  
  if(doParallel) {if (length(cores)>0) registerDoParallel(cores=cores) else registerDoParallel()
                  boot=foreach(i=1:bootstrap.samples, .packages=c("MASS","EMcorrProbit","tmvtnorm")) %dopar% {
                    
                    random.int=mvrnorm(n=l, mu=rep(0,ncol(sigma.rand.estimate)), Sigma=sigma.rand.estimate)
                    
                    pred.fixed=sapply(1:mult.obs,function(obs) as.matrix(x$predictors.fixed[,,obs])%*%beta.estimate, simplify=TRUE)
                    if(ncol(sigma.rand.estimate)==1) pred.rand=sapply(1:mult.obs,function(obs) x$predictors.random[,,obs]*random.int,simplify=TRUE) else
                      pred.rand=sapply(1:mult.obs,function(obs) apply(x$predictors.random[,,obs]*random.int,1,sum),simplify=TRUE)
                    
                    y1=pred.fixed+pred.rand+matrix(rnorm(l*mult.obs),ncol=mult.obs)
                    
                    data.ordinal.new=matrix(cut(y1,c(min(y1)-1,x$thresholds,max(y1)+1), labels=FALSE),ncol=mult.obs)
                    data.ordinal.new=ifelse(miss, NA, data.ordinal.new)
                    
                    ecm.one.ordinal(data.ordinal.new,x$predictors.fixed,
                                    x$predictors.random,start.values.beta,
                                    start.values.delta,start.values.sigma.rand,
                                    exact=x$exact,montecarlo=x$montecarlo,epsilon=epsilon, additional=FALSE)
                  } #foreach
  } else {boot=list(NA)
          
          for(i in 1:bootstrap.samples) {
            random.int=mvrnorm(n=l, mu=rep(0,ncol(sigma.rand.estimate)), Sigma=sigma.rand.estimate)
            
            pred.fixed=sapply(1:mult.obs,function(obs) as.matrix(x$predictors.fixed[,,obs])%*%beta.estimate, simplify=TRUE)
            if(ncol(sigma.rand.estimate)==1) pred.rand=sapply(1:mult.obs,function(obs) x$predictors.random[,,obs]*random.int, simplify=TRUE) else
              pred.rand=sapply(1:mult.obs,function(obs) apply(x$predictors.random[,,obs]*random.int,1,sum),simplify="array")
            
            y1=pred.fixed+pred.rand+matrix(rnorm(l*mult.obs),ncol=mult.obs)
            
            data.ordinal.new=matrix(cut(y1,c(min(y1)-1,x$thresholds,max(y1)+1), labels=FALSE),ncol=mult.obs)
            data.ordinal.new=ifelse(miss, NA, data.ordinal.new)
            
            boot[[i]]=ecm.one.ordinal(data.ordinal.new,x$predictors.fixed,
                                      x$predictors.random,start.values.beta,
                                      start.values.delta,start.values.sigma.rand,
                                      exact=x$exact,montecarlo=x$montecarlo,epsilon=epsilon, additional=FALSE)
            
          }  #for
  } #else  
  est=t(sapply(1:bootstrap.samples, function(i) 
    c(boot[[i]][[1]][lower.tri(boot[[i]][[1]],diag=TRUE)], boot[[i]][[2]],boot[[i]][[3]]), simplify=TRUE))
  res=var(est)
  nam <- matrix(paste("sigma ", outer(1:ncol(sigma.rand.estimate),1:ncol(sigma.rand.estimate),
                                      function(x,y) paste(x,y,sep="")),sep=""),ncol=ncol(sigma.rand.estimate))
  nam <- nam[lower.tri(nam,diag=TRUE)]
  colnames(res)=c(nam, paste("Predictor",1:length(beta.estimate)), 
                  paste("Diff",1:length(delta.estimates)))
  rownames(res)=c(nam, paste("Predictor",1:length(beta.estimate)), 
                  paste("Diff",1:length(delta.estimates)))
  res
} # function standard error



# ecm.ord.plus.cont <- function(data.ordinal,
#                               data.continuous,
#                               predictors.fixed.ordinal,
#                               predictors.fixed.continous,
#                               predictors.random.ordinal,
#                               predictors.random.continuous,
#                               start.values.beta.ordinal,
#                               start.values.beta.continuous,
#                               start.values.delta,
#                               start.values.sigma.rand,
#                               start.values.sigma22, 
#                               start.values.lambda,
#                               exact=F, montecarlo=100, epsilon=0.001,
#                               additional=FALSE)
# {
  
#} #function ecm.ord.plus.cont

#################################################################
condNormal <- function(x.given, mu, sigma, given.ind, req.ind){
  # Returns conditional mean and variance of x[req.ind] 
  # Given x[given.ind] = x.given
  # where X is multivariate Normal with
  # mean = mu and covariance = sigma
  # 
  B <- sigma[req.ind, req.ind]
  C <- sigma[req.ind, given.ind, drop=FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% solve(D)
  cMu <- c(mu[req.ind] + CDinv %*% (x.given - mu[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean=cMu, condVar=cVar)
}




#################################################
######################   ecm algorithm
########################################################
ecm.ord.plus.cont <- function(data.ordinal,
                              data.continuous,
                              predictors.fixed.ordinal,
                              predictors.fixed.continous,
                              predictors.random.ordinal,
                              predictors.random.continuous,
                              start.values.beta.ordinal,
                              start.values.beta.continuous,
                              start.values.delta,
                              start.values.sigma.rand,
                              start.values.sigma22, 
                              start.values.lambda,
                              exact=FALSE, montecarlo=100, epsilon=0.001,
                              additional=FALSE) {  
  ########################################
  ### sigma.rand is the covaraince matrix of the random effects
  ##################################
  
 
  ################################################################3
  
  
  
  ### definition of the new estimates
  betanew.ordinal <- start.values.beta.ordinal
  betanew.continuous <- start.values.beta.continuous
  deltanew <- start.values.delta
  sigma.rand.new <- start.values.sigma.rand
  sigma22new <- start.values.sigma22
  lambdanew <- start.values.lambda
  
  num.categories <- max(data.ordinal,na.rm=TRUE)
  if(num.categories==2) deltanew <- NULL
  
  new.est <- c(start.values.sigma.rand,start.values.beta.ordinal,
               start.values.beta.continuous,start.values.delta,
               start.values.lambda,start.values.sigma22)
  old.est <- new.est+1.5*epsilon
  
  number.it <- 0
  
  #############################
  ## this is necessary
  ##############################
  data.continuousA=data.continuous
  data.continuousA[is.na(data.continuous)]=0
  
  # definition of lower and upper boundary of the truncated normal distr. 
  # (the new variable given observed data)
  trunc.lower <- ifelse(data.ordinal==1,-Inf,0)
  trunc.upper <- ifelse(data.ordinal==1,0,ifelse(data.ordinal<num.categories,1,Inf))
  
  # number of observations in each category
  n <- table(data.ordinal)
  individuals <- length(data.ordinal[,1])
  mult.obs <- dim(predictors.fixed.ordinal)[3]
  
  dimension.sigma <- ncol(sigma.rand.new)
  p <- 0
  p[1] <- dim(predictors.random.ordinal)[2]
  #the dimension of the random effects for the ordinal response
  p[2] <- dim(predictors.random.continuous)[2]
  #the dimension of the random effects for the continuous response
  
  # definition of the matrices for the predictors for the random effects
  z.ordinal <- sapply(1:individuals, 
                      function(i) matrix(t(predictors.random.ordinal[i,,]), 
                                         ncol=p[1],byrow=FALSE), 
                      simplify="array")
  ## mult.obs x p[1] x individuals
  z.continuous <- sapply(1:individuals, 
                         function(i) matrix(t(predictors.random.continuous[i,,]), 
                                            ncol=p[2],byrow=FALSE), 
                         simplify="array")
  ## mult.obs x p[2] x individuals
  z <- sapply(1:individuals, 
              function(i) rbind(cbind(z.ordinal[,,i],matrix(0,ncol=p[2],nrow=mult.obs)),
                                cbind(matrix(0,ncol=p[1],nrow=mult.obs),z.continuous[,,i])), 
              simplify="array")
  ### (2xmult.obs) x dimension.sigma x individuals
  
  ### definitions of k, knot i kk
  # missing values
  k.ordinal <- lapply(1:individuals, function(i) which(is.na(data.ordinal[i,])))
  # observed values
  knot.ordinal <- lapply(1:individuals, function(i) which(!is.na(data.ordinal[i,])))
  # how much are the missing values
  kk.ordinal <- sapply(1:individuals, function(i) length(k.ordinal[[i]]), simplify=TRUE)
  # analogous for the continuous response
  k.continuous <- lapply(1:individuals, function(i) which(is.na(data.continuous[i,])))
  knot.continuous <- lapply(1:individuals, function(i) which(!is.na(data.continuous[i,])))
  kk.continuous <- sapply(1:individuals, function(i) length(k.continuous[[i]]), simplify=TRUE)
  
  ######################################################################
  ############             while loop           ##########################
  ##########################################################################
  while(any(abs(new.est-old.est)>epsilon)) {
    
    number.it <- number.it+1
    old.est <- new.est
    
    #######################################################################
    delta.exp <- c(1,deltanew,1)
    delta <- matrix(delta.exp[data.ordinal], ncol=mult.obs)
    alpha.exp <- c(0,0,cumsum(deltanew))
    alpha <- matrix(alpha.exp[data.ordinal], ncol=mult.obs)
    ## NA where data.ordinal is NA
    
    #############################################
    ### marginal expectation of ynew #####
    ####################################
    if(length(start.values.beta.ordinal)>1) {
      predyneword <- sapply(1:mult.obs, function(i) 
        (predictors.fixed.ordinal[,,i]%*%betanew.ordinal), 
        simplify=TRUE)}  else 
          predyneword <- betanew.ordinal
    ##  matrix: individuals x mult.obs
    
    mean.yneword <- (predyneword-alpha)/delta
    ## NA where data.ordinal is NA
    
    #### needed - for the expectation of b_1ib_2i
    mean.ynewordA <- ifelse(is.na(data.ordinal),0,mean.yneword)
    ## 0 where data.ordinal is NA
    
    if(length(start.values.beta.continuous)>1) {
      mean.ycontinuous <- sapply(1:mult.obs, function(i) 
        (predictors.fixed.continuous[,,i]%*%betanew.continuous), 
        simplify=TRUE)}  else 
          mean.ycontinuous <- array(betanew.continuous, dim=c(individuals, mult.obs))
    
    mean.ycontinuous[is.na(data.continuous)] <- NA
    ##  matrix: observations x mult.obs
    
    ######## needed 
    mean.ycontinuousA <- ifelse(is.na(data.continuous),0,mean.ycontinuous)
    ## 0 where data.continuous is NA
    
    mean.ynew <- cbind(mean.yneword,mean.ycontinuous)
    ## observations x 2*mult.obs
    ######## needed 
    mean.ynewA <- cbind(mean.ynewordA,mean.ycontinuousA)
    
    
    varord <- sapply(1:individuals, function(i)
      (z.ordinal[,,i]%*%
         matrix(sigma.rand.new[1:p[1],1:p[1]],ncol=p[1])%*%
         t(z.ordinal[,,i])+
         diag(mult.obs)*(1+lambdanew^2*sigma22new))/(delta[i,]%*%t(delta[i,])),
      simplify="array")
    ## if delta is NA all calculations with it are NA
    # mult.obs x mult.obs x individuals
    
    varordcont <- sapply(1:individuals, function(i)
      (z.ordinal[,,i]%*%
         matrix(sigma.rand.new[1:p[1],(p[1]+1):dimension.sigma],ncol=p[2])%*%
         t(z.continuous[,,i])+
         diag(mult.obs)*(lambdanew*sigma22new))/(delta[i,]),
      simplify="array")
    # NA for the rows related to NA values in the ordinal response (because of delta)
    ### NA for the columns related to NA values in the continuous variable
    for(i in 1:individuals) varordcont[,k.continuous[[i]],i] <- NA
    
    varcont=sapply(1:individuals, function(i)
      (z.continuous[,,i]%*%
         matrix(sigma.rand.new[(p[1]+1):dimension.sigma,(p[1]+1):dimension.sigma],ncol=p[2])%*%
         t(z.continuous[,,i])+
         diag(mult.obs)*sigma22new),
      simplify="array")
    # mult.obs x mult.obs x individuals
    # NA rows and columns for the missing values
    for(i in 1:individuals) {
      varcont[,k.continuous[[i]],i] <- NA
      varcont[k.continuous[[i]],,i] <- NA
    }
    
    variance.ynew <- sapply(1:individuals, function(i)
      rbind(cbind(varord[,,i],varordcont[,,i]),
            cbind(t(varordcont[,,i]),varcont[,,i])),
      simplify="array")
    ## 2*mult.obs x 2*mult.obs x observations
    
    ##################### conditional distribution of y latent given observed normal
    ################################################################3
    latent <- lapply(1:individuals, function(i) 
      condNormal(x=data.continuous[i,knot.continuous[[i]]], 
                 mu=mean.ynew[i,!is.na(mean.ynew[i,])], 
                 sigma=matrix(variance.ynew[,,i][!is.na(variance.ynew[,,i])], 
                              ncol=length(c(knot.ordinal[[i]],knot.continuous[[i]]))), 
                 req=1:length(knot.ordinal[[i]]), 
                 given=(length(knot.ordinal[[i]])+1):
                   (2*mult.obs-kk.ordinal[i]-kk.continuous[i])))
    
    
    if(exact) { ### computation of exact moments
      simulations.exact <- lapply(1:individuals, function(i) 
        mtmvnorm(mean=latent[[i]]$condMean, 
                 sigma=matrix(as.vector(latent[[i]]$condVar),ncol=mult.obs-kk.ordinal[i]),
                 ##########################################
                 ## explicitly to define sigma as matrix
                 ##########################################
                 lower=trunc.lower[i,knot.ordinal[[i]]],
                 upper=trunc.upper[i,knot.ordinal[[i]]]))
      ##  mult.obs x individuals
      moments1 <- lapply(1:individuals, function(i) simulations.exact[[i]]$tmean)
      moments2 <- lapply(1:individuals, function(i) simulations.exact[[i]]$tvar)
    } else     {### Monte Carlo computation of moments
      simulations <- lapply(1:individuals, function(i) 
        rtmvnorm(montecarlo, mean=latent[[i]]$condMean, 
                 sigma=matrix(as.vector(latent[[i]]$condVar),ncol=mult.obs-kk.ordinal[i]),
                 ##########################################
                 ## explicitly to define sigma as matrix
                 ########################################
                 lower=trunc.lower[i,knot.ordinal[[i]]],
                 upper=trunc.upper[i,knot.ordinal[[i]]], algorithm="gibbs"))
      ##### Sometimes problems with the generation of random numbers
      lipsi=numeric(0)
      for(i in 1:individuals) if(any(is.nan(simulations[[i]]))) lipsi=c(lipsi,i)
      for(i in lipsi)
        for(j in 1:length(knot.ordinal[[i]])) 
          simulations[[i]][,j] <- rtmvnorm(montecarlo, mean=latent[[i]]$condMean[j], 
                                           sigma=matrix(as.vector(latent[[i]]$condVar[j,j]),ncol=1),
                                           lower=trunc.lower[i,knot.ordinal[[i]]][j],
                                           upper=trunc.upper[i,knot.ordinal[[i]]][j], 
                                           algorithm="gibbs")
      ## montecarlo x mult.obs x observations
      
      moments1 <- lapply(1:individuals, function(i) if(kk.ordinal[i]<mult.obs-1) 
        apply(simulations[[i]],2,mean) else mean(simulations[[i]]))
      moments2 <- lapply(simulations, var)
    }
    
    
    ## I can not work with NA values only
    #moments1A <- list(0)
    moments1Aarray <- array(0, dim=c(individuals, mult.obs))
    for(i in 1:individuals) {#moments1A[[i]] <- rep(0, mult.obs) 
      #moments1A[[i]][knot.ordinal[[i]]] <- moments1[[i]]
      moments1Aarray[i,knot.ordinal[[i]]] <- moments1[[i]]
    }
    
    
    moments2A <- list(0)
    for(i in 1:individuals) {
      moments2A[[i]] <- array(0, dim=c(mult.obs, mult.obs))
      moments2A[[i]][knot.ordinal[[i]],knot.ordinal[[i]]] <- moments2[[i]]
    }
    ### list - to make it an array - not needed
    
    ########## to calculate the expectation of the random effects
    covyordinalb <- sapply(1:individuals, function(i) 
      matrix(z.ordinal[,,i],ncol=p[1])%*%sigma.rand.new[1:p[1],]/delta[i,],simplify="array")
    ## mult.obs x dimension.sigma x individuals
    ## NA rows for unobserved
    
    covycontb <- sapply(1:individuals, function(i) 
      matrix(z.continuous[,,i], ncol=p[2])%*%sigma.rand.new[(p[1]+1):dimension.sigma,],simplify="array")
    # mult.obs x dim.sigma x individuals
    for(i in 1:individuals) covycontb[k.continuous[[i]],,i]=NA
    ## NA rows for unobserved
    
    covynewb <- sapply(1:individuals, function(i)
      rbind(covyordinalb[,,i],covycontb[,,i]), simplify="array")
    #2*mult.obs x dimension.sigma x individuals
    
    ######## needed - for sigma !!!!!!!!
    covynewbA <- covynewb
    for(i in 1:individuals) covynewbA[,,i][is.na(covynewb[,,i])]=0
    ## 0 rows for unobserved
    
    
    sigmab <- lapply(1:individuals, 
                     function(i) t(matrix(covynewb[,,i][!is.na(covynewb[,,i])], ncol=dimension.sigma))%*%
                       solve(matrix(variance.ynew[,,i][!is.na(variance.ynew[,,i])],
                                    ncol=length(c(knot.ordinal[[i]],knot.continuous[[i]]))))
    )
    # dimension.sigma x obs.ord+obs.cont x individuals
    
    sigmabA=list(NA)
    for(i in 1:individuals) {sigmabA[[i]]=array(0, dim=c(dimension.sigma, 2*mult.obs))
    sigmabA[[i]][,c(knot.ordinal[[i]],mult.obs+knot.continuous[[i]])]=sigmab[[i]]
    }
    ## 0 columns for NA values
    
    firstmomentb <- t(sapply(1:individuals, 
                             function(i) sigmab[[i]]%*%
                               (c(moments1[[i]],data.continuous[i,knot.continuous[[i]]])-
                                  mean.ynew[i,!is.na(mean.ynew[i,])]), simplify=TRUE))
    ####  individuals x dimension.sigma 
    ################################################
    if (p[1]>1)  randompart.ordinal <- t(sapply(1:individuals, function(i) 
      z.ordinal[,,i]%*%firstmomentb[i,1:p[1]], simplify=TRUE))   else 
        randompart.ordinal <- t(sapply(1:individuals, function(i) 
          z.ordinal[,,i]%*%t(firstmomentb[i,1:p[1]]),simplify=TRUE))  
    ### individuals x mult.obs
    
    if (p[2]>1)  randompart.continuous <- t(sapply(1:individuals, function(i) 
      z.continuous[,,i]%*%firstmomentb[i,(p[1]+1):dimension.sigma], simplify=TRUE))   else 
        randompart.continuous <- t(sapply(1:individuals, function(i) 
          z.continuous[,,i]%*%t(firstmomentb[i,(p[1]+1):dimension.sigma]),simplify=TRUE))  
    ### individuals x mult.obs
    
    #   conditional expectation of the mean of y observed normal on random effects given observed data
    mu2 <- mean.ycontinuous+randompart.continuous
    # matrix: individuals x obs.cont 
    ###? to define ???? mu2A
    
    
    
    #################################################################################################
    ###############################                        sigma.rand estimate
    #################################################################################################
    # the variance of ynew given observed data
    
    zero <- lapply(1:individuals, function(i)  {
      a <- matrix(0,ncol=2*mult.obs,nrow=2*mult.obs) 
      a[knot.ordinal[[i]],knot.ordinal[[i]]] <- moments2[[i]]
      a}  )
    
    sigma.rand <- sapply(1:individuals, function(i)
      sigma.rand.new-sigmabA[[i]]%*%covynewbA[,,i]+sigmabA[[i]]%*%
        (zero[[i]]+c(moments1Aarray[i,],data.continuousA[i,])%*%
           t(c(moments1Aarray[i,],data.continuousA[i,]))-
           c(moments1Aarray[i,],data.continuousA[i,])%*%t(mean.ynewA[i,])-
           mean.ynewA[i,]%*%
           t(c(moments1Aarray[i,],data.continuousA[i,]))+
           (mean.ynewA[i,])%*%t(mean.ynewA[i,]))%*%
        t(sigmabA[[i]]),
      simplify=TRUE)
    #print(which(is.na(sigma.rand)))
    
    sigma.rand.matrix <- sapply(1:individuals,
                                function(i) matrix(sigma.rand[,i], ncol=dimension.sigma),simplify="array")
    ## dimension.sigma x dimension.sigma x individuals
    
    
    #if(dimension.sigma==1) sigma.rand.new=mean(sigma.rand) else 
    sigma.rand.new <- matrix(apply(sigma.rand,1,mean),ncol=dimension.sigma)
    sigma.rand.new[lower.tri(sigma.rand.new)] <- t(sigma.rand.new)[lower.tri(sigma.rand.new)]
    
    
    ###update
    
    
    
    ##################################################################################3
    ################################ beta.ordinal estimates
    ########################################################################################
    
    ######### some 
    continuous.part <- lambdanew*(data.continuous-mu2)
    continuous.part <- ifelse(is.na(data.continuous),0,continuous.part)
    
    ywave <- delta*moments1Aarray-randompart.ordinal+alpha-continuous.part 
    
    if(length(start.values.beta.ordinal)==1) betanew.ordinal <- lm(c(ywave)~1)$coef else {
      pred.beta <- predictors.fixed.ordinal[,,1]
      for(i in 2:mult.obs) pred.beta <- rbind(pred.beta,predictors.fixed.ordinal[,,i])
      betanew.ordinal <- lm(c(ywave)~pred.beta-1)$coef 
      rm(pred.beta)}
    #print(betanew.ordinal)
    
    
    #### update
    
    ####################################################################
    ##############################        sigma22 estimate
    ########################################################################
    ### expectation of mu2 squared
    if(p[2]>1) randompart.continuous2 <- t(sapply(1:individuals, function(i) 
      diag(z.continuous[,,i]%*%
             sigma.rand.matrix[(p[1]+1):dimension.sigma,(p[1]+1):dimension.sigma,i]%*%
             t(z.continuous[,,i])), simplify=TRUE)) else
               randompart.continuous2 <- t(sapply(1:individuals, function(i) 
                 diag(z.continuous[,,i]%*%
                        matrix(sigma.rand.matrix[(p[1]+1):dimension.sigma,(p[1]+1):dimension.sigma,i])%*%
                        t(z.continuous[,,i])), simplify=TRUE))
    mu22<- mean.ycontinuous*(mean.ycontinuous+2*randompart.continuous)+randompart.continuous2
    
    sigma22new <- sum(data.continuous*(data.continuous-2*mu2)+mu22, 
                      na.rm=T)/(mult.obs*individuals-sum(kk.continuous))
    #print(sigma22new)
    
    
    ### update
    
    
    
    
    #######################################################################
    ################################################          lambda estimate
    ######################################################################
    
    meanylatent <- (predyneword+randompart.ordinal-alpha)/delta
    ## mean of conditional expectation of y latent on random effects given observed data
    
    
    ### lambda is related to the correlation between osebravations at the same time point. 
    ####### I need only paired individuals, not only on one of the variables
    
    expectation1 <- lapply(1:individuals, function(i) 
      moments2A[[i]]+moments1Aarray[i,]%*%
        t(moments1Aarray[i,])-mean.ynewordA[i,]%*%t(moments1Aarray[i,]))
    ## mult.obs x mult.obs x individuals
    
    expectation2 <- lapply(1:individuals, function(i) 
      (data.continuousA[i,]-mean.ycontinuousA[i,])%*%t(moments1Aarray[i,]))
    ## mult.obs x mult.obs x individuals
    expectation <- lapply(1:individuals, function(i)
      rbind(expectation1[[i]],expectation2[[i]]))
    ## 2mult.obs x mult.obs x individuals
    
    ## the expectation of z_i*b_i*y_1i_new'
    expectationbig <- lapply(1:individuals, function(i) 
      z[,,i]%*%sigmabA[[i]]%*%expectation[[i]])
    ## 2xmult.obs x obs.ord x individuals
    ## zero columns for NAs
    
    expectationbigb1 <- t(sapply(1:individuals, function(i)
      diag(expectationbig[[i]]), simplify=TRUE))
    
    expectationbigb2 <- t(sapply(1:individuals, function(i) 
      diag(expectationbig[[i]][(mult.obs+1):(2*mult.obs),]), simplify=TRUE))
    
    
    expectationb1ib2i <- t(sapply(1:individuals, function(i)
      diag(z.ordinal[,,i]%*%
             matrix(sigma.rand.matrix[1:p[1],(p[1]+1):dimension.sigma,i],ncol=p[2])%*%
             t(z.continuous[,,i])), simplify="array"))
    ## mult.obs x mult.obs
    ## the first column is z[,,i]*e(b1ib2i)*z[1,,i]
    ## i need only the diagonal element !!!!!!! ###
    
    lambdanew <- sum(delta*(moments1Aarray-meanylatent)*(data.continuous-mean.ycontinuous)-
                       randompart.continuous*(alpha-predyneword)+
                       expectationb1ib2i-delta*expectationbigb2, na.rm=T)/
      sum(data.continuous*(data.continuous-2*mu2)+mu22, na.rm=T)
    #print(lambdanew)
    
    ####################################################
    ##update
    ##########################################################3
    
    ##################################################################
    ###################                 delta estimates
    ################################################################## 
    randompart.continuousA <- ifelse(is.na(data.continuous),0,randompart.continuous)
    
    if(num.categories>2) {
      
      first <- t(sapply(1:individuals, function(i) 
        diag(moments2A[[i]])+(moments1Aarray[i,])^2, 
        simplify="array"))
      ## individuals x mult.obs
      ## zero for NAs
      for (k in 1:(num.categories-2)) {
        second <- moments1Aarray*(predyneword-alpha+lambdanew*
                                    (data.continuousA-mean.ycontinuousA))+
          expectationbigb1-lambdanew*expectationbigb2
        
        alpha.exp.term <- c(0,0)
        for(i in 1:length(deltanew)) alpha.exp.term[i+2] <- ifelse(k<=i, alpha.exp[i+2]-deltanew[k], 0)
        alpha.term <- matrix(alpha.exp.term[data.ordinal],ncol=mult.obs)
        
        third <- moments1Aarray*delta-predyneword-randompart.ordinal+alpha.term-
          lambdanew*(data.continuousA-mean.ycontinuousA-randompart.continuousA)
        
        
        a <- sum(first[data.ordinal==(k+1)], na.rm=T)+sum(n[(k+2):num.categories])
        
        b <- sum(second[data.ordinal==k+1], na.rm=T)-sum(third[data.ordinal>k+1], na.rm=T)
        
        c <- -n[(k+1)]
        
        deltanew[k] <- (b+sqrt(b^2-4*a*c))/(2*a) 
        #print(deltanew)
        
        ## update
        ##############################
        ########3 do tuk
        ##############################
        
      } ## delta update
      
      
      
    }
    #############################################################################
    ######################                 beta continuous
    ############################################################################
    term <- 1+lambdanew^2*sigma22new
    
    ordinal.part <- delta*moments1Aarray+alpha-predyneword-randompart.ordinal
    ordinal.part <- ifelse(is.na(data.ordinal),0,ordinal.part)
    
    sumsum <- term*(data.continuous-randompart.continuous)-
      lambdanew*sigma22new*ordinal.part
    
    ## individuals x mult.obs
    ## NA - if data.continuous is NA
    
    if (length(start.values.beta.ordinal)>1) {
      vect <- apply(sapply(1:mult.obs, function(j) { 
        apply(sumsum[,j]*predictors.fixed.continuous[,,j], 2, sum, na.rm=TRUE)},
        simplify=TRUE), 1, sum, na.rm=TRUE)  
      
      #### I need to remove some elements of m!! 
      ### those with NAs for data.continuous
      m <- apply(sapply(1:mult.obs, function(j) 
        apply(apply(predictors.fixed.continuous[,,j][!is.na(data.continuous[,j]),], 1, 
                    function(i) i%*%t(i)),1,sum), simplify="array"),1,sum)
      
      m <- term*matrix(m, ncol=length(betanew.continuous), byrow=FALSE)
      #ncol=dim(predictors.fixed.continuous)[2],byrow=FALSE)
      #m=term*m
    } else {
      ### to check what follows !!!!!!!!!!!!!!!!
      vect <- sum(sapply(1:mult.obs, function(j) { 
        sum(sumsum[,j]*predictors.fixed.continuous[,,j], na.rm=TRUE)},
        simplify=TRUE))
      
      m <-term*sum( sapply(1:mult.obs, function(j) 
        sum(predictors.fixed.continuous[,,j][!is.na(data.continuous[,j])]^2), 
        simplify=TRUE))
      
    }   
    betanew.continuous <- solve(m,vect)
    #print(betanew.continuous)
    
    new.est <- c(sigma.rand.new,betanew.ordinal,
                 betanew.continuous,deltanew,
                 lambdanew,sigma22new)
    ########################
    #cat(number.it," ")
    #cat("number.it: ",number.it,"\n")
    #cat("sigma.rand.new: ","\n")
    #print(sigma.rand.new)
    #cat(" betanew.ordinal: ",betanew.ordinal,"\n",
    #"betanew.continuous: ",betanew.continuous,"\n","deltanew: ",deltanew,"\n",
    #"lambdanew: ",lambdanew,"\n",
    #"sigma22new: ",sigma22new,"\n")
    #print(Sys.time())
    #      cat(number.it,sigma.rand.new,betanew.ordinal,betanew.continuous,deltanew,lambdanew,sigma22new, "\n")
    #write.table(c(number.it,sigma.rand.new,betanew.ordinal,betanew.continuous,deltanew,lambdanew,sigma22new), 
    #            "last-estimates-ordinal-plus-continuous.txt",row.names=F, col.names=F)
    #cat(number.it," ")
    
  } # while
  #  )## time elapsed
  
  
  
  
  if(additional) {
    ###########################################################################
    ################### random effects
    ###########################################################################
    delta.exp <- c(1,deltanew,1)
    delta <- matrix(delta.exp[data.ordinal], ncol=mult.obs)
    alpha.exp <- c(0,0,cumsum(deltanew))
    alpha <- matrix(alpha.exp[data.ordinal], ncol=mult.obs)
    ## NA where data.ordinal is NA
    
    #############################################
    ### marginal expectation of ynew #####
    ####################################
    if(length(start.values.beta.ordinal)>1) {
      predyneword <- sapply(1:mult.obs, function(i) 
        (predictors.fixed.ordinal[,,i]%*%betanew.ordinal), 
        simplify=TRUE)}  else 
          predyneword <- betanew.ordinal
    ##  matrix: individuals x mult.obs
    
    mean.yneword <- (predyneword-alpha)/delta
    ## NA where data.ordinal is NA
    
    #### needed - for the expectation of b_1ib_2i
    mean.ynewordA <- ifelse(is.na(data.ordinal),0,mean.yneword)
    ## 0 where data.ordinal is NA
    
    if(length(start.values.beta.continuous)>1) {
      mean.ycontinuous <- sapply(1:mult.obs, function(i) 
        (predictors.fixed.continuous[,,i]%*%betanew.continuous), 
        simplify=TRUE)}  else 
          mean.ycontinuous <- array(betanew.continuous, dim=c(individuals, mult.obs))
    
    mean.ycontinuous[is.na(data.continuous)] <- NA
    ##  matrix: observations x mult.obs
    
    ######## needed 
    mean.ycontinuousA <- ifelse(is.na(data.continuous),0,mean.ycontinuous)
    ## 0 where data.continuous is NA
    
    mean.ynew <- cbind(mean.yneword,mean.ycontinuous)
    ## observations x 2*mult.obs
    ######## needed 
    mean.ynewA <- cbind(mean.ynewordA,mean.ycontinuousA)
    
    
    varord <- sapply(1:individuals, function(i)
      (z.ordinal[,,i]%*%
         matrix(sigma.rand.new[1:p[1],1:p[1]],ncol=p[1])%*%
         t(z.ordinal[,,i])+
         diag(mult.obs)*(1+lambdanew^2*sigma22new))/(delta[i,]%*%t(delta[i,])),
      simplify="array")
    ## if delta is NA all calculations with it are NA
    # mult.obs x mult.obs x individuals
    
    varordcont <- sapply(1:individuals, function(i)
      (z.ordinal[,,i]%*%
         matrix(sigma.rand.new[1:p[1],(p[1]+1):dimension.sigma],ncol=p[2])%*%
         t(z.continuous[,,i])+
         diag(mult.obs)*(lambdanew*sigma22new))/(delta[i,]),
      simplify="array")
    # NA for the rows related to NA values in the ordinal response (because of delta)
    ### NA for the columns related to NA values in the continuous variable
    for(i in 1:individuals) varordcont[,k.continuous[[i]],i] <- NA
    
    varcont=sapply(1:individuals, function(i)
      (z.continuous[,,i]%*%
         matrix(sigma.rand.new[(p[1]+1):dimension.sigma,(p[1]+1):dimension.sigma],ncol=p[2])%*%
         t(z.continuous[,,i])+
         diag(mult.obs)*sigma22new),
      simplify="array")
    # mult.obs x mult.obs x individuals
    # NA rows and columns for the missing values
    for(i in 1:individuals) {
      varcont[,k.continuous[[i]],i] <- NA
      varcont[k.continuous[[i]],,i] <- NA
    }
    
    variance.ynew <- sapply(1:individuals, function(i)
      rbind(cbind(varord[,,i],varordcont[,,i]),
            cbind(t(varordcont[,,i]),varcont[,,i])),
      simplify="array")
    ## 2*mult.obs x 2*mult.obs x observations
    
    ##################### conditional distribution of y latent given observed normal
    ################################################################3
    latent <- lapply(1:individuals, function(i) 
      condNormal(x=data.continuous[i,knot.continuous[[i]]], 
                 mu=mean.ynew[i,!is.na(mean.ynew[i,])], 
                 sigma=matrix(variance.ynew[,,i][!is.na(variance.ynew[,,i])], 
                              ncol=length(c(knot.ordinal[[i]],knot.continuous[[i]]))), 
                 req=1:length(knot.ordinal[[i]]), 
                 given=(length(knot.ordinal[[i]])+1):
                   (2*mult.obs-kk.ordinal[i]-kk.continuous[i])))
    
    
    if(exact) { ### computation of exact moments
      simulations.exact <- lapply(1:individuals, function(i) 
        mtmvnorm(mean=latent[[i]]$condMean, 
                 sigma=matrix(as.vector(latent[[i]]$condVar),ncol=mult.obs-kk.ordinal[i]),
                 ##########################################
                 ## explicitly to define sigma as matrix
                 ##########################################
                 lower=trunc.lower[i,knot.ordinal[[i]]],
                 upper=trunc.upper[i,knot.ordinal[[i]]]))
      ##  mult.obs x individuals
      moments1 <- lapply(1:individuals, function(i) simulations.exact[[i]]$tmean)
      moments2 <- lapply(1:individuals, function(i) simulations.exact[[i]]$tvar)
    } else     {### Monte Carlo computation of moments
      simulations <- lapply(1:individuals, function(i) 
        rtmvnorm(montecarlo, mean=latent[[i]]$condMean, 
                 sigma=matrix(as.vector(latent[[i]]$condVar),ncol=mult.obs-kk.ordinal[i]),
                 ##########################################
                 ## explicitly to define sigma as matrix
                 ########################################
                 lower=trunc.lower[i,knot.ordinal[[i]]],
                 upper=trunc.upper[i,knot.ordinal[[i]]], algorithm="gibbs"))
      ##### Sometimes problems with the generation of random numbers
      lipsi=numeric(0)
      for(i in 1:individuals) if(any(is.nan(simulations[[i]]))) lipsi=c(lipsi,i)
      for(i in lipsi)
        for(j in 1:length(knot.ordinal[[i]])) 
          simulations[[i]][,j] <- rtmvnorm(montecarlo, mean=latent[[i]]$condMean[j], 
                                           sigma=matrix(as.vector(latent[[i]]$condVar[j,j]),ncol=1),
                                           lower=trunc.lower[i,knot.ordinal[[i]]][j],
                                           upper=trunc.upper[i,knot.ordinal[[i]]][j], 
                                           algorithm="gibbs")
      ## montecarlo x mult.obs x observations
      
      moments1 <- lapply(1:individuals, function(i) if(kk.ordinal[i]<mult.obs-1) 
        apply(simulations[[i]],2,mean) else mean(simulations[[i]]))
      moments2 <- lapply(simulations, var)
    }
    
    
    ## I can not work with NA values only
    #moments1A <- list(0)
    moments1Aarray <- array(0, dim=c(individuals, mult.obs))
    for(i in 1:individuals) {#moments1A[[i]] <- rep(0, mult.obs) 
      #moments1A[[i]][knot.ordinal[[i]]] <- moments1[[i]]
      moments1Aarray[i,knot.ordinal[[i]]] <- moments1[[i]]
    }
    
    
    moments2A <- list(0)
    for(i in 1:individuals) {
      moments2A[[i]] <- array(0, dim=c(mult.obs, mult.obs))
      moments2A[[i]][knot.ordinal[[i]],knot.ordinal[[i]]] <- moments2[[i]]
    }
    ### list - to make it an array - not needed
    
    ########## to calculate the expectation of the random effects
    covyordinalb <- sapply(1:individuals, function(i) 
      matrix(z.ordinal[,,i],ncol=p[1])%*%sigma.rand.new[1:p[1],]/delta[i,],simplify="array")
    ## mult.obs x dimension.sigma x individuals
    ## NA rows for unobserved
    
    covycontb <- sapply(1:individuals, function(i) 
      matrix(z.continuous[,,i], ncol=p[2])%*%sigma.rand.new[(p[1]+1):dimension.sigma,],simplify="array")
    # mult.obs x dim.sigma x individuals
    for(i in 1:individuals) covycontb[k.continuous[[i]],,i]=NA
    ## NA rows for unobserved
    
    covynewb <- sapply(1:individuals, function(i)
      rbind(covyordinalb[,,i],covycontb[,,i]), simplify="array")
    #2*mult.obs x dimension.sigma x individuals
    
    ######## needed - for sigma !!!!!!!!
    covynewbA <- covynewb
    for(i in 1:individuals) covynewbA[,,i][is.na(covynewb[,,i])]=0
    ## 0 rows for unobserved
    
    
    sigmab <- lapply(1:individuals, 
                     function(i) t(matrix(covynewb[,,i][!is.na(covynewb[,,i])], ncol=dimension.sigma))%*%
                       solve(matrix(variance.ynew[,,i][!is.na(variance.ynew[,,i])],
                                    ncol=length(c(knot.ordinal[[i]],knot.continuous[[i]]))))
    )
    # dimension.sigma x obs.ord+obs.cont x individuals
    
    sigmabA=list(NA)
    for(i in 1:individuals) {sigmabA[[i]]=array(0, dim=c(dimension.sigma, 2*mult.obs))
    sigmabA[[i]][,c(knot.ordinal[[i]],mult.obs+knot.continuous[[i]])]=sigmab[[i]]
    }
    ## 0 columns for NA values
    
    firstmomentb <- t(sapply(1:individuals, 
                             function(i) sigmab[[i]]%*%
                               (c(moments1[[i]],data.continuous[i,knot.continuous[[i]]])-
                                  mean.ynew[i,!is.na(mean.ynew[i,])]), simplify=TRUE))
    ####  individuals x dimension.sigma 
    #########################################################################################
    #########################################################################################
    #########################################################################################
    ################ log-likelihood #########################################################
    #########################################################################################
    #########################################################################################
    continuous.part=sapply(1:individuals, function(i) 
      log(dmvnorm(data.continuous[i,knot.continuous[[i]]], 
                  mean=mean.ycontinuous[i,knot.continuous[[i]]], 
                  sigma=matrix(varcont[,,i][!is.na(varcont[,,i])],ncol=mult.obs-kk.continuous[i])
      )), simplify=TRUE)
    
    ordinal.part=sapply(1:individuals, function(i) log(pmvnorm(mean=latent[[i]]$condMean,
                                                               sigma=matrix(latent[[i]]$condVar, ncol=mult.obs-kk.ordinal[i]),
                                                               lower=trunc.lower[i,knot.ordinal[[i]]], upper=trunc.upper[i,knot.ordinal[[i]]])),
                        simplify="array")
    
    
    loglikelihood <- sum(continuous.part)+sum(ordinal.part)
    
    
    ###################################################################################
    ############## AIC and BIC ########################################################
    ###################################################################################
    AIC=-2*loglikelihood+2*(length(start.values.beta.ordinal)+length(start.values.beta.continuous)+
                              length(start.values.delta)+2+choose(dimension.sigma+1,2))
    ###############I think that this should be not +1+3, but only +2 (lambda and sigma)
    
    BIC=-2*loglikelihood+log(length(which(!is.na(data.ordinal)))+length(which(!is.na(data.continuous))))*
      (length(start.values.beta.ordinal)+length(start.values.beta.continuous)
       +length(start.values.delta)+2+choose(dimension.sigma+1,2))
    
    
  } else {
    firstmomentb=NULL
    loglikelihood=NULL
    AIC=NULL
    BIC=NULL 
  }   ### if (additional) 
#  print(head(firstmomentb))
#  print(loglikelihood)
#  print(AIC)
#  print(BIC)
  
  
  #print(c(number.it,sigma.rand.new,betanew.ordinal,betanew.continuous,deltanew,lambdanew,sigma22new))
  #  cat("Final estimates with accuracy=",epsilon," :","\n",
  #      "number.it: ",number.it,"\n",
  #      "sigma.rand.estimate: ",sigma.rand.new,"\n",
  #      "betanew.ordinal.estimate: ",betanew.ordinal,"\n",
  #      "betanew.continuous.estimate: ",betanew.continuous,"\n",
  #      "deltanew.estimate: ",deltanew,"\n",
  #      "lambdanew.estimate: ",lambdanew,"sigma22new.estimate: ",sigma22new,"\n")
  
  #  c(number.it,
  c(sigma.rand.new,betanew.ordinal,betanew.continuous,deltanew,lambdanew,sigma22new,
    cumsum(c(0,deltanew)),1+lambdanew^2*sigma22new ,lambdanew*sigma22new)
  
   list(number.iterations=number.it,
        Sigma.rand.effects=sigma.rand.new,
        beta.ordinal=betanew.ordinal,
        beta.continuous=betanew.continuous,
        regression.coefficients=list(ordinal=betanew.ordinal,continuous=betanew.continuous),
        differences.in.thresholds=deltanew,
        thresholds=cumsum(c(0,deltanew)),
        lambda=lambdanew,
        sigma22=sigma22new,
        cov.errors=lambdanew*sigma22new,
        sigma11=1+lambdanew^2*sigma22new,  
        random.effects=firstmomentb,
        loglikelihood=loglikelihood,
        AIC=AIC,
        BIC=BIC,
        number.iterations=number.it,
        data.ordinal=data.ordinal,
        data.continuous=data.continuous,
        predictors.fixed.ordinal=predictors.fixed.ordinal,
        predictors.fixed.continous=predictors.fixed.continous,
        predictors.random.ordinal=predictors.random.ordinal,
        predictors.random.continuous=predictors.random.continuous,
        exact=exact,
        montecarlo=montecarlo,
        epsilon=epsilon,
        additional=additional
   )
   
  
} #function ecm.ord.plus.cont
