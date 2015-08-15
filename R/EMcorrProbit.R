
emcorrprobit <- function(model, y, xfixed, xrand, start.values.beta, 
                         start.values.delta=NULL,  start.values.sigma.rand, 
                         exact, montecarlo, epsilon, ...)
{
  obj <- list(model = model, y=y, xfixed=xfixed, xrand=xrand, 
              start.values.beta=start.values.beta, 
              start.values.delta=start.values.delta,  
              start.values.sigma.rand=start.values.sigma.rand, 
              exact=exact, montecarlo=montecarlo, epsilon=epsilon, ...)
  obj$call <-match.call()
  class(obj) <- c(model)
  emcorrprobitFit(obj)
}

emcorrprobitFit <- function(obj, ...) UseMethod("emcorrprobitFit")

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
                         exact,montecarlo,epsilon, additional=TRUE)
  
  
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

ecm.one.ordinal <- function(data.ordinal,predictors.fixed,predictors.random,start.values.beta,start.values.delta,start.values.sigma.rand,
                         exact,montecarlo,epsilon, additional) 
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




