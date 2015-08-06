#' Estimates via EM.
#' 
#' @param xfixed Fixed predictors.
#' @param xrand Random predictors.
#' @param y Ordinal data.
#' @return The estimation of correlated probit model.
#' @examples
#' EMcorrProbit(xfixed, xrand, y, start.values.beta, start.values.delta=NULL,  start.values.sigma.rand, exact, montecarlo=100, epsilon=.001, ...)

emcorrprobit <- function(xfixed, xrand, ...) UseMethod("emcorrprobit")

emcorrprobit.default <- function(y, xfixed, xrand, start.values.beta, 
                                 start.values.delta=NULL,  start.values.sigma.rand, 
                                 exact, montecarlo=100, epsilon=.001, ...)
{
  xfixed <- as.matrix(xfixed)
  xrand <- as.matrix(xrand)
  y <- as.numeric(y)
  
  est <- ecm.one.ordinal(start.values.beta,start.values.delta=NULL,
                         start.values.sigma.rand,
                         exact,montecarlo=100,y,xfixed,xrandom,epsilon=.001)
  
  
  est$call <-match.call()
  
  class(est) <- "emcorrprobit"
  est
}

print.emcorrprobit <- function(x, ...)
{
  print(x)
}
summary.emcorrprobit <- function(x, ...)
{
  print(x)
}

print.summary.emcorrprobit <- function(x, ...)
{
  print(x)
}

formula.emcorrprobit <- function(formula, data=list(), ...)
{
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  
  est <- emcorrprobit.default(x, y, ...)
  est$call <-match.call
  est$formula <-formula
  est
}

ecm.one.ordinal=function(start.values.beta,start.values.delta=NULL,start.values.sigma.rand,
                         exact,montecarlo=100,data.ordinal,predictors.fixed,predictors.random,epsilon=.001) {
  
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
  
  
  #######################
  ## check
  ############################
  if(class(data.ordinal)!="matrix") print("Warning message: data.ordinal must be a matrix")
  if(class(start.values.sigma.rand)!="matrix") print("Warning message: start.values.sigma.rand must be a matrix")
  if(any(sort(unique(c(data.ordinal)))!=1:max(data.ordinal,na.rm=T))) print("Warning message: missing levels in the response variable")
  if(length(start.values.delta)!=max(data.ordinal,na.rm=T)-2) print("Warning message: incorrect dimension of start.values.delta")
  #if(exact==F & is.na(montecarlo)) {print("Warning message: montecarlo parameter isn't defined, taken by default 100")
  #montecarlo=100
  #}
  
  
  ##############################################################
  
  betanew=start.values.beta
  deltanew=start.values.delta
  sigma.rand.new=start.values.sigma.rand
  
  num.categories=max(data.ordinal,na.rm=T)
  if(num.categories==2) deltanew=NULL
  
  
  new.est=c(start.values.sigma.rand,start.values.beta,start.values.delta)
  old.est=new.est+1.5*epsilon
  
  number.it=0
  
  # definition of lower and upper boundary of the truncated normal distr. 
  # (the new variable given observed data)
  trunc.lower=ifelse(data.ordinal==1,-Inf,0)
  trunc.upper=ifelse(data.ordinal==1,0,ifelse(data.ordinal<num.categories,1,Inf))
  
  #number of observations in each category
  n=table(data.ordinal)
  individuals=length(data.ordinal[,1])
  mult.obs=dim(predictors.fixed)[3]
  dimension.sigma=ncol(sigma.rand.new)
  
  z=sapply(1:individuals, function(i) matrix(t(predictors.random[i,,]), ncol=dimension.sigma,byrow=F), 
           simplify="array")
  ## mult.obs x dimension.sigma x individuals
  
  ### definition na k, knot i kk
  k=sapply(1:individuals, function(i) which(is.na(data.ordinal[i,])), simplify=F)
  ### list: missings for each individual
  knot=sapply(1:individuals, function(i) which(!is.na(data.ordinal[i,])),simplify=F)
  ### list: observed for each individual
  kk=sapply(1:individuals, function(i) length(k[[i]]), simplify=T)
  ### individuals: length of missings for ecah individual
  
  ######################################################################
  ############             while loop           ##########################
  ##########################################################################
  while(any(abs(new.est-old.est)>epsilon)) {
    
    number.it=number.it+1
    
    old.est=new.est
    
    #######################################################################
    delta.exp=c(1,deltanew,1)
    delta=matrix(delta.exp[data.ordinal],ncol=mult.obs)
    alpha.exp=c(0,0,cumsum(deltanew))
    alpha=matrix(alpha.exp[data.ordinal],ncol=mult.obs)
    
    ###expectation of ynew
    if(length(start.values.beta)>1) pred=sapply(1:mult.obs, function(i) predictors.fixed[,,i]%*%betanew, simplify=TRUE)  else pred=betanew
    ## makes pred a matrix: observations x mult.obs
    mean.ynew=(pred-alpha)/delta
    
    variance.ynew=sapply(1:individuals, function(i)
      (z[,,i]%*%sigma.rand.new%*%t(z[,,i])+diag(mult.obs))/(delta[i,]%*%t(delta[i,])),
      simplify="array")
    ###
    
    #######################################################
    ############### if there is a problem for some individuals - second version!!! independent generation
    ##################################################################################
    if(exact==F) {
      simulations=sapply(1:individuals, function(i) { 
        if(kk[i]>0)  rtmvnorm(montecarlo, mean=mean.ynew[i,-k[[i]]], 
                              sigma=matrix(as.vector(variance.ynew[-k[[i]],-k[[i]],i]),ncol=mult.obs-kk[i]),
                              lower=trunc.lower[i,-k[[i]]],upper=trunc.upper[i,-k[[i]]], algorithm="gibbs") else  rtmvnorm(
                                montecarlo, mean=mean.ynew[i,], 
                                sigma=matrix(as.vector(variance.ynew[,,i]),ncol=mult.obs),
                                lower=trunc.lower[i,],upper=trunc.upper[i,], algorithm="gibbs")
      }, # function
      simplify=F)  
      ## simulations is a list
      moments1=sapply(1:individuals, function(i) apply(as.matrix(simulations[[i]]),2,mean),simplify=F)
      ## list
      moments2=sapply(1:individuals, function(i) var(simulations[[i]]),simplify="array")
      ## list
      ## montecarlo x mult.obs x observations
    } else {
      simulations.exact=sapply(1:individuals, function(i) {
        if(kk[i]>0) mtmvnorm(mean=mean.ynew[i,-k[[i]]], 
                             sigma=matrix(as.vector(variance.ynew[-k[[i]],-k[[i]],i]),ncol=mult.obs-kk[i]),
                             lower=trunc.lower[i,-k[[i]]],upper=trunc.upper[i,-k[[i]]]) else 	mtmvnorm(mean=mean.ynew[i,], 
                                                                                                       sigma=matrix(as.vector(variance.ynew[,,i]),ncol=mult.obs),
                                                                                                       lower=trunc.lower[i,],upper=trunc.upper[i,])
      } #function
      ,simplify=F)
      ##  mult.obs x observations
      moments1=sapply(1:individuals, function(i) simulations.exact[[i]]$tmean,simplify=F)
      ## list
      moments2=sapply(1:individuals, function(i) simulations.exact[[i]]$tvar,simplify=F)
      ## the variance of y latent
    }
    
    moments1A=array(NA,dim=c(individuals,mult.obs))
    for(i in 1:individuals) moments1A[i,knot[[i]]]=moments1[[i]]
    
    moments2A=list(0)
    for(i in 1:individuals) {moments2A[[i]]=array(NA, dim=c(mult.obs, mult.obs))
                             moments2A[[i]][knot[[i]],knot[[i]]]=moments2[[i]]}
    
    
    ######### to calculate the expectation of the random effects
    
    sigmab=sapply(1:individuals, function(i) {
      if(kk[i]>0 & kk[i]<(mult.obs-1)) t(t(sigma.rand.new%*%t(z[-k[[i]],,i]))/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]) else
        if(kk[i]==(mult.obs-1)) t(t(sigma.rand.new%*%z[-k[[i]],,i])/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]) else 
          t(t(sigma.rand.new%*%t(z[,,i]))/delta[i,])%*%solve(variance.ynew[,,i])},simplify=F)
    ### list: dimension.sigma x mult.obs x individuals
    
    firstmomentb=sapply(1:individuals, function(i) {
      if(kk[i]>0 & kk[i]<(mult.obs-1)) sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,-k[[i]]]) else 
        if(kk[i]==(mult.obs-1))    sigmab[[i]]%*%as.matrix(moments1[[i]]-mean.ynew[i,-k[[i]]])  else 
          sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,])}, simplify="array")
    ### array: dimension.sigma x 1 x individuals
    ### when dimension.sigma=1: vector : individuals 
    
    
    if (dimension.sigma>1)  pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
      pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
    ### individuals x mult.obs
    
    
    
    ##################################################
    #########################  beta estimate
    ##################################################
    
    
    ywave=delta*moments1A-pred.random+alpha
    if(length(start.values.beta)==1) betanew=lm(c(ywave)~1)$coef else {
      pred.beta=predictors.fixed[,,1]
      for(i in 2:mult.obs) pred.beta=rbind(pred.beta,predictors.fixed[,,i])
      betanew=lm(c(ywave)~pred.beta-1)$coef 
      rm(pred.beta)}
    #print(betanew)
    
    ############ update
    
    ###expectation of ynew
    if(length(start.values.beta)>1) pred=sapply(1:mult.obs, function(i) predictors.fixed[,,i]%*%betanew, simplify=TRUE)  else pred=betanew
    ## makes pred a matrix: observations x mult.obs
    mean.ynew=(pred-alpha)/delta
    
    
    
    firstmomentb=sapply(1:individuals, function(i) {
      if(kk[i]>0 & kk[i]<(mult.obs-1)) sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,-k[[i]]]) else 
        if(kk[i]==(mult.obs-1))    sigmab[[i]]%*%as.matrix(moments1[[i]]-mean.ynew[i,-k[[i]]])  else 
          sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,])}, simplify="array")
    ### array: dimension.sigma x 1 x individuals
    ### when dimension.sigma=1: vector : individuals 
    
    
    if (dimension.sigma>1)  pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
      pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
    ### individuals x mult.obs
    
    
    
    #######################################
    ################### delta estimates
    ##########################################
    if(num.categories>2) {
      
      ### where NA - with 0 - for first and second1
      first=t(sapply(1:individuals, function(i) diag(moments2A[[i]])+(moments1A[i,])^2, simplify="array"))
      ### individuals x mult.obs
      
      for (delta.index in 1:(num.categories-2)) {
        
        second1=moments1A*(pred-alpha)
        ### individuals x mult.lbs
        
        alpha.exp.term=c(0,0)
        for(i in 1:length(deltanew)) alpha.exp.term[i+2]=ifelse(delta.index<=i, alpha.exp[i+2]-deltanew[delta.index], 0)
        alpha.term=matrix(alpha.exp.term[data.ordinal],ncol=mult.obs)
        
        
        #### to check it carefully - if I have only one observation for some individual
        second2=sapply(1:individuals, function(i)
          if(kk[i]>0) sigmab[[i]]%*%(moments2[[i]]+moments1[[i]]%*%t(moments1[[i]])-mean.ynew[i,-k[[i]]]%*%t(moments1[[i]])
          ) else sigmab[[i]]%*%(moments2[[i]]+moments1[[i]]%*%t(moments1[[i]])-mean.ynew[i,]%*%t(moments1[[i]])),
          simplify=F)
        ## for each i: q x n_i matrix 
        
        if(dimension.sigma==1)  second2=sapply(1:individuals, function(i) {if(kk[i]>0) 
          z[-k[[i]],,i]*as.vector(second2[[i]]) else z[,,i]*as.vector(second2[[i]]) },simplify=F)  else
            second2=sapply(1:individuals, function(i) {if(kk[i]>0) diag(z[-k[[i]],,i]%*%second2[[i]]) else diag(z[,,i]%*%second2[[i]])},
                           simplify=F)
        
        second2new=array(NA, dim=c(l,mult.obs))
        for(i in 1:individuals) {
          if(kk[i]>0) second2new[i,knot[[i]]]=second2[[i]] else second2new[i,]=second2[[i]]
        }
        
        third=moments1A*delta-pred-pred.random+alpha.term
        ## observations x mult.obs
        
        a=sum(first[data.ordinal==(delta.index+1)],na.rm=T)+sum(n[(delta.index+2):num.categories])
        b=sum((second1+second2new)[data.ordinal==delta.index+1], na.rm=T)-sum(third[data.ordinal>delta.index+1], na.rm=T)
        c=-n[(delta.index+1)]
        
        deltanew[delta.index]=(b+sqrt(b^2-4*a*c))/(2*a) 
        #print(deltanew[delta.index])
        
        
        #################### update
        
        delta.exp=c(1,deltanew,1)
        delta=matrix(delta.exp[data.ordinal],ncol=mult.obs)
        alpha.exp=c(0,0,cumsum(deltanew))
        alpha=matrix(alpha.exp[data.ordinal],ncol=mult.obs)
        
        mean.ynew=(pred-alpha)/delta
        
        variance.ynew=sapply(1:individuals, function(i)
          (z[,,i]%*%sigma.rand.new%*%t(z[,,i])+diag(mult.obs))/(delta[i,]%*%t(delta[i,])),
          simplify="array")
        
        sigmab=sapply(1:individuals, function(i) {
          if(kk[i]>0 & kk[i]<(mult.obs-1)) t(t(sigma.rand.new%*%t(z[-k[[i]],,i]))/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]) else
            if(kk[i]==(mult.obs-1)) t(t(sigma.rand.new%*%z[-k[[i]],,i])/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]) else 
              t(t(sigma.rand.new%*%t(z[,,i]))/delta[i,])%*%solve(variance.ynew[,,i])},simplify=F)
        ### list: dimension.sigma x mult.obs x individuals
        
        firstmomentb=sapply(1:individuals, function(i) {
          if(kk[i]>0 & kk[i]<(mult.obs-1)) sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,-k[[i]]]) else 
            if(kk[i]==(mult.obs-1))    sigmab[[i]]%*%as.matrix(moments1[[i]]-mean.ynew[i,-k[[i]]])  else 
              sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,])}, simplify="array")
        ### array: dimension.sigma x 1 x individuals
        ### when dimension.sigma=1: vector : individuals 
        
        if (dimension.sigma>1)  pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
          pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
        ### individuals x mult.obs
      } 
      
    } #if
    #############################################################################################
    ############################### sigma.rand estimate
    #############################################################################################
    
    sigma.rand=sapply(1:individuals, function(i) {
      if(kk[i]>0 & kk[i]<(mult.obs-1)) sigma.rand.new-sigmab[[i]]%*%((z[-k[[i]],,i]%*%sigma.rand.new)/delta[i,-k[[i]]])+
        sigmab[[i]]%*%(moments2[[i]]+moments1[[i]]%*%t(moments1[[i]])
                       -(moments1[[i]])%*%t(mean.ynew[i,-k[[i]]])
                       -(mean.ynew[i,-k[[i]]])%*%t(moments1[[i]])
                       +(mean.ynew[i,-k[[i]]])%*%t(mean.ynew[i,-k[[i]]])
        )%*%t(sigmab[[i]]) else   if(kk[i]==(mult.obs-1)) sigma.rand.new-sigmab[[i]]%*%((z[-k[[i]],,i]%*%sigma.rand.new)/delta[i,-k[[i]]])+
        sigmab[[i]]%*%(moments2[[i]]+moments1[[i]]%*%t(moments1[[i]])-moments1[[i]]%*%
                         t(mean.ynew[i,-k[[i]]])-(mean.ynew[i,-k[[i]]])%*%t(moments1[[i]])+(mean.ynew[i,-k[[i]]])%*%t(mean.ynew[i,-k[[i]]])
        )%*%t(sigmab[[i]])    else sigma.rand.new-sigmab[[i]]%*%((z[,,i]%*%sigma.rand.new)/delta[i,])+
        sigmab[[i]]%*%(moments2[[i]]+(moments1[[i]])%*%t(moments1[[i]])-(moments1[[i]])%*%t(mean.ynew[i,])
                       -(mean.ynew[i,])%*%t(moments1[[i]])+(mean.ynew[i,])%*%t(mean.ynew[i,])
        )%*%t(sigmab[[i]])
    }, simplify=T)
    
    
    if(dimension.sigma==1) sigma.rand.new=matrix(mean(sigma.rand)) else {sigma.rand.new=matrix(apply(sigma.rand,1,mean), ncol=dimension.sigma)
                                                                         sigma.rand.new[lower.tri(sigma.rand.new)] = t(sigma.rand.new)[lower.tri(sigma.rand.new)] }
    #print(sigma.rand.new)
    
    ########################
    new.est=c(sigma.rand.new,betanew,deltanew)
    ########################
    cat(c(number.it,sigma.rand.new,betanew,deltanew),"\n")
  } # while
  
  
  cat("Final estimates with accuracy=",epsilon," :","\n",
      "number of iterations: ",number.it,"\n",
      "sigma estimate: ",sigma.rand.new,"\n",
      "betanew.estimate: ",betanew,"\n")
  if(num.categories>2) cat("delta estimate: ",deltanew,"\n")
  
  
  ###########################################################################
  ################### random effects
  ###########################################################################
  delta.exp=c(1,deltanew,1)
  delta=matrix(delta.exp[data.ordinal],ncol=mult.obs)
  alpha.exp=c(0,0,cumsum(deltanew))
  alpha=matrix(alpha.exp[data.ordinal],ncol=mult.obs)
  
  ###expectation of ynew
  if(length(start.values.beta)>1) pred=sapply(1:mult.obs, function(i) predictors.fixed[,,i]%*%betanew, simplify=TRUE)  else pred=betanew
  ## makes pred a matrix: observations x mult.obs
  mean.ynew=(pred-alpha)/delta
  
  variance.ynew=sapply(1:individuals, function(i)
    (z[,,i]%*%sigma.rand.new%*%t(z[,,i])+diag(mult.obs))/(delta[i,]%*%t(delta[i,])),
    simplify="array")
  ###
  
  #######################################################
  ############### if there is a problem for some individuals - second version!!! independent generation
  ##################################################################################
  if(exact==F) {
    simulations=sapply(1:individuals, function(i) { 
      if(kk[i]>0)  rtmvnorm(montecarlo, mean=mean.ynew[i,-k[[i]]], 
                            sigma=matrix(as.vector(variance.ynew[-k[[i]],-k[[i]],i]),ncol=mult.obs-kk[i]),
                            lower=trunc.lower[i,-k[[i]]],upper=trunc.upper[i,-k[[i]]], algorithm="gibbs") else  rtmvnorm(
                              montecarlo, mean=mean.ynew[i,], 
                              sigma=matrix(as.vector(variance.ynew[,,i]),ncol=mult.obs),
                              lower=trunc.lower[i,],upper=trunc.upper[i,], algorithm="gibbs")
    }, # function
    simplify=F)  
    ## simulations is a list
    moments1=sapply(1:individuals, function(i) apply(as.matrix(simulations[[i]]),2,mean),simplify=F)
    ## list
    moments2=sapply(1:individuals, function(i) var(simulations[[i]]),simplify="array")
    ## list
    ## montecarlo x mult.obs x observations
  } else {
    simulations.exact=sapply(1:individuals, function(i) {
      if(kk[i]>0) mtmvnorm(mean=mean.ynew[i,-k[[i]]], 
                           sigma=matrix(as.vector(variance.ynew[-k[[i]],-k[[i]],i]),ncol=mult.obs-kk[i]),
                           lower=trunc.lower[i,-k[[i]]],upper=trunc.upper[i,-k[[i]]]) else 	mtmvnorm(mean=mean.ynew[i,], 
                                                                                                     sigma=matrix(as.vector(variance.ynew[,,i]),ncol=mult.obs),
                                                                                                     lower=trunc.lower[i,],upper=trunc.upper[i,])
    } #function
    ,simplify=F)
    ##  mult.obs x observations
    moments1=sapply(1:individuals, function(i) simulations.exact[[i]]$tmean,simplify=F)
    ## list
    moments2=sapply(1:individuals, function(i) simulations.exact[[i]]$tvar,simplify=F)
    ## the variance of y latent
  }
  
  moments1A=array(NA,dim=c(individuals,mult.obs))
  for(i in 1:individuals) moments1A[i,knot[[i]]]=moments1[[i]]
  
  moments2A=list(0)
  for(i in 1:individuals) {moments2A[[i]]=array(NA, dim=c(mult.obs, mult.obs))
                           moments2A[[i]][knot[[i]],knot[[i]]]=moments2[[i]]}
  
  
  ######### to calculate the expectation of the random effects
  
  sigmab=sapply(1:individuals, function(i) {
    if(kk[i]>0 & kk[i]<(mult.obs-1)) t(t(sigma.rand.new%*%t(z[-k[[i]],,i]))/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]) else
      if(kk[i]==(mult.obs-1)) t(t(sigma.rand.new%*%z[-k[[i]],,i])/delta[i,-k[[i]]])%*%solve(variance.ynew[-k[[i]],-k[[i]],i]) else 
        t(t(sigma.rand.new%*%t(z[,,i]))/delta[i,])%*%solve(variance.ynew[,,i])},simplify=F)
  ### list: dimension.sigma x mult.obs x individuals
  
  firstmomentb=sapply(1:individuals, function(i) {
    if(kk[i]>0 & kk[i]<(mult.obs-1)) sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,-k[[i]]]) else 
      if(kk[i]==(mult.obs-1))    sigmab[[i]]%*%as.matrix(moments1[[i]]-mean.ynew[i,-k[[i]]])  else 
        sigmab[[i]]%*%(moments1[[i]]-mean.ynew[i,])}, simplify="array")
  ### array: dimension.sigma x 1 x individuals
  ### when dimension.sigma=1: vector : individuals 
  
  
  #########################################################################################
  #########################################################################################
  #########################################################################################
  ################ log-likelihood #########################################################
  #########################################################################################
  #########################################################################################
  
  ordinal.part=sapply(1:individuals, function(i) log(pmvnorm(mean=mean.ynew[i,knot[[i]]],
                                                             sigma=matrix(variance.ynew[knot[[i]],knot[[i]],i], ncol=mult.obs-kk[i]),
                                                             lower=trunc.lower[i,knot[[i]]], upper=trunc.upper[i,knot[[i]]])),
                      simplify="array")
  
  loglikelihood=sum(ordinal.part)
  
  cat("Log-likelihood: ",loglikelihood,"\n")
  
  ###################################################################################
  ############## AIC and BIC ########################################################
  ###################################################################################
  AIC=-2*loglikelihood+2*(length(betanew)+length(deltanew)+choose(dimension.sigma+1,2))
  
  BIC=-2*loglikelihood+log(length(which(!is.na(data.ordinal))))*(length(betanew)+length(deltanew)+choose(dimension.sigma+1,2))
  
  
  cat("AIC: ",AIC,"\n")
  cat("BIC: ",BIC,"\n")
  
  #########################################################################################
  #########################################################################################
  list(Sigma.rand.effects=sigma.rand.new,
       regression.coefficients=betanew,
       differences.in.thresholds=deltanew,
       thresholds=c(0,cumsum(deltanew)),
       random.effects=firstmomentb,
       loglikelihood=loglikelihood,
       AIC=AIC,
       BIC=BIC)
  
} #function

#################################################
######################   ecm algorithm
########################################################
ecm.one.ordinal.complete.cases=function(start.values.beta,start.values.delta=NULL,start.values.sigma.rand,
                                        exact,montecarlo=100,data.ordinal,predictors.fixed,predictors.random,epsilon=.001) {
  
  ########################################
  ### wide format of longitudinal data: first level denoted by 1, second - 2 and so on
  ### predictors.fixed is a 3-dimentional array: individuals x dimentions predictros fixed x time points 
  ### predictors.random is a 3-dimentional array: individuals x dimentions predictros fixed x time points  
  ### sigma.rand is the covaraince matrix of the random effects
  ##################################
  #library(MASS)
  #library(tmvtnorm)
  
  #######################
  ## check
  ############################
  if(class(data.ordinal)!="matrix") print("Warning message: data.ordinal must be a matrix")
  if(class(start.values.sigma.rand)!="matrix") print("Warning message: start.values.sigma.rand must be a matrix")
  if(any(sort(unique(c(data.ordinal)))!=1:max(data.ordinal,na.rm=T))) print("Warning message: missing levels in the response variable")
  if(length(start.values.delta)!=max(data.ordinal,na.rm=T)-2) print("Warning message: incorrect dimension of start.values.delta")
  #if(exact==F & !exists("montecarlo")) {print("Warning message: montecarlo parameter isn't defined, taken by default 100")
  #montecarlo=100
  #}
  
  ################################################################3
  
  
  
  betanew=start.values.beta
  deltanew=start.values.delta
  sigma.rand.new=start.values.sigma.rand
  
  num.categories=max(data.ordinal)
  if(num.categories==2) deltanew=NULL
  
  new.est=c(start.values.sigma.rand,start.values.beta,start.values.delta)
  old.est=new.est+1.5*epsilon
  
  number.it=0
  
  # definition of lower and upper boundary of the truncated normal distr. 
  # (the new variable given observed data)
  trunc.lower=ifelse(data.ordinal==1,-Inf,0)
  trunc.upper=ifelse(data.ordinal==1,0,ifelse(data.ordinal<num.categories,1,Inf))
  #number of observations in each category
  n=0
  for (i in 1:num.categories) n[i]=sum(data.ordinal==i)
  individuals=length(data.ordinal[,1])
  mult.obs=dim(predictors.fixed)[3]
  dimension.sigma=ncol(sigma.rand.new)
  
  z=sapply(1:individuals, function(i) matrix(t(predictors.random[i,,]), ncol=dimension.sigma,byrow=F), 
           simplify="array")
  ## mult.obs x dimension.sigma x individuals
  
  ######################################################################
  ############             while loop           ##########################
  ##########################################################################
  
  #while(any(abs(c(betanew-betaold,deltanew-deltaold,sigma.rand.new-sigma.rand.old))>epsilon)) {
  while(any(abs(new.est-old.est)>epsilon)) {
    
    number.it=number.it+1
    
    old.est=new.est
    
    #######################################################################
    delta.exp=c(1,deltanew,1)
    delta=matrix(delta.exp[data.ordinal],ncol=mult.obs)
    alpha.exp=c(0,0,cumsum(deltanew))
    alpha=matrix(alpha.exp[data.ordinal],ncol=mult.obs)
    
    
    ###expectation of ynew
    if(length(start.values.beta)>1) pred=sapply(1:mult.obs, function(i) predictors.fixed[,,i]%*%betanew, simplify=TRUE)  else pred=betanew
    ## makes pred a mtarix: individuals x mult.obs
    mean.ynew=(pred-alpha)/delta
    
    variance.ynew=sapply(1:individuals, function(i)
      (z[,,i]%*%sigma.rand.new%*%t(z[,,i])+diag(mult.obs))/(delta[i,]%*%t(delta[i,])),
      simplify="array")
    
    ### simulations
    if(exact==F) {simulations=sapply(1:individuals, function(i) rtmvnorm(montecarlo, mean=mean.ynew[i,], 
                                                                         sigma=matrix(as.vector(variance.ynew[,,i]),ncol=mult.obs),
                                                                         lower=trunc.lower[i,],upper=trunc.upper[i,], algorithm="gibbs"),
                                     simplify="array")
                  ## montecarlo x mult.obs x individuals
                  moments1=t(sapply(1:individuals, function(i) apply(simulations[,,i],2,mean),simplify="array"))
                  moments2=sapply(1:individuals, function(i) var(simulations[,,i]),simplify="array")} else
                  {simulations.exact=sapply(1:individuals, function(i) mtmvnorm(mean=mean.ynew[i,], 
                                                                                sigma=matrix(as.vector(variance.ynew[,,i]),ncol=mult.obs),
                                                                                lower=trunc.lower[i,],upper=trunc.upper[i,]),simplify="array")
                   ##  mult.obs x individuals
                   moments1=t(sapply(1:individuals, function(i) simulations.exact[,i]$tmean,simplify="array"))
                   moments2=sapply(1:individuals, function(i) simulations.exact[,i]$tvar,simplify="array")
                   ## the variance of y latent
                  }
    
    
    ######### to calculate the expectation of the random effects
    sigmab=sapply(1:individuals, function(i) t(t(sigma.rand.new%*%t(z[,,i]))/delta[i,])%*%
                    solve(variance.ynew[,,i]),simplify="array")
    
    firstmomentb=sapply(1:individuals, function(i) sigmab[,,i]%*%((moments1[i,])-
                                                                    mean.ynew[i,]), simplify="array")
    
    if (dimension.sigma>1)  
      pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
        pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
    ### individuals x mult.obs
    
    
    
    ###################################################################
    ########################  beta estimate
    ###################################################################
    
    ywave=delta*moments1-pred.random+alpha
    if(length(start.values.beta)==1) betanew=lm(c(ywave)~1)$coef else {
      pred.beta=predictors.fixed[,,1]
      for(i in 2:mult.obs) pred.beta=rbind(pred.beta,predictors.fixed[,,i])
      betanew=lm(c(ywave)~pred.beta-1)$coef 
      rm(pred.beta)}
    #print(betanew)
    
    
    #### update  
    
    ###expectation of ynew
    if(length(start.values.beta)>1) pred=sapply(1:mult.obs, function(i) predictors.fixed[,,i]%*%betanew, simplify=TRUE)  else pred=betanew
    ## makes pred a mtarix: individuals x mult.obs
    mean.ynew=(pred-alpha)/delta
    
    firstmomentb=sapply(1:individuals, function(i) sigmab[,,i]%*%((moments1[i,])-
                                                                    mean.ynew[i,]), simplify="array")
    
    if (dimension.sigma>1)  
      pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
        pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
    ### individuals x mult.obs
    
    
    
    #######################################
    ################### delta estimates
    ##########################################
    if(num.categories>2) {
      
      first=t(sapply(1:individuals, function(i) diag(moments2[,,i])+(moments1[i,])^2, simplify="array"))
      
      for (k in 1:(num.categories-2)) {
        
        second1=moments1*(pred-alpha)
        
        alpha.exp.term=c(0,0)
        for(i in 1:length(deltanew)) alpha.exp.term[i+2]=ifelse(k<=i, alpha.exp[i+2]-deltanew[k], 0)
        alpha.term=matrix(alpha.exp.term[data.ordinal],ncol=mult.obs)
        ## checked
        
        second2=sapply(1:individuals, function(i)
          sigmab[,,i]%*%(moments2[,,i]+moments1[i,]%*%t(moments1[i,])-mean.ynew[i,]%*%t(moments1[i,])),
          simplify="array")
        ## for each individual:  q x n_i matrix 
        
        if (dimension.sigma==1) second2=t(sapply(1:individuals, function(i) diag(z[,,i]%*%t(second2[,,i])),
                                                 simplify="array")) else second2=t(sapply(1:individuals, function(i) diag(z[,,i]%*%second2[,,i]),simplify="array"))
        
        third=moments1*delta-pred-pred.random+alpha.term
        
        a=sum(first[data.ordinal==(k+1)])+sum(n[(k+2):num.categories])
        b=sum((second1+second2)[data.ordinal==k+1])-sum(third[data.ordinal>k+1])
        c=-n[(k+1)]
        
        deltanew[k]=(b+sqrt(b^2-4*a*c))/(2*a) 
        #print(deltanew[k])
        
        ######### update
        delta.exp=c(1,deltanew,1)
        delta=matrix(delta.exp[data.ordinal],ncol=mult.obs)
        alpha.exp=c(0,0,cumsum(deltanew))
        alpha=matrix(alpha.exp[data.ordinal],ncol=mult.obs)
        
        mean.ynew=(pred-alpha)/delta
        
        variance.ynew=sapply(1:individuals, function(i)
          (z[,,i]%*%sigma.rand.new%*%t(z[,,i])+diag(mult.obs))/(delta[i,]%*%t(delta[i,])),
          simplify="array")
        
        sigmab=sapply(1:individuals, function(i) t(t(sigma.rand.new%*%t(z[,,i]))/delta[i,])%*%
                        solve(variance.ynew[,,i]),simplify="array")
        
        firstmomentb=sapply(1:individuals, function(i) sigmab[,,i]%*%((moments1[i,])-
                                                                        mean.ynew[i,]), simplify="array")
        
        if (dimension.sigma>1)  
          pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%firstmomentb[,,i],simplify=T))       else 
            pred.random=t(sapply(1:individuals, function(i) z[,,i]%*%t(firstmomentb[i]),simplify=T))  
        ### individuals x mult.obs
        
      }
      
    } #if
    ############################### sigma.rand estimate
    if(dimension.sigma==1)
      sigma.rand=sapply(1:individuals, function(i)
        sigma.rand.new-sigmab[,,i]%*%((z[,,i]%*%sigma.rand.new)/delta[i,])+
          sigmab[,,i]%*%(moments2[,,i]+(moments1[i,])%*%t(moments1[i,])-(moments1[i,])%*%t(mean.ynew[i,]) 
                         -(mean.ynew[i,])%*%t(moments1[i,])+(mean.ynew[i,])%*%t(mean.ynew[i,])
          )%*%sigmab[,,i], simplify=T) else
            sigma.rand=sapply(1:individuals, function(i)
              sigma.rand.new-sigmab[,,i]%*%((z[,,i]%*%sigma.rand.new)/delta[i,])+
                sigmab[,,i]%*%(moments2[,,i]+(moments1[i,])%*%t(moments1[i,])-(moments1[i,])%*%t(mean.ynew[i,])
                               -(mean.ynew[i,])%*%t(moments1[i,])+(mean.ynew[i,])%*%t(mean.ynew[i,])
                )%*%t(sigmab[,,i]), simplify=T)
    
    
    if(dimension.sigma==1) sigma.rand.new=matrix(mean(sigma.rand)) else {sigma.rand.new=matrix(apply(sigma.rand,1,mean), ncol=dimension.sigma)
                                                                         sigma.rand.new[lower.tri(sigma.rand.new)] = t(sigma.rand.new)[lower.tri(sigma.rand.new)] }
    #print(sigma.rand.new)
    
    ########################
    new.est=c(sigma.rand.new,betanew,deltanew)
    ########################
    cat(c(number.it,sigma.rand.new,betanew,deltanew,"\n"))
  } # while
  
  
  cat("Final estimates with accuracy=",epsilon," :","\n",
      "number of iterations: ",number.it,"\n",
      "sigma estimate: ",sigma.rand.new,"\n",
      "betanew.estimate: ",betanew,"\n")
  if(num.categories>2) cat("delta estimate: ",deltanew,"\n")
  
  #############################################################################################
  ############################################################################################
  
  
  
  
  
  
  
  
  #############################################################################################
  ############################################################################################
  #########################################################################################
  #########################################################################################
  #########################################################################################
  ################ log-likelihood #########################################################
  #########################################################################################
  #########################################################################################
  
  ordinal.part=sapply(1:individuals, function(i) log(pmvnorm(mean=mean.ynew[i,],
                                                             sigma=variance.ynew[,,i], lower=trunc.lower[i,], upper=trunc.upper[i,])),
                      simplify="array")
  
  loglikelihood=sum(ordinal.part)
  
  cat("Log-likelihood: ",loglikelihood,"\n")
  
  ###################################################################################
  ############## AIC and BIC ########################################################
  ###################################################################################
  AIC=-2*loglikelihood+2*(length(betanew)+length(deltanew)+choose(dimension.sigma+1,2))
  
  BIC=-2*loglikelihood+log(mult.obs*individuals)*(length(betanew)+length(deltanew)+choose(dimension.sigma+1,2))
  
  cat("AIC: ",AIC,"\n")
  cat("BIC: ",BIC,"\n")
  
  
  #########################################################################################
  #########################################################################################
  list(Sigma.rand.effects=sigma.rand.new,
       regression.coefficients=betanew,
       differences.in.thresholds=deltanew,
       thresholds=c(0,cumsum(deltanew)),
       random.effects=firstmomentb,
       loglikelihood=loglikelihood,
       AIC=AIC,
       BIC=BIC)
  
} #function
