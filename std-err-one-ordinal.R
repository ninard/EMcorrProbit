x=fffnew
epsilon=0.002
exact=T
montecarlo=100

library(MASS)
library(doParallel)
registerDoParallel()

standard.error.bootstrap.one.ordinal=function(data.ordinal,predictors.fixed,
                                              predictors.random, x, exact, montecarlo, 
                                              epsilon, bootstrap.samples) {
  
  #boot=foreach(i=1:bootstrap.samples, .packages=c("MASS","EMcorrProbit")) %dopar% {
  
  beta.estimate=x$regression.coefficients
  delta.estimates=x$differences.in.thresholds
  sigma.rand.estimate=x$Sigma.rand.effects
  
  start.values.beta=beta.estimate
  start.values.delta=delta.estimates
  start.values.sigma.rand=sigma.rand.estimate
  
  l=length(data.ordinal[,1])
  mult.obs=ncol(data.ordinal)
  miss=is.na(data.ordinal)
  
  boot=list(NA)
  
  for(i in 1:bootstrap.samples) {
  
  
    random.int=mvrnorm(n=l, mu=rep(0,ncol(sigma.rand.estimate)), Sigma=sigma.rand.estimate)
    
    pred.fixed=sapply(1:mult.obs,function(obs) as.matrix(predictors.fixed[,,obs])%*%beta.estimate, simplify=T)
    if(ncol(sigma.rand.estimate==1)) pred.rand=sapply(1:mult.obs,function(obs) predictors.random[,,obs]*random.int,simplify=T) else
      pred.rand=sapply(1:mult.obs,function(obs) apply(predictors.random[,,obs]*random.int,1,sum),simplify=T)
    
    y1=pred.fixed+pred.rand+matrix(rnorm(l*mult.obs),ncol=mult.obs)
    
    data.ordinal.new=matrix(cut(y1,c(min(y1)-1,x$thresholds,max(y1)+1), labels=F),ncol=mult.obs)
    data.ordinal.new=ifelse(miss==T, NA, data.ordinal.new)
    
    boot[[i]]=
    ecm.one.ordinal(data.ordinal.new,
                    predictors.fixed,
                    predictors.random,
                    start.values.beta,
                    start.values.delta,
                    start.values.sigma.rand,
                    exact=exact,
                    montecarlo,epsilon, additional=F)
  } #forech
  
  est=t(sapply(1:bootstrap.samples, function(i) c(boot[[i]][[1]], boot[[i]][[2]],boot[[i]][[3]]), simplify=T))
  var(est)
  
} # function standard error

dada=standard.error.bootstrap.one.ordinal(data.ordinal,predictors.fixed,
                                     predictors.random, 
                                     x, exact, montecarlo, epsilon, 
                                     bootstrap.samples=5) 