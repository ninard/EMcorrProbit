# 578 - regression.coefficients=list(ordinal=betanew) ??


standard.error.bootstrap.one.ordinal.one.continuous=function(x, bootstrap.samples = 50, epsilon=NULL,
                                              doParallel = FALSE, cores=NULL) {
    
  beta.ordinal.estimate=(x$regression.coefficients)$ordinal
  beta.continuous.estimate=(x$regression.coefficients)$continuous
  delta.estimates=x$differences.in.thresholds
  sigma.rand.estimate=x$Sigma.rand.effects
  sigma22.estimate <- x$sigma22
  lambda.estimate <- x$lambda
     
  l=nrow(x$data.ordinal)
  mult.obs=ncol(x$data.ordinal)
  miss.ordinal=is.na(x$data.ordinal)
  miss.continuous=is.na(x$data.continuous)
  ### to check
  re.ord=dim(x$predictors.random.ordinal)[2]
  re.cont=dim(x$predictors.random.continuous)[2]
  
  if (length(epsilon)==0) epsilon=x$epsilon*10
  
  if(doParallel) {if (length(cores)>0) registerDoParallel(cores=cores) else registerDoParallel()
                  boot=foreach(i=1:bootstrap.samples, .packages=c("MASS","EMcorrProbit","tmvtnorm")) %dopar% {
                    
                    random.int=mvrnorm(n=l, mu=rep(0,ncol(sigma.rand.estimate)), Sigma=sigma.rand.estimate)
                    
                    pred.fixed.ordinal=sapply(1:mult.obs,function(obs) 
                      as.matrix(x$predictors.fixed.ordinal[,,obs])%*%beta.ordinal.estimate, simplify=TRUE)
                    pred.fixed.continuous=sapply(1:mult.obs,function(obs) 
                      as.matrix(x$predictors.fixed.continuous[,,obs])%*%beta.continuous.estimate, simplify=TRUE)
       ### to check this !!!!!!
       
                    if(re.ord==1) 
                      pred.rand.ord=sapply(1:mult.obs,function(obs) x$predictors.random.ordinal[,,obs]*random.int[,1:re.ord],simplify=TRUE) else
                      pred.rand.ord=sapply(1:mult.obs,function(obs) apply(x$predictors.random.ordinal[,,obs]*random.int,1,sum),simplify=TRUE)
                   if(re.cont==1) 
                      pred.rand.cont=sapply(1:mult.obs,function(obs) 
                        x$predictors.random.continuous[,,obs]*random.int[,(re.cont+1):ncol(sigma.rand.estimate)],simplify=TRUE) else
                      pred.rand.cont=sapply(1:mult.obs,function(obs) 
                        apply(x$predictors.random.continuous[,,obs]*random.int[,(re.cont+1):ncol(sigma.rand.estimate)],1,sum),simplify=TRUE)
                    
                    y1=pred.fixed.ordinal+pred.rand.ord+matrix(rnorm(l*mult.obs),ncol=mult.obs)
                    y2=pred.fixed.continuous+pred.rand.cont+matrix(rnorm(l*mult.obs),ncol=mult.obs)
                    
                    data.ordinal.new=matrix(cut(y1,c(min(y1)-1,x$thresholds,max(y1)+1), labels=FALSE),ncol=mult.obs)
                    data.ordinal.new=ifelse(miss.ordinal, NA, data.ordinal.new)
                    data.continuous.new=ifelse(miss.continuous, NA, y2)
                    
       
       ecm.ord.plus.cont(data.ordinal.new,data.continuous.new,x$predictors.fixed.ordinal,x$predictors.fixed.continous,
                         x$predictors.random.ordinal,x$predictors.random.continuous,
                         beta.ordinal.estimate,beta.continuous.estimate,delta.estimates,
                         start.values.sigma.rand.estimate,sigma22.estimate,lambda.estimate,
                         exact=x$exact,montecarlo=x$montecarlo, epsilon=epsilon,
                                     additional=FALSE) {  
                  } #foreach
  } else {boot=list(NA)
          
          for(i in 1:bootstrap.samples) {
#             random.int=mvrnorm(n=l, mu=rep(0,ncol(sigma.rand.estimate)), Sigma=sigma.rand.estimate)
#             
#             pred.fixed=sapply(1:mult.obs,function(obs) as.matrix(x$predictors.fixed[,,obs])%*%beta.estimate, simplify=TRUE)
#             if(ncol(sigma.rand.estimate)==1) pred.rand=sapply(1:mult.obs,function(obs) x$predictors.random[,,obs]*random.int, simplify=TRUE) else
#               pred.rand=sapply(1:mult.obs,function(obs) apply(x$predictors.random[,,obs]*random.int,1,sum),simplify="array")
#             
#             y1=pred.fixed+pred.rand+matrix(rnorm(l*mult.obs),ncol=mult.obs)
#             
#             data.ordinal.new=matrix(cut(y1,c(min(y1)-1,x$thresholds,max(y1)+1), labels=FALSE),ncol=mult.obs)
#             data.ordinal.new=ifelse(miss, NA, data.ordinal.new)
#             
            boot[[i]] <- ecm.ord.plus.cont(data.ordinal.new,data.continuous.new,x$predictors.fixed.ordinal,x$predictors.fixed.continous,
                                           x$predictors.random.ordinal,x$predictors.random.continuous,
                                           beta.ordinal.estimate,beta.continuous.estimate,delta.estimates,
                                           start.values.sigma.rand.estimate,sigma22.estimate,lambda.estimate,
                                           exact=x$exact,montecarlo=x$montecarlo, epsilon=epsilon,
                                           additional=FALSE)
          }  #for
  } #else  
  est=t(sapply(1:bootstrap.samples, function(i) 
    c(boot[[i]][[2]][lower.tri(boot[[i]][[2]],diag=TRUE)], boot[[i]][[3]][[1]], boot[[i]][[3]][[2]],
      boot[[i]][[4]], boot[[i]][[5]], boot[[i]][[6]], simplify=TRUE))
  res=var(est)
  nam <- matrix(paste("sigma ", outer(1:ncol(sigma.rand.estimate),1:ncol(sigma.rand.estimate),
                                      function(x,y) paste(x,y,sep="")),sep=""),ncol=ncol(sigma.rand.estimate))
  nam <- nam[lower.tri(nam,diag=TRUE)]
  colnames(res) <- c(nam, paste("Predictor random",1:length(beta.ordinal.estimate)), 
                  paste("Predictor continuous",1:length(beta.continuous.estimate)), 
                  paste("Diff",1:length(delta.estimates)), "Lambda", "Sigma 22")
  rownames(res) <- c(nam, paste("Predictor random",1:length(beta.ordinal.estimate)), 
                  paste("Predictor continuous",1:length(beta.continuous.estimate)), 
                  paste("Diff",1:length(delta.estimates)), "Lambda", "Sigma 22")
  
  res
} # function standard error
