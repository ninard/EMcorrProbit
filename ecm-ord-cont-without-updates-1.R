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
ecm.ord.plus.cont <- function(start.values.beta.ordinal,
                              start.values.beta.continuous,
                              start.values.delta,
                              start.values.sigma.rand,
                              start.values.sigma22, 
                              start.values.lambda,
                              data.ordinal,
                              data.continuous,
                              predictors.fixed.ordinal,
                              predictors.fixed.continous,
                              predictors.random.ordinal,
                              predictors.random.continuous, 
                              exact=F, montecarlo=100, epsilon=0.001,
                              additional=FALSE) {
  
  ########################################
  ### sigma.rand is the covaraince matrix of the random effects
  ##################################
  
  ### not needed in the package
  library(MASS)
  library(tmvtnorm)
  
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
 print(head(firstmomentb))
print(loglikelihood)
print(AIC)
print(BIC)


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
  
#  list(number.of.iterations=number.it,
#       sigma.matrix=sigma.rand.new,
#       beta.ordinal=betanew.ordinal,
#       beta.continuous=betanew.continuous,
#       delta=deltanew,
#       thresholds=cumsum(c(0,deltanew)),
#       lambda=lambdanew,
#       sigma22=sigma22new,
#       cov.errors=lambdanew*sigma22new,
#       sigma11=1+lambdanew^2*sigma22new       
#  )
  
  
  
} #function



### data simulation
## RIAS models

library(MASS)
random.int=mvrnorm(n = 5000, mu=c(0,0,0,0), Sigma=matrix(c(1,0,0.8,0,
0,.1,0,0,
0.8,0,1,-0.05,
0,0,-0.05,.1),ncol=4))
l=length(random.int[,1])
mult.obs=7
random.error=sapply(1:mult.obs, function(i) mvrnorm(n = l, mu=c(0,0), Sigma=matrix(c(4/3,2/sqrt(3),2/sqrt(3),4),ncol=2)),
simplify="array")
#l x 2 x mult.obs
int1=-0.5
int2=1
slope1=1
slope2=-.5
time=matrix(rep(1:mult.obs,l),byrow=T,ncol=mult.obs)
# l x mult.obs
y1=int1+slope1*time+sapply(1:mult.obs, function(i) random.int[,1:2]%*%c(1,i), simplify=T)+random.error[,1,]
y2=int2+slope2*time+sapply(1:mult.obs, function(i) random.int[,3:4]%*%c(1,i), simplify=T)+random.error[,2,]

data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.2,2,ifelse(y1<3,3,ifelse(y1<4.5,4,5))))
data.continuous=y2
table(data.ordinal)
data.ordinal[2,4]=NA
data.continuous[3,3]=NA

#predictors.random.ordinal=matrix(rep(1,l*mult.obs),ncol=mult.obs)
#predictors.random.continuous=predictors.random.ordinal
predictors.fixed.ordinal=sapply(1:mult.obs, function(i) cbind(1,time[,i]), simplify="array")
predictors.fixed.continuous=predictors.fixed.ordinal
predictors.random.ordinal=sapply(1:mult.obs, function(i) cbind(1,time[,i]), simplify="array")
predictors.random.continuous=predictors.random.ordinal

start.values.beta.ordinal=c(-0.5,1)
start.values.beta.continuous=c(1,-0.5)
start.values.delta=c(1.2,1.8,1.5)

start.values.sigma.rand=matrix(c(1,0,0.8,0,
0,.1,0,0,
0.8,0,1,-0.05,
0,0,-0.05,.1),ncol=4)

start.values.sigma22=4
start.values.lambda=1/(2*sqrt(3))
epsilon=0.001
exact=FALSE
montecarlo=100

ecm.ord.plus.cont(start.values.beta.ordinal,start.values.beta.continuous,start.values.delta,
start.values.sigma.rand,start.values.sigma22, start.values.lambda,
data.ordinal,data.continuous,predictors.fixed.ordinal,predictors.fixed.continous,
predictors.random.ordinal,predictors.random.continuous, exact=F,montecarlo=75, epsilon=0.002)->proba0002

#################################################
######## RI model
rm(list=ls())
library(MASS)
random.int=mvrnorm(n = 5000, mu=c(0,0), Sigma=matrix(c(1,0.8,
0.8,1),ncol=2))
l=length(random.int[,1])
mult.obs=7
random.error=sapply(1:mult.obs, function(i) mvrnorm(n = l, mu=c(0,0), Sigma=matrix(c(4/3,2/sqrt(3),2/sqrt(3),4),ncol=2)),
simplify="array")
#l x 2 x mult.obs
int1=-0.5
int2=1
slope1=1
slope2=-.5
time=matrix(rep(1:mult.obs,l),byrow=T,ncol=mult.obs)
# l x mult.obs
y1=int1+slope1*time+random.int[,1]+random.error[,1,]
y2=int2+slope2*time+random.int[,2]+random.error[,2,]

data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.2,2,ifelse(y1<3,3,ifelse(y1<4.5,4,5))))
data.continuous=y2
table(data.ordinal)
data.ordinal[2,4]=NA
data.continuous[3,3]=NA

#predictors.random.ordinal=matrix(rep(1,l*mult.obs),ncol=mult.obs)
#predictors.random.continuous=predictors.random.ordinal
predictors.fixed.ordinal=sapply(1:mult.obs, function(i) cbind(1,time[,i]), simplify="array")
predictors.fixed.continuous=predictors.fixed.ordinal
predictors.random.ordinal=sapply(1:mult.obs, function(i) matrix(1,ncol=1,nrow=l), simplify="array")
predictors.random.continuous=predictors.random.ordinal

start.values.beta.ordinal=c(-0.5,1)
start.values.beta.continuous=c(1,-0.5)
start.values.delta=c(1.2,1.8,1.5)

start.values.sigma.rand=matrix(c(1,0.8,
0.8,1),ncol=2)

start.values.sigma22=4
start.values.lambda=1/(2*sqrt(3))
epsilon=0.001
exact=FALSE
montecarlo=100


ecm.ord.plus.cont (start.values.beta.ordinal,start.values.beta.continuous,start.values.delta,
start.values.sigma.rand,start.values.sigma22, start.values.lambda,
data.ordinal,data.continuous,predictors.fixed.ordinal,predictors.fixed.continous,
predictors.random.ordinal,predictors.random.continuous, exact=F, montecarlo=75, epsilon=0.001)->proba0002


#################################################
######## RI model
rm(list=ls())
library(MASS)
random.int=mvrnorm(n = 1000, mu=c(0,0), Sigma=matrix(c(1,0.8,
                                                       0.8,1),ncol=2))
l=length(random.int[,1])
mult.obs=7
random.error=sapply(1:mult.obs, function(i) mvrnorm(n = l, mu=c(0,0), Sigma=matrix(c(4/3,2/sqrt(3),2/sqrt(3),4),ncol=2)),
                    simplify="array")
#l x 2 x mult.obs
int1=-0.5
int2=1
slope1=1
slope2=-.5
time=matrix(rep(1:mult.obs,l),byrow=T,ncol=mult.obs)
# l x mult.obs
y1=int1+random.int[,1]+random.error[,1,]
y2=int2+random.int[,2]+random.error[,2,]

data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.2,2,ifelse(y1<3,3,ifelse(y1<4.5,4,5))))
data.continuous=y2
table(data.ordinal)
data.ordinal[2,2:7]=NA
data.continuous[3,2:7]=NA

#predictors.random.ordinal=matrix(rep(1,l*mult.obs),ncol=mult.obs)
#predictors.random.continuous=predictors.random.ordinal
predictors.random.ordinal=sapply(1:mult.obs, function(i) matrix(1,ncol=1,nrow=l), simplify="array")
predictors.random.continuous=predictors.random.ordinal
predictors.fixed.ordinal=predictors.random.ordinal
predictors.fixed.continuous=predictors.fixed.ordinal

start.values.beta.ordinal=c(-0.5)
start.values.beta.continuous=c(1)
start.values.delta=c(1.2,1.8,1.5)

start.values.sigma.rand=matrix(c(1,0.8,
                                 0.8,1),ncol=2)

start.values.sigma22=4
start.values.lambda=1/(2*sqrt(3))
epsilon=0.001
exact=FALSE
montecarlo=100


ecm.ord.plus.cont (start.values.beta.ordinal,start.values.beta.continuous,start.values.delta,
                   start.values.sigma.rand,start.values.sigma22, start.values.lambda,
                   data.ordinal,data.continuous,predictors.fixed.ordinal,predictors.fixed.continous,
                   predictors.random.ordinal,predictors.random.continuous, exact=F, montecarlo=75, 
                   epsilon=0.001,additional=TRUE)->proba0002


## What if for some individual I have only 1 observation either on the ordinal or continuous response
# and all others are missing - not very logical to 

