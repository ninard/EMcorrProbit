############################################################
### RI model, 2 predictors, 3-level ordinal variable #######
### 750 individuals - 2 observations per subject ###########
############################################################
rm(list=ls())
random.int=rnorm(750,0,0.1)
l=length(random.int)
int=-0.5
mult.obs=2
y1=sapply(1:mult.obs,function(i) random.int+int+i+rnorm(l), simplify="array")
data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.5,2,3))
#table(data.ordinal)
#head(data.ordinal)
time=sapply(1:mult.obs, function(i) rep(i,l), simplify="array")
predictors.fixed=sapply(1:mult.obs, function(i) cbind(1,time[,i]), simplify="array")
predictors.random=sapply(1:mult.obs, function(i) matrix(rep(1,l),ncol=1), simplify="array")

sigma.rand=matrix(.01)
beta=c(-0.55,0.95)
delta=c(1.5)
mc=500
e=T

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta=beta,start.values.delta=delta,
                    start.values.sigma.rand=sigma.rand,
                    exact=e,montecarlo=mc,epsilon=.0002)
fffnew
cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])


fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta=beta,start.values.delta=delta,
                    start.values.sigma.rand=sigma.rand,
                    exact=F,montecarlo=mc,epsilon=.0002)

cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])



fffcm=ecm.one.ordinal.complete.cases(data.ordinal,predictors.fixed,
                    predictors.random,
                    start.values.beta=beta,start.values.delta=delta,
                    start.values.sigma.rand=sigma.rand,
                    exact=e,montecarlo=mc,epsilon=.0002)

cat(fffcm[[1]],fffcm[[2]],fffcm[[3]])

fffcm=ecm.one.ordinal.complete.cases(data.ordinal,predictors.fixed,
                                     predictors.random,
                                     start.values.beta=beta,start.values.delta=delta,
                                     start.values.sigma.rand=sigma.rand,
                                     exact=F,montecarlo=mc,epsilon=.0002)

cat(fffcm[[1]],fffcm[[2]],fffcm[[3]])


data.ordinal[1,2]=NA
#head(data.ordinal)

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta=beta,start.values.delta=delta,
                    start.values.sigma.rand=sigma.rand,
                    exact=F,montecarlo=mc,epsilon=.0002)
cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])


#######################################################################
############################################################
### RI model, 2 predictors, 4-level ordinal variable #######
### 500 individuals - 7 observations per subject ###########
############################################################
rm(list=ls())
l=500
mult.obs=7
random.int=rnorm(l,0,0.1)
int=-2

y1=sapply(1:mult.obs,function(i) random.int+i+int+rnorm(l), simplify="array")
data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.5,2,ifelse(y1<=3,3,4)))

#head(data.ordinal)
#data.ordinal
#table(data.ordinal)
time=sapply(1:mult.obs, function(i) rep(i,l), simplify="array")
predictors.fixed=sapply(1:mult.obs, function(i) cbind(1,time[,i]), simplify="array")
predictors.random=sapply(1:mult.obs, function(i) matrix(rep(1,l),ncol=1), simplify="array")

start.values.beta=c(-2,1)
start.values.delta=c(1.5,1.5)
start.values.sigma.rand=matrix(.01)
montecarlo=200
exact=F

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
               xrand=predictors.random,
               start.values.beta,start.values.delta,start.values.sigma.rand,
               exact=F,montecarlo=200,epsilon=.0005)
cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])

fffcm=ecm.one.ordinal.complete.cases(data.ordinal,predictors.fixed,
                                     predictors.random,
                                     start.values.beta,start.values.delta,
                                     start.values.sigma.rand,
                                     exact=F,montecarlo=200,epsilon=.0005)

cat(fffcm[[1]],fffcm[[2]],fffcm[[3]])


data.ordinal[1,2]=NA
#head(data.ordinal)

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta,start.values.delta,
                    start.values.sigma.rand,
                    exact=F,montecarlo=200,epsilon=.0005)
cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])






#####################################################################
########### RIAS model, 4 predictors, 5-level ordinal variable ######
### 2500 individuals - 10 observations per subject ##################
#####################################################################
rm(list=ls())
library(MASS)
l=2500
mult.obs=10
random.effects=mvrnorm(n=l, mu=c(0,0),Sigma=matrix(c(0.01,0.0005,0.0005,
                                                        0.0001),ncol=2))
time=1:mult.obs

beta=c(0.5,1,2,-2.5)
pred1=rnorm(l)
pred2=rexp(l)
predictors.fixed=sapply(1:mult.obs, function(i) cbind(1,rep(i,l),pred1[i],pred2[i]),
                        simplify="array")
y1=t(sapply(1:l, function(j) sapply(1:mult.obs, function(i) c(1,time[i])%*%
                                      random.effects[i,]+predictors.fixed[j,,i]%*%beta+rnorm(1)), simplify="array"))
data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.5,2,ifelse(y1<=3,3,ifelse(y1<=4,4,5))))
#table(data.ordinal)
predictors.random=sapply(1:mult.obs,function(i) cbind(1,rep(i,l)),simplify="array")

start.values.beta=c(0.5,1,2,-2.5)
start.values.delta=c(1.5,1.5,1)
start.values.sigma.rand=matrix(c(0.01,0.0005,0.0005,0.0001),ncol=2)

montecarlo=150
exact=F


a=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
               xrand=predictors.random,
               start.values.beta,start.values.delta,start.values.sigma.rand,
               exact=exact,montecarlo=montecarlo,epsilon=.005)

cat(a[[1]],a[[2]],a[[3]])

fffcm=ecm.one.ordinal.complete.cases(data.ordinal,predictors.fixed,
                                     predictors.random,
                                     start.values.beta,start.values.delta,
                                     start.values.sigma.rand,
                                     exact=F,montecarlo=100,epsilon=.005)

cat(fffcm[[1]],fffcm[[2]],fffcm[[3]])


data.ordinal[1,2]=NA
#head(data.ordinal)

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta,start.values.delta,
                    start.values.sigma.rand,
                    exact=F,montecarlo=100,epsilon=.005)
cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])



#######################################################################
############################################################
########## only intercept ##################################
### RI model, 1 predictor, 5-level ordinal variable #######
### 500 individuals - 7 observations per subject ###########
############################################################
l=500
mult.obs=7
random.int=rnorm(l,0,0.1)
int=2

y1=sapply(1:mult.obs,function(i) random.int+int+rnorm(l), simplify="array")
data.ordinal=ifelse(y1<=0,1,ifelse(y1<=1.5,2,ifelse(y1<=3,3,ifelse(y1<=4,4,5))))

head(data.ordinal)
#data.ordinal
table(data.ordinal)

predictors.random=sapply(1:mult.obs, function(i) matrix(rep(1,l),ncol=1), simplify="array")
predictors.fixed=predictors.random
start.values.beta=c(2.05)
start.values.delta=c(1.5,1.5,1)
start.values.sigma.rand=matrix(.01)
montecarlo=200
exact=F

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
               xrand=predictors.random,
               start.values.beta,start.values.delta,start.values.sigma.rand,
               exact=F,montecarlo,epsilon=.0005)

cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])

fffcm=ecm.one.ordinal.complete.cases(data.ordinal,predictors.fixed,
                                     predictors.random,
                                     start.values.beta,start.values.delta,
                                     start.values.sigma.rand,
                                     exact=F,montecarlo,epsilon=.0005)

cat(fffcm[[1]],fffcm[[2]],fffcm[[3]])


data.ordinal[1,2]=NA
#head(data.ordinal)

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta,start.values.delta,
                    start.values.sigma.rand,
                    exact=F,montecarlo=100,epsilon=.0005)
cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])


#######################################################################
######## binary data ###############################
### RI model, 2 predictors, 2-level ordinal variable #######
### 750 individuals - 5 observations per subject ###########
############################################################


random.int=rnorm(750,0,0.1)
l=length(random.int)
int=-1.5
mult.obs=5
y1=sapply(1:mult.obs,function(i) random.int+int+i+rnorm(l), simplify="array")
data.ordinal=ifelse(y1<=0,1,2)
table(data.ordinal)
head(data.ordinal)
time=sapply(1:mult.obs, function(i) rep(i,l), simplify="array")
predictors.fixed=sapply(1:mult.obs, function(i) cbind(1,time[,i]), simplify="array")
predictors.random=sapply(1:mult.obs, function(i) matrix(rep(1,l),ncol=1), simplify="array")
start.values.sigma.rand=matrix(.01)

start.values.beta=c(-1.55,0.995)
start.values.delta=NULL
montecarlo=200
exact=F

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
               xrand=predictors.random,
               start.values.beta,start.values.delta,start.values.sigma.rand,
               exact=F,montecarlo=200,epsilon=.0005)
cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])

fffcm=ecm.one.ordinal.complete.cases(data.ordinal,predictors.fixed,
                                     predictors.random,
                                     start.values.beta,start.values.delta,
                                     start.values.sigma.rand,
                                     exact=F,montecarlo,epsilon=.0005)

cat(fffcm[[1]],fffcm[[2]],fffcm[[3]])


data.ordinal[1,2]=NA
#head(data.ordinal)

fffnew=emcorrprobit(y=data.ordinal,xfixed=predictors.fixed,
                    xrand=predictors.random,
                    start.values.beta,start.values.delta,
                    start.values.sigma.rand,
                    exact=F,montecarlo=100,epsilon=.0005)
cat(fffnew[[1]],fffnew[[2]],fffnew[[3]])
