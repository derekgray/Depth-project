#1 The random walk with noise

randomwalk <- function(x) dlmModPoly(1, dV = x[1], dW = x[2])
StructTS(response, "level") #level is W, epsilon is V

randomwalk.fit <- dlmMLE(response, parm = c(1, sd(response)),build = randomwalk)
mod <- randomwalk(randomwalk.fit$par)
unlist(randomwalk(randomwalk.fit$par)[c("V", "W")])

bob<-dlmFilter(response,mod)
par(mfrow=c(2,2))
plot(dropFirst(bob$m), main="filtered state") #filtered values of state vector
plot(dropFirst(bob$f),main="one step ahead forecast") #one step ahead forecast
plot(dropFirst(bob$a), main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input


AICdlm(filteredobject="bob",MLEfitobject="randomwalk.fit",responsedata="response",burn_in=10)

#2 random walk with seasonal component
library(dlm)
randomwalk.seasonal <- function(x) dlmModPoly(order=1, dV = exp(x[1]), dW = exp(x[2])) +dlmModSeas(frequency=12, dV=sd(response))
randomwalk.seasonal.fit<-dlmMLE(response,par=rep(0,2), build=randomwalk.seasonal, lower=c(1e-6,0))
mod.randomwalk.seasonal <- randomwalk.seasonal(randomwalk.seasonal.fit$par)

unlist(randomwalk.seasonal(randomwalk.seasonal.fit$par)[c("V", "W")])
randomwalk.seasonal.filtered<-dlmFilter(response,mod.randomwalk.seasonal)
plot(randomwalk.seasonal.filtered$f) #plot one step ahead forecasts
plot(randomwalk.seasonal.filtered$m[,1:3]) #filtered values of state vectors
lines(response,col="red")
AICdlm(filteredobject="randomwalk.seasonal.filtered",MLEfitobject="randomwalk.seasonal.fit",responsedata="response",burn_in=10)

#for comparison, try the StructTS function
bob2<-StructTS(response,type="trend") #local trend model, V=epsilon, W=level
plot(bob2$fitted)
bob3<-StructTS(response,type="BSM") #local trend plus seasonal component
plot(bob3$fitted)


#3 local trend model with seasonal component
library(dlm)
localtrend.seasonal <- function(x) dlmModPoly(order=2, dV = exp(x[1]), dW = c(rep(exp(x[2]),2))) +dlmModSeas(frequency=12, dV=sd(response))
localtrend.seasonal.fit<-dlmMLE(response,par=rep(0,2), build=localtrend.seasonal, lower=c(1e-6,0))
mod <- localtrend.seasonal(localtrend.seasonal.fit$par)

unlist(localtrend.seasonal(localtrend.seasonal.fit$par)[c("V", "W")])
localtrend.seasonal.filtered<-dlmFilter(response,mod)
plot(localtrend.seasonal.filtered$f) #plot one step ahead forecasts
plot(localtrend.seasonal.filtered$m[,1:3]) #filtered values of state vectors
lines(response,col="red")
AICdlm(filteredobject="localtrend.seasonal.filtered",MLEfitobject="localtrend.seasonal.fit",responsedata="response",burn_in=10)

#for comparison, try the StructTS function
bob2<-StructTS(response,type="trend") #local trend model, V=epsilon, W=level
plot(bob2$fitted)
bob3<-StructTS(response,type="BSM") #local trend plus seasonal component
plot(bob3$fitted)



#4 Just temperature

datas<-ts.intersect(response,temppred)
response<-datas[,1]
temppred<-datas[,2]

temp.function <- function(x) dlmModReg(temppred,dW=sd(temppred),dV=exp(x[2]),addInt=FALSE)
temp.fit <- dlmMLE(response,par=c(sd(temppred), 1), build=temp.function)
mod.temp <- temp.function(temp.fit$par)

bob.temp<-dlmFilter(response,mod.temp) 
unlist(temp.function(temp.fit$par)[c("V", "W")])
par(mfrow=c(2,2))
plot(bob$m, main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a, main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob.temp",MLEfitobject="temp.fit",responsedata="response",burn_in=10)

outS <- dlmSmooth(response, mod)
plot(dropFirst(outS$s))


#5 Temperature plus seasonal component

datas<-ts.intersect(response,temppred)
response<-datas[,1]
temppred<-datas[,2]
library(vegan)
tempstand<-decostand(temppred, method="standardize")

temp.seasonal.function <- function(x) dlmModReg(tempstand,dW=sd(tempstand),dV=exp(x[2]),addInt=FALSE) +dlmModSeas(frequency=12,dV=sd(response))
temp.seasonal.fit <- dlmMLE(response,par=c(sd(tempstand), 1), build=temp.seasonal.function)
mod.temp.seasonal <- temp.seasonal.function(temp.seasonal.fit$par)

bob.temp.seasonal<-dlmFilter(response,mod.temp.seasonal) 
unlist(temp.seasonal.function(temp.seasonal.fit$par)[c("V", "W")])
par(mfrow=c(2,2))
plot(bob$m[,1:3], main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a[,1:3], main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob.temp.seasonal",MLEfitobject="temp.seasonal.fit",responsedata="response",burn_in=10)


#6 stratification plus seasonal component

datas<-ts.intersect(response,stratpred)
response<-datas[,1]
stratpred<-datas[,2]
library(vegan)
stratstand<-decostand(stratpred, method="standardize")

temp.seasonal.function <- function(x) dlmModReg(stratstand,dW=sd(stratstand),dV=exp(x[2]),addInt=FALSE) +dlmModSeas(frequency=12,dV=sd(response))
temp.seasonal.fit <- dlmMLE(response,par=c(sd(stratstand), 1), build=temp.seasonal.function)
mod.temp.seasonal <- temp.seasonal.function(temp.seasonal.fit$par)

bob.temp.seasonal<-dlmFilter(response,mod.temp.seasonal) 
unlist(temp.seasonal.function(temp.seasonal.fit$par)[c("V", "W")])
par(mfrow=c(2,2))
plot(bob$m[,1:3], main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a[,1:3], main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob.temp.seasonal",MLEfitobject="temp.seasonal.fit",responsedata="response",burn_in=10)


#2 random walk plus avg temperature in top 250m----------------------------------

datas<-ts.intersect(response,temppred)
temppred<-datas[,2]
response<-datas[,1]

#combine models
random.temp.function <- function(x) dlmModPoly(1, dV = x[1], dW = x[2])+ dlmModReg(temppred,dW=sd(temppred),addInt=FALSE)

random.temp.fit <- dlmMLE(response,par=c(1,sd(response),sd(temppred)), build=random.temp.function)
mod <- random.temp.function(random.temp.fit$par)

#try the model with discount factors between 0.9 and 0.99
res<-list()
for (DFA in c(seq(from=0.9,to=1,by=0.01))){
  bob<-dlmFilterDF(response,mod,DF=DFA)
  res[[as.character(DFA)]]<-AICdlm(filteredobject="bob",MLEfitobject="random.temp.fit",responsedata="response",burn_in=10)
}   

bob<-dlmFilterDF(response,mod, DF=0.98) #0.98 best discount

par(mfrow=c(2,2))
plot(bob$m, main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a, main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob1",MLEfitobject="random.temp.fit",responsedata="response",burn_in=10)

#3 random walk plus avg temperature in top 250m plus day length----------------------------------

random.temp.function <- function(x) dlmModPoly(order=1, dV = x[1], dW = x[2])+ dlmModReg(temppred,dW=sd(temppred), addInt=FALSE)

random.temp.fit <- dlmMLE(response,par=c(1,sd(response),sd(temppred)), build=random.temp.function)
mod <- random.temp.function(random.temp.fit$par)
bob<-dlmFilter(response,mod)

par(mfrow=c(2,2))
plot(bob$m, main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a, main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob",MLEfitobject="random.temp.fit",responsedata="response",burn_in=10)

#4a random walk plus seasonality plus statification intensity----------------------------------
datas<-ts.intersect(response,stratpred)
response<-datas[,1]; stratpred<-datas[,2]

random.temp.function <- function(x) dlmModPoly(order=1, dV = x[1], dW = x[2])+ dlmModSeas(frequency=12)

random.temp.fit <- dlmMLE(response,par=c(1,sd(response)), build=random.temp.function)
mod <- random.temp.function(random.temp.fit$par)
bob<-dlmFilter(response,mod)

par(mfrow=c(2,2))
plot(bob$m, main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a, main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob",MLEfitobject="random.temp.fit",responsedata="response",burn_in=10)



#4c random walk plus seasonality plus phytoplankton DWA----------------------------------
datas<-ts.intersect(response,NAfillts(allphyt))
response<-datas[,1]; phytpred<-datas[,2]

random.temp.function <- function(x) dlmModReg(log10(phytpred+1),dW=1, addInt=FALSE)

random.temp.fit <- dlmMLE(response,par=0, build=random.temp.function)
mod <- random.temp.function(random.temp.fit$par)
bob<-dlmFilter(response,mod)

par(mfrow=c(2,2))
plot(bob$m, main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a, main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob",MLEfitobject="random.temp.fit",responsedata="response",burn_in=10)

#4b random walk plus seasonality plus statification intensity----------------------------------
datas<-ts.intersect(response,stratpred)
response<-datas[,1]; stratpred<-datas[,2]

random.temp.function <- function(x) dlmModPoly(order=1, dV = x[1], dW = x[2])+ dlmModSeas(frequency=12) + dlmModReg(stratpred,dW=sd(stratpred), addInt=FALSE)

random.temp.fit <- dlmMLE(response,par=c(1,sd(response),sd(stratpred)), build=random.temp.function)
mod <- random.temp.function(random.temp.fit$par)
bob<-dlmFilter(response,mod)

par(mfrow=c(2,2))
plot(bob$m, main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a, main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob",MLEfitobject="random.temp.fit",responsedata="response",burn_in=10)



#function to calculate statistics for model selection
AICdlm<-function(filteredobject,responsedata,MLEfitobject,burn_in=10){
MSE<-c(); MAD<-c(); MAPE<-c(); U<-c(); logLik<-c(); k<-c()
linearTrend_resid <- get(filteredobject)$y-get(filteredobject)$f
linearTrend_resid <- tail(linearTrend_resid, -burn_in)
MSE[as.character(MLEfitobject)] <- mean(linearTrend_resid^2) #mean squared error
MAD[as.character(MLEfitobject)] <- mean(abs(linearTrend_resid)) #mean absolute deviation (measure of dispersion)
MAPE[as.character(MLEfitobject)] <- mean(abs(linearTrend_resid) / tail(get(responsedata), -burn_in)) #Mean absolute percentage error, 0 is a perfect fit
U[as.character(MLEfitobject)] <- sqrt(mean(linearTrend_resid^2) / mean(tail(diff(get(responsedata))^2, -(burn_in-1)))) #Theil's U, < 1 better than guessing, 1 means same as guessing, >1 worse than guessing 
logLik[as.character(MLEfitobject)] <- -get(MLEfitobject)$value
k[as.character(MLEfitobject)] <- length(get(MLEfitobject)$par)
AIC <- -2 * (logLik - k) # vectorized, for all models at once
result<-data.frame(MSE,MAD,MAPE,U,logLik,k,AIC)
return(result)
}

#Function to use the dlmFilter with varying discount factors.-------------------- 
dlmFilterDF <- function (y, mod, simplify = FALSE, DF) 
{
  ## storage.mode(y) <- "double"
  mod1 <- mod
  yAttr <- attributes(y)
  ytsp <- tsp(y)
  y <- as.matrix(y)
  timeNames <- dimnames(y)[[1]]
  stateNames <- names(mod$m0)
  m <- rbind(mod$m0, matrix(0, nr = nrow(y), nc = length(mod$m0)))
  a <- matrix(0, nr = nrow(y), nc = length(mod$m0))
  f <- matrix(0, nr = nrow(y), nc = ncol(y))
  U.C <- vector(1 + nrow(y), mode = "list")
  D.C <- matrix(0, 1 + nrow(y), length(mod$m0))
  U.R <- vector(nrow(y), mode = "list")
  D.R <- matrix(0, nrow(y), length(mod$m0))
  U.W <- vector(nrow(y), mode = "list")
  D.W <- matrix(0, nrow(y), length(mod$m0))
  Wliste <- vector(nrow(y), mode = "list")
  P <- vector(nrow(y), mode = "list")
  tmp <- La.svd(mod$V, nu = 0)
  Uv <- t(tmp$vt)
  Dv <- sqrt(tmp$d)
  Dv.inv <- 1/Dv
  Dv.inv[abs(Dv.inv) == Inf] <- 0
  sqrtVinv <- Dv.inv * t(Uv)
  sqrtV <- Dv * Uv
  tmp <- La.svd(mod$C0, nu = 0)
  U.C[[1]] <- t(tmp$vt)
  D.C[1, ] <- sqrt(tmp$d)
  for (i in seq(length = nrow(y))) {
    tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
    a[i, ] <- mod$GG %*% m[i, ]
    P[[i]] <- mod$GG %*% crossprod(D.C[i,] * t(U.C[[i]])) %*% t(mod$GG)
    Wliste[[i]] <- P[[i]]* ((1-DF)/DF)
    svdW <- La.svd( Wliste[[i]] , nu = 0)
    sqrtW <- sqrt(svdW$d) * svdW$vt
    U.W[[i]] <- t(svdW$vt)
    D.W[i, ] <- sqrt(svdW$d)
    tmp <- La.svd(rbind(D.C[i, ] * t(mod$GG %*% U.C[[i]]), 
                        sqrtW), nu = 0)
    U.R[[i]] <- t(tmp$vt)
    D.R[i, ] <- tmp$d
    f[i, ] <- mod$FF %*% a[i, ]
    D.Rinv <- 1/D.R[i, ]
    D.Rinv[abs(D.Rinv) == Inf] <- 0
    tmp <- La.svd(rbind(sqrtVinv %*% mod$FF %*% U.R[[i]], 
                        diag(x = D.Rinv, nrow = length(D.Rinv))), nu = 0)
    U.C[[i + 1]] <- U.R[[i]] %*% t(tmp$vt)
    foo <- 1/tmp$d
    foo[abs(foo) == Inf] <- 0
    D.C[i + 1, ] <- foo
    m[i + 1, ] <- a[i, ] + crossprod(D.C[i + 1, ] * t(U.C[[i
                                                           + 1]])) %*% tF.Vinv %*% as.matrix(y[i, ] - f[i,])
  }        
  m <- drop(m)
  a <- drop(a)
  f <- drop(f)
  attributes(f) <- yAttr
  ans <- list(m = m, U.C = U.C, D.C = D.C, a = a, U.R = U.R, 
              D.R = D.R, f = f, U.W=U.W, D.W=D.W)
  ans$m <- drop(ans$m)
  ans$a <- drop(ans$a)
  ans$f <- drop(ans$f)
  attributes(ans$f) <- yAttr
  if (!is.null(ytsp)) {
    tsp(ans$a) <- ytsp
    tsp(ans$m) <- c(ytsp[1] - 1/ytsp[3], ytsp[2:3])
    class(ans$a) <- class(ans$m) <- if (length(mod$m0) > 
                                          1) 
      c("mts", "ts")
    else "ts"
  }
  if (!(is.null(timeNames) && is.null(stateNames))) {
    dimnames(ans$a) <- list(timeNames, stateNames)
    dimnames(ans$m) <- list(if (is.null(timeNames)) NULL else c("", 
                                                                timeNames), stateNames)
  }
  if (simplify) 
    return(c(mod = list(mod1), ans))
  else {
    attributes(y) <- yAttr
    ans <- c(y = list(y), mod = list(mod1), ans)
    class(ans) <- "dlmFiltered"
    return(ans)
  }
}




bob<-c()
set.seed(1234)
r <- rnorm(100)
X <- r
u <- -1*X + 0.5*rnorm(100)
MyModel <- function(x)  dlmModReg(X, FALSE, dV = x[1]^2)
fit <- dlmMLE(u, parm = c(0.3), build = MyModel)
mod <- MyModel(fit$par)
bob<-dlmFilter(u,mod)
plot(bob$a)
points(bob$f,col="red")


StructTS(response, "level")
#level is W, epsilon is V




myMod <- dlmModReg(ts.intersect(response,temppred),addInt=F, dV = 0.3228)
                   bob<-dlmFilter(response, myMod)

#p is number of parameters
p=1
Akaike<-(-2)*(-dlmLL(response, myMod))+2*p

2*k-(2*log(-dlmLL(response, myMod)))

buildFun <- function(x) {dlmModReg(response, dV = 0.3228)}



mod <- dlmModPoly(1, dV = 16300, C0 = 1e8)
X(mod) <- matrix(0, nrow = length(Nile))
X(mod)[time(Nile) == 1899] <- 60580
JW(mod) <- 1

library(dlm)
?stdata<-BJsales.lead
BJmod <- dlmModReg(X = cbind(BJsales.lead, log(BJsales.lead)))
X(BJmod)
JFF(BJmod)
FF(BJmod)
bob<-dlmFilter(data, BJmod)
plot(data)
lines(bob$a)

http://definetti.uark.edu/~gpetris/UseR-2011/SSMwR-useR2011handout.pdf

FF(BJmod)

nileBuild <- function(par) {
  dlmModPoly(1, dV = exp(par[1]), dW = exp(par[2]))
}
nileMLE <- dlmMLE(Nile, rep(0,2), nileBuild); nileMLE$conv
nileMod <- nileBuild(nileMLE$par)
V(nileMod)
W(nileMod)
nileFilt <- dlmFilter(Nile, nileMod)
nileSmooth <- dlmSmooth(nileFilt)
plot(cbind(Nile, nileFilt$m[-1], nileSmooth$s[-1]), plot.type='s',
     col=c("black","red","blue"), ylab="Level", main="Nile river", lwd=c(1,2,2))

library(sspir)

phistart <- StructTS(UKgas)$coef
  
  c(3.7e-4,0,1.7e-5,7.1e-4)
StructTS(UKgas)

gasmodel <- ssm(log10(UKgas) ~ -1 + tvar(polytime(time,degree=1)) + tvar(sumseason(time,period=4)), phi=NULL)
fit <- getFit(gasmodel)
plot(fit$m[,1:3])

a<-ssm(response~ tvar(polytime(index(response),degree=1)))
fit<-getFit(a)
plot(fit$m[,2])
lines(response,col="red")

library(dlmodeler)


data(vandrivers)
vandrivers$y <- ts(vandrivers$y,start=1969,frequency=12)
vd.time <- time(vandrivers$y)
vd <- ssm( y ~ tvar(1) + seatbelt + sumseason(vd.time,12),
           family=poisson(link="log"),
           data=vandrivers,
           phi = c(1,0.0004),
           C0=diag(13)*100,
           fit=FALSE
)
phi(vd)["(Intercept)"] <- exp(- 2*3.703307 )
C0(vd) <- diag(13)*1000
vd.res <- kfs(vd)
plot( vd.res$m[,1:3] )


dwa.time <- time(response)
dwa.mod <- ssm( response ~ tvar(1) + tvar(temppred),family=gaussian,phi =NULL, fit=T)

dwa.mod.f<-kfilter(dwa.mod)
plot(dwa.mod$ss$y)
lines(dwa.mod.f$m[,2], col="red")

phi(vd)["(Intercept)"] <- exp(- 2*3.703307 )
C0(vd) <- diag(13)*1000
vd.res <- kfs(vd)
plot( vd.res$m[,1:3] )


data(kurit)
m1 <- SS(kurit)
phi(m1) <- c(100,5)
m0(m1) <- matrix(130)
C0(m1) <- matrix(400)
m1.f <- kfilter(m1)
plot(m1$y)
lines(m1.f$m,lty=2)


# generate some data
N <- 365*5
t <- c(1:N,rep(NA,365))
a <- rnorm(N+365,0,.5)
y <- pi + cos(2*pi*t/365.25) + .25*sin(2*pi*t/365.25*3) +
  exp(1)*a + rnorm(N+365,0,.5)
# build a model for this data

library(dlmodeler)
datas<-ts.intersect(response,stratpred)
resp<-datas[,1]
temp<-datas[,2]
m1 <- dlmodeler.build.polynomial(0,sigmaH=NA,sigmaQ=NA,name="level")
m4 <- dlmodeler.build.regression(temp,sigmaH=NA,name="reg")

m5 <- dlmodeler.add(m1,m4,name="mymodel")
fit<-dlmodeler.fit(resp,model=m1, method="MLE")
AIC(fit)
f <- dlmodeler.filter(resp, fit$model, raw.result=TRUE)
# extract all the components
lines(as.vector(f$f))

m.state.mean <- dlmodeler.extract(f,m4,type="state", value="mean")
m.state.cov <- dlmodeler.extract(f,m4,type="state", value="covariance")
m.obs.mean <- dlmodeler.extract(f,m4,type="observation",value="mean")
m.obs.cov <- dlmodeler.extract(f,m4,type="observation",value="covariance")
m.obs.int <- dlmodeler.extract(f,m4,type="observation",value="interval",prob=.99)


plot(as.vector(resp))
lines(m.obs.int$level$upper[1,],col="red")
lines(m.obs.int$level$lower[1,],col="red")
lines(m.obs.int$level$mean[1,],col=2)



#create test dataset
bob<-c()
set.seed(1234)
r <- rnorm(100)
u <- -5*X + 0.5*rnorm(100) #known regression parameters

#Fit regression model with predictor (r)
MyModel <- function(x)  dlmModReg(r, FALSE, dV = x[1]^2,dW=x[2]^2) #r is explantory variable
fit <- dlmMLE(u, parm = c(0.3,0.3), build = MyModel) #u is response variable
mod <- MyModel(fit$par)
unlist(MyModel(fit$par)[c("V", "W")])
bob<-dlmFilter(u,mod)
plot(bob$a) #plot regression coefficients through time
plot(bob$f)
AICdlm(filteredobject="bob",MLEfitobject="fit",responsedata="u",burn_in=10)
lines(u,col="red")

#random walk model
Mymodel2 <- function(x) dlmModPoly(order=1, dV = exp(x[1]), dW = exp(x[2]))

result<-list()
for (i in 1:10){
  fit2<-dlmMLE(u,par=rep(i-1,2), build=Mymodel2, lower=c(1e-6,0))
  result[[i]] <-unlist(Mymodel2(fit2$par)[c("V", "W")])
}
result

mod <- Mymodel2(fit2$par)

unlist(Mymodel2(fit2$par)[c("V", "W")])
bob<-dlmFilter(u,mod)
plot(bob$f) #plot regression coefficients through time
lines(u,col="red")
AICdlm(filteredobject="bob",MLEfitobject="fit2",responsedata="u",burn_in=10)

StructTS(u,type="level") #first order polynomial, V=epsilon, W=level
