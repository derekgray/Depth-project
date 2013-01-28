#Decomposition with Kalman Smoother

response<-read.csv(file.choose())

library(dlm)
Depth <- ts(response, frequency=12, start=c(1948,1))

#designate model
dlmdepth <- dlmModPoly(order=1) + dlmModSeas(frequency=12)
buildFun <- function(x) {
  diag(W(dlmdepth))[2:3] <- exp(x[1:2])
  V(dlmdepth)<-exp(x[3])
  return(dlmdepth)
}

#determine parameters for model using maximum likelihood
(fit <- dlmMLE(Depth, parm = rep(0, 3), build = buildFun))$conv

#construct model using estimated parameters
dlmdepth <- buildFun(fit$par)

#Use Kalman smoother to estimate state variables
DepthSmooth <- dlmSmooth(Depth, mod = dlmdepth)

#Bind data together and plot
x <- cbind(Depth, dropFirst(DepthSmooth$s[,c(1,3)]),(Depth-rowSums(dropFirst(DepthSmooth$s[,c(1,3)]))))
colnames(x) <- c("DWA", "Trend", "Seasonal","Residual")
plot(x, type = 'o', main = "Density-weighted depth")

#another method (seasonal decomposition by Loess)
plot(stl(response, s.window="periodic"))