#dlm to explain depth changes through time

#need to run "DWA updated.R" to generate response variables (density-weighted depth for species of interest)
DWAmonthly$DWA #response variable

#Need to run "phytoplankton DWA.R" to get phytoplankton predictor variables
allphyt<-ts(allDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(allDWA$Date[1], 1,4)),as.numeric(substr(allDWA$Date[1],6,7)))) # all groups predictor variable
diatom<-ts(diatomDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(diatomDWA$Date[1], 1,4)),as.numeric(substr(diatomDWA$Date[1],6,7)))) # all groups predictor variable
cyano<-ts(cyanoDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(cyanoDWA$Date[1], 1,4)),as.numeric(substr(cyanoDWA$Date[1],6,7)))) # all groups predictor variable
green<-ts(greenDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(greenDWA$Date[1], 1,4)),as.numeric(substr(greenDWA$Date[1],6,7)))) # all groups predictor variable
euglenophytes<-ts(euglenophytesDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(euglenophytesDWA$Date[1], 1,4)),as.numeric(substr(euglenophytesDWA$Date[1],6,7)))) # all groups predictor variable
chrysophyte<-ts(chrysophtyeDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(chrysophtyeDWA$Date[1], 1,4)),as.numeric(substr(chrysophtyeDWA$Date[1],6,7)))) # all groups predictor variable
cryptophyte<-ts(cryptophyteDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(cryptophyteDWA$Date[1], 1,4)),as.numeric(substr(cryptophyteDWA$Date[1],6,7)))) # all groups predictor variable
dinoflagellate<-ts(dinoflagellateDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(dinoflagellateDWA$Date[1], 1,4)),as.numeric(substr(dinoflagellateDWA$Date[1],6,7)))) 
bacteria<-ts(BacteriaDWA$DWA.mean,frequency=12, start=c(as.numeric(substr(BacteriaDWA$Date[1], 1,4)),as.numeric(substr(BacteriaDWA$Date[1],6,7)))) 

response<-ts(DWAmonthly$DWA, frequency=12, start=c(1955,1))
temppred<-ts(temperature$Temp,frequency=12, start=c(1948,1))
stratpred<-ts(stratexp$strat,frequency=12, start=c(1948,1))
chlpred<-ts(chlDWA$DWA.Mean, frequency=12, start=c(1990,6)) #not a good predictor
sechpred<-ts(secchi$Value.mean, frequency=12, start=c(1979,1)) 
spec.pgram(response,spans=c(3,3)) #seasonal oscillations?

#cross correlations to define lags
predictors<-c("temppred","stratpred","chlpred","sechpred","dayl","allphyt","diatom","cyano","green","euglenophytes","chrysophyte","cryptophyte","dinoflagellate","bacteria")
lagg<-c()
for (i in 1:length(predictors)){
tte<-ccf(response, get(predictors[i]),na.action=na.pass)
lagg[i]<-tte$lag[which(tte$lag>0)][which.max(tte$acf[which(tte$lag>0)])]
}
lags<-data.frame("predictor"=predictors, "lag"=lagg)


library(dynlm)
library (forecast)
tsdisplay(response)
seasonplot(response)
decompose(response)

#dlm 
mod1 <- dynlm(response ~ trend(response)+harmon(response,2))
summary(mod1)
  plot(response, type = "p")
  lines(fitted(mod1), col = 2)

a<-StructTS(response, type="trend")


mod2<-dynlm(response ~ temppred+stratpred)
summary(mod2)
plot(response, type = "p")
lines(fitted(mod2), col = 2)
vif(mod2)

mod3<-dynlm(response ~ L(stratpred, 2)*L(temppred,2)) #maybe green, def chrysophyte, def dinoflagellate
summary(mod3)
plot(response, type = "l")
lines(fitted(mod3), col = 2)
vif(mod3)

mod4<-dynlm(response ~ temppred+stratpred+sechpred) #maybe green, def chrysophyte, def dinoflagellate
summary(mod4)
plot(response, type = "p")
lines(fitted(mod3), col = 2)
vif(mod3)

anova(mod3,mod4)
AIC(mod3)
AIC(mod4)
edata("USDistLag", package = "lmtest")
dfm1 <- dynlm(consumption ~ gnp + L(consumption), data = USDistLag)
dfm2 <- dynlm(consumption ~ gnp + L(gnp), data = USDistLag)
plot(USDistLag[, "consumption"])
lines(fitted(dfm1), col = 2)
lines(fitted(dfm2), col = 4)
encomptest(dfm1, dfm2)










StructTS(response, type="BSM")







mod<-dlmModSeas(12, dV = 8.265)

modFilt <- dlmFilter(response, mod)
plot(response, type = 'o', col = "seagreen")
attach(modFilt)
v <- unlist(dlmSvd2var(U.C, D.C))
pl <- dropFirst(m) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(m) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown")
lines(pu, lty = 2, col = "brown")



buildFun <- function(x) { dlmModPoly(1, dV = exp(x[1]), dW = exp(x[2])) }
fit <- dlmMLE(response, parm = c(0,0), build = buildFun)
fit$conv
dlmresponse <- buildFun(fit$par)
V(dlmresponse)

#work needed to create predictor variables:

#Function to place NAs in missing dates for monthly time series. x is a dataframe with a "Date" column
NAdates<-function(x){
  years<-substr(x$Date,1,4)
  start = as.Date(x$Date[which.min(x$Date)])
  full <- seq(start, by="1 month", length=((max(as.numeric(years))+1)-min(as.numeric(years)))*12)
  allsp<-data.frame("Date"=full)
  for (i in 1:ncol(x)){
    t1<-data.frame("Count"=with(x, x[,i][match(full, as.Date(x$Date))]))
    allsp<-cbind(allsp, t1)}
  allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(x); allsp$Date<-full
  return(allsp)
}


#identify NAs and then fill with monthly average for whole time series
NAfill<-function(x){
  years<-substr(x$Date,1,4)
  start = as.Date(x$Date[which.min(x$Date)])
  full <- seq(start, by="1 month", length=((max(as.numeric(years))+1)-min(as.numeric(years)))*12)
  allsp<-data.frame("Date"=full)
  for (i in 1:ncol(x)){
    t1<-data.frame("Count"=with(x, x[,i][match(full, as.Date(x$Date))]))
    allsp<-cbind(allsp, t1)}
  allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(x); allsp$Date<-full
  #names(allsp)<-c("Date","strat")
  missing<-allsp$Date[which((is.na(allsp[,2]))==TRUE)] 
  for (i in 1:length(missing)){
    allsp[,2][which(allsp$Date==as.Date(missing[i]))]<-mean(allsp[,2][which(as.numeric(substr(allsp$Date,6,7))==as.numeric(substr(missing[i],6,7)))],na.rm=T)
  }
  return(allsp)
}

#NA fill for monthly time series objects

NAfillts<-function(x){
  years<-seq(from=start(x)[1], to=end(x)[1],by=1)
  startt = start(x)
  endd=end(x)
  full <- seq(from=as.Date(paste(startt[1], startt[2],"01",sep="-")), to=as.Date(paste(endd[1], endd[2],"01",sep="-")),by="1 month")
  allsp<-data.frame("Date"=full)
  temp<-data.frame("Date"=full, x)
  for (i in 1:ncol(temp)){
    t1<-data.frame("Count"=with(temp, temp[,i][match(full, as.Date(temp$Date))]))
    allsp<-cbind(allsp, t1) }
  allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(x); allsp$Date<-full
  #names(allsp)<-c("Date","strat")
  missing<-allsp$Date[which(is.na(allsp[,2])==T)] 
  for (i in 1:length(missing)){
    allsp[,2][which(allsp$Date==as.Date(missing[i]))]<-mean(allsp[,2][which(as.numeric(substr(allsp$Date,6,7))==as.numeric(substr(missing[i],6,7)))],na.rm=T)
  }
  return(ts(allsp[,2],start=start(x),frequency=12))
}

#function to return proper dates from a time series object
tsdates<-function(x){
  years<-seq(from=start(x)[1], to=end(x)[1],by=1)
  startt = start(x)
  endd=end(x)
  full <- seq(from=as.Date(paste(startt[1], startt[2],"01",sep="-")), to=as.Date(paste(endd[1], endd[2],"01",sep="-")),by="1 month")
  return(full)}


#Secchi depth
physical<-read.csv(file="chlorophyll_and_secchi_date2.csv");physical$Date<-as.Date(physical$Date,format="%d/%m/%y")
secchi<-physical[which(physical$The.code.parameter==1),]
secchi$Date<-paste(substr(secchi$Date,1,7),"-01",sep="")
library(doBy)
secchi<-summaryBy(Value~Date,FUN=mean,data=secchi)
secchi$Date<-(as.Date(secchi$Date))
secchi<-NAfill(secchi)

# average monthly daylength
library(geosphere)
dl <- daylength(51.9018, 1:365)
dayl<-tapply(dl, rep(1:12, c(31,28,31,30,31,30,31,31,30,31,30,31)), mean)
dayl<-ts(rep(dayl, 55), frequency=12, start=c(1948,1))

#temperature
temperature<-read.csv(file="temperature19482002.csv", stringsAsFactors=F)
temperature$DATE<-paste(substr(temperature$DATE,1,8),"01",sep="")
temperature$DEPTH<-as.numeric(temperature$DEPTH)

par(mfrow=c(3,3))
for(i in c(0,5,10, 20,25,30,50)){
temperature2<-temperature[which(temperature$DEPTH==i),]
library(doBy)
meanfun<-function(x) {mean(x,na.rm=T)}
temperature2<-summaryBy(TEMPERATURE~DATE,FUN=mean,data=temperature2)
temperature2$DATE<-(as.Date(temperature2$DATE))
names(temperature2)<-c("Date","Temp")
temperature2<-NAfill(temperature2)
plot(temperature2$Temp~temperature2$Date)
temperature2$Depth<-rep(i, nrow(temperature2))
if(i==0){templist<-data.frame("Date"=temperature2$Date, "Temp"=temperature2$Temp, "Depth"=temperature2$Depth)}
if(i>0){templist<-rbind(templist,temperature2)}
}
#templist<-templist[which(as.numeric(substr(templist$Date,6,7))==7|as.numeric(substr(templist$Date,6,7))==8|as.numeric(substr(templist$Date,6,7))==9),]
templist<-templist[which(as.numeric(substr(templist$Date,6,7))==8),]

#templist$Date<-paste(substr(templist$Date,1,5),"01-01",sep="")
templist<-summaryBy(Temp~Date+Depth,data=templist, FUN=mean)

record<-summaryBy(Temp~Date,data=templist,FUN=length)

dates<-unique(templist$Date)
for (i in 1:length(dates)){
  templist$Temp[which(templist$Date==dates[i])]<-templist$Temp[which(templist$Date==dates[i])]-(templist$Temp[which(templist$Date==dates[i])][which(templist$Depth[which(templist$Date==dates[i])]==0)])
}

colrs<-c(rep("red",10),rep("blue",10),rep("green",10),rep("orange",10),rep("purple",10))
plot(templist$Depth[which(templist$Date==dates[1])]~templist$Temp[which(templist$Date==dates[1])],type="l",ylim=c(50,0), xlim=c(-15,0),col="red")
for(i in 2:50){
lines(templist$Depth[which(templist$Date==dates[i])]~templist$Temp[which(templist$Date==dates[i])], col=colrs[i])
}

#quarterly temperatures
temperature<-read.csv(file="temperature19482002.csv", stringsAsFactors=F)
temperature$DATE<-paste(substr(temperature$DATE,1,8),"01",sep="")
temperature$DEPTH<-as.numeric(temperature$DEPTH)

tempdepth=0
temperature2<-temperature[which(temperature$DEPTH==tempdepth),]
library(doBy)
meanfun<-function(x) {mean(x,na.rm=T)}
temperature2<-summaryBy(TEMPERATURE~DATE,FUN=mean,data=temperature2)
temperature2$DATE<-(as.Date(temperature2$DATE))
names(temperature2)<-c("Date","Temp")
temperature2<-NAfill(temperature2)
temperature2$Month<-as.numeric(substr(temperature2$Date,6,7))
temperature2$Year<-as.numeric(substr(temperature2$Date,1,4))
temperature3<-temperature2[which(temperature2$Month==7|temperature2$Month==8|temperature2$Month==9),]
temperature4<-summaryBy(Temp~Year,FUN=mean,data=temperature3)

if(tempdepth==0){temp0<-temperature4}
if(tempdepth==50){temp50<-temperature4}

temps<-merge(temp0,temp50,by="Year"); names(temps)<-c("Year","temp0","temp50")
temps$diff<-temps$temp0-temps$temp50
plot(temps$diff~temps$Year, ylab="Difference surface to 50m (ºC)", pch=16, xlab="Year", cex=1.5,cex.axis=1.5, cex.lab=1.5) 
mod<-lm(temps$diff~temps$Year, na.action=na.exclude)
Years<-temps$Year[which(temps$Year!=1961)]
lines(fitted(mod)~temps$Year, lwd=7, col="darkgrey")
summary(mod)

#stratification intensity
temperature<-read.csv(file="temperature19482002.csv", stringsAsFactors=F)
temperature$DATE<-paste(substr(temperature$DATE,1,8),"01",sep="")
library(doBy)
temperature<-summaryBy(TEMPERATURE~DATE+DEPTH,FUN=mean,data=temperature)
temperature$DATE<-(as.Date(temperature$DATE))
temp0<-temperature[which(temperature$DEPTH==0),]
temp10<-temperature[which(temperature$DEPTH==10),]
temp25<-temperature[which(temperature$DEPTH==25),]
temp50<-temperature[which(temperature$DEPTH==50),]

names(temp0)<-c("Date","Depth","Temp")
names(temp50)<-c("Date","Depth","Temp")
temp0$Date<-as.Date(temp0$Date)
temp50$Date<-as.Date(temp50$Date)
stratexp<-merge(temp0,temp50,by="Date")
depths0.50<-data.frame("dep0"=stratexp$Temp.x,"dep50"=stratexp$Temp.y)
stratexp$strat<-apply(depths0.50,FUN=function(x){sd(x,na.rm=F)},MARGIN=1)
stratexp<-data.frame("Date"=stratexp$Date, "strat"=stratexp$strat)
stratexp<-NAfill(stratexp)


#concentration weighted chl-a depths
physical<-read.csv(file="chlorophyll_and_secchi_date2.csv");physical$Date<-as.Date(physical$Date,format="%d/%m/%y")
chl<-physical[which(physical$The.code.parameter==3),]
chl$Horizon[which(chl$Horizon==0)]<-1
chl<-chl[which(chl$Horizon==1|chl$Horizon==5|chl$Horizon==10|chl$Horizon==25|chl$Horizon==50|chl$Horizon==100|chl$Horizon==150|chl$Horizon==200),]
chl$DWA<-chl$Horizon*chl$Value
chl$Date<-paste(substr(chl$Date,1,7),"-01",sep="")

#figure out which dates have all depths
temp<-summaryBy(DWA~Date+Horizon,FUN=mean,data=chl)
plot(temp$DWA.mean[which(temp$Horizon==250)]~as.Date(temp$Date[which(temp$Horizon==250)]))
barplot(table(temp$Horizon))
temp2<-summaryBy(Horizon~Date,FUN=length,data=temp)
plot(temp2[,2]~as.Date(temp2[,1]))
temp3<-temp2[which(temp2$Horizon.length==8),]
chl<-chl[which(chl$Date%in%temp3$Date==T),]

library(doBy)
Mean<-function(x){mean(x, na.rm=T)}
chlDWA<-summaryBy(DWA~Date, data=chl, FUN=c(Mean,length), order=TRUE)
chlDWA$Date<-as.Date(chlDWA$Date)
chlDWA<-chlDWA[which(chlDWA$Date>as.Date("1980-01-01")),]

#get summer only
chlDWA$Month<-as.numeric(substr(chlDWA$Date,6,7))
chlDWA$Year<-as.numeric(substr(chlDWA$Date,1,4))
chlDWA<-chlDWA[which(chlDWA$Month==7|chlDWA$Month==8|chlDWA$Month==9),]
chlDWA<-summaryBy(DWA.Mean~Year,FUN=Mean,data=chlDWA)

plot(chlDWA[,2]~chlDWA[,1])

chlDWA<-NAfill(chlDWA)


physical<-read.csv(file="chlorophyll_and_secchi_date2.csv");physical$Date<-as.Date(physical$Date,format="%d/%m/%y")
chl<-physical[which(physical$The.code.parameter==3),]
chl$DWA<-chl$Horizon*chl$Value


library(MARSS)