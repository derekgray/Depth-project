################################################################################
###################Data processing for depth analyses###########################
########This code starts with Lizzie's zooplankton files in long################
########form and pulls the requested data out by depth, group, #################
########start date, species of interest, and/or lifestage ######################
########as specified at the top of the code. Caution should be##################
########used when looking at a particular species since species#################
########names are matched using the "grep" function that might #################
########return locations for species with similar names. In other###############
########words, it is necessary to inspect the species lists to make#############
########sure that the data summary is providing the desired output##############
################################################################################

######Output consists of these objects:
######1. "zooplankton": A time series by month (i.e. abundances are monthly averages)
######2. "species": a list of the unique species from the processed data
######Code requires the following packages: doBy, xts, zoo


########Before starting two files must be read in to the R environment: 
combined<-list()
combined.dates<-list()
for (depth in c("0to10","10to25","50to100","100to150","150to250","lessthan250")){
######## 1. The main data file containing all observations
# Check if main file has been loaded, if not read data into R
if(exists("zoodata")==FALSE){zoodata<-read.csv(file="zoopzeros250key1.csv", stringsAsFactors=F)}       #"c:/Users/gray/Dropbox/Postdoc/R data/zoopzeros250key1.csv")

######## 2. The zooplankton key
#Check if the key has been loaded yet, if not read it in 
if(exists("key")==FALSE){key<-read.csv(file="key.csv")}                            #"c:/Users/gray/Dropbox/Postdoc/R data/key.csv")

#######
#Enter depth of interest in meters (can be <, >, or range, e.g. "0to10")
#Depths available: 0to10, 10to25, 50to100, 100to150, 150to250, lessthan250 (previous 5 combined), "<251"=all depths less than 251m
#depth<-"0to10"

#start date of data extracted
startdate<-"1955-01-01"

#group of interest ("all", "Copepod", "Cladoceran", "Rotifer")
taxa<-"all"

#Enter species of interest ("all", or e.g. "Epischura baicalensis")
sppec<-"Epischura baicalensis"

#Enter lifestage of interest("all", "adult", "copep","naup") #This has not been tested for anything except for "all" or "adult" although many other options are in the key
stage<-"all"

#remove doublecounts? ("yes","no")
dbl<-"yes"

######################################################
#################start data processing################
######################################################

allzoo<-zoodata

if(taxa!="all"){
  allzoo<-allzoo[which(allzoo$Group_General==taxa),] }#take records for group of interest

if(stage!="all"){
  allzoo<-allzoo[which(allzoo$Lifestage_Cop==stage),] }#take records for lifestage of interest

#Take data only from species of interest (or take all data if an individual species wasn't specified)
if (sppec!="all"){
  allzoo<-allzoo[grep(sppec, paste(allzoo$Genus,allzoo$Species)),]}

#Take after requested date
allzoo<-allzoo[which(allzoo$Date1>=as.Date(startdate)),]

#create a volume column and total number of individuals column for later calculations
allzoo$Volume<-(as.numeric(allzoo$Nig_Gr)-as.numeric(allzoo$Ver_Gr))*(pi*(0.375/2)^2)*1000 # 0.375 diameter net, multiply by 1000 to convert to litres. Note NAs will be introduced due to some records having no depth information
allzoo<-allzoo[which(allzoo$Volume!="NA"),] #removes records with no volume
allzoo$Tot<-allzoo$Volume*allzoo$Count #Add a column that shows total number of individuals (volume sampled multiplied by Lyubov's densities (#/m3))

#create single depth column by merging Ver_Gr and Nig_Gr
allzoo$Depth<-paste(allzoo$Ver_Gr,allzoo$Nig_Gr, sep="-")

#create separate data sets for each depth category (0-10, 10-25, 50-100, 100-150, 150-250)
allzoo0to10<-allzoo[which(allzoo$Depth=="0-10"),]; allzoo0to10$Depth<-NULL
allzoo10to25<-allzoo[which(allzoo$Depth=="10-25"),]; allzoo10to25$Depth<-NULL
allzoo50to100<-allzoo[which(allzoo$Depth=="50-100"),]; allzoo50to100$Depth<-NULL
allzoo100to150<-allzoo[which(allzoo$Depth=="100-150"),]; allzoo100to150$Depth<-NULL
allzoo150to250<-allzoo[which(allzoo$Depth=="150-250"),]; allzoo150to250$Depth<-NULL
allzoolessthan250<-allzoo[which(allzoo$Depth=="0-10"|allzoo$Depth=="10-25"|allzoo$Depth=="50-100"|allzoo$Depth=="100-150"|allzoo$Depth=="150-250"),]

#check which depth levels were requested
if (charmatch(">", depth, nomatch=0)>0|charmatch("<", depth, nomatch=0)>0){
  allzoocommondepths<-allzoo[which(as.numeric(allzoo$Nig_Gr)<as.numeric(gsub("\\D", "", depth))),] 
}

if (charmatch(">", depth, nomatch=0)>0|charmatch("<", depth, nomatch=0)==0){
  allzoocommondepths<-get(paste("allzoo", depth, sep=""))}

#remove doublecounts
if(dbl=="yes"){
  KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")] #which codes are not doubles (according to Lizzie)
  allzoocommondepths<-allzoocommondepths[which(allzoocommondepths$KOD%in% KODsnotdouble),]
}

#create a list of species that entered analysis (useful for double-checking in the case where you only one species included)
species<-unique(paste(allzoo$Genus,allzoo$Species))

#create a column of genus +species
allzoocommondepths$genspp<-paste(allzoocommondepths$Genus,allzoocommondepths$Species,sep=" ")

#make all sample dates first of the month (makes it easier to process with the summaryBy function into monthly means)
allzoocommondepths$Date1<-paste(substr(allzoocommondepths$Date1,1,7),"-01",sep="") #convert all samples to 1st of month to make a regular time series

library(doBy)
allzoo3<-summaryBy(Tot+Volume~Date1+genspp, data=allzoocommondepths, FUN=sum, order=TRUE)
allzoo3$Density<-allzoo3$Tot.sum/allzoo3$Volume.sum #Figure out density by adding all individuals and dividing by total volume sampled
allzoo3<-data.frame("Date"=allzoo3$Date1, "Density"=allzoo3$Density, "species"=allzoo3$genspp) #Make a nice dataframe with Date, Density, and KOD
clustall<-reshape(allzoo3,idvar="Date",v.names="Density",timevar="species",direction="wide") #Convert from long form to wide form (species names as column names instead of row names)
names(clustall)<-(gsub("Density.", "", names(clustall)))  #take only numbers from column names. Necessary because text is introduced by the reshape command
row.names(clustall)<-unique(allzoo3$Date)[1:nrow(clustall)]; clustall[,1]<-NULL #Make rownames sample dates

if (ncol(clustall)!=1){
clustall<-clustall[,which(colSums(clustall,na.rm=T)>0)] }#get rid of columns that are all zeroes

zooplankton<-clustall
zooplankton$Date<-as.Date(row.names(zooplankton))

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

NAdates(zooplankton)

combined[[as.character(depth)]]<-zooplankton[,1]
combined.dates[[as.character(depth)]]<-zooplankton$Date}
#write.csv(zooplankton, file=paste("zooplankton",depth,".csv"),row.names=F) #write the result to a csv if desired


#plots
plot(unlist(combined[1])~as.Date(unlist(combined.dates[1])), type="l")
lines(unlist(combined[2])~as.Date(unlist(combined.dates[2])), col="red")
lines(unlist(combined[3])~as.Date(unlist(combined.dates[3])), col="green")
lines(unlist(combined[4])~as.Date(unlist(combined.dates[4])), col="blue")
lines(unlist(combined[4])~as.Date(unlist(combined.dates[4])), col="orange")

par(mfrow=c(2,3))
for (Depth.int in c("0to10","10to25","50to100","100to150","150to250","lessthan250")){
#convert to Julian date
dts<-as.Date(unlist(combined.dates[[Depth.int]]))
juldat<-as.POSIXlt(dts)$yday+1

#zoocol=c("red","blue","green","yellow","black","purple","orange","pink","cyan","brown")
years<-substr(as.Date(unlist(combined.dates[[Depth.int]])),1,4)
zoocol<-c(rep("red",9), rep("blue",9),rep("green",9),rep("yellow",9),rep("cyan",9))

x=1
plot(unlist(log10(combined[[Depth.int]])[which(years==1970)]+1)~juldat[which(years==1970)],col="red", ylim=c(0,max(log10(unlist(combined[[Depth.int]])+1),na.rm=T)), ylab="log(Abundance+1)", xlab="Julian Day", main=Depth.int)
for (i in 1955:2002){
  x=x+1
  try(points(log10(unlist(combined[[Depth.int]])[which(years==i)]+1)~juldat[which(years==i)],col="black"))
}

#Smooth plot of peaks by day
dt<-as.Date(unlist(combined.dates[[Depth.int]]))

#par(oma=c(1,3,1,1))
new<-data.frame("Count"=unlist(combined[[Depth.int]]),"Day"=juldat)
res<-summaryBy(Count~Day,data=new, FUN=mean)
te<-loess.smooth(res$Day,y=log10(res$Count.mean+1), span=2/6, evaluation=25)
#plot(log10(res$Count.mean+1)~res$Day, xlab="Day of the year", ylab="log (x+1) abundance per litre", cex.lab=1.5, cex.axis=1.5)
lines(te$y~te$x, type="l", col="red")
}


#time series plots
par(mfrow=c(2,3))

library(xts)
density.months<-list()
for (Depth.int in c("0to10","10to25","50to100","100to150","150to250","lessthan250")){
temp<-xts(log10(unlist(combined[[Depth.int]])+1), order.by=unlist(combined.dates[[Depth.int]]))  
locations<-endpoints(temp, "months")
density.months[[Depth.int]]<-period.apply(temp, INDEX=locations, FUN=function(x) mean(x,na.rm=T))
plot(density.months[[Depth.int]], type="l", main=Depth.int)}

#cross correlations
par(mfrow=c(2,2))
ccf(drop(density.months[["0to10"]]), drop(density.months[["10to25"]]), main="0to10 vs 10to25")
ccf(drop(density.months[["10to25"]]), drop(density.months[["50to100"]]), main="10to25 vs 50to100")
ccf(drop(density.months[["50to100"]]), drop(density.months[["100to150"]]), main="50to100 vs 100to150")
ccf(drop(density.months[["100to150"]]), drop(density.months[["150to250"]]), main="100to150 vs 150to250")

ccf(drop(density.months[["0to10"]]), drop(density.months[["50to100"]]), main="0to10 vs 50to100")



#########Experiments with DLM

#1 The random walk with noise

response<-ts(as.vector(density.months[["0to10"]]), frequency=12, start=c(1955,1))
plot(response)

library(dlm)
randomwalk <- function(x) dlmModPoly(1, dV = x[1], dW = x[2])
mod<-StructTS(response, "level") #level is W, epsilon is V
plot(response)
lines(fitted(mod), col="red")

randomwalk.fit <- dlmMLE(response, parm = c(1, sd(response)),build = randomwalk)
mod <- randomwalk(randomwalk.fit$par)
unlist(randomwalk(randomwalk.fit$par)[c("V", "W")])

bob<-dlmFilter(response,mod)
par(mfrow=c(2,2))
plot(dropFirst(bob$m), main="filtered state") #filtered values of state vector
plot(dropFirst(bob$f),main="one step ahead forecast") #one step ahead forecast
plot(dropFirst(bob$a), main="predicted") #predicted values given data up to and including previous time unit
points(bob$y, main="raw data", col="red") #data input


AICdlm(filteredobject="bob",MLEfitobject="randomwalk.fit",responsedata="response",burn_in=10)


#2 Just temperature

#temperature
temperature<-read.csv(file="temperature19482002.csv", stringsAsFactors=F)
temperature$DATE<-paste(substr(temperature$DATE,1,8),"01",sep="")
temperature<-temperature[which(temperature$DEPTH<151),]
library(doBy)
temperature<-summaryBy(TEMPERATURE~DATE,FUN=mean,data=temperature)
temperature$DATE<-(as.Date(temperature$DATE))
names(temperature)<-c("Date","Temp")
temperature<-NAdates(temperature)
plot(temperature$Temp~temperature$Date)
temppred<-ts(temperature$Temp,frequency=12, start=c(1955,1))

datas<-ts.intersect(response,temppred)
response<-datas[,1]
temppred<-datas[,2]

temp.function <- function(x) dlmModReg(temppred,dW=sd(temppred),dV=exp(x[2]),addInt=FALSE)
temp.fit <- dlmMLE(response,par=c(sd(temppred,na.rm=T), 1), build=temp.function)
mod.temp <- temp.function(temp.fit$par)

bob.temp<-dlmFilter(response,mod.temp) 
unlist(temp.function(temp.fit$par)[c("V", "W")])
par(mfrow=c(2,2))
plot(bob$m, main="filtered state") #filtered values of state vector
plot(bob$f,main="one step ahead forecast") #one step ahead forecast
plot(bob$a, main="predicted") #predicted values given data up to and including previous time unit
plot(bob$y, main="raw data") #data input

AICdlm(filteredobject="bob.temp",MLEfitobject="temp.fit",responsedata="response",burn_in=10)

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
