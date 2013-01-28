##Depth overlap analyses

#Need 4 files for this:
#1. Lizzie's zooplankton file
#2. Latest phytoplankton file (one checked by Stephanie)
#3. File called "edible.csv" that categorizes the common phytoplankton genera
#4. Zooplankton key

#Need to run the custom functions at the end of this file first

#Need the following packages:
#1. doBy, xts
require(doBy); require(xts); require(dyn)

#phyto group ("all", "diatom","cyano","green","euglenophytes","chrysophtye","cryptophyte","dinoflagellate","Bacteria","Flagellates","PicoAlgae unID","Flagellate","chrysoCyst")
phytgrp<-"all"

#Edible phytoplankton? ("all", "Yes","No") yes=only edible; no=only inedible; all= all species regardless of edibility
edi<-"all"

#zooplankton group ("Copepod", "Cladoceran", "Rotifer")
zoogrp<-"Copepod"

#zooplankton species? ("all" for both if you just want all copepods, OR zoogen="Epischura" zoosp="baicalensis")
zoogen<-"all"
zoosp<-"all"

#remove double counts? ("yes", "no")
dbl<-"yes"

#fill NA values for overlap with monthly means from time series? ("yes","no")
fillNA<-"no"

#Enter lifestage of interest("all", "adult", "naup", "copep") #This has not been tested for anything except for "all" or "adult" although many other options are in the key
stage<-"copep"

#read in zooplankton and phytoplankton files, make sure date column is in correct format, make better column names for phyto
if(exists("phyto")==FALSE){phyto<-read.csv(file="phytoplankton.csv", stringsAsFactors=F)
                           phyto$DATE<-as.Date(phyto$DATE,"%m/%d/%Y")
                           names(phyto)<-c("Date", "Month","Year","Depth","Code","Count","Group","Genus","Species")
                           phyto<-orderBy(Date~.,data=phyto)
                           phyto$Group[which(phyto$Group=="chrysophtye")]<-"chrysophyte"
                           phyto$Group[which(phyto$Group=="Bacteria")]<-"cyano"
                           phyto$Group[which(phyto$Group=="Flagellates")]<-"Flagellate"
                           phyto$Group[which(phyto$PhytoGenus=="Romeria")]<-"PicoAlgae unID"
                           phyto$Group[which(phyto$PhytoGenus=="Synechococcus")]<-"PicoAlgae unID"
                           }

if(exists("zooplankton")==FALSE){zooplankton<-read.csv(file="zoopzeros250key1.csv", stringsAsFactors=F)
zooplankton$Date1<-as.Date(zooplankton$Date1)} 


if(exists("edible")==FALSE){edible<-read.csv(file="edible.csv",stringsAsFactors=F)}


#Start date
stdate<-"1975-01-01" #most consistent for depths and after preservation change

####PREY
#take select depths and common genera and check start date
phyto2<-phyto[which(phyto$Depth==0|phyto$Depth==5|phyto$Depth==10|phyto$Depth==50|phyto$Depth==100|phyto$Depth==150|phyto$Depth==200|phyto$Depth==250),]
phyto2<-phyto2[which(phyto2$Date>=stdate),]

if (phytgrp!="all"){phyto2<-phyto2[which(phyto2$Group==phytgrp),]} #which phyto group was requested?
if (phytgrp=="all"){phyto2<-phyto2[which(phyto2$Group=="diatom"|phyto2$Group=="cyano"|phyto2$Group=="green"|phyto2$Group=="chrysophtye"|phyto2$Group=="cryptophyte"|phyto2$Group=="dinoflagellate"|phyto2$Group=="PicoAlgae unID"),]} #which phyto group was requested?


if (edi!="all"){
  phyed<-edible$Genus[which(edible$Edible==edi)]
  phyto2<-phyto2[which(phyto2$Genus%in%phyed),]} #seperate out edible species?


#summarize data by density according to date and depth
library(doBy)
smfun<-function(x){mean(x,na.rm=T)}
phyto3<-summaryBy(Count~Date+Depth, data=phyto2, FUN=smfun, order=TRUE)

#classify depths into categories by depth range
phyto3$Depth[which(phyto3$Depth==0|phyto3$Depth==5|phyto3$Depth==10)]<-"0-10"
phyto3$Depth[which(phyto3$Depth==50|phyto3$Depth==100)]<-"50-100"
phyto3$Depth[which(phyto3$Depth==150|phyto3$Depth==200|phyto3$Depth==250)]<-"150-250"

phyto4<-summaryBy(Count.smfun~Date+Depth, data=phyto3, FUN=smfun, order=TRUE)
names(phyto4)<-c("Date","Depth","Phytcount")

#####PREDATORS

#create single depth column by merging Ver_Gr and Nig_Gr
zooplankton$Depth<-paste(zooplankton$Ver_Gr,zooplankton$Nig_Gr, sep="-")
zoo2<-zooplankton[which(zooplankton$Depth=="0-10"|zooplankton$Depth=="50-100"|zooplankton$Depth=="150-250"),]
zoo2<-zoo2[which(zoo2$Date1>=stdate),]
if (zoogrp!="all"){zoo2<-zoo2[which(zoo2$Group_General==zoogrp),]} #which group was requested?
if (zoogen!="all"){zoo2<-zoo2[which(zoo2$Genus==zoogen),]} #which genus was requested?
if (zoosp!="all"){zoo2<-zoo2[which(zoo2$Species==zoosp),]} #which species was requested?

#which stage was requested (for copepods)
if(stage!="all"){
  zoo2<-zoo2[which(zoo2$Lifestage_Cop==stage),] }#take records for stage of interest

#remove doublecounts?
if(dbl=="yes"){
  #Check if the key has been loaded yet, if not read it in 
  if(exists("key")==FALSE){key<-read.csv(file="key.csv")}                            #"c:/Users/gray/Dropbox/Postdoc/R data/key.csv")
  KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")] #which codes are not doubles (according to Lizzie)
  zoo2<-zoo2[which(zoo2$KOD%in% KODsnotdouble),]
}
smfun<-function(x){mean(x,na.rm=T)}
zoo3<-summaryBy(Count~Date1+Depth, data=zoo2, FUN=smfun, order=TRUE)
names(zoo3)<-c("Date", "Depth","Zoocount")

#Merge phytoplankton and zooplankton into one dataframe by date
merged<-merge(zoo3,phyto4,all.x=T)

#bin into months
merged$Date<-paste(substr(merged$Date,1,7),"-01",sep="")
merged<-summaryBy(Zoocount+Phytcount~Date+Depth, data=merged,FUN=smfun,order=T)
merged$Date<-as.Date(merged$Date)
names(merged)<-c("Date","Depth","Zoocount","Phytcount")

#calculate overlap  
merged$calc1<-merged$Zoocount*merged$Phytcount #step1
calc2fun<-function(x){c("ct"=sum(x,na.rm=T)*(length(na.omit(x))), "ln"=length(na.omit(x)))}

step2<-summaryBy(calc1~Date,data=merged,FUN=calc2fun)
step2<-step2[which(step2$calc1.ln>=2),1:2] #take only dates when all 3 depths were available

sumfun<-function(x){sum(x,na.rm=T)}

step3<-summaryBy(Zoocount~Date,data=merged,FUN=sumfun)
step4<-summaryBy(Phytcount~Date,data=merged,FUN=sumfun)
#step4$Date[which(step4[,2]==0)]

#merge portions of the calculation into one dataframe
allsteps<-merge(step2,step3)
allsteps<-merge(allsteps,step4)
names(allsteps)<-c("Date","step2","step3","step4")

#do final calculations
overlap<-allsteps$step2/(allsteps$step3*allsteps$step4)
overlap2<-xts(overlap, order.by=allsteps$Date)
locations<-endpoints(overlap2,on="quarters")
overlap3<-period.apply(overlap2,INDEX=locations,FUN=function(x){mean(x,na.rm=T)})
index(overlap3)<-as.Date(paste(substr(index(overlap3),1,7),"01",sep="-"))

years<-substr(index(overlap3),1,4)
start = as.Date(index(overlap3[which.min(index(overlap3))]))
full <- seq(start, by="3 months", length=length(unique(years))*12)
dta<-as.vector(overlap3)
overlap3<-xts(dta[match(full, as.Date(index(overlap3)))],order.by=full)
overlap3<-na.trim(overlap3)

if(fillNA=="yes"){
overlap4<-na.aggregate(overlap3, quarters)
#overlap4<-na.trim(na.approx(overlap3, x=as.Date(full)))
}
if(fillNA=="no"){
  overlap4<-overlap3}

dates<-as.Date(index(overlap4))
overlap4<-as.vector(overlap4)


#Plot over whole time span
par(mfrow=c(1,1))
plot(overlap4~as.Date(dates), ylab="Overlap", xlab="Date", cex.axis=1.5, cex=1.5, cex.lab=1.5,pch=16, main=paste("Overlap:","Phytoplankton=",phytgrp,"/","Zooplankton=",zoogrp, "/ Edible:",edi, "/ Species:",zoogen,zoosp))
mod<-lm(overlap4~dates,na.action="na.exclude")
b<-summary(mod);  p<-coefficients(b)[2,4]
lines(fitted(mod)~dates)
legend("topright",lty=1,legend=paste("linear trend,","p=", round(p,4)),bty="n", cex=1.5)


#Summer plot
Summer=c(7,8,9)

months<-substr(dates,6,7)
  subs<-overlap4[which(as.numeric(months)%in%Summer)]
  subdate<-as.Date(dates[which(as.numeric(months)%in%Summer)])
datas<-na.omit(data.frame(subs,subdate))
  plot(subs~subdate, main=paste("Summer Overlap:","Phyto group=",phytgrp,"/","Zooplankton=",zoogrp,"/","Stage=",stage),xlab="Date",ylab="Overlap",xlim=c(as.Date("1975/1/1"),as.Date("2000/1/1")),cex.main=0.75)
  mod<-dyn$lm(subs~subdate,na.action="na.exclude")
  library(lmtest)
  auto<-dwtest(mod) #test for autocorrelation
  DW<-auto$p.value
  library(nlme)
  if (DW<0.05){mod<-gls(subs~subdate, correlation=corAR1(form=~1),na.action=na.exclude)}
  if(DW>0.05){b<-summary(mod);  p<-coefficients(b)[2,4]}
  if(DW<0.05){b<-summary(mod)$tTable; p<-b[2,4]}
  lines(fitted(mod)~subdate,col="red",lwd=3)
  legend("topright",lty=1,legend=paste("linear trend,","p=", round(p,4)),bty="n", cex=1)
  
#this is for making a dataframe for the linear model
if(zoogrp=="Copepod"&stage=="copep"){Copepod.copep<-data.frame("Overlap"=subs, "Year"=subdate, "Taxa"=rep("Copep", length(subs)))} #copepodites
if(zoogrp=="Copepod"&stage=="naup"){Copepod.naup<-data.frame("Overlap"=subs, "Year"=subdate, "Taxa"=rep("Naup", length(subs)))} #nauplii
if(zoogrp=="Copepod"&stage=="adult"){Copepod.adult<-data.frame("Overlap"=subs, "Year"=subdate, "Taxa"=rep("Adult", length(subs)))} #adult
if(zoogrp=="Rotifer"){Rotifer<-data.frame("Overlap"=subs, "Year"=subdate, "Taxa"=rep("Rotifer", length(subs)))} #rotifers
if(zoogrp=="Cladoceran"){Cladoceran<-data.frame("Overlap"=subs, "Year"=subdate, "Taxa"=rep("Cladoceran", length(subs)))} #cladoceran
#allDWA<-rbind(Copepod.copep, Copepod.naup,Copepod.adult,Rotifer,Cladoceran)
allDWA$Year<-as.numeric(substr(allDWA$Year,1,4))
allDWA$Taxa<-as.factor(allDWA$Taxa)

#linear model with taxa nested within year
model1<-lm(Overlap~Taxa/Year, na.action=na.exclude, data=allDWA) #no temporal autocorrelation
library(nlme)
model2<-gls(Overlap~Taxa/Year,na.action=na.exclude, data=allDWA)

model3<-gls(Overlap~Taxa/Year, correlation=corAR1(form=~1),na.action=na.exclude, data=allDWA) #temporal autocorrelation
anova(model2,model3) #which is better? temporal autocorrelation or no?

summary(model2)


### Test if trends are different from zero
library(multcomp)
results<-glht(model2, linfct = c("TaxaAdult:Year=0","TaxaNaup:Year = 0","TaxaRotifer:Year = 0","TaxaCladoceran:Year = 0"))
summary(results)

#Model fit
par(mar=c(4,5,2,2))
pdata <- expand.grid(Year=seq(1975, 1999, by=1), Taxa=c("Copep","Naup","Adult","Rotifer","Cladoceran"))
pdata$Temp <- allDWA$Overlap
plot(pdata$Year,pdata$Temp, type="n",xlab="Year", ylab="Overlap", cex.lab=1.5,cex.axis=1.5)
points(pdata$Year[1:25], pdata$Temp[1:25], type="p", pch=15)
points(pdata$Year[26:50], pdata$Temp[26:50], type="p", pch=16, col="red")
points(pdata$Year[51:75], pdata$Temp[51:75], type="p", pch=17, col="green")
points(pdata$Year[76:100], pdata$Temp[76:100], type="p", pch=18, col="blue")
points(pdata$Year[101:125], pdata$Temp[101:125], type="p", pch=25, col="orange")
lines(fitted(model2)[1:25]~pdata$Year[1:25])
lines(fitted(model2)[26:50]~pdata$Year[26:50], col="red")
lines(fitted(model2)[51:75]~pdata$Year[51:75], col="green")
lines(fitted(model2)[76:100]~pdata$Year[76:100], col="blue")
lines(fitted(model2)[101:125]~pdata$Year[101:125], col="orange")

legend("topleft", c("Copepodites","Nauplii","Adult Copepods","Rotifers","Cladocerans"), pch=c(15, 16, 17,18,25),col=c("black","red","green","blue","orange"), lty=1, lwd=2, cex=0.8)

                                 
#Testing if overlap index works

#create dataframe of randomly distributed zoo and phyto
randzoo<-rnorm(600,mean=5,sd=1)
randphyt<-rnorm(600,mean=10,sd=1)
dt<-rep(seq(as.Date("1975/1/1"), by="2 weeks", length.out=200),times=3)
depth<-rep(c("0-10","50-150","150-250"), times=200)
test<-data.frame("Zoocount"=randzoo,"Phytcount"=randphyt,"Date"=dt, "Depth"=depth)
cor(test$Zoocount,test$Phytcount) #correlation coefficient

#create a dataframe of zoo and phyto highly correlated by depth
randzoo<-rnorm(600,mean=c(10,1,2.5),sd=c(1,0.2,0.6))
randphyt<-rnorm(600,mean=c(50,5,25),sd=c(10,1,3)) #change sd to get different levels of correlation
dt<-rep(seq(as.Date("1975/1/1"), by="2 weeks", length.out=200),times=3)
depth<-rep(c("0-10","50-150","150-250"), times=200)
test<-data.frame("Zoocount"=randzoo,"Phytcount"=randphyt,"Date"=dt,"Depth"=depth)

#create a dataframe of zoo and phyto with signficant depth differences
randzoo<-rnorm(600,mean=c(5,50,25),sd=c(1,10,3))
randphyt<-rnorm(600,mean=c(50,5,25),sd=c(10,1,3)) #change sd to get different levels of correlation
dt<-rep(seq(as.Date("1975/1/1"), by="2 weeks", length.out=200),times=3)
depth<-rep(c("0-10","50-150","150-250"), times=200)
test<-data.frame("Zoocount"=randzoo,"Phytcount"=randphyt,"Date"=dt,"Depth"=depth)

#order by date for easy manual calcs
test<-orderBy(Date~.,data=test)

#calculate overlap  
test$calc1<-test$Zoocount*test$Phytcount #step1
calc2fun<-function(x){sum(x,na.rm=T)*(length(x))}
step2<-summaryBy(calc1~Date,data=test,FUN=calc2fun)
sumfun<-function(x){sum(x,na.rm=T)}
step3<-summaryBy(Zoocount~Date,data=test,FUN=sumfun)
step4<-summaryBy(Phytcount~Date,data=test,FUN=sumfun)

#merge portions of the calculation into one dataframe
allsteps<-merge(step2,step3)
allsteps<-merge(allsteps,step4)
names(allsteps)<-c("Date","step2","step3","step4")

#do final calculations
overlap2<-allsteps$step2/(allsteps$step3*allsteps$step4)
plot(overlap2~allsteps$Date, ylab="Overlap", xlab="Date")


##########Custom functions needed for code above------------------------------

##Identify missing dates in a monthly time series and insert an "NA"
##Input is a dataframe with a "Date" column and one other column of observations
NAdates<-function(x){
  years<-substr(x$Date,1,4)
  start = as.Date(x$Date[which.min(x$Date)])
  full <- seq(start, by="1 month", length=length(unique(years))*12)
  allsp<-data.frame("Date"=full)
  for (i in 1:ncol(x)){
    t1<-data.frame("Count"=with(x, x[,i][match(full, as.Date(x$Date))]))
    allsp<-cbind(allsp, t1)}
  allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(x); allsp$Date<-full
  #names(allsp)<-c("Date","strat")
  return(allsp)
}

#Function to identify NAs and then fill with monthly average for whole time series
#Requires a dataframe with a "Date" column and one other column of observations
NAfill<-function(x){
  years<-substr(x$Date,1,4)
  start = as.Date(x$Date[which.min(x$Date)])
  full <- seq(start, by="1 month", length=length(unique(years))*12)
  allsp<-data.frame("Date"=full)
  for (i in 1:ncol(x)){
    t1<-data.frame("Count"=with(x, x[,i][match(full, as.Date(x$Date))]))
    allsp<-cbind(allsp, t1)}
  allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(x); allsp$Date<-full
  #names(allsp)<-c("Date","strat")
  missing<-allsp$Date[which(is.na(allsp[,2])==T)] 
  for (i in 1:length(missing)){
    allsp[,2][which(allsp$Date==as.Date(missing[i]))]<-mean(allsp[,2][which(as.numeric(substr(allsp$Date,6,7))==as.numeric(substr(missing[i],6,7)))],na.rm=T)
  }
  return(allsp)
}