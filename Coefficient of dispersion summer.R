##############################################################
################coefficient of dispersion#####################
##############################################################

######## 1. Need to run code for custom functions at bottom of this file

######## 2. The main data file containing all observations
# Check if main file has been loaded, if not read data into R
if(exists("zoodata")==FALSE){zoodata<-read.csv(file="zoopzeros250key1.csv", stringsAsFactors=F)}       #"c:/Users/gray/Dropbox/Postdoc/R data/zoopzeros250key1.csv")

######## 3. The zooplankton key
#Check if the key has been loaded yet, if not read it in 
if(exists("key")==FALSE){key<-read.csv(file="key.csv")}                            #"c:/Users/gray/Dropbox/Postdoc/R data/key.csv")

#lists to be filled with DWA data for species of interest for each quarter
majspp.summer<-list()
majspp.density<-list()
Abundance.monthly<-list()
missingdata<-list() #percent of months with zero density or missing values e.g. 72 months of 720 = 10%

#group of interest ("all", "Copepod", "Cladoceran", "Rotifer")
taxa<-"all" #check stage... needs to be "all" if this is not copepod

if (taxa=="all"){spplist<-c("all")} #stage must be "all" for "all" taxa
if (taxa=="Copepod"){splist<-c("all","Epischura","Cyclops")}
if (taxa=="Cladoceran"){splist<-c("all","Bosmina")} #note stage must be "all" for cladocerans
if (taxa=="Rotifer"){splist<-c("all","Conochilus","Filinia","Kellicottia","Keratella")} #note stage must be "all" for cladocerans

#Enter lifestage of interest("all", "adult", "naup", "copep") #This has not been tested for anything except for "all" or "adult" although many other options are in the key
stage<-"all"

#species of interest
for (hg in 1: length(splist)){

sppec<-splist[hg]

#Type of analysis? ("DWA"=Density weighted depth, COD"=coefficient of dispersion)
analysis<-"DWA"

#start date of data extracted
startdate<-"1955-01-01"

#replace NAs with monthly averages?
replaceNA<-"no"
  
#Quarter type ("DJF" or "JFM") quarters DJF as in Hampton et al. or standard quarters
quarter.type<-"JFM"



#remove doublecounts? ("yes","no")
dbl<-"yes"

######################################################
#################start data processing################
######################################################

allzoo<-zoodata

if(taxa!="all"){
  allzoo<-allzoo[which(allzoo$Group_General==taxa),] }#take records for group of interest

if (taxa=="Copepod"){
if(stage!="all"){
  allzoo<-allzoo[which(allzoo$Lifestage_Cop==stage),] }#take records for stage of interest
}

#Take data only from species of interest (or take all data if an individual species wasn't specified)
if (sppec!="all"){
  if (length(strsplit(sppec,' ')[[1]])==1){
  allzoo<-allzoo[grep(sppec, allzoo$Genus),]}
  if (length(strsplit(sppec,' ')[[1]])==2){
  allzoo<-allzoo[grep(sppec, paste(allzoo$Genus,allzoo$Species)),]}
}

#Take after requested date
allzoo<-allzoo[which(allzoo$Date1>=as.Date(startdate)),]

#create a volume column and total number of individuals column for later calculations
allzoo$Volume<-(as.numeric(allzoo$Nig_Gr)-as.numeric(allzoo$Ver_Gr))*(pi*(0.375/2)^2)*1000 # 0.375 diameter net, multiply by 1000 to convert to litres. Note NAs will be introduced due to some records having no depth information
allzoo<-allzoo[which(allzoo$Volume!="NA"),] #removes records with no volume
allzoo$Tot<-allzoo$Volume*allzoo$Count #Add a column that shows total number of individuals (volume sampled multiplied by Lyubov's densities (#/m3))

#create single depth column by merging Ver_Gr and Nig_Gr
allzoo$Depth<-paste(allzoo$Ver_Gr,allzoo$Nig_Gr, sep="-")

#For each date multiply the depth of the sample by the density
allzoo$DWA<-((as.numeric(allzoo$Nig_Gr)+as.numeric(allzoo$Ver_Gr))/2)*allzoo$Count #create a column for DWA calculations

allzoocommondepths<-allzoo[which(allzoo$Depth=="0-10"|allzoo$Depth=="10-25"|allzoo$Depth=="50-100"|allzoo$Depth=="100-150"|allzoo$Depth=="150-250"),]

#remove doublecounts
if(dbl=="yes"){
KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")] #which codes are not doubles (according to Lizzie)
allzoocommondepths<-allzoocommondepths[which(allzoocommondepths$KOD%in% KODsnotdouble),]
}

#Check if all depths were sampled on a particular date
exc<-c()
da<-c()
cat<-list()
undates<-unique(allzoocommondepths$Date1)
for (i in 1:length(undates)){
  this<-allzoocommondepths[which(allzoocommondepths$Date1==undates[i]),]
  if(nrow(this)>0){
  cat[[i]]<-(unique(this$Depth))
  da[i]<-length(unique(this$Depth))
  if(da[i]==5){exc[i]<-"N"}; if(da[i]<5){exc[i]="Y"}}
  }
gooddates<-undates[which(exc=="N")]
allzoocommondepths<-allzoocommondepths[which(allzoocommondepths$Date1%in%gooddates),]

#Calculations for coefficient of dispersion
if(analysis=="COD"){
analysis<-"Coefficient of dispersion"
library(doBy)
DWAdispstep1<-summaryBy(DWA+Count~Date1+Depth, FUN=sum, order=TRUE,data=allzoocommondepths)
custom<-function(x){sd(x, na.rm=T)/mean(x, na.rm=T)} #sd over mean
DWAstep2<-summaryBy(DWA.sum~Date1,FUN=custom,data=DWAdispstep1)
names(DWAstep2)<-c("Date","DWA")
}

#Calculations for Density weighted depth
if(analysis=="DWA"){
library(doBy)
analysis<-"Density weighted depth"
sum.noNA<-function(x) sum(x,na.rm=T)
DWAstep1<-summaryBy(DWA+Count~Date1,data=allzoocommondepths,FUN=sum.noNA,order=TRUE) #sum count*depth and count by date
DWAstep2<-data.frame("DWA"=DWAstep1$DWA.sum.noNA/DWAstep1$Count.sum.noNA, "Date"=DWAstep1$Date1)
}

#########################
#summarize by month######
#########################

mean.noNA<-function(x) mean(x,na.rm=T)

#for density weighted depth
DWAstep2$Date<-paste(substr(DWAstep2$Date,1,7),"-01",sep="")
DWAmonthly<-summaryBy(DWA~Date, data=DWAstep2, FUN=mean.noNA, order=TRUE)
DWAmonthly$Date<-as.Date(DWAmonthly$Date)
missingdata[as.character(sppec)]<-(length(which(is.na(NAdates(DWAmonthly)[,2]==T)))/length(DWAmonthly[,2]))*100
if (replaceNA=="no"){DWAmonthly<-NAdates(DWAmonthly)}
if (replaceNA=="yes"){DWAmonthly<-NAfill(DWAmonthly)}
names(DWAmonthly)<-c("Date","DWA")
  
#for abundance
Density<-summaryBy(Count+Volume~Date1+Depth,data=allzoocommondepths,FUN=c(sum,mean),order=TRUE) 
Density$Tot<-Density$Count.sum*Density$Volume.mean
Density$Date1<-paste(substr(Density$Date1,1,7),"-01",sep="")
Density.monthly<-summaryBy(Tot+Volume.mean~Date1, data=Density, FUN=sum, order=TRUE)
Density.monthly$density<-Density.monthly$Tot.sum/Density.monthly$Volume.mean.sum
Density.monthly<-data.frame("Date"=Density.monthly$Date1,"Density"=Density.monthly$density)
Density.monthly$Date<-as.Date(Density.monthly$Date)
if (replaceNA=="no"){Density.monthly<-NAdates(Density.monthly)}
if (replaceNA=="yes"){Density.monthly<-NAfill(Density.monthly)}
Abundance.monthly[[sppec]]<-Density.monthly

allsp<-DWAmonthly

##########################
#summarize data by quarters
##########################

if (quarter.type=="DJF"){
  library(xts)
  months<-as.numeric(substr(allsp$Date,6,7))+1
  months[which(months==13)]<-1
  years<-as.numeric(substr(allsp$Date,1,4))
  years[which(months==1)]<-years[which(months==1)]+1
  months[which(nchar(months)==1)]<-paste("0",months[which(nchar(months)==1)],sep="")
  days<-rep("01",length(years))
  allsp$Date<-as.Date(paste(years,months,days,sep="-"))
  dates<-as.Date(allsp$Date)
}

#For DWA
library(xts)
qtnb<-xts(x=allsp$DWA, order.by=as.Date(allsp$Date))
locations<-endpoints(qtnb, "quarters")
qtrnb<-as.data.frame(period.apply(qtnb, INDEX=locations, FUN=function(x) mean(x,na.rm=T)))

quarterly<-data.frame(qtrnb)
quarters<-(as.numeric(substr(row.names(qtrnb),6,7)))
quarters[which(quarters==1|quarters==2|quarters==3)]<-1
quarters[which(quarters==4|quarters==5|quarters==6)]<-2
quarters[which(quarters==7|quarters==8|quarters==9)]<-3
quarters[which(quarters==10|quarters==11|quarters==12)]<-4
qtrnb$quarter<-quarters

summer<-qtrnb[which(qtrnb$quarter==3),1]


#For Density

qtnb2<-xts(x=Density.monthly$Density, order.by=as.Date(Density.monthly$Date))
locations<-endpoints(qtnb2, "quarters")
qtrnb2<-as.data.frame(period.apply(qtnb2, INDEX=locations, FUN=function(x) mean(x,na.rm=T)))

quarterly<-data.frame(qtrnb2)
quarters<-(as.numeric(substr(row.names(qtrnb2),6,7)))
quarters[which(quarters==1|quarters==2|quarters==3)]<-1
quarters[which(quarters==4|quarters==5|quarters==6)]<-2
quarters[which(quarters==7|quarters==8|quarters==9)]<-3
quarters[which(quarters==10|quarters==11|quarters==12)]<-4
qtrnb2$quarter<-quarters

Density<-qtrnb2[which(qtrnb2$quarter==3),1]



years<-seq(min(as.numeric(unique(substr(DWAstep2$Date,1,4)))), max(as.numeric(unique(substr(DWAstep2$Date,1,4)))), by=1)


majspp.summer[[as.character(sppec)]]<-summer
majspp.density[[as.character(sppec)]]<-Density

}

#this is for making a dataframe for the linear model
if(taxa=="Copepod"&stage=="copep"){Copepod.copep<-data.frame("DWA"=majspp.summer$all, "Year"=years[1:length(majspp.summer$all)], "Taxa"=rep("Copep", length(majspp.summer$all)))} #copepodites
if(taxa=="Copepod"&stage=="naup"){Copepod.naup<-data.frame("DWA"=majspp.summer$all, "Year"=years[1:length(majspp.summer$all)], "Taxa"=rep("Naup", length(majspp.summer$all)))} #nauplii
if(taxa=="Copepod"&stage=="adult"){Copepod.adult<-data.frame("DWA"=majspp.summer$all, "Year"=years[1:length(majspp.summer$all)], "Taxa"=rep("Adult", length(majspp.summer$all)))} #adult
if(taxa=="Rotifer"){Rotifer<-data.frame("DWA"=majspp.summer$all, "Year"=years[1:length(majspp.summer$all)], "Taxa"=rep("Rotifer", length(majspp.summer$all)))} #rotifers
if(taxa=="Cladoceran"){Cladoceran<-data.frame("DWA"=majspp.summer$all, "Year"=years[1:length(majspp.summer$all)], "Taxa"=rep("Cladoceran", length(majspp.summer$all)))} #cladoceran
if(taxa=="all"){all<-data.frame("DWA"=majspp.summer$all, "Year"=years[1:length(majspp.summer$all)], "Taxa"=rep("all", length(majspp.summer$all)))} #all

#allDWA<-rbind(Copepod.copep, Copepod.naup,Copepod.adult,Rotifer,Cladoceran,all)
allDWA$Year<-as.numeric(allDWA$Year)
allDWA$Taxa<-as.factor(allDWA$Taxa)



######Plots###########----------------------------------------------------

#Summer plot
library(ggplot2)
library(dyn)
dta<-data.frame()
dta<-cbind(majspp.summer$all,as.numeric(years[1:length(majspp.summer$all)]))
mod<-dyn$lm(dta[,1]~dta[,2],na.action=na.exclude)

library(lmtest)
auto<-dwtest(mod) #test for autocorrelation
DW<-auto$p.value
library(nlme)
if (DW<0.05){mod<-gls(dta[,1]~dta[,2], correlation=corAR1(form=~1),na.action=na.exclude)}

if(DW>0.05){b<-summary(mod);  p<-coefficients(b)[2,4]}
if(DW<0.05){b<-summary(mod)$tTable; p<-b[2,4]}

results<-data.frame("Year"=as.numeric(years)[1:length(fitted(mod))], "Fitted"=fitted(mod))

groupkey<-data.frame(inkey=c("all", "Copepod","Cladoceran","Rotifer"), real=c("All Zooplankton", "Copepod","Cladoceran","Rotifer"))
currentgroup<-groupkey$real[which(groupkey$inkey==taxa)]
datas<-na.omit(data.frame("Date"=years[1:length(majspp.summer$all)],"DWA"=majspp.summer$all,"Density"=majspp.density$all))

ggplot(datas, aes(x=Date, y=DWA, size=Density),legend=FALSE)+
  theme_bw()+
  geom_point(colour="black", shape=16)+ 
  opts(axis.text.x = theme_text(colour = "black", size = 20))+
  opts(axis.text.y = theme_text(colour = "black", size = 20))+
  opts(axis.title.x = theme_text(size = 20, colour = "black"))+
  opts(axis.title.y = theme_text(size = 20, colour = "black", angle = 90), hjust=-0.3)+ 
  scale_y_reverse(name=paste(currentgroup,"Coefficient of dipsersion"))+
  geom_line(aes(x = Year, y = Fitted),size=1,data=results, colour="red") +
  annotate("text", label = paste("linear trend,","p=",round(p,4),",","DW=",round(DW,2)), x = 1990, y = (1.02*min(datas$DWA)), size = 7, colour = "red")

paste(taxa,"lifestage=",stage,"Quarters=",quarter.type,sep=" ")

####for major species------------------------------------------------------

majspp2<-names(majspp.summer)[2:length(majspp.summer)]


for (jj in 1:length(majspp2)){
library(dyn)
dta<-data.frame()
dta<-cbind(majspp.summer[[jj+1]],as.numeric(years[1:length(majspp.summer[[jj+1]])]))

mod<-dyn$lm(dta[,1]~dta[,2],na.action=na.exclude)

library(lmtest)
auto<-dwtest(mod) #test for autocorrelation
DW<-auto$p.value
library(nlme)
if (DW<0.05){mod<-gls(dta[,1]~dta[,2], correlation=corAR1(form=~1),na.action=na.exclude)}

if(DW>0.05){b<-summary(mod);  p<-coefficients(b)[2,4]}
if(DW<0.05){b<-summary(mod)$tTable; p<-b[2,4]}

results<-data.frame("Year"=as.numeric(years)[1:length(fitted(mod))], "Fitted"=fitted(mod))

groupkey<-data.frame(inkey=c("all", "Copepod","Cladoceran","Rotifer"), real=c("All Zooplankton", "Copepod","Cladoceran","Rotifer"))
currentgroup<-groupkey$real[which(groupkey$inkey==taxa)]
datas<-na.omit(data.frame("Date"=years[1:length(majspp.summer[[jj+1]])],"DWA"=majspp.summer[[jj+1]],"Density"=majspp.density[[jj+1]]))

pl<-ggplot(datas, aes(x=Date, y=DWA, size=Density),legend=FALSE)+
  theme_bw()+
  geom_point(colour="black", shape=16)+ 
  opts(axis.text.x = theme_text(colour = "black", size = 12))+
  opts(axis.text.y = theme_text(colour = "black", size = 12))+
  opts(axis.title.x = theme_text(size = 12, colour = "black"))+
  opts(axis.title.y = theme_text(size = 12, colour = "black", angle = 90), hjust=-0.3)+ 
  scale_y_reverse(name=paste(majspp2[[jj]],"Coefficient of dipsersion"))+
  geom_line(aes(x = Year, y = Fitted),size=1,data=results, colour="red") +
  annotate("text", label = paste("linear trend,","p=",round(p,4),",","DW=",round(DW,2)), x = 1980, y =(0.98*max(datas$DWA)), size = 4, colour = "red")

assign(paste("pl",jj,sep=""),pl)
}

library(gridExtra)
if(taxa=="Rotifer"){grid.arrange(pl1, pl2, pl3,pl4,ncol=2)}
if(taxa=="Copepod"){grid.arrange(pl1, pl2, ncol=1)}
if(taxa=="Cladoceran"){pl1}

paste(taxa,"major species","lifestage=",stage,"Quarters=",quarter.type,sep=" ")




#Regressions------------------------------------------------------------------------------------

model2<-gls(DWA~Taxa/Year, correlation=corAR1(form=~1),na.action=na.exclude, data=allDWA) #temporal autocorrelation
model2<-gls(DWA~Taxa/Year, na.action=na.exclude, data=allDWA) #no temporal autocorrelation
model3<-gls(DWA~Year*Taxa, na.action=na.exclude, data=allDWA) #no temporal autocorrelation
anova(model1,model2)

model2<-lm(DWA~Taxa/Year, na.action=na.exclude, data=allDWA) #no temporal autocorrelation
summary(model2)

#decided to go with model 4 for simplicity
model4<-lm(DWA~Year*Taxa, na.action=na.exclude, data=allDWA) #no temporal autocorrelation
summary(model4)

### Test if trends are different from zero
results<-glht(model4, linfct = c("Year=0","Year:TaxaNaup = 0",
                                 "Year:TaxaAdult = 0",
                                 "Year:TaxaRotifer = 0",
                                 "Year:TaxaCladoceran = 0"))
summary(results)

#Model fit
par(mar=c(4,5,2,2))
pdata <- expand.grid(Year=seq(1955, 2002, by=1), Taxa=c("Copep","Naup","Adult","Rotifer","Cladoceran"))
pdata$Temp <- allDWA$DWA
plot(pdata$Year,pdata$Temp, type="n",xlab="Year", ylim=range(pdata$Temp, na.rm=T),ylab=analysis, cex.lab=1.5,cex.axis=1.5)
points(pdata$Year[which(pdata$Taxa=="Copep")], pdata$Temp[which(pdata$Taxa=="Copep")], type="p", pch=15)
points(pdata$Year[which(pdata$Taxa=="Naup")], pdata$Temp[which(pdata$Taxa=="Naup")], type="p", pch=16, col="red")
points(pdata$Year[which(pdata$Taxa=="Adult")], pdata$Temp[which(pdata$Taxa=="Adult")], type="p", pch=17, col="green")
points(pdata$Year[which(pdata$Taxa=="Rotifer")], pdata$Temp[which(pdata$Taxa=="Rotifer")], type="p", pch=18, col="blue")
points(pdata$Year[which(pdata$Taxa=="Cladoceran")], pdata$Temp[which(pdata$Taxa=="Cladoceran")], type="p", pch=25, col="orange")
lines(fitted(model4)[which(pdata$Taxa=="Copep")]~pdata$Year[which(pdata$Taxa=="Copep")])
lines(fitted(model4)[which(pdata$Taxa=="Naup")]~pdata$Year[which(pdata$Taxa=="Naup")], col="red")
lines(fitted(model4)[which(pdata$Taxa=="Adult")]~pdata$Year[which(pdata$Taxa=="Adult")], col="green")
lines(fitted(model4)[which(pdata$Taxa=="Rotifer")]~pdata$Year[which(pdata$Taxa=="Rotifer")], col="blue")
lines(fitted(model4)[which(pdata$Taxa=="Cladoceran")]~pdata$Year[which(pdata$Taxa=="Cladoceran")], col="orange")
#if (class(model4)=="lm"){temp<-round(summary(model2)$coefficients[6:10,4],4)} #for lm model
#if (class(model4)=="gls"){temp<-round(summary(model2)$tTable[6:10,4],4)} #for gls model

temp<-round(summary(results)[9]$test$pvalues[1:5],4)
temp[which(temp<0.0001)]<-"<0.0001"
legend("topright", c(paste("Copepodites, p=",temp[1]),paste("Nauplii, p=",temp[2]),paste("Adult Copepods, p=",temp[3]),paste("Rotifers, p=",temp[4]),paste("Cladocerans, p=",temp[5])), pch=c(15, 16, 17,18,25),col=c("black","red","green","blue","orange"), lty=1, lwd=2, cex=0.8)

#random intercept and slope
model<-lme(DWA~Year*Taxa,data=allDWA,method="ML",na.action=na.exclude)

library(nlme)
model3<-lme(DWA~Year*Taxa,random=~1|Taxa,data=allDWA,method="ML",na.action=na.exclude)
summary(model3)$tTable

### set up general linear hypothesis
results1<-glht(model3, linfct = K)
results<-glht(model3, linfct = c("Year=0","Year:TaxaNaup = 0",
                                 "Year:TaxaAdult = 0",
                                 "Year:TaxaRotifer = 0",
                                 "Year:TaxaCladoceran = 0"))

#Predictions from hierarchical model
par(mar=c(4,5,2,2))
pdata <- expand.grid(Year=seq(1955, 2002, by=1), Taxa=c("Copep","Naup","Adult","Rotifer","Cladoceran"))
pdata$Temp <- allDWA$DWA
plot(pdata$Year,pdata$Temp, type="n",xlab="Year", ylim=range(pdata$Temp, na.rm=T),ylab=analysis, cex.lab=1.5,cex.axis=1.5)
points(pdata$Year[which(pdata$Taxa=="Copep")], pdata$Temp[which(pdata$Taxa=="Copep")], type="p", pch=15)
points(pdata$Year[which(pdata$Taxa=="Naup")], pdata$Temp[which(pdata$Taxa=="Naup")], type="p", pch=16, col="red")
points(pdata$Year[which(pdata$Taxa=="Adult")], pdata$Temp[which(pdata$Taxa=="Adult")], type="p", pch=17, col="green")
points(pdata$Year[which(pdata$Taxa=="Rotifer")], pdata$Temp[which(pdata$Taxa=="Rotifer")], type="p", pch=18, col="blue")
points(pdata$Year[which(pdata$Taxa=="Cladoceran")], pdata$Temp[which(pdata$Taxa=="Cladoceran")], type="p", pch=25, col="orange")
lines(fitted(model3)[which(pdata$Taxa=="Copep")]~pdata$Year[which(pdata$Taxa=="Copep")])
lines(fitted(model3)[which(pdata$Taxa=="Naup")]~pdata$Year[which(pdata$Taxa=="Naup")], col="red")
lines(fitted(model3)[which(pdata$Taxa=="Adult")]~pdata$Year[which(pdata$Taxa=="Adult")], col="green")
lines(fitted(model3)[which(pdata$Taxa=="Rotifer")]~pdata$Year[which(pdata$Taxa=="Rotifer")], col="blue")
lines(fitted(model3)[which(pdata$Taxa=="Cladoceran")]~pdata$Year[which(pdata$Taxa=="Cladoceran")], col="orange")
if (class(model2)=="lm"){temp<-round(summary(model3)$coefficients[6:10,4],4)} #for lm model
if (class(model2)=="gls"){temp<-round(summary(model3)$tTable[6:10,4],4)} #for gls model
temp[which(temp<0.0001)]<-"<0.0001"
legend("topright", c(paste("Copepodites, p=",temp[1]),paste("Nauplii, p=",temp[2]),paste("Adult Copepods, p=",temp[3]),paste("Rotifers, p=",temp[4]),paste("Cladocerans, p=",temp[5])), pch=c(15, 16, 17,18,25),col=c("black","red","green","blue","orange"), lty=1, lwd=2, cex=0.8)


#Individual regressions
model1<-dyn$lm(Copepod.copep$DWA~Copepod.copep$Year)
library(lmtest)
dwtest(model1)
summary(model1)

model2<-dyn$lm(Copepod.naup$DWA~Copepod.naup$Year)
dwtest(model2)
summary(model2)

model3<-dyn$lm(Copepod.adult$DWA~Copepod.adult$Year)
dwtest(model3)
summary(model3)

model4<-dyn$lm(Rotifer$DWA~Rotifer$Year)
dwtest(model4)
summary(model4)

model5<-dyn$lm(Cladoceran$DWA~Cladoceran$Year)
dwtest(model5)
summary(model5)

#correct p-values for number of tests
ps<-c(1.82*10^-6, 0.000704, 0.866, 1.24*10^-7, 0.00926)
p.adjust(ps, method = "bonferroni", n = length(ps))

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