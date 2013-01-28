##############################################################
#############Density weighted average depth###################
######################## OR ##################################
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
majspp.winter<-list()
majspp.spring<-list()
majspp.summer<-list()
majspp.fall<-list()
Abundance.monthly<-list()
missingdata<-list() #percent of months with zero density or missing values e.g. 72 months of 720 = 10%

#species of interest
for (sppec in c("all", "Epischura")){

#group of interest ("all", "Copepod", "Cladoceran", "Rotifer")
taxa<-"Copepod"

#Type of analysis? ("DWA"=Density weighted depth, COD"=coefficient of dispersion)
analysis<-"DWA"

#start date of data extracted
startdate<-"1955-01-01"

#replace NAs with monthly averages?
replaceNA<-"yes"
  
#Quarter type ("DJF" or "JFM") quarters DJF as in Hampton et al. or standard quarters
quarter.type<-"JFM"

#Enter lifestage of interest("all", "adult", "naup", "copep") #This has not been tested for anything except for "all" or "adult" although many other options are in the key
stage<-"naup"

#remove doublecounts? ("yes","no")
dbl<-"yes"

######################################################
#################start data processing################
######################################################

allzoo<-zoodata

if(taxa!="all"){
  allzoo<-allzoo[which(allzoo$Group_General==taxa),] }#take records for group of interest

if(stage!="all"){
  allzoo<-allzoo[which(allzoo$Lifestage_Cop==stage),] }#take records for group of interest

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
custom<-function(x){sd(x, na.rm=T)/mean(x, na.rm=T)}
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

winter<-qtrnb[which(qtrnb$quarter==1),1]
spring<-qtrnb[which(qtrnb$quarter==2),1]
summer<-qtrnb[which(qtrnb$quarter==3),1]
fall<-qtrnb[which(qtrnb$quarter==4),1]

years<-unique(substr(DWAstep2$Date,1,4))

majspp.winter[[as.character(sppec)]]<-winter
majspp.spring[[as.character(sppec)]]<-spring
majspp.summer[[as.character(sppec)]]<-summer
majspp.fall[[as.character(sppec)]]<-fall
}

if (stage=="all"){all.abundance<-Abundance.monthly$Cyclops}
if (stage=="adult"){adult.abundance<-Abundance.monthly$Cyclops}
if (stage=="copep"){copep.abundance<-Abundance.monthly$Cyclops}
if (stage=="naup"){naup.abundance<-Abundance.monthly$Cyclops}


######Plots###########----------------------------------------------------
#plot by season
par(mfrow=c(2,2),mar=c(1,2,1,1),oma = c(2,3,1,0))

for (ssn in c("majspp.winter","majspp.spring","majspp.summer","majspp.fall")){
  plot(get(ssn)$all~years[1:(length(get(ssn)$all))],xaxt="n",cex.axis=1.5, main=substr(ssn,8,nchar(ssn)))
  if(ssn=="majspp.summer"|ssn=="majspp.fall"){axis(side=1, cex.axis=1.5)}
  
  library(dyn)
  dta<-data.frame()
  dta<-cbind(get(ssn)$all,as.numeric(years[1:length(get(ssn)$all)]))
  mod<-dyn$lm(dta[,1]~dta[,2],na.action=na.exclude)
  
  library(lmtest)
  auto<-dwtest(mod) #test for autocorrelation
  DW<-auto$p.value
  library(nlme)
  if (DW<0.05){mod<-gls(dta[,1]~dta[,2], correlation=corAR1(form=~1),na.action=na.exclude)}
  
  lines(fitted(mod)~years[1:(length(get(ssn)$all))], col="red", lwd=3)
  if(DW>0.05){b<-summary(mod);  p<-coefficients(b)[2,4]}
  if(DW<0.05){b<-summary(mod)$tTable; p<-b[2,4]}
  library(lmtest)
  
  legend("topleft",lty=1,lwd=3,cex=1,col="red",legend=paste("linear trend,","p=",round(p,4),",","DW=",round(DW,2)),bty="n")}

mtext("Year", side=1, line=1, outer=TRUE, cex=1.5)
mtext(analysis, side=2, line=1, outer=TRUE, cex=1.5)


####for major species------------------------------------------------------

majspp2<-names(majspp.winter)[2:length(majspp.winter)]

par(mfrow=c(2,2),mar=c(2,2,2,2),oma = c(2,3,1,2))

u=0

for (ssn in c("majspp.winter", "majspp.spring","majspp.summer","majspp.fall")){
  u=u+1
  fix(majspp2)
  spp<-majspp2
  clrs<-c("black","red","blue", "green")
  points<-c(15,16,17,18)
  dta<-list()
  mod<-list()
  p<-c()
  
  for (i in 1:length(spp)){
    
    specofint<-spp[i]
    
    kt<-grep(specofint,names(get(ssn)))
    if (length(kt)>1){dta[[i]]<-unlist(get(ssn)[kt])
                      mod[[i]]<-lm(unlist(get(ssn)[kt])~as.numeric(unique(years)[1:length(unlist(get(ssn)[kt]))]), na.action=na.exclude)}
    
    if (length(kt)==1){dta[[i]]<-unlist(get(ssn)[kt])
                       mod[[i]]<-lm(unlist(get(ssn)[kt])~as.numeric(unique(years)[1:length(unlist(get(ssn)[kt]))]), na.action=na.exclude)}
  }
  
  #plot
    allyears<-unique(years); allyears<-allyears[which(allyears!="NA")]
  plot(dta[[1]]~as.numeric(unique(years)[1:length(dta[[1]])]), xlab="Year", ylab=analysis, main=substr(ssn,8,nchar(ssn)), pch=points[1], col=clrs[1], cex.lab=1.5, cex.axis=1.5, ylim=c(0.95*min(unlist(dta),na.rm=T),1.2*(max(unlist(dta),na.rm=T))))
  lines(fitted(mod[[1]])~allyears[1:length(fitted(mod[[1]]))],col=clrs[1]) 
  b<-summary(mod[[1]])
  p[1]<-coefficients(b)[2,4]
  
  if (length(spp)>1){
    for (i in 2:length(spp)){
      points(dta[[i]]~as.numeric(unique(years)[1:length(dta[[i]])]),pch=points[i],col=clrs[i])
      lines(fitted(mod[[i]])~allyears[1:length(fitted(mod[[i]]))],col=clrs[i]) 
      b<-summary(mod[[i]])
      p[i]<-coefficients(b)[2,4]
    }}
  
  legend("topright",pch=points, col=clrs,legend=paste(spp,"p=",round(p,4)), cex=0.7)
}
mtext("Year", side=1, line=1, outer=TRUE, cex=1.5)
mtext(analysis, side=2, line=1, outer=TRUE, cex=1.5)

#paste("Cyclops",taxa,"lifestage=",stage,"Quarters=",quarter.type,sep=" ")


#Plots for abundance through time for epischura life stages. Must run code above for each life stage first

plot(all.abundance, ylim=c(0,350), ylab="Density (Number/m3)")
plot(naup.abundance, col="red", type="l", ylim=c(0,350), ylab="Density (/m3)")
lines(copep.abundance, col="blue")
lines(adult.abundance,col="green")
legend("topleft",legend=c("nauplii","copepodite","adult"), lty=1, col=c("red","blue","green"))


plot(adult.abundance[,2]+copep.abundance[,2]+naup.abundance[,2]~adult.abundance[,1])
points(all.abundance,col="red")



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