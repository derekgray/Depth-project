#calculate DWA for phytoplankton by group

for (depth in c(0,10,50,100,200)){
#for (phytgrp in c("diatom","cyano","green","chrysophyte","cryptophyte","dinoflagellate","PicoAlgae unID")){
library(doBy)

missing.dates<-c() #when fill NA is "no" this gives the missing dates
percent.missing<-c() #when fill NA is "no" this gives the proportion of missing months

#phyto group ("all", "diatom","cyano","green","euglenophytes","chrysophyte","cryptophyte","dinoflagellate","Bacteria", "PicoAlgae unID","Flagellate","chrysoCyst")
phytgrp<-"diatom"

#Type of quarter ("DJF" or "JFM")
quarter.type<-"JFM"

#should NA dates be filled?
fillNAs<-"no" 

#Start date
stdate<-as.Date("1974-01-01") #most consistent for depths and after preservation change

#End date
enddate<-as.Date("1999-12-31")

#Need to read in phytoplankton file
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

#take select depths that have been sampled regularly throughout the program
#phyto2<-phyto[which(phyto$Depth==0|phyto$Depth==10|phyto$Depth==50|phyto$Depth==100|phyto$Depth==200),] #5, 250? |phyto$Depth==150
phyto2<-phyto[which(phyto$Depth==depth),] #250?

#Take data after the specified start date
phyto2<-phyto2[which(phyto2$Date>=stdate&phyto2$Date<=enddate),]

#Some functions needed for processing----------------------

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
  missing<-allsp$Date[which(is.na(allsp[,2])==T)] 
  for (i in 1:length(missing)){
    allsp[,2][which(allsp$Date==as.Date(missing[i]))]<-mean(allsp[,2][which(as.numeric(substr(allsp$Date,6,7))==as.numeric(substr(missing[i],6,7)))],na.rm=T)
  }
  return(allsp)
}

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

#---------------------------------------------------------

#Check which phyto group was requested
if (phytgrp!="all"){phyto2<-phyto2[which(phyto2$Group==phytgrp),]} #which phyto group was requested?
if (phytgrp=="all"){phyto2<-phyto2[which(phyto2$Group=="diatom"|phyto2$Group=="green"|phyto2$Group=="chrysophtye"|phyto2$Group=="cryptophyte"|phyto2$Group=="dinoflagellate"|phyto2$Group=="PicoAlgae unID"|phyto2$Group=="cyano"),]} #which phyto group was requested? 

#summarize data by density according to date and depth
phyto2$Date<-paste(substr(phyto2$Date,1,7),"-01",sep="") #make all dates 1st of month for easy monthly averages using summaryBy

#check if all depths were sampled on a particular date and exclude dates that have missing depths
exc<-c()
da<-c()
cat<-list()
undates<-unique(phyto2$Date)
for (i in 1:length(undates)){
  this<-phyto2[which(phyto2$Date==undates[i]),]
  if(nrow(this)>0){
    cat[[i]]<-(unique(this$Depth))
    da[i]<-length(unique(this$Depth))
    if(da[i]==length(unique(phyto2$Depth))){exc[i]<-"N"}; if(da[i]<length(unique(phyto2$Depth))){exc[i]="Y"}}
}
gooddates<-undates[which(exc=="N")]
phyto2<-phyto2[which(phyto2$Date%in%gooddates),]

#Do the calculations for DWA monthly
library(doBy)
smfun<-function(x){sum(x,na.rm=T)}
phyto3<-summaryBy(Count~Date+Depth, data=phyto2, FUN=smfun, order=TRUE)
phyto3$Depth[which(phyto3$Depth==0)]<-1
phyto3$DWA<-phyto3$Depth*phyto3$Count.smfun
phyto4<-summaryBy(DWA+Count.smfun~Date, FUN=smfun,data=phyto3)
phyto4$DWA<-phyto4$DWA/phyto4$Count.smfun.smfun
phyto4$Date<-as.Date(phyto4$Date)
phyto5<-data.frame("Date"=phyto4$Date, "DWA"=phyto4$DWA,stringsAsFactors=F)

#Figure out densities throught ime
mnfun<-function(x){mean(x,na.rm=T)}
phytdens<-summaryBy(Count~Date,data=phyto2, FUN=mnfun)
phytdens$Date<-as.Date(phytdens$Date)
names(phytdens)<-c("Date","Density(#/m3)")

#identify missing data
temp<-NAdates(data.frame("Date"=phyto4$Date, "Count"=phyto4$Count.smfun.smfun))
missing.dates<-temp$Date[which(is.na(temp$Count)==T)]
percent.missing<-length(temp$Date[which(is.na(temp$Count)==T)])/nrow(temp)
missing.months<-substr(missing.dates,6,7)

if (quarter.type=="JFM"){
missing.winter<-length(which(missing.months=="01"|missing.months=="02"|missing.months=="03"))/length(which(substr(phyto5$Date,6,7)=="01"|substr(phyto5$Date,6,7)=="02"|substr(phyto5$Date,6,7)=="03"))
missing.spring<-length(which(missing.months=="04"|missing.months=="05"|missing.months=="06"))/length(which(substr(phyto5$Date,6,7)=="04"|substr(phyto5$Date,6,7)=="05"|substr(phyto5$Date,6,7)=="06"))
missing.summer<-length(which(missing.months=="07"|missing.months=="08"|missing.months=="09"))/length(which(substr(phyto5$Date,6,7)=="07"|substr(phyto5$Date,6,7)=="08"|substr(phyto5$Date,6,7)=="09"))
missing.fall<-length(which(missing.months=="10"|missing.months=="11"|missing.months=="12"))/length(which(substr(phyto5$Date,6,7)=="10"|substr(phyto5$Date,6,7)=="11"|substr(phyto5$Date,6,7)=="12"))
}
if (quarter.type=="DJF"){
  missing.winter<-length(which(missing.months=="01"|missing.months=="02"|missing.months=="12"))/length(which(substr(phyto5$Date,6,7)=="01"|substr(phyto5$Date,6,7)=="02"|substr(phyto5$Date,6,7)=="12"))
  missing.spring<-length(which(missing.months=="04"|missing.months=="05"|missing.months=="03"))/length(which(substr(phyto5$Date,6,7)=="04"|substr(phyto5$Date,6,7)=="05"|substr(phyto5$Date,6,7)=="03"))
  missing.summer<-length(which(missing.months=="07"|missing.months=="08"|missing.months=="06"))/length(which(substr(phyto5$Date,6,7)=="07"|substr(phyto5$Date,6,7)=="08"|substr(phyto5$Date,6,7)=="06"))
  missing.fall<-length(which(missing.months=="10"|missing.months=="11"|missing.months=="09"))/length(which(substr(phyto5$Date,6,7)=="10"|substr(phyto5$Date,6,7)=="11"|substr(phyto5$Date,6,7)=="09"))
}

#Fill NAs?
if (fillNAs=="yes"){
  phyto5<-NAfill(phyto5)
  
  #NAs with zero abundances should be NA, while sample dates missed should be filled with averages
  phyto5$DWA[which(temp$Count==0)]<-NA
  phytdens<-NAfill(phytdens)}

if (fillNAs=="no"){
phyto5<-NAdates(phyto5)
phytdens<-NAdates(phytdens)
}


#Change to quarterly time series

if (quarter.type=="DJF"){
  library(xts)
  months<-as.numeric(substr(phyto5$Date,6,7))+1
  months[which(months==13)]<-1
  years<-as.numeric(substr(phyto5$Date,1,4))
  years[which(months==1)]<-years[which(months==1)]+1
  months[which(nchar(months)==1)]<-paste("0",months[which(nchar(months)==1)],sep="")
  days<-rep("01",length(years))
  phyto5$Date<-as.Date(paste(years,months,days,sep="-"))
  dates<-as.Date(phyto5$Date)
  phytdens$Date<-phyto5$Date
}

#first for DWA 
library(xts)
qtnb<-xts(x=phyto5$DWA, order.by=as.Date(phyto5$Date))
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

years<-as.numeric(unique(substr(phyto5$Date,1,4)))

#second for abundance (needed for bubble sizes)
qtnb2<-xts(x=phytdens[,2], order.by=as.Date(phytdens$Date))
locations2<-endpoints(qtnb2, "quarters")
qtrnb2<-as.data.frame(period.apply(qtnb2, INDEX=locations, FUN=function(x) mean(x,na.rm=T)))

quarterly2<-data.frame(qtrnb2)
quarters2<-(as.numeric(substr(row.names(qtrnb2),6,7)))
quarters2[which(quarters2==1|quarters2==2|quarters2==3)]<-1
quarters2[which(quarters2==4|quarters2==5|quarters2==6)]<-2
quarters2[which(quarters2==7|quarters2==8|quarters2==9)]<-3
quarters2[which(quarters2==10|quarters2==11|quarters2==12)]<-4
qtrnb2$quarter<-quarters2

Density<-qtrnb2[which(qtrnb2$quarter==3),1]
}


#this is for making a dataframe for the linear model
#phyto group ("all", "diatom","cyano","green","euglenophytes","chrysophyte","cryptophyte","dinoflagellate","Bacteria", "PicoAlgae unID","Flagellate","chrysoCyst")

if(phytgrp=="diatom"){Diatom<-data.frame("DWA"=summer, "Year"=years, "Taxa"=rep("diatom", length(summer)))} #diatoms
if(phytgrp=="cyano"){Cyano<-data.frame("DWA"=summer, "Year"=years, "Taxa"=rep("cyanobacteria", length(summer)))} #cyanobacteria
if(phytgrp=="green"){Green<-data.frame("DWA"=summer, "Year"=years, "Taxa"=rep("green", length(summer)))} #green algae
if(phytgrp=="chrysophyte"){Chrysophyte<-data.frame("DWA"=summer, "Year"=years, "Taxa"=rep("chrysophyte", length(summer)))} #chrysophyte
if(phytgrp=="cryptophyte"){Cryptophyte<-data.frame("DWA"=summer, "Year"=years, "Taxa"=rep("cryptophyte", length(summer)))} #cryptophyte
if(phytgrp=="dinoflagellate"){Dinoflagellate<-data.frame("DWA"=summer, "Year"=years, "Taxa"=rep("dinoflagellate", length(summer)))} #dinoflagellate
if(phytgrp=="PicoAlgae unID"){Picoplankton<-data.frame("DWA"=summer, "Year"=years, "Taxa"=rep("picoplankton", length(summer)))} #picoplankton

phyto.depth<-summer #for overlap analysis
#}

allDWA<-rbind(Diatom,Cyano,Green,Chrysophyte,Cryptophyte,Dinoflagellate,Picoplankton)
allDWA$Year<-as.numeric(allDWA$Year)
allDWA$Taxa<-as.factor(allDWA$Taxa)


#Summer plot

library(dyn)
library(ggplot2)
dta<-data.frame()
dta<-cbind(summer,as.numeric(years[1:length(summer)]))
mod<-dyn$lm(dta[,1]~dta[,2],na.action=na.exclude)

library(lmtest)
auto<-dwtest(mod) #test for autocorrelation
#DW<-auto$p.value
DW=1
library(nlme)
#if (DW<0.05){mod<-gls(dta[,1]~dta[,2], correlation=corAR1(form=~1),na.action=na.exclude)}

lines(fitted(mod)~years[1:(length(summer))], col="red", lwd=3)
if(DW>0.05){b<-summary(mod);  p<-coefficients(b)[2,4]}
#if(DW<0.05){b<-summary(mod)$tTable; p<-b[2,4]}

results<-data.frame("Year"=as.numeric(years), "Fitted"=fitted(mod))
years<-as.numeric(years)

groupkey<-data.frame(inkey=c("all", "diatom","cyano","green","euglenophytes","chrysophyte","cryptophyte","dinoflagellate","Bacteria", "PicoAlgae unID","Flagellate","chrysoCyst"), real=c("All Groups", "Diatom","Cyanobacteria","Green Algae","Euglenophytes","Chrysophyte","Cryptophyte","Dinoflagellate","Bacteria", "Unidentified Picoplankton","Flagellate","ChrysoCyst"))
currentgroup<-groupkey$real[which(groupkey$inkey==phytgrp)]
datas<-data.frame("Date"=years,"DWA"=summer,"Density"=Density)

ggplot(datas, aes(x=Date, y=DWA, size=Density),legend=FALSE)+
  theme_bw()+
  geom_point(colour="black", shape=16)+ 
   opts(axis.text.x = theme_text(colour = "black", size = 20))+
  opts(axis.text.y = theme_text(colour = "black", size = 20))+
  opts(axis.title.x = theme_text(size = 20, colour = "black"))+
  opts(axis.title.y = theme_text(size = 20, colour = "black", angle = 90), hjust=-0.3)+ 
  scale_y_reverse(name=paste(currentgroup,"Density-weighted Depth"))+
  geom_line(aes(x = Year, y = Fitted),size=1,data=results, colour="red") +
  annotate("text", label = paste("linear trend,","p=",round(p,4),",","DW=",round(DW,2)), x = 1990, y = 20, size = 7, colour = "red")




#linear model with taxa nested within year
model1<-lm(DWA~Taxa/Year, na.action=na.exclude, data=allDWA) #no temporal autocorrelation

model2<-gls(DWA~Taxa/Year, correlation=corAR1(form=~1),na.action=na.exclude, data=allDWA) #temporal autocorrelation
model3<-gls(DWA~Taxa/Year,na.action=na.exclude, data=allDWA)
anova(model2,model3) #which is better? temporal autocorrelation or no?

summary(model2)

model4<-lm(DWA~Year*Taxa,na.action=na.exclude, data=allDWA)

### Test if trends are different from zero
results<-glht(model4, linfct = c("Year = 0","Year:Taxacyanobacteria = 0",
                                 "Year:Taxagreen = 0",
                                 "Year:Taxachrysophyte = 0",
                                 "Year:Taxacryptophyte = 0",
                                 "Year:Taxadinoflagellate = 0",
                                 "Year:Taxapicoplankton = 0"))
summary(results)

#Model fit
par(mar=c(4,5,2,2))
pdata <- expand.grid(Year=seq(1975, 1999, by=1), Taxa=c("Diatom","Cyano","Green","Chrysophyte","Cryptophyte","Dinoflagellate","Picoplankton"))
pdata$Temp <- allDWA$DWA
plot(pdata$Year,pdata$Temp, type="n",xlab="Year", ylab="Density-weighted Depth (m)", cex.lab=1.5,cex.axis=1.5)
points(pdata$Year[1:25], pdata$Temp[1:25], type="p", pch=15)
points(pdata$Year[26:50], pdata$Temp[26:50], type="p", pch=16, col="red")
points(pdata$Year[51:75], pdata$Temp[51:75], type="p", pch=17, col="green")
points(pdata$Year[76:100], pdata$Temp[76:100], type="p", pch=18, col="blue")
points(pdata$Year[101:125], pdata$Temp[101:125], type="p", pch=25, col="orange")
points(pdata$Year[126:150],pdata$Temp[126:150],type="p",pch=4,col="purple")
points(pdata$Year[151:175],pdata$Temp[151:175], type="p",pch=3,col="cyan")
lines(fitted(model4)[1:25]~pdata$Year[1:25])
lines(fitted(model4)[26:50]~pdata$Year[26:50], col="red")
lines(fitted(model4)[51:75]~pdata$Year[51:75], col="green")
lines(fitted(model4)[76:100]~pdata$Year[76:100], col="blue")
lines(fitted(model4)[101:125]~pdata$Year[101:125], col="orange")
lines(fitted(model4)[126:150]~pdata$Year[126:150], col="purple")
lines(fitted(model4)[151:175]~pdata$Year[151:175], col="cyan")
legend("topleft", c("Diatom","Cyano","Green","Chrysophyte","Cryptophyte","Dinoflagellate","Picoplankton"), pch=c(15, 16, 17,18,25),col=c("black","red","green","blue","orange","purple","cyan"), lty=1, lwd=2, cex=0.8)






#regression results and correct p-values







summary(mod)
ps<-c(0.4620,0.0035,0.2211,0.8200,0.9983,0.7230,0.6503)
p.adjust(ps, method = "bonferroni", n = length(ps))