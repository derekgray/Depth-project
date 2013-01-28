#calculate DWA for phytoplankton by group

library(doBy)

missing.dates<-c() #when fill NA is "no" this gives the missing dates
percent.missing<-c() #when fill NA is "no" this gives the proportion of missing months

#phyto group ("all", "diatom","cyano","green","euglenophytes","chrysophyte","cryptophyte","dinoflagellate","Bacteria", "PicoAlgae unID","Flagellate","chrysoCyst")
phytgrp<-"PicoAlgae unID"

#Type of quarter ("DJF" or "JFM")
quarter.type<-"JFM"

#should NA dates be filled?
fillNAs<-"yes" 

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
                           phyto$Group[which(phyto$Group=="Flagellates")]<-"Flagellate"
                           phyto$Group[which(phyto$PhytoGenus=="Romeria")]<-"PicoAlgae unID"
                           phyto$Group[which(phyto$PhytoGenus=="Synechococcus")]<-"PicoAlgae unID"
                           }

#take select depths that have been sampled regularly throughout the program
phyto2<-phyto[which(phyto$Depth==0|phyto$Depth==10|phyto$Depth==50|phyto$Depth==100|phyto$Depth==200),] #5, 250? |phyto$Depth==150
#phyto2<-phyto[which(phyto$Depth==0|phyto$Depth==10|phyto$Depth==50),] #250?

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
  phyto5$DWA[which(temp$Count==0)]<-NA}

if (fillNAs=="no"){
phyto5<-NAdates(phyto5)}


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
}

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

years<-unique(substr(phyto5$Date,1,4))



#plot by season
par(mfrow=c(2,2),mar=c(1,2,1,1),oma = c(2,3,1,0))

for (ssn in c("winter","spring","summer","fall")){
  plot(get(ssn)~years[1:(length(get(ssn)))],xaxt="n",ylim=rev(range(get(ssn),na.rm=T)),cex.axis=1.5, main=paste(paste(ssn,":",sep=""), "% months missing=",(round(get(paste("missing.",ssn,sep="")),2)*100)))
  if(ssn=="summer"|ssn=="fall"){axis(side=1, cex.axis=1.5)}
  
  library(dyn)
  dta<-data.frame()
  dta<-cbind(get(ssn),as.numeric(years[1:length(get(ssn))]))
  mod<-dyn$lm(dta[,1]~dta[,2],na.action=na.exclude)
  
  library(lmtest)
  auto<-dwtest(mod) #test for autocorrelation
  DW<-auto$p.value
  library(nlme)
  if (DW<0.05){mod<-gls(dta[,1]~dta[,2], correlation=corAR1(form=~1),na.action=na.exclude)}
  
  lines(fitted(mod)~years[1:(length(get(ssn)))], col="red", lwd=3)
  if(DW>0.05){b<-summary(mod);  p<-coefficients(b)[2,4]}
  if(DW<0.05){b<-summary(mod)$tTable; p<-b[2,4]}
  
  legend("bottomleft",lty=1,lwd=3,cex=1,col="red",legend=paste("linear trend,","p=",round(p,4),",","DW=",round(DW,2)),bty="n")}

mtext("Year", side=1, line=1, outer=TRUE, cex=1.5)
mtext(paste(phytgrp, "DWA"), side=2, line=1, outer=TRUE, cex=1.5)
percent.missing


#make a continuous plot through time------------------------------
par(mfrow=c(1,1))
plot(phyto5, ylim=rev(range(phyto5$DWA,na.rm=T)), xlab="Year",ylab=paste(phytgrp,"DWA"))
library(dyn)
dta<-data.frame()
dta<-phyto5
mod<-dyn$lm(dta[,2]~dta[,1],na.action=na.exclude)

library(lmtest)
auto<-dwtest(mod) #test for autocorrelation
DW<-auto$p.value
library(nlme)
dts<-as.vector(dta[,1])
dpths<-as.vector(dta[,2])

if (DW<0.05){mod<-gls(dpths~dts, correlation=corAR1(form=~1),na.action=na.exclude)}

lines(fitted(mod)~phyto5$Date, col="red", lwd=3)
if(DW>0.05){b<-summary(mod);  p<-coefficients(b)[2,4]}
if(DW<0.05){b<-summary(mod)$tTable; p<-b[2,4]}

legend("bottomleft",lty=1,lwd=3,cex=1,col="red",legend=paste("linear trend,","p=",round(p,4),",","DW=",round(DW,2)),bty="n")

