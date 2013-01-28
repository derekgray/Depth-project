#calculate DWA for phytoplankton by group

#phyto group ("all", "diatom","cyano","green","euglenophytes","chrysophtye","cryptophyte","dinoflagellate","Bacteria","Flagellates","PicoAlgae unID","Flagellate","chrysoCyst")
#phytgrp<-"all"

for (phytgrp in c("all", "diatom","cyano","green","euglenophytes","chrysophtye","cryptophyte","dinoflagellate","Bacteria")){

#Start date
stdate<-"1975-01-01" #most consistent for depths and after preservation change


if(exists("phyto")==FALSE){phyto<-read.csv(file="phytoplankton.csv", stringsAsFactors=F)
                           phyto$DATE<-as.Date(phyto$DATE,"%m/%d/%Y")
                           names(phyto)<-c("Date", "Month","Year","Depth","Code","Count","Group","Genus","Species")
                           phyto<-orderBy(Date~.,data=phyto)}


#take select depths and common genera and check start date
phyto2<-phyto[which(phyto$Depth==0|phyto$Depth==5|phyto$Depth==10|phyto$Depth==50|phyto$Depth==100|phyto$Depth==150|phyto$Depth==200|phyto$Depth==250),]
phyto2<-phyto2[which(phyto2$Date>=stdate),]

if (phytgrp!="all"){phyto2<-phyto2[which(phyto2$Group==phytgrp),]} #which phyto group was requested?

#summarize data by density according to date and depth
phyto2$Date<-paste(substr(phyto2$Date,1,7),"-01",sep="") #make all dates 1st of month for easy monthly averages using summaryBy

library(doBy)
smfun<-function(x){mean(x,na.rm=T)}
phyto3<-summaryBy(Count~Date+Depth, data=phyto2, FUN=smfun, order=TRUE)
phyto3$Depth[which(phyto3$Depth==0)]<-1
phyto3$DWA<-phyto3$Depth*phyto3$Count.smfun
#phyto4<-summaryBy(DWA~Date, FUN=c(mean,length),data=phyto3)
phyto4<-summaryBy(DWA~Date, FUN=mean,data=phyto3)
phyto4$Date<-as.Date(phyto4$Date)
phyto5<-NAdates(phyto4)
assign(paste(phytgrp, "DWA", sep=""),phyto5)}
#phyto3[which(as.Date(phyto3$Date)>as.Date("1995-01-01")),]