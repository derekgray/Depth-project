##############################################################
#############Density weighted average depth###################
######################## OR ##################################
################coefficient of dispersion#####################
##############################################################

##Note: these analyses start with the allzoocommondepths object from the "depth analysis.R" script (must run that first)
######## 1. The main data file containing all observations
# Check if main file has been loaded, if not read data into R
if(exists("zoodata")==FALSE){zoodata<-read.csv(file="zoopzeros250key1.csv", stringsAsFactors=F)}       #"c:/Users/gray/Dropbox/Postdoc/R data/zoopzeros250key1.csv")

######## 2. The zooplankton key
#Check if the key has been loaded yet, if not read it in 
if(exists("key")==FALSE){key<-read.csv(file="key.csv")}                            #"c:/Users/gray/Dropbox/Postdoc/R data/key.csv")

#group of interest ("Copepod", "Cladoceran", "Rotifer")
taxa<-c("Cladoceran")

#species of interest
majspp.stratified<-list()
majspp.unstratified<-list()
missingdata<-list() #percent of months missing e.g. 72 months of 720 = 10%

for (sppec in c("all","Bosmina")){#Conochilus","Filinia","Kellicottia","Keratella","Notholca","Synchaeta")){
  
  
  #######
  
  #Type of analysis? ("DWA"=Density weighted depth, COD"=coefficient of dispersion)
  analysis<-"DWA"
  
  #start date of data extracted
  startdate=c("1955-01-01")
  
  
  #replace NAs with monthly averages?
  replaceNA<-"no"
  
  #Enter species of interest ("all", or e.g. "Epischura baicalensis")
  #sppec<-c("all")
  
  #Enter lifestage of interest("all", "adult", "naup", "copep") #This has not been tested for anything except for "all" or "adult" although many other options are in the key
  stage=c("all")
  
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
    allzoo<-allzoo[grep(sppec, paste(allzoo$Genus,allzoo$Species)),]}
  
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
  
  
  #Check if all depths were sampled on
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
    #plot(DWAstep2$DWA~as.Date(DWAstep2$Date))
  }
  
  #Calculations for Density weighted depth
  if(analysis=="DWA"){
    library(doBy)
    analysis<-"Density weighted depth"
    DWAstep1<-summaryBy(DWA+Count~Date1,data=allzoocommondepths,FUN=sum,order=TRUE) #sum count*depth and count by date
    DWAstep2<-data.frame("DWA"=DWAstep1$DWA.sum/DWAstep1$Count.sum, "Date"=DWAstep1$Date1)
    #plot(DWAstep2$DWA~DWAstep2$Date)
  }
  
  #########################
  #summarize by month######
  #########################
  DWAstep2$Date<-paste(substr(DWAstep2$Date,1,7),"-01",sep="")
  DWAmonthly<-summaryBy(DWA~Date, data=DWAstep2, FUN=mean, order=TRUE)
  DWAmonthly$Date<-as.Date(DWAmonthly$Date)
  
  #Identify months with no data and put an NA where data should be (analagous to joining a table with all dates to identify missing samples)
  years<-substr(DWAmonthly$Date,1,4)
  start = as.Date(DWAmonthly$Date[which.min(DWAmonthly$Date)])
  full <- seq(start, by="1 month", length=length(unique(years))*12)
  allsp<-data.frame("Date"=full)
  for (i in 1:ncol(DWAmonthly)){
    ebaicalensis<-data.frame("Count"=with(DWAmonthly, DWAmonthly[,i][match(full, as.Date(DWAmonthly$Date))]))
    allsp<-cbind(allsp, ebaicalensis)}
  allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(DWAmonthly); allsp$Date<-full
  names(allsp)<-c("Date","DWA")
  
  missingdata[as.character(sppec)]<-(length(which(is.na(allsp$DWA)))/length(allsp$DWA))*100
  
  #replace NAs with monthly average for whole time series?
  if (replaceNA=="yes"){
    missing<-allsp$Date[which(is.na(allsp$DWA)==T)] 
    for (i in 1:length(missing)){
      allsp$DWA[which(allsp$Date==as.Date(missing[i]))]<-mean(allsp$DWA[which(as.numeric(substr(allsp$Date,6,7))==as.numeric(substr(missing[i],6,7)))],na.rm=T)
    }}
  
  DWAmonthly<-allsp
 
  #Divide into two seasons: stratified and unstratified
  years<-as.numeric(substr(DWAmonthly$Date,1,4))
  
  JASO<-c()
  FMAM<-c()
  for (i in 1:length(unique(years))){
    thisyear<-DWAmonthly[which(years==unique(years)[i]),]
    months<-as.numeric(substr(thisyear$Date,6,7))
    
    JASO[i]<-mean(thisyear$DWA[which(months==7|months==8|months==9|months==10)],na.rm=T)
    FMAM[i]<-mean(thisyear$DWA[which(months==2|months==3|months==4|months==5)],na.rm=T)
  }
  
  majspp.stratified[[as.character(sppec)]]<-JASO
  majspp.unstratified[[as.character(sppec)]]<-FMAM
    
}



#plot by season
years<-unique(as.numeric(substr(DWAmonthly$Date,1,4)))

par(mfrow=c(1,2),mar=c(1,2,1,1),oma = c(2,3,1,0))

for (ssn in c("majspp.stratified","majspp.unstratified")){
  plot(get(ssn)$all~years[1:(length(get(ssn)$all))],xaxt="n",cex.axis=1.5, main=substr(ssn,8,nchar(ssn)))
  axis(side=1, cex.axis=1.5)
  
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


#for major species

majspp2<-names(majspp.stratified)[2:length(majspp.stratified)]

par(mfrow=c(1,2)) #,mar=c(2,2,2,2),oma = c(2,3,1,2))

u=0
for (ssn in c("majspp.stratified","majspp.unstratified")){
  u=u+1
  fix(majspp2)
  spp<-majspp2
  clrs<-c("black","red","blue", "green")
  points<-c(15,16,17,18)
  dta<-list()
  mod<-list()
  p<-c()
  library(Hmisc)
  
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
  plot(dta[[1]]~as.numeric(unique(years)[1:length(dta[[1]])]), xlab="Year", ylab=analysis, main=substr(ssn,8,nchar(ssn)), pch=points[1], col=clrs[1], cex.lab=1, cex.axis=1, ylim=c(0.95*min(unlist(dta),na.rm=T),1.2*(max(unlist(dta),na.rm=T))))
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