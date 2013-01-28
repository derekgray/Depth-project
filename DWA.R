##############################################################
#############Density weighted average depth###################
##############################################################

##Note: these analyses start with the allzoocommondepths object from the "depth analysis.R" script (must run that first)

#Run this section for analyses of coefficient of dispersion (do not run the section immediately following)
analysis<-"Coefficient of dispersion" 
DWAdispstep1<-summaryBy(DWA+Count~Date1+Depth+KOD, FUN=sum, order=TRUE,data=allzoocommondepths)
DWAdispstep2<-data.frame("DWA"=DWAdispstep1$DWA.sum/DWAdispstep1$Count.sum, "Depth"=DWAdispstep1$Depth,"Date"=DWAdispstep1$Date1, "KOD"=DWAdispstep1$KOD)
custom<-function(x){sd(x, na.rm=T)/mean(x, na.rm=T)}
DWAdispstep3<-summaryBy(DWA~Date+KOD,FUN=custom,data=DWAdispstep2)
DWAdispstep4<-reshape(DWAdispstep3,idvar="Date",v.names="DWA.custom",timevar="KOD",direction="wide") #Convert from long form to wide form (species names as column names instead of row names)
names(DWAdispstep4)<-as.numeric(gsub("\\D", "", names(DWAdispstep4)))  #take only numbers from column names. Necessary because text is introduced by the reshape command
row.names(DWAdispstep4)<-unique(DWAdispstep3$Date)[1:nrow(DWAdispstep4)]; DWAdispstep4[,1]<-NULL #Make rownames sample dates
DWAdispstep4<-DWAdispstep4[,which(colSums(DWAdispstep4,na.rm=T)>0)]
DWAstep4<-DWAdispstep4; DWAstep4<-DWAstep4[which(rownames(DWAstep4)>="1955-01-01"),]

#Run this section for analyses of density weighted mean depth (preceeding section not necessary)
analysis<-"Density weighted depth"
DWAstep1<-summaryBy(DWA+Count~Date1+KOD,data=allzoocommondepths,FUN=sum,order=TRUE) #sum count*depth and count by date
DWAstep2<-data.frame("DWA"=DWAstep1$DWA.sum/DWAstep1$Count.sum, "Date"=DWAstep1$Date1, "KOD"=DWAstep1$KOD)
DWAstep3<-reshape(DWAstep2,idvar="Date",v.names="DWA",timevar="KOD",direction="wide") #Convert from long form to wide form (species names as column names instead of row names)
names(DWAstep3)<-as.numeric(gsub("\\D", "", names(DWAstep3)))  #take only numbers from column names. Necessary because text is introduced by the reshape command
row.names(DWAstep3)<-unique(DWAstep2$Date)[1:nrow(clustall)]; DWAstep3[,1]<-NULL #Make rownames sample dates
DWAstep4<-DWAstep3[,which(colSums(DWAstep3,na.rm=T)>0)]


#Start subsequent analysis common to both of the above sections

#Figure out which KODs represent individual species (i.e. species may be represented by multiple KODs)
matchlocations<-which(names(DWAstep4) %in% key$KOD)
result<-data.frame("Genus"=key$Genus[matchlocations], "Species"=key$Species[matchlocations], "Lifestage"=key$Lifestage_Gen[matchlocations], "Copepod stage"=key$Lifestage_Cop[matchlocations])
wholesplist<-paste(key$Genus,key$Species, sep=" ")
unspecies<-unique(wholesplist)

#This section combines data from different KODs if they represent the same species (but it eliminates double-counts as classified in Lizzie's key)
anydoubles<-c() #For quality control, see in loop below
if(stage=="all"){
  #combine adults and juveniles (combine columns with different KODs)
  for (i in 1:length(unspecies)){
    group<-key$KOD[which(wholesplist %in% unspecies[i])]
    KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")] #which codes are not doubles (according to Lizzie)
    group<-group[which(group%in%KODsnotdouble==TRUE)] #include only codes that are not doubles
    
    #looks for doubles that could indicate KODs that refer to the same count
    if (length(which(duplicated(t(as.data.frame(DWAstep4[,which(names(DWAstep4)%in% group)])),MARGIN=1)==TRUE))>0){
      anydoubles[i]<-names(DWAstep4)[which(duplicated(t(as.data.frame(DWAstep4[,which(names(DWAstep4)%in% group)])),MARGIN=1)==TRUE)]}
    
    if(is.data.frame(DWAstep4[,which(names(DWAstep4)%in% group)])==TRUE){
      if (i==1){sumdataframe<-data.frame(rowMeans(DWAstep4[,which(names(DWAstep4)%in% group)],na.rm=T)); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if(i>1){sumdataframe[,i]<-rowMeans(DWAstep4[,which(names(DWAstep4)%in% group)],na.rm=T); names(sumdataframe)[i]<-as.character(unspecies[i])}}
    if(is.data.frame(DWAstep4[,which(names(DWAstep4)%in% group)])==FALSE){
      if (i==1){sumdataframe<-data.frame(DWAstep4[,which(names(DWAstep4)%in% group)]); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if (i>1){sumdataframe[,i]<-DWAstep4[,which(names(DWAstep4)%in% group)]; names(sumdataframe)[i]<-as.character(unspecies[i])}
    }}
  sumdataframe<-sumdataframe[,which(colSums(sumdataframe, na.rm=T)!=0)] #REMOVE species with no observations
  #sumdataframe<-sumdataframe[,setdiff(1:ncol(sumdataframe), grep("unknown", names(sumdataframe)))] #get rid of species names with "unknown"
  sumdataframe[,which(names(sumdataframe)=="Golomyanka ")]<-NULL
  sumdataframe[,which(names(sumdataframe)=="Macrohectopus branickii")]<-NULL
}

if(stage!="all"){
  #Get data for a select species and separate by life stage
  for (i in 1:length(unspecies)){
    group<-key$KOD[which(wholesplist %in% unspecies[i])]
    KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")] #which codes are not doubles (according to Lizzie's key)
    group<-group[which(group%in%KODsnotdouble==TRUE)] #include only codes that are not doubles
    lifestage<-as.vector(key$Lifestage_Cop[which(wholesplist %in% unspecies[i])])
    newgroup<-grep(stage, lifestage)
    if(is.data.frame(DWAstep4[,which(names(DWAstep4)%in% newgroup)])==TRUE){
      if (i==1){sumdataframe<-data.frame(rowMeans(DWAstep4[,which(names(DWAstep4)%in% newgroup)],na.rm=T)); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if(i>1){sumdataframe[,i]<-rowMeans(DWAstep4[,which(names(DWAstep4)%in% newgroup)],na.rm=T); names(sumdataframe)[i]<-as.character(unspecies[i])}}
    if(is.data.frame(DWAstep4[,which(names(DWAstep4)%in% newgroup)])==FALSE){
      if (i==1){sumdataframe<-data.frame(DWAstep4[,which(names(DWAstep4)%in% group)]); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if (i>1){sumdataframe[,i]<-DWAstep4[,which(names(DWAstep4)%in% group)]; names(sumdataframe)[i]<-as.character(unspecies[i])}
    }}}
row.names(sumdataframe)<-row.names(DWAstep4)

#Take data only from species of interest (or take all data if an individual species wasn't specified)
if (sppec=="all"){series<-sumdataframe
                  series$Date<-as.Date(row.names(sumdataframe))}
if (sppec!="all"){
  series<-data.frame("Count"=sumdataframe[,grep(spp, names(sumdataframe))], "Date"=as.Date(row.names(sumdataframe)))}
newnb<-series

#Identify months with no data and put an NA where data should be (analagous to joining a table with all dates to identify missing samples)
years<-substr(newnb$Date,1,4)
start = as.Date(newnb$Date[which.min(newnb$Date)])
full <- seq(start, by="1 month", length=length(unique(years))*12)
allsp<-data.frame("Date"=full)
for (i in 1:ncol(newnb)){
  ebaicalensis<-data.frame("Count"=with(newnb, newnb[,i][match(full, as.Date(newnb$Date))]))
  allsp<-cbind(allsp, ebaicalensis)}
allsp<-allsp[,2:ncol(allsp)]; names(allsp)<-names(newnb); allsp$Date<-full

nb<-allsp[which(allsp$Date>=startdate),] #take data only after startdate
years<-substr(nb$Date,1,4) #just to use for graphs below

#log transform data?
if(lg=="yes"){nb<-log10(nb+1)}

zooplankton<-nb; row.names(zooplankton)<-zooplankton$Date; zooplankton$Date<-NULL; zooplankton<-zooplankton[,which(colSums(zooplankton,na.rm=T)>0)]

#This removes species with observations below the "minn" threshold set at the top of the code
belowmin<-c()
for (i in 1:ncol(zooplankton)){
  belowmin[i]<-length(which(zooplankton[,i]>0))}
zooplankton<-zooplankton[,which(belowmin>minn)]

##########################
#summarize data by seasons?
##########################

if (seasons=="seasons"){
  for(i in 1:(ncol(zooplankton))){
    library(xts)
    qtnb<-xts(x=zooplankton[,i], order.by=as.Date(row.names(zooplankton)), frequency=12)
    locations<-endpoints(qtnb, "quarters")-1
    locations<-locations[which(locations>0)]
    if(i==1){qtrnb<-data.frame("test"=rep(0,times=(length(locations)-1)))}
    qtrnb<-as.data.frame(period.apply(qtnb, INDEX=locations, FUN=function(x) mean(x,na.rm=T)))
    
    names(qtrnb)<-names(zooplankton)[i]
    if (i==1){quarterly<-data.frame(qtrnb)}
    if (i>1){quarterly[,i]<-qtrnb}
    
    f<-c(1,2,3,4); h<-rep(f, times=1000)
    startmnth<-as.numeric(substr(nb$Date[locations[2]],6,7))
    if (startmnth==12|startmnth==1|startmnth==2){startqtr<-1}
    if (startmnth==3|startmnth==4|startmnth==5){startqtr<-2}
    if (startmnth==6|startmnth==7|startmnth==8){startqtr<-3}
    if (startmnth==9|startmnth==10|startmnth==11){startqtr<-4}
    qtrnb$quarter<-h[(min(which(h==startqtr))):(nrow(qtrnb)+min(which(h==startqtr))-1)]
    
    if(i==1){winter<-data.frame("test"=rep(0,times=(length(which(qtrnb$quarter==1)))))}
    if(i==1){spring<-data.frame("test"=rep(0,times=(length(which(qtrnb$quarter==2)))))}
    if(i==1){summer<-data.frame("test"=rep(0,times=(length(which(qtrnb$quarter==3)))))}
    if(i==1){fall<-data.frame("test"=rep(0,times=(length(which(qtrnb$quarter==4)))))}
    
    winter[,i]<-qtrnb[which(qtrnb$quarter==1),1]; winter$quarter<-NULL
    spring[,i]<-qtrnb[which(qtrnb$quarter==2),1]; spring$quarter<-NULL
    summer[,i]<-qtrnb[which(qtrnb$quarter==3),1]; summer$quarter<-NULL
    fall[,i]<-qtrnb[which(qtrnb$quarter==4),1]; fall$quarter<-NULL
  }}

#Remove species with no observations and add species names as column headings
nonzerolocs<-which(colSums(winter,na.rm=T)>0)
winter<-winter[,nonzerolocs]
names(winter)<-names(zooplankton)[nonzerolocs]; winter2<-winter[,order(names(winter))]
winter<-as.data.frame(winter2); names(winter)<-names(winter2)
nonzerolocs<-which(colSums(spring,na.rm=T)>0)
spring<-spring[,nonzerolocs]; 
names(spring)<-names(zooplankton)[nonzerolocs]; spring2<-spring[,order(names(spring))]
spring<-as.data.frame(spring2); names(spring)<-names(spring2)
nonzerolocs<-which(colSums(summer,na.rm=T)>0)
summer<-summer[,nonzerolocs]
names(summer)<-names(zooplankton)[nonzerolocs]; summer2<-summer[,order(names(summer))] 
summer<-as.data.frame(summer2); names(summer)<-names(summer2)
nonzerolocs<-which(colSums(fall,na.rm=T)>0)
fall<-fall[,nonzerolocs]
names(fall)<-names(zooplankton)[nonzerolocs]; fall2<-fall[,order(names(fall))]
fall<-as.data.frame(fall2); names(fall)<-names(fall2)

#Choose 50 random cells for quality control
if (seasons=="seasons"){
  spp<-list(); num<-list(); qrt<-list()
  
  for (i in 1:1000){
    rowloc<-round(runif(1,1, nrow(winter)),0)
    colloc<-round(runif(1,1,ncol(winter)),0)
    spp[[i]]<-(names(winter)[colloc])
    num[[i]]<-winter[rowloc,colloc]
    qrt[[i]]<-row.names(winter)[rowloc]
  }
  
  QC<-data.frame("Species"=unlist(spp), "Abundance"=unlist(num), "Quarter"=unlist(qrt))
  QC<-orderBy(Abundance~., data=QC)  
}

KODsnotdouble<-key$KOD[which(key$DoubleCount=="N")]
cells<-zoodata[which(zoodata$Species=="baicalensis"&zoodata$Date1<"1948-03-01"&zoodata$Date1>"1947-11-01"),]
cells<-cells[which(as.numeric(cells$KOD)%in%KODsnotdouble),]
write.csv(cells,file="QC2.csv")

if (seasons=="monthly"){
  
  spp<-list(); num<-list(); date<-list()
  
  for (i in 1:1000){
    rowloc<-round(runif(1,1, nrow(winter)),0)
    colloc<-round(runif(1,1,ncol(winter)),0)
    spp[[i]]<-(names(zooplankton)[colloc])
    num[[i]]<-zooplankton[rowloc,colloc]
    date[[i]]<-row.names(zooplankton)[rowloc]
  }
  QC<-data.frame("Species"=unlist(spp), "Abundance"=unlist(num), "Date"=unlist(date))
  QC<-orderBy(Abundance~., data=QC)  }



######Plots###########


#for all members of group (rotifers, copepods, or cladocerans)
par(mfrow=(c(2,2)))
for (ssn in c("winter","spring","summer","fall")){
DA<-rowMeans(get(ssn), na.rm=T)
plot(DA~as.numeric(unique(years)[1:nrow(get(ssn))]), pch=15, xlab="Date", ylab=analysis,main=paste("All",taxa,"Combined,",ssn))
model<-lm(DA~as.numeric(unique(years)[1:nrow(get(ssn))]), na.action=na.exclude)
b<-summary(model); p<-coefficients(b)[2,4]
lines(fitted(model)~as.numeric(unique(years)[1:nrow(get(ssn))]))
legend("topright",lty=1, legend=paste("Linear model,","p=",round(p,4)))
}

#for major species, note uses majspp2 object from depth analyses.R, so need to run that first

fix(majspp2) #Edit the list of major species

par(mfrow=(c(2,2)))
u=0
for (ssn in c("winter", "spring","summer","fall")){
  u=u+1
  spp<-majspp2[[u]] #c("Epischura baicalensis", "Cyclops") #Bosmina", "Daphnia") #Conochilus","Filinia","Kellicottia","Keratella")
  clrs<-c("black","red","blue", "green")
  points<-c(15,16,17,18)
  dta<-list()
  mod<-list()
  p<-c()
  library(Hmisc)
  
  for (i in 1:length(spp)){
    
    specofint<-spp[i]
    
    kt<-grep(first.word(specofint),names(get(ssn)))
    if (length(kt)>1){dta[[i]]<-(rowMeans(get(ssn)[,kt],na.rm=T))
                      mod[[i]]<-lm(rowMeans(get(ssn)[,kt],na.rm=T)~as.numeric(unique(years)[1:nrow(get(ssn))]), na.action=na.exclude)}
    
    if (length(kt)==1){dta[[i]]<-get(ssn)[,kt]
                       mod[[i]]<-lm(get(ssn)[,kt]~as.numeric(unique(years)[1:length(get(ssn)[,kt])]), na.action=na.exclude)}
  }
  
  #plot
  plot(dta[[1]]~as.numeric(unique(years)[1:length(dta[[1]])]), xlab="Year", ylab=analysis, main=ssn, pch=points[1], col=clrs[1], cex.lab=1.5, cex.axis=1.5, ylim=c(0.95*min(unlist(dta),na.rm=T),1.2*(max(unlist(dta),na.rm=T))))
  lines(fitted(mod[[1]])~allyears[1:length(fitted(mod[[1]]))],col=clrs[1]) 
  b<-summary(mod[[1]])
  p[1]<-coefficients(b)[2,4]
  allyears<-unique(years); allyears<-allyears[which(allyears!="NA")]
  
  if (length(spp)>1){
    for (i in 2:length(spp)){
      points(dta[[i]]~as.numeric(unique(years)[1:length(dta[[i]])]),pch=points[i],col=clrs[i])
      lines(fitted(mod[[i]])~allyears[1:length(fitted(mod[[i]]))],col=clrs[i]) 
      b<-summary(mod[[i]])
      p[i]<-coefficients(b)[2,4]
    }}
  
  legend("topright",pch=points, col=clrs,legend=paste(spp,"p=",round(p,4)), cex=0.7)
  
}
