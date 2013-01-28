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
########sure that the data summary is providing the correct output##############
################################################################################

########Output consists of these objects:
######1. "zooplankton": If monthly time series is requested, "zooplankton" is the monthly time series for the data requested (abundance over time)
######2. "winter","spring","summer","fall": If seasonal time series is requested, there will be four objects "winter", "spring", "summer", "fall" that summarize abundance over time by seasons for the species requested
######3. "species": list of species found in the data
######4. "Rabundance": % of total individuals by species. For seasonal data the list is composed of 4 elements in this order: winter, spring, summer, fall
######5. "Major": list of the major species (>10% abundance)
######6. "QC": 50 randomly chosen cells from the dataframe for comparison with Lyubov's original data
######7. "anydoubles": For each species, were any codes doubles? Uses the anyDuplicated function to look for duplicated columns in dataframe. If all zeroes, then no duplicates detected
######8. "quarterly": Average abundance by quarter for all species
########Code requires the following packages: doBy, xts, zoo


########Before starting two files must be read in to the R environment: 
######## 1. The main data file containing all observations
# Check if main file has been loaded, if not read data into R
if(exists("zoodata")==FALSE){zoodata<-read.csv(file="zoopzeros250key1.csv", stringsAsFactors=F)}       #"c:/Users/gray/Dropbox/Postdoc/R data/zoopzeros250key1.csv")
                            
######## 2. The zooplankton key
#Check if the key has been loaded yet, if not read it in 
if(exists("key")==FALSE){key<-read.csv(file="key.csv")}                            #"c:/Users/gray/Dropbox/Postdoc/R data/key.csv")


#######
#Enter depth of interest in meters (can be <, >, or range, e.g. "0to10")
#Depths available: 0to10, 10to25, 50to100, 100to150, 150to250, lessthan250
interest<-c("lessthan250")

#start date of data extracted
startdate=c("1955-01-01")

#monthly, or summarize by seasons (options= "seasons" or "monthly")
seasons<-c("seasons")

#group of interest ("all", "Copepod", "Cladoceran", "Rotifer")
taxa<-c("Rotifer")

#Enter species of interest ("all", or e.g. "Epischura baicalensis")
sppec<-c("all")

#Enter lifestage of interest("all", "adult")
stage=c("all")

#min number of observations per species by month (approx 720 months in dataset) species with fewer observations discarded)
minn=1 #5%

#log transform the data?
lg<-c("no")


#start data processing
allzoo<-zoodata

allzoo<-allzoo[which(allzoo$Group_General==taxa),] #take records for group of interest

#create a volume column and total number of individuals column for later calculations
allzoo$Volume<-(as.numeric(allzoo$Nig_Gr)-as.numeric(allzoo$Ver_Gr))*(pi*(0.375/2)^2)*1000 # 0.375 diameter net, multiply by 1000 to convert to litres
allzoo<-allzoo[which(allzoo$Volume!="NA"),] #removes records with no volume
allzoo$Tot<-allzoo$Volume*allzoo$Count #Add a column that shows total number of individuals (volume sampled multiplied by Lyubov's densities (#/m3))
#create single depth column by merging Ver_Gr and Nig_Gr
allzoo$Depth<-paste(allzoo$Ver_Gr,allzoo$Nig_Gr, sep="-")
allzoo$DWA<-((as.numeric(allzoo$Nig_Gr)+as.numeric(allzoo$Ver_Gr))/2)*allzoo$Count #create a column for DWA calculations
  
#create separate data sets for each depth category (0-10, 10-25, 50-100, 100-150, 150-250)
allzoo0to10<-allzoo[which(allzoo$Depth=="0-10"),]; allzoo0to10$Depth<-NULL
allzoo10to25<-allzoo[which(allzoo$Depth=="10-25"),]; allzoo10to25$Depth<-NULL
allzoo50to100<-allzoo[which(allzoo$Depth=="50-100"),]; allzoo50to100$Depth<-NULL
allzoo100to150<-allzoo[which(allzoo$Depth=="100-150"),]; allzoo100to150$Depth<-NULL
allzoo150to250<-allzoo[which(allzoo$Depth=="150-250"),]; allzoo150to250$Depth<-NULL
allzoolessthan250<-allzoo[which(allzoo$Depth=="0-10"|allzoo$Depth=="10-25"|allzoo$Depth=="50-100"|allzoo$Depth=="100-150"|allzoo$Depth=="150-250"),]

#Check if a depth range was requested
allzoocommondepths<-get(paste("allzoo", interest, sep=""))

#summarize by month
allzoocommondepths$Date1<-paste(substr(allzoocommondepths$Date1,1,7),"-01",sep="") #convert all samples to 1st of month to make a regular time series

library(doBy)
allzoo3<-summaryBy(Tot+Volume~Date1+KOD, data=allzoocommondepths, FUN=sum, order=TRUE)
allzoo3$Density<-allzoo3$Tot.sum/allzoo3$Volume.sum #Figure out density by adding all individuals and dividing by total volume sampled
allzoo3<-data.frame("Date"=allzoo3$Date1, "Density"=allzoo3$Density, "KOD"=allzoo3$KOD) #Make a nice dataframe with Date, Density, and KOD

clustall<-reshape(allzoo3,idvar="Date",v.names="Density",timevar="KOD",direction="wide") #Convert from long form to wide form (species names as column names instead of row names)
names(clustall)<-as.numeric(gsub("\\D", "", names(clustall)))  #take only numbers from column names. Necessary because text is introduced by the reshape command
row.names(clustall)<-unique(allzoo3$Date)[1:nrow(clustall)]; clustall[,1]<-NULL #Make rownames sample dates
clustall<-clustall[,which(colSums(clustall,na.rm=T)>0)]

#Figure out which KODs represent individual species (i.e. species may be represented by multiple KODs)
matchlocations<-which(names(clustall) %in% key$KOD)
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
    if (length(which(duplicated(t(as.data.frame(clustall[,which(names(clustall)%in% group)])),MARGIN=1)==TRUE))>0){
    anydoubles[i]<-names(clustall)[which(duplicated(t(as.data.frame(clustall[,which(names(clustall)%in% group)])),MARGIN=1)==TRUE)]}
   
    if(is.data.frame(clustall[,which(names(clustall)%in% group)])==TRUE){
      if (i==1){sumdataframe<-data.frame(rowMeans(clustall[,which(names(clustall)%in% group)],na.rm=T)); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if(i>1){sumdataframe[,i]<-rowMeans(clustall[,which(names(clustall)%in% group)],na.rm=T); names(sumdataframe)[i]<-as.character(unspecies[i])}}
    if(is.data.frame(clustall[,which(names(clustall)%in% group)])==FALSE){
      if (i==1){sumdataframe<-data.frame(clustall[,which(names(clustall)%in% group)]); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if (i>1){sumdataframe[,i]<-clustall[,which(names(clustall)%in% group)]; names(sumdataframe)[i]<-as.character(unspecies[i])}
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
    if(is.data.frame(clustall[,which(names(clustall)%in% newgroup)])==TRUE){
      if (i==1){sumdataframe<-data.frame(rowMeans(clustall[,which(names(clustall)%in% newgroup)],na.rm=T)); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if(i>1){sumdataframe[,i]<-rowMeans(clustall[,which(names(clustall)%in% newgroup)],na.rm=T); names(sumdataframe)[i]<-as.character(unspecies[i])}}
    if(is.data.frame(clustall[,which(names(clustall)%in% newgroup)])==FALSE){
      if (i==1){sumdataframe<-data.frame(clustall[,which(names(clustall)%in% group)]); names(sumdataframe)[1]<-as.character(unspecies[1])}
      if (i>1){sumdataframe[,i]<-clustall[,which(names(clustall)%in% group)]; names(sumdataframe)[i]<-as.character(unspecies[i])}
      }}}
row.names(sumdataframe)<-row.names(clustall)

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

if(i==1){winter<-data.frame("test"=rep(0,times=length(which(qtrnb$quarter==1))))}
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


#Identify "major" species representing >10% of individuals
if (seasons=="seasons"){
  b=0
  Rabundance<-list()
  Major<-list()
  for (i in c("winter","spring","summer","fall")){  
  b=b+1
  nc<-0
  try(nc<-ncol(get(i)))
  if(length(nc)>0){Rabundance[[b]]<-colSums(get(i),na.rm=T)/sum(rowSums(get(i),na.rm=T))#percent of individuals over entire season
  Major[[b]]<-Rabundance[[b]][which(Rabundance[[b]]>0.1)]}
  
  if(length(nc)==0)  {Rabundance[[b]]<-1; names(Rabundance[[b]])<-names(get(i))[1]}
  } 
  #identify names of major species
  majspp<-list()
  for (i in 1:length(Major)){
    majspp[[i]]<-names(unlist(Major[i]))
  library(stringr)
  majspp2<-list()
  }
  for (i in 1:4){
    majspp2[[i]]<-unique(str_split_fixed(unlist(majspp[[i]]), " ", n = 2)[, 1])}
  }
        
if (seasons=="monthly"){  Rabundance<-colSums(zooplankton,na.rm=T)/sum(rowSums(zooplankton,na.rm=T))
                           Major<-Rabundance[which(Rabundance>0.1)]}


###write data to file

write.csv(summer[1:48,],file=paste("Summer",taxa,"abundance",".csv"), row.names=unique(years)[1:48])
write.csv(fall[1:48,],file=paste("Fall",taxa,"abundance",".csv"), row.names=unique(years)[1:48])

######Plots for major species (>10% abundance) by season

#Manual edits to major species, if needed
#fix(majspp2)

pdf(file=paste("taxa_",taxa,"_","abundance",".pdf",sep=""), width=8.5, height=11) #open graphics device for pdf

par(mfrow=(c(2,1)))
u=2 #seasons are numbered 1-4 (winter,spring,summer,fall)
for (ssn in c("summer","fall")){
u=u+1
spp<-majspp2[[u]] #c("Epischura baicalensis", "Cyclops") #Bosmina", "Daphnia") #Conochilus","Filinia","Kellicottia","Keratella")
clrs<-c("black","red","blue", "green")
points<-c(15,16,17,18) #plot characters for species
dta<-list()
mod<-list()
p<-c()
library(Hmisc)

for (i in 1:length(spp)){

specofint<-spp[i]

kt<-grep(first.word(specofint),names(get(ssn)))
if (length(kt)>1){dta[[i]]<-rowSums(get(ssn)[1:48,kt],na.rm=T)
                  mod[[i]]<-lm(rowSums(get(ssn)[1:48,kt])~as.numeric(unique(years)[1:48]), na.action=na.omit)}

if (length(kt)==1){dta[[i]]<-get(ssn)[1:48,kt]
mod[[i]]<-lm(get(ssn)[1:48,kt]~as.numeric(unique(years)[1:length(get(ssn)[1:48,kt])]), na.action=na.omit)}
}

#plot
allyears<-unique(years); allyears<-allyears[which(allyears!="NA")]

plot(dta[[1]]~as.numeric(unique(years)[1:length(dta[[1]])]), xlab="Year", ylab="Number per liter", main=ssn, pch=points[1], col=clrs[1], cex.lab=1.5, cex.axis=1.5, ylim=c(0,1.05*(max(unlist(dta),na.rm=T))))
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

legend("topleft",pch=points, col=clrs,legend=paste(spp,"p=",round(p,4)), cex=1.5)
}
dev.off() #close graphics device to flush figure to pdf
