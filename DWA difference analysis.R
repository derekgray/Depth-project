#Difference in DWA through time


allDWA<-allDWA[which(allDWA$Year>=1975&allDWA$Year<=1999),]
allphyt<-rep(phyto.depth,6)
allDWA$DWA<-allDWA$DWA-allphyt #Difference between phytoplankton DWA and zooplankton DWA

library(dyn)
model2<-dyn$lm(DWA~Year*Taxa, na.action=na.exclude, data=allDWA) #no temporal autocorrelation
summary(model2)

library(multcomp)
results<-glht(model2, linfct = c("Year=0","Year:TaxaNaup = 0",
                                 "Year:TaxaAdult = 0",
                                 "Year:TaxaRotifer = 0",
                                 "Year:TaxaCladoceran = 0",
                                 "Year:Taxaall=0"))

summary(results)

#Model fit
pdata <- expand.grid(Year=seq(1975, 1999, by=1), Taxa=c("Copep","Naup","Adult","Rotifer","Cladoceran","all"))
pdata$Temp <- allDWA$DWA     
#plots
par(mfrow=c(2,3))
taxas<-unique(allDWA$Taxa)
for (i in 1:length(taxas)){
plot(pdata$Year[which(pdata$Taxa==taxas[i])],pdata$Temp[which(pdata$Taxa==taxas[i])], xlab="Year", ylab="Difference in DWA", cex.lab=1.5,cex.axis=1.5, pch=16, main=taxas[i], cex=1.5)
lines(fitted(model2)[which(allDWA$Taxa==taxas[i])]~pdata$Year[which(allDWA$Taxa==taxas[i])], col="black", lty=1)
}

