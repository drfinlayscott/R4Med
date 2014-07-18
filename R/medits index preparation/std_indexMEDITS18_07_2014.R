#################################################################################
#    Copyright (C) <2014>  <Chato Osio, Finlay Scott, Tristan Rouyer>

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/
###################################################################################


###################################
library(RODBC) # odbcConnectAccess
library(doBy)  # summaryBy
library(mapdata)  	# 'worldHires' map
library(maps)     	# mapping 'worldHires'
library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)


# FIRST STEP IS TO SET THE WORK SPACE TO MAKE TO CODE RUN, chose the directory where you have copied the "script" folder

#setwd("S:/data coverage report 2013")
setwd("C:/EWG14_09")

# select area (GSA) and species according to MEDIST standard species coding
sp <- "MERL MER"
gsa <- "17"
len.unit <- "cm"  # choose between "mm"=millimeters and "cm"=centimeters


# connection and data upload
# output: dataset by age 'TB' 

source("db_connection_16_07_2014FINAL.R")

# response variable, choose between weight (Ptot) or number (NBTot)
y <- "CPUE"
TB$CPUE <- TB$NBTot/TB$swept
#TB$CPUE <- TB$PTot/TB$swept

#TB$MONTH<-as.numeric(as.character(TB$MONTH))
#TB$YEAR<-as.numeric(as.character(TB$YEAR))

# Check the data
lines(ddply(TC, "YEAR", summarise, ntot = sum(NBLEN.raised)), type="l")
lines(ddply(TB, "YEAR", summarise, nbtot = sum(NBTot)), type="b", col="red")

# plot CPUE
plot(ddply(TB, "YEAR", summarise, cpue = sum(CPUE)), type="l")

###################################
# For the calculation of a standardized index, chose either the Cochrane stratified means



# these are already defined in the std_indexMEDITS
gsa<- "7"
species<-unique(TC$genspe)
# Need to select the 
size<-seq(min(TC$LENGTH_CLASS), max(TC$LENGTH_CLASS), 5) # size threshold for data aggregation (below/above). If 0 no account for size

year<-seq(min(TC$YEAR), max(TC$YEAR)) # number of years for which stratified numbers at length should be calculated, as is will select all the year range, it can be overridden by manual selection of years by using:   year<-seq(1994,2013,1)

source("stratifiedmeans16_07_2014.r")

#this returns a stratified_N_atlength data frame

# OR 

# data exploration and model selection
# open the script "model_selection.R" and select the appropriate model



###################################
# Examples of Regression models that could be used for CPUE standardization

library(mgcv)


# model formulation
mod <- glm(log(CPUE)~factor(YEAR)+MONTH.x+Latitude+Longitude, family=gaussian, data=subset(TB, CPUE>0))

modqp <- glm(CPUE~factor(YEAR)+MONTH.x+Latitude+Longitude+DEPTH, family=quasipoisson, data=TB)
modqp <- glm(CPUE~factor(YEAR)+factor(MONTH.x)+Latitude*Longitude+DEPTH, family=quasipoisson, data=TB)

mod<-gam(log(CPUE)~s(YEAR, k=5)+s(MONTH.x, k=3)+s(DAY.x, k=4)+s(Latitude,Longitude, k=4)+s(DEPTH, k=4),data=subset(TB, CPUE>0))

mod1<-gam(log(PTot+1)~s(YEAR)+s(MONTH.x, k=5)+s(Latitude,Longitude)+s(DEPTH)+s(swept),data=subset(TB, swept>0))
mod2<-gam(log(NBTot+0.001)~s(YEAR, k=5)+s(MONTH.x, k=2)+s(Latitude,Longitude)+s(DEPTH, k=4)+offset(log(swept)),data=subset(TB, swept>0))

#################################################################
#################################################################
#################################################################
# model summary and CPUE predictions

### save model summary statistics
sink(paste("mod_summary_",substring(sp,1,4),substring(sp,6,8),".txt",sep=""))
summary(mod)
sink()

# if LateX is used the table can be 
   
  print(xtable(summary(mod)), type="latex", file="../output/mod_summary_",floating=TRUE, floating.environment="table", table.placement = "ht", align="p{3cm}")

### diagnostic plot of the fitted model
tiff(file=paste("../output/mod_residuals_",substring(sp,1,4),substring(sp,6,8),".tiff",sep=""), width=15, height=15, bg="white", units="cm", compression="none", res=200)

par(mfrow=c(2,2))
plot(mod)


dev.off()


### compile and save annual CPUE observed and predicted
yr.ref <- as.numeric(names(which.max(table(TB$YEAR))))

pred.grid <- TB[TB$YEAR==yr.ref,]
tmp <- rep(unique(TB$YEAR),rep(dim(pred.grid)[1],length(unique(TB$YEAR))))

for(i in 1:(length(unique(TB$YEAR))-1)){
  pred.grid <- rbind(pred.grid,TB[TB$YEAR==yr.ref,])
}
pred.grid$YEAR <- tmp

pred <- predict(mod, pred.grid, type="response", se.fit=TRUE) 
obs.yr <- aggregate(TB[,y], list(TB$YEAR), mean)
pred.yr <- aggregate(pred$fit, list(pred.grid$YEAR), mean)
dat.out <- data.frame(YEAR=obs.yr[,1], obs=round(obs.yr$x,3), pred=round(pred.yr$x,3))


write.table(dat.out,paste("CPUE_predictions_",substring(sp,1,4),substring(sp,6,8),".txt",sep=""), col.names=TRUE, row.names=FALSE, sep="\t")





#################################################################
#################################################################
#################################################################
# annual maps of CPUE observed "...
# diagnostic plots of model fitting "...
# annual CPUE observed and predicted "CPUEobsANDpred_MERLMER.tiff"

### mapping observed mean CPUE by annual
library(mapdata)  	# 'worldHires' map
library(maps)     	# mapping 'worldHires'

# function tiff will produce a tiff image, however uncompressed, as in R 2.8 there is a bug when using compress="jpeg", in R 2.10 this should be fixed
tiff(file=paste("map_logCPUEobs_",substring(sp,1,4),substring(sp,6,8),Sys.Date(),".tif",sep=""), width=20, height=20, bg="white", units="cm", compression="none", res=200)
par(mfrow=c(4,5), mai=c(0.3,0.4,0,0), omi=c(0.6,0.6,0.1,0.1))

yr.vec <- unique(TB$YEAR)

for(i in 1:length(yr.vec)){
  plot(1,1,type="n",xlim=c(min(TB$Longitude)-0.1, max(TB$Longitude)+0.1), ylim=c(min(TB$Latitude)-0.1, max(TB$Latitude)+0.1), xlab="Longitude", ylab="Latitude")

  map("worldHires", fill=T, col="black",add=T)

  points(TB$Longitude[TB$YEAR==yr.vec[i] & TB[TB$YEAR==yr.vec[i],y]>0], TB$Latitude[TB$YEAR==yr.vec[i] & TB[TB$YEAR==yr.vec[i],y]>0], cex=3*(log(TB[TB$YEAR==yr.vec[i] & TB[TB$YEAR==yr.vec[i],y]>0,y]+1)/max(log(TB[,y]+1))), col="red")
  points(TB$Longitude[TB$YEAR==yr.vec[i] & TB[TB$YEAR==yr.vec[i],y]==0], TB$Latitude[TB$YEAR==yr.vec[i] & TB[TB$YEAR==yr.vec[i],y]==0], pch=4, col="red")

  legend("topleft",paste(yr.vec[i]), text.col="red", bty="n")
}

mtext(expression(paste("Longitude (",degree,"E)")), side=1, outer=T, line=2)
mtext(expression(paste("Latitude (",degree,"N)")), side=2, outer=T, line=2)

dev.off()






### annual CPUE observed and predicted 
# THis is now set up to predict for a glm or a gam where the CPUE is modeled on the Log scale
yr.ref <- as.numeric(names(which.max(table(TB$YEAR))))

pred.grid <- TB[TB$YEAR==yr.ref,]
tmp <- rep(unique(TB$YEAR),rep(dim(pred.grid)[1],length(unique(TB$YEAR))))

for(i in 1:(length(unique(TB$YEAR))-1)){
  pred.grid <- rbind(pred.grid,TB[TB$YEAR==yr.ref,])
}
pred.grid$YEAR <- tmp

pred <- predict(mod, pred.grid, type="response", se.fit=TRUE)
pred.yr <- aggregate(exp(pred$fit), list(pred.grid$YEAR), mean)
pred.yr.se.up <- aggregate(exp(pred$fit+1.96*pred$se.fit), list(pred.grid$YEAR), mean)
pred.yr.se.lo <- aggregate(exp(pred$fit-1.96*pred$se.fit), list(pred.grid$YEAR), mean)

obs.yr <- aggregate((TB[,y]), list(TB$YEAR), mean)


#tiff(file=paste("../output/CPUEobsANDpred_",substring(sp,1,4),substring(sp,6,8),".tiff",sep=""), width=20, height=15, bg="white", units="cm", compression="none", res=200)

plot(obs.yr[,1], obs.yr$x, type="p", ylim=c(min(pred.yr.se.lo$x, obs.yr$x),max(pred.yr.se.up$x, obs.yr$x)*1.1), xlab="Year", ylab="mean CPUE", cex=0.8, pch=19)
lines(pred.yr[,1], pred.yr$x, lty=1)
lines(pred.yr.se.up[,1], pred.yr.se.up$x, lty=2)
lines(pred.yr.se.lo[,1], pred.yr.se.lo$x, lty=2)

legend("topleft", paste(c("observed","predicted","95% CI")), pch=c(1,NA,NA), lty=c(NA,1,2), bty="n", cex=0.8)

dev.off()




