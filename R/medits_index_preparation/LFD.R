###========================================================================================####
#    Title:   LFD_4_EWG.R
#    Release: 0.1
#                            ______><((((?>______
# 
#    Description: FUNCTION THAT PRODUCES STANDARDIZED LFD FROM MEDITS DATA
#    Authors:  Matteo MURENU, Alessandro MANNINI, Tristan ROUYER, Chato OSIO, Finlay SCOTT
#    Date: created on June 2015, ISPRA, EWG 15-06 meeting
#    Updates:
#      on August 2015 during EWG 15-11  by Matteo Murenu
#      on December 2015 by Finlay Scott
#
###========================================================================================####
#
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
#
###==================================================================================###

library(reshape2)
library(ggplot2)
library(data.table)
source("LFD_func.R")

#---------------------------------------
# Parameters

datadir <- "./data/"
tabdir <- "./out_test"
plotdir <- "./out_test"
# do you like to split sex?
sex.by <-"n"  # sex combined
#sex.by <-"y"  # split by sex 

# select area (GSA) and species according to MEDIST standard species coding
# which species?
spp <-c("MERL MER","MULL BAR") 
gen <- unlist(lapply(strsplit(spp, " "), function(x) x[1]))
spec <- unlist(lapply(strsplit(spp, " "), function(x) x[2]))
# which GSAs?
gsa <- c("9","11") 
# which length unit? 
# Is this not stored somewhere?
len.unit <- "cm"  # choose between "mm"=millimeters and "cm"=centimeters

# LOAD raw data (MEDITS dataset tables) from csv and subset the GSAs and species we want
TAn <- fread(paste0(datadir, "medits_ta", ".csv"), sep=",", header=T)
TAn <- subset(TAn, area %in% gsa)
TBn <- fread(paste0(datadir, "medits_tb", ".csv"), sep=",", header=T)
TBn <- subset(TBn, (area %in% gsa) & (genus %in% gen) & (species %in% spec))
TCn <- fread(paste0(datadir, "medits_tc", ".csv"), sep=",", header=T)
TCn <- subset(TCn, (area %in% gsa) & (genus %in% gen) & (species %in% spec))

# Convert to data.frames so existing code works (update later)
TAn <- as.data.frame(TAn)
TBn <- as.data.frame(TBn)
TCn <- as.data.frame(TCn)

# Editing
names(TAn)<-toupper(names(TAn))
names(TBn)<-toupper(names(TBn))
names(TCn)<-toupper(names(TCn))

# Drop columns we don't need, e.g. Rigging
# Or better, just keep the ones we do need

# Remove invalid hauls from TA 
TAn<-subset(TAn, VALIDITY=="V")

# Checks
# Check if there are duplicated rows in the files, if there are, running unique will remove duplicates, was a TEMP FIX FOR MERL in AREA 7
if (dim(unique(TAn))[1] != dim(TAn)[1]){
    stop("TA table contains duplicates\n")
}
if (dim(unique(TBn))[1] != dim(TBn)[1]){
    stop("TB table contains duplicates\n")
}
if (dim(unique(TCn))[1] != dim(TCn)[1]){
    stop("TC table contains duplicates\n")
}
# Check Wing opening
if (!all(TAn$WING_OPENING > 0)){
    stop("WING_OPENING less than 0 in TA table\n")
}

# Sometimes NUMBER_OF_THE_STRATUM is not reported (-1) or misreported (not congruent with mean depth)
# What is Number of the stratum (NSTRATE)?

# Need to set CODESTATA to match with medstrata
# Assign the depth strata code (CODESTRATA) using mean depth values (DEPTH)
# i.e "1-50m","51-100m","101-200m","201-500m","501-800m"
TAn$HAULING_DEPTH<-as.numeric(as.character(TAn$HAULING_DEPTH))
TAn$SHOOTING_DEPTH<-as.numeric(as.character(TAn$SHOOTING_DEPTH))
TAn$DEPTH<-(TAn$HAULING_DEPTH+TAn$SHOOTING_DEPTH)/2
TAn$CODESTRATA<-'-'
TAn[which(TAn$DEPTH>0   & TAn$DEPTH<=50), ]$CODESTRATA<-"A"
TAn[which(TAn$DEPTH>50  & TAn$DEPTH<=100),]$CODESTRATA<-"B"
TAn[which(TAn$DEPTH>100 & TAn$DEPTH<=200),]$CODESTRATA<-"C"
TAn[which(TAn$DEPTH>200 & TAn$DEPTH<=500),]$CODESTRATA<-"D"
TAn[which(TAn$DEPTH>500),                 ]$CODESTRATA<-"E"

#------------------------------------------------------
# Merging tables

# Calculate swept area, adjust units of measure to km 
# Wing opening is measured in decimeteres and distance in meters, then:
TAn$DISTANCE<-as.numeric(as.character(TAn$DISTANCE))
TAn$WING_OPENING<-as.numeric(as.character(TAn$WING_OPENING))
TAn$SWEPT<-(TAn$DISTANCE *((TAn$WING_OPENING)/10000))/1000 #Swept area in Km^2

# Merge TA and strata info 
# Load medstra dataframe
medstra<-read.csv2(paste0(datadir,"MEDITS_Strata.csv"), sep=",", header=T) # medits strata
names(medstra)[3:5]<-c("NSTRATE", "CODZONE","DEPTHSTRATA")
# Merge and add columns DEPTHSTRATA, AREASTRATA and CODZONE to TAn
TAn <- merge(medstra,TAn,  by=c("COUNTRY","AREA" ,"CODESTRATA","NSTRATE"), all.y=T)

# merge GENUS and SPECIE in TB and TC
TBn$GENSPE<-paste(TBn$GENUS, TBn$SPECIES,sep = " ")
TCn$GENSPE<-paste(TCn$GENUS, TCn$SPECIES,sep = " ")

# Set columns to be used to merge tables
clmn.TA<- c("COUNTRY","VESSEL","YEAR","HAUL_NUMBER","MONTH","DAY","HAULING_DURATION",
            "VALIDITY","DISTANCE","VERTICAL_OPENING","WING_OPENING", "AREA","DEPTH",
            "DEPTHSTRATA","NSTRATE", "CODESTRATA","AREASTRATA","SWEPT")
clmn.TB<- c("ID_MEDITS_TB","COUNTRY","AREA","VESSEL","YEAR","HAUL_NUMBER","PTOT","NBTOT","GENSPE") 
clmn.TC<- c("COUNTRY","AREA","VESSEL","YEAR","HAUL_NUMBER","CODLON","PFRAC","PECHAN","SEX","NBSEX",
            "LENGTH_CLASS","MATURITY","NBLON","GENSPE","MATSUB","CATFAU")


# Why are we merging here? 
# merge tables - is this right way round - lose some rows here
# If using data.table need to set keys
# SO converted to data.frames above
# Why merge this - do we need TB?
TBn <- merge(TBn[clmn.TB], TAn[clmn.TA], by=c("COUNTRY","YEAR","VESSEL","AREA" ,"HAUL_NUMBER"), all.y=TRUE) # UPDATE 15/7/2014 now also with by= "MONTH", "DAY"
# Add columns to TCn - which ones?
TCn <- merge(TAn[clmn.TA], TCn[clmn.TC], by=c("COUNTRY","YEAR","VESSEL","AREA" ,"HAUL_NUMBER"), all.y=TRUE)

# rename column names
colnames(TBn)[colnames(TBn)=="NBTOT"] <- "NBTot" #"TOTAL_NUMBER_IN_HAUL"
colnames(TBn)[colnames(TBn)=="PTOT"] <- "PTot" #"TOTAL_WEIGHT_IN_HAUL"
colnames(TCn)[colnames(TCn)=="NBLON"] <- "NBLEN"
colnames(TCn)[colnames(TCn)=="PECHAN"] <- "WEIGHT_OF_THE_SAMPLE_MEASURED"
colnames(TCn)[colnames(TCn)=="PFRAC"] <-   "WEIGHT_OF_THE_FRACTION"
colnames(TCn)[colnames(TCn)=="CODLON"] <-   "LENGTH_CLASSES_CODE"
colnames(TCn)[colnames(TCn)=="MATSUB"] <-   "MATURITY_SUB_STAGING"

# some factors to numbers
TCn$NBLEN<-as.numeric(as.character(TCn$NBLEN))
TCn$NBSEX<-as.numeric(as.character(TCn$NBSEX))
TCn$WEIGHT_OF_THE_FRACTION<-as.numeric(as.character(TCn$WEIGHT_OF_THE_FRACTION))
TCn$WEIGHT_OF_THE_SAMPLE_MEASURED<-as.numeric(as.character(TCn$WEIGHT_OF_THE_SAMPLE_MEASURED))
TCn$LENGTH_CLASS<-as.numeric(as.character(TCn$LENGTH_CLASS))
TBn$NBTot<-as.numeric(as.character(TBn$NBTot))
TBn$PTot<-as.numeric(as.character(TBn$PTot))


# Adjust lengths based on unit
# What is LENGTH_CLASSES_CODE?
if(len.unit=="cm"){
  TCn$length <- ifelse(TCn$LENGTH_CLASSES_CODE=="M", TCn$LENGTH_CLASS/10, TCn$LENGTH_CLASS)
} else {
  TCn$length <- ifelse(TCn$LENGTH_CLASSES_CODE!="M", TCn$LENGTH_CLASS*10, TCn$LENGTH_CLASS)
}

# Raise the measured sample to the whole catch 
TCn$NBLEN.raised <- TCn$NBLEN*TCn$WEIGHT_OF_THE_FRACTION/TCn$WEIGHT_OF_THE_SAMPLE_MEASURED

#--------------------------------------------------------------
#  PRODUCE INDICES FROM MEDITS DATA 

# Loop over GSA
for (j in 1:length(gsa)) {
    # Subset the GSA
    TA<-subset(TAn, AREA==gsa[j])
    #TBgsa<-subset(TBn, AREA==gsa[j])
    #TCgsa<-subset(TCn, AREA==gsa[j])
    cat("Processing gsa: ", gsa[j], "\n")
    # Loop over species
    for (k in 1:length(spp)) {
        TB <-subset(TBn, GENSPE==spp[k] & AREA == gsa[j])
        TC <-subset(TCn, GENSPE==spp[k] & AREA == gsa[j])
        TC <-droplevels(TC)
        cat("Processing species: ", spp[k], "\n")
        if (sex.by != "y") {
            sex <-"all"
            cat("Processing sex combined: ", sex, "\n")
            #source(paste0(scriptdir,"LFD_fun.R"))
            lfd(TA, TB, TC, medstra, sex, len.unit, tabdir, plotdir, plots=TRUE)
        }
        else {
            # Loop over sex 
            for (sex in c("F","M","I")) {  # N= not determinated give problems 
                TC.sx<-TC
                TC<-TC.sx[as.character(TC.sx$SEX)==sex,]
                cat("Processing sex: ", sex, "\n")
                #source(paste0(scriptdir,"LFD_fun.R"))
                lfd(TA, TB, TC, medstra, sex, len.unit, tabdir, plotdir, plots=TRUE)
                TC<-TC.sx
            }
        }
    }
} 

# Remove some not needed staff  ####
#rm(lat,lat_end,lat_end2,LatEndDeg,LatEndMin,LatEndSec,
#   lat_start,lat_start2,LatStartDeg,LatStartMin,LatStartSec,
#   lon,lon_end,lon_end2,LonEndDeg,LonEndMin,LonEndSec,
#   lon_start,lon_start2,LonStartDeg,LonStartMin,LonStartSec)
#rm(clmn.TA,clmn.TB,clmn.TC,i,j,k,m,year,tbl,tmp,dbtype,database,
#   stra.fac,stra.sur,stra.sur.t,stratified.N,swept.y)
#rm(lfd.std,lfd.std.t,lfd.y,lfd.y2)
#rm(TA,TB,TC,TB1,TC1)

# rm(len.unit,codspe,EWG,gsa,hm.gsa,hm.sp,p,sex.by,species,spp,nspp)
# rm(cat04,catches,lan09,land,tor5,tor5.ls,effort,disc,medstra,nstra11)
# rm(TA.m,TAn,TBn,TCn)


# Lat and Lon stuff

#Estimate mid point in Haul and convert to Decimal Degrees for plotting
#lat_start<-as.numeric(as.character(TAn$SHOOTING_LATITUDE))
#lon_start<-as.numeric(as.character(TAn$SHOOTING_LONGITUDE))
## lat_start<-TAn$SHOOTING_LATITUDE
## lon_start<-TAn$SHOOTING_LONGITUDE
#LatStartSec = (lat_start - floor(lat_start))/100;
#LonStartSec = (lon_start - floor(lon_start))/100;
#LatStartDeg = floor(floor(lat_start)/100);
#LonStartDeg = floor(floor(lon_start)/100);
#LatStartMin = (floor(lat_start)/100 - LatStartDeg)/60*100;
#LonStartMin = (floor(lon_start)/100 - LonStartDeg)/60*100;
#
#lat_end<-as.numeric(as.character(TAn$HAULING_LATITUDE))
#lon_end<-as.numeric(as.character(TAn$HAULING_LONGITUDE))
## lat_end<-TAn$HAULING_LATITUDE
## lon_end<-TAn$HAULING_LONGITUDE
#LatEndSec = (lat_end - floor(lat_end))/100;
#LonEndSec = (lon_end - floor(lon_end))/100;
#LatEndDeg = floor(floor(lat_end)/100);
#LonEndDeg = floor(floor(lon_end)/100);
#LatEndMin = (floor(lat_end)/100 - LatEndDeg)/60*100;
#LonEndMin = (floor(lon_end)/100 - LonEndDeg)/60*100;
#
#lat_start2 = LatStartDeg + LatStartMin + LatStartSec;
#lon_start2 = LonStartDeg + LonStartMin + LonStartSec;
#
#lat_end2 = LatEndDeg + LatEndMin + LatEndSec;
#lon_end2 = LonEndDeg + LonEndMin + LonEndSec;
#
#lat_start = lat_start2;
#lon_start = lon_start2;
## use the quadrant to identify negative longitudes
#lon_start<-ifelse(TAn$SHOOTING_QUADRANT==7, lon_start*-1, lon_start)
#lat_end = lat_end2;
#lon_end = lon_end2;
## use the quadrant to identify negative longitudes
#lon_end<-ifelse(TAn$HAULING_QUADRANT==7, lon_end*-1, lon_end)
#
##FIXED MID HAUL POSITION
#lat = (lat_start+lat_end)/2 
#lon = (lon_start+lon_end)/2 
#
#TAn$Latitude<-lat
#TAn$Longitude<-lon

