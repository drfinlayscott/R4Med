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
library(dplyr)
library(plyr)

#---------------------------------------
# Parameters

# Set these for your computer
# Location of TA, TB, TC and Medits Strata tables
datadir <- "../../../tables/medits/"
# Where do you want the results to go
#tabdir <- "../../../../medits_test/tables"
#plotdir <- "../../../../medits_test/figures"
# Do you want to split sex?
sex.split <- FALSE # TRUE  or FALSE

# select area (GSA) and species according to MEDIST standard species coding
# which species?
#spp <-c("MERL MER","MULL BAR") 
#spp <-c("MERL MER")
#spp <-c("NEPR NOR")
spp <-c("PAPE LON")
gen <- unlist(lapply(strsplit(spp, " "), function(x) x[1]))
spec <- unlist(lapply(strsplit(spp, " "), function(x) x[2]))
# which GSAs?
#gsa <- c("9","11") 
#gsa <- 9
gsa <- 17 
# which length unit? 
# Is this not stored somewhere?
len.unit <- "cm"  # choose between "mm"=millimeters and "cm"=centimeters

#---------------------------------------------------
# Data loading and preparation

# Load raw data (MEDITS dataset tables) from csv and subset the GSAs and species we want
# Check the filenames of the tables and the serator (here ;)
TAn <- fread(paste0(datadir, "ta", ".csv"), sep=";", header=T)
TAn <- subset(TAn, area %in% gsa)
TCn <- fread(paste0(datadir, "tc", ".csv"), sep=";", header=T)
TCn <- subset(TCn, (area %in% gsa) & (genus %in% gen) & (species %in% spec))
# You may get warnings - may just be character conversions for maturity column

# For consistency
names(TAn)<-toupper(names(TAn))
names(TCn)<-toupper(names(TCn))

# Just pull out one year
# 2008 has lots of NAs
#TAn <- subset(TAn, YEAR==2005) 
#TCn <- subset(TCn, YEAR==2005) 

# Convert to data.frames until I understand data.tables properly
#TAn <- as.data.frame(TAn)
#TBn <- as.data.frame(TBn)
#TCn <- as.data.frame(TCn)

## Fix problem in TCn length class for NERP NOR
#cn <- names(TCn)
##TCntemp <- TCn
#TCn$LENGTH_CLASS <- round(TCn$LENGTH_CLASS)
## Sum indivs in same length class
## (no need)
##(subset(TCn, YEAR==2004))
#TCn <- ddply(TCn, .(ID_MEDITS_TC, COUNTRY, AREA, VESSEL, YEAR, HAUL_NUMBER, CODEND_CLOSING, PARTIT, GENUS, SPECIES, CODLON, PFRAC,PECHAN,SEX,NBSEX,LENGTH_CLASS, MATURITY, UPLOAD_DATE, USER_ID, MATSUB, TF,  MONTH, DAY, CATFAU,MD5,CORRECT), summarise, NBLON = sum(NBLON))

#sort(unique(subset(TAn, YEAR == 1995 & COUNTRY=="ITA")$HAUL_NUMBER))
#sort(unique(subset(TCn, YEAR == 1995 & COUNTRY=="ITA")$HAUL_NUMBER))

# rename column names
colnames(TCn)[colnames(TCn)=="NBLON"] <- "NBLEN"
colnames(TCn)[colnames(TCn)=="PECHAN"] <- "WEIGHT_OF_THE_SAMPLE_MEASURED"
colnames(TCn)[colnames(TCn)=="PFRAC"] <-   "WEIGHT_OF_THE_FRACTION"
colnames(TCn)[colnames(TCn)=="CODLON"] <-   "LENGTH_CLASSES_CODE"
colnames(TCn)[colnames(TCn)=="MATSUB"] <-   "MATURITY_SUB_STAGING"

# Add other conversions (e.g. to numeric) here
# some factors to numbers
TCn$NBLEN<-as.numeric(as.character(TCn$NBLEN))
TCn$NBSEX<-as.numeric(as.character(TCn$NBSEX))
TCn$WEIGHT_OF_THE_FRACTION<-as.numeric(as.character(TCn$WEIGHT_OF_THE_FRACTION))
TCn$WEIGHT_OF_THE_SAMPLE_MEASURED<-as.numeric(as.character(TCn$WEIGHT_OF_THE_SAMPLE_MEASURED))
TCn$LENGTH_CLASS<-as.numeric(as.character(TCn$LENGTH_CLASS))
#TBn$NBTOT<-as.numeric(as.character(TBn$NBTOT))
#TBn$PTOT<-as.numeric(as.character(TBn$PTOT))

# Remove invalid hauls from TA 
TAn<-subset(TAn, VALIDITY=="V")

# Checks
# Check if there are duplicated rows in the files, if there are, running unique will remove duplicates, was a TEMP FIX FOR MERL in AREA 7
if (dim(unique(TAn))[1] != dim(TAn)[1]){
    stop("TA table contains duplicates\n")
}
#if (dim(unique(TBn))[1] != dim(TBn)[1]){
#    stop("TB table contains duplicates\n")
#}
if (dim(unique(TCn))[1] != dim(TCn)[1]){
    stop("TC table contains duplicates\n")
}
# Check Wing opening
if (!all(TAn$WING_OPENING > 0)){
    stop("WING_OPENING less than 0 in TA table\n")
}

# Are all hauls in TCn in TAn - if not join will not work properly

haul_check_fail <- rep(FALSE, length(unique(TCn$YEAR)))
names(haul_check_fail) <- unique(TCn$YEAR)
for (i in unique(TCn$YEAR)){
    if(!all(unique(subset(TCn, YEAR==i)$HAUL_NUMBER) %in% unique(subset(TAn, YEAR==i)$HAUL_NUMBER))){
        haul_check_fail[as.character(i)] <- TRUE
    }
}
if(any(haul_check_fail)){
    stop("Not all hauls in TC are in TA\n")
}

# Load MEDITS strata dataframe
medstrata<-fread(paste0(datadir,"MEDITS_Strata.csv"), sep=",", header=T)
names(medstrata)[3:5]<-c("NSTRATE", "CODEZONE","DEPTHSTRATA")
# Check NUMBER_OF_STRATA in TA has values in medstra only - else merge will fail
if(!all(TAn$NSTRATE %in% medstrata$NSTRATE)){
    stop("Not all strata numbers in TA are in the MEDITS strata table (some missing values or -1?).\n")
}

#--------------------------------------------------------------
# Correct Latitude and Longitude for mapping
#Estimate mid point in Haul and convert to Decimal Degrees for plotting
lat_start<-as.numeric(as.character(TAn$SHOOTING_LATITUDE))
lon_start<-as.numeric(as.character(TAn$SHOOTING_LONGITUDE))
LatStartSec = (lat_start - floor(lat_start))/100;
LonStartSec = (lon_start - floor(lon_start))/100;
LatStartDeg = floor(floor(lat_start)/100);
LonStartDeg = floor(floor(lon_start)/100);
LatStartMin = (floor(lat_start)/100 - LatStartDeg)/60*100;
LonStartMin = (floor(lon_start)/100 - LonStartDeg)/60*100;
lat_end<-as.numeric(as.character(TAn$HAULING_LATITUDE))
lon_end<-as.numeric(as.character(TAn$HAULING_LONGITUDE))
LatEndSec = (lat_end - floor(lat_end))/100;
LonEndSec = (lon_end - floor(lon_end))/100;
LatEndDeg = floor(floor(lat_end)/100);
LonEndDeg = floor(floor(lon_end)/100);
LatEndMin = (floor(lat_end)/100 - LatEndDeg)/60*100;
LonEndMin = (floor(lon_end)/100 - LonEndDeg)/60*100;
lat_start2 = LatStartDeg + LatStartMin + LatStartSec;
lon_start2 = LonStartDeg + LonStartMin + LonStartSec;
lat_end2 = LatEndDeg + LatEndMin + LatEndSec;
lon_end2 = LonEndDeg + LonEndMin + LonEndSec;
lat_start = lat_start2;
lon_start = lon_start2;
# use the quadrant to identify negative longitudes
lon_start<-ifelse(TAn$SHOOTING_QUADRANT==7, lon_start*-1, lon_start)
lat_end = lat_end2;
lon_end = lon_end2;
# use the quadrant to identify negative longitudes
lon_end<-ifelse(TAn$HAULING_QUADRANT==7, lon_end*-1, lon_end)
#FIXED MID HAUL POSITION
lat = (lat_start+lat_end)/2 
lon = (lon_start+lon_end)/2 
TAn$Latitude<-lat
TAn$Longitude<-lon

#----------------------------------------------------

# Add sanity check for DEPTH

# spp <-c("PAPE LON") in area 17
# 2005
# TAn has CODESTRATA E (no NSTRATE 21110)
# Wrong coding of STRATA or no hauls that year?


# Need to set CODESTATA to match with medstrata
# Assign the depth strata code (CODESTRATA) using mean depth values (DEPTH)
# i.e "1-50m","51-100m","101-200m","201-500m","501-800m"
#TAn$HAULING_DEPTH<-as.numeric(as.character(TAn$HAULING_DEPTH))
#TAn$SHOOTING_DEPTH<-as.numeric(as.character(TAn$SHOOTING_DEPTH))
#TAn$DEPTH<-(TAn$HAULING_DEPTH+TAn$SHOOTING_DEPTH)/2
#TAn$CODESTRATA<-'-'
#TAn[which(TAn$DEPTH>0   & TAn$DEPTH<=50), ]$CODESTRATA<-"A"
#TAn[which(TAn$DEPTH>50  & TAn$DEPTH<=100),]$CODESTRATA<-"B"
#TAn[which(TAn$DEPTH>100 & TAn$DEPTH<=200),]$CODESTRATA<-"C"
#TAn[which(TAn$DEPTH>200 & TAn$DEPTH<=500),]$CODESTRATA<-"D"
#TAn[which(TAn$DEPTH>500),                 ]$CODESTRATA<-"E"

# Actually no samples at that Depth


#----------------------------------------------------
# Adjust lengths based on unit
# What is LENGTH_CLASSES_CODE?
if(len.unit=="cm"){
  TCn$length <- ifelse(TCn$LENGTH_CLASSES_CODE=="M", TCn$LENGTH_CLASS/10, TCn$LENGTH_CLASS)
} else {
  TCn$length <- ifelse(TCn$LENGTH_CLASSES_CODE!="M", TCn$LENGTH_CLASS*10, TCn$LENGTH_CLASS)
}

# Standardise
# See page 56 of MEDITS manual 2012

# GSA is split into several zones (CODEZONE) in medstrata
# Each zone has 5 depths (the STRATA)
# Need to combine across zones, keeping strata separate
# Total area of each strata, summing across zones
strata_weight <- ddply(medstrata, .(AREA, CODESTRATA), summarise, TOTAL_STRATA_AREA = sum(AREASTRATA))

#temp <- medstrata[,sum(AREASTRATA), by=list(AREA,CODESTRATA)]

# Get weight of each strata in each area: Wk
strata_weight <- ddply(strata_weight, .(AREA), transform, STRATA_WEIGHT = TOTAL_STRATA_AREA / sum(TOTAL_STRATA_AREA))

# Get total swept area of hauls in each strata: Ak
# First join TA with MEDITS strata info 
# Joins by AREA, COUNTRY, NSTRATE
#TAn <- join(TAn, medstrata)
#temp <- join(TAn, medstrata) # mismatch between COUNTRY+NSTRATE in TAn and medstrata - TAn has SLOVENIA and CROATIA
#temp <- join(TAn, medstrata, by="NSTRATE") # Now we have two Country
TAn <- join(TAn, medstrata[,c("NSTRATE", "AREASTRATA", "CODESTRATA"), with=FALSE], by="NSTRATE") # Now we have two Country


#TAnd <- as.data.frame(TAn)
#medstratad <- as.data.frame(medstrata)
#temp <- merge(TAnd, medstratad)
# Now TAn has CODESTRATA
# Calc SWEPT AREA of each haul in TAn
# Calculate swept area, adjust units of measure to km 
# Wing opening is measured in decimeteres and distance in meters, then:
TAn$DISTANCE<-as.numeric(as.character(TAn$DISTANCE))
TAn$WING_OPENING<-as.numeric(as.character(TAn$WING_OPENING))
TAn$SWEPT<-(TAn$DISTANCE *((TAn$WING_OPENING)/10000))/1000 #Swept area in Km^2
# Get total swept area of hauls in each stratum
total_swept <- ddply(TAn, .(YEAR, AREA, CODESTRATA), summarise, TOTAL_AREA_SWEPT=sum(SWEPT))

# Add SPECIES
# Sum N across hauls for each strata and length
# Join GENUS and SPECIES columns
TCn$SP <- paste(TCn$GENUS, TCn$SPECIES)
# Raise the measured sample to the whole catch in the haul
TCn$NBLEN.sample.raised <- TCn$NBLEN*TCn$WEIGHT_OF_THE_FRACTION/TCn$WEIGHT_OF_THE_SAMPLE_MEASURED
# Need CODESTRATA in TC table
# Join TC with TA to bring CODESTRATA into TC
# Joins by COUNTRY, AREA, YEAR, HAUL_NUMBER
TCn <- join(TCn, TAn[,c("YEAR", "COUNTRY", "AREA", "HAUL_NUMBER", "CODESTRATA"), with=FALSE])

# This should be a sanity check at the beginning
# Join is bad as not all hauls in TC are in TA - this is a data problem
# Anyway of determining Haul info (codestrata) from TCn?
temp <- join(TCn, TAn[,c("YEAR", "COUNTRY", "AREA", "HAUL_NUMBER", "CODESTRATA"), with=FALSE])
unique(TAn$HAUL_NUMBER) # No Haul 83
all(unique(TCn$HAUL_NUMBER) %in% unique(TAn$HAUL_NUMBER))
unique(TCn$HAUL_NUMBER)[which(!(unique(TCn$HAUL_NUMBER) %in% unique(TAn$HAUL_NUMBER)))]
subset(TAn, HAUL_NUMBER==83)
subset(TCn, HAUL_NUMBER==83)
# Hence join is bad - NA in CODESTRATA
# What do we do in this case?
# Flag it up - can we correct? Depends on how many hauls are missing?
# Assign them to another strata?


# TCn now has CODESTRATA

# Calculate the mean N per area (by length and stratum)
mean_n <- ddply(TCn, .(YEAR, AREA, CODESTRATA, SP, LENGTH_CLASS), summarise, TOTAL_N = sum(NBLEN.sample.raised))
# Bring in the total swept area in each stratum
mean_n <- join(mean_n, total_swept)
mean_n$MEAN_N <- mean_n$TOTAL_N / mean_n$TOTAL_AREA_SWEPT

# Multiply the mean n by the weight of the strata by the 
mean_n <- join(mean_n, strata_weight)
mean_n$WEIGHTED_MEAN_N <- mean_n$MEAN_N * mean_n$STRATA_WEIGHT

# Sum over strata
# By SEX?
index <- ddply(mean_n, .(AREA, YEAR, SP, LENGTH_CLASS), summarise, data = sum(WEIGHTED_MEAN_N))





#mb11 <- subset(index, SP == "MULL BAR" & AREA == 11)



ggplot(index, aes(x=LENGTH_CLASS, y=data)) +geom_bar(stat="identity") + facet_wrap(~YEAR)






