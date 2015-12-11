# R code update July 2014 STECF EWG 14-09 

# Authors Chato Osio, Finlay Scott and Tristan Rouyer
#################################################################################

# R code update 12/12/2013 STECF EWG 13-19 

# Authors Chato Osio, Finlay Scott and Tristan Rouyer


###########################################################################################################################
#                   R code developed for SGMED assessment routines on MEDITS data                                     #
#                 Authors: Valerio Bartolino, Chato Osio, Graham Pilling and Finlay Scott                                 #
#                         March 2010                                                                                      #
#                         Minor revisions September 2011                                                                  #one line to give the program's name and a brief idea of what it does.>
#    Copyright (C) <2012>  <Valerio Bartolino, Chato Osio, Graham Pilling and Finlay Scott>

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
###################################
# TO DO LIST

###################################
# connection and data upload
require(RODBC) # odbcConnectAccess
require(doBy)  # summaryBy

spsplit<-strsplit(sp," ")[[1]]
genus<-spsplit[1]
species<-spsplit[[2]]

# queries
# to get only table TA
#qry0 <- paste("SELECT * FROM medits_ta")

qry1 <- paste("SELECT * FROM medits_ta WHERE area='", gsa,"'",sep="")

qry2 <- paste("SELECT * FROM medits_tb WHERE genus='", genus,"' AND species='", species, "' AND area='",   gsa,"'",sep="")

qry3 <- paste("SELECT * FROM medits_tc WHERE genus='", genus,"' AND species='", species, "' AND area='", gsa,"'",sep="" )

# Need to specify the path and the database name of the ACCESS Database, make sure you point to a .mdb file, will not work if it is .accdb

ch <- odbcConnectAccess("C:/Users/osioogi/Desktop/med_data_access/Surveys.mdb")

#sqlTables(ch)

#temp <- sqlQuery(ch, qry0) # run if you want to download all the TA table, not only for one area

TAn <- sqlQuery(ch, qry1)

TBn <- sqlQuery(ch, qry2)

TCn <- sqlQuery(ch, qry3)


odbcClose(ch)

###################################

#merge GENUS and SPECIE
TBn$genspe<-paste(TBn$genus, TBn$species,sep = " ")
TCn$genspe<-paste(TCn$genus, TCn$species,sep = " ")

# In 2014 database columns have changed names, so need to revert to original file names, convert to upper case here and later rename some columns

names(TAn)<-toupper(names(TAn))
names(TBn)<-toupper(names(TBn))
names(TCn)<-toupper(names(TCn))

# remove invalid hauls from TA right away
#TAn<-subset(TAn, VALIDITY=="V")

# For medit 2014 in gsa there are intercalibration tows, need to remove from dataset before producing index
table(TAn$YEAR)



## Find mean haul dept
TAn$DEPTH<-(TAn$HAULING_DEPTH+TAn$SHOOTING_DEPTH)/2

# Check if there are duplicated rows in the files, if there are, running unique will remove duplicates, was a TEMP FIX FOR MERL in AREA 7

#table(duplicated(TBn))
#table(duplicated(TCn))
#TBn<-unique(TBn)
#TCn<-unique(TCn)

#Estimate mid point in Haul and convert to Decimal Degrees for plotting

lat_start<-TAn$SHOOTING_LATITUDE
lon_start<-TAn$SHOOTING_LONGITUDE

LatStartSec = (lat_start - floor(lat_start))/100;
LonStartSec = (lon_start - floor(lon_start))/100;
LatStartDeg = floor(floor(lat_start)/100);
LonStartDeg = floor(floor(lon_start)/100);
LatStartMin = (floor(lat_start)/100 - LatStartDeg)/60*100;
LonStartMin = (floor(lon_start)/100 - LonStartDeg)/60*100;
lat_end<-TAn$HAULING_LATITUDE
lon_end<-TAn$HAULING_LONGITUDE
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
################################3
# drop WING_OPENING=0 in TAn
#TAn <- TAn[TAn$WING_OPENING>0,]  # This is a hack to avoid script crash, there should be no 0 wing opening with valid tows

#################################
# MERGING TAnBLES
colnames(TAn)[colnames(TAn)=="NSTRATE"] <- "NUMBER_OF_THE_STRATUM"


#tempTBn <- merge(TAn, TBn, by=c("YEAR","VESSEL","AREA" ,"HAUL_NUMBER"), all.x=TRUE)
TBn <- merge(TAn, TBn, by=c("YEAR","VESSEL","AREA" ,"HAUL_NUMBER", "COUNTRY"), all.x=TRUE) # UPDATE 15/7/2014 now also with by= "MONTH", "DAY"

TCn <- merge(TAn, TCn, by=c("YEAR","VESSEL","AREA" ,"HAUL_NUMBER", "COUNTRY"), all.x=TRUE)

# change too long column names (why not keeping MEDITS labels?)
colnames(TBn)[colnames(TBn)=="NBTOT"] <- "TOTAL_NUMBER_IN_HAUL"
colnames(TBn)[colnames(TBn)=="PTOT"] <- "TOTAL_WEIGHT_IN_HAUL"

# pruning columns, keep only what we need and overwrigth the object
TBn <- TBn[,c("COUNTRY","VESSEL","YEAR","HAUL_NUMBER","MONTH.x","DAY.x","SHOOTING_TIME","HAULING_DURATION","VALIDITY","DISTANCE","VERTICAL_OPENING","WING_OPENING", "GEOMETRICAL_PRECISION","BRIDLES_LENGTH","WARP_LENGTH","WARP_DIAMETER","AREA","TOTAL_WEIGHT_IN_HAUL","TOTAL_NUMBER_IN_HAUL","GENSPE", "Latitude", "Longitude", "DEPTH", "NUMBER_OF_THE_STRATUM","UPLOAD_DATE.x" )]


#colnames(TCn)[colnames(TCn)=="NO_OF_INDIVIDUAL_OF_THE_ABOVE_SEX_MEASURED"] <- "NBSEX"
colnames(TCn)[colnames(TCn)=="NBLON"] <- "NBLEN"
#colnames(TCn)[colnames(TCn)=="NSTRATE"] <- "NUMBER_OF_THE_STRATUM"
colnames(TCn)[colnames(TCn)=="PECHAN"] <- "WEIGHT_OF_THE_SAMPLE_MEASURED"
colnames(TCn)[colnames(TCn)=="PFRAC"] <-   "WEIGHT_OF_THE_FRACTION"
colnames(TCn)[colnames(TCn)=="CODLON"] <-   "LENGTH_CLASSES_CODE"
colnames(TCn)[colnames(TCn)=="MATSUB"] <-   "MATURITY_SUB_STAGING"
#colnames(TCn)[colnames(TCn)=="CODLON"] <-   "LENGTH_CLASSES_CODE"


TCn <- TCn[,c("COUNTRY","VESSEL","YEAR","HAUL_NUMBER","MONTH.x","DAY.x","SHOOTING_TIME","HAULING_DURATION","VALIDITY","DISTANCE","VERTICAL_OPENING","WING_OPENING","GEOMETRICAL_PRECISION","BRIDLES_LENGTH","WARP_LENGTH","WARP_DIAMETER","AREA","LENGTH_CLASSES_CODE","WEIGHT_OF_THE_FRACTION","WEIGHT_OF_THE_SAMPLE_MEASURED","SEX","NBSEX","LENGTH_CLASS","MATURITY","NBLEN","GENSPE", "Latitude", "Longitude", "DEPTH", "MATURITY_SUB_STAGING", "NUMBER_OF_THE_STRATUM","UPLOAD_DATE.x" )]  

# distinguish between real "zero-hauls" and "no-measure-hauls"
# comparing TBn and TCn


TBn.0 <- TBn[which(is.na(TBn$GENSPE)==TRUE),]
TCn.0 <- TCn[which(is.na(TCn$GENSPE)==TRUE),]

TBn.0$IDX <- paste(TBn.0$YEAR,TBn.0$HAUL_NUMBER, sep="_")
TCn.0$IDX <- paste(TCn.0$YEAR,TCn.0$HAUL_NUMBER, sep="_")

truezero.haul <- which(is.na(match(TCn.0$IDX,TBn.0$IDX))==FALSE)
truezero.haul <- TCn.0[truezero.haul,]

TCn <- TCn[is.na(TCn$GENSPE)==FALSE,] #exclude truezero and no-measure-hauls


# true "zero-hauls" and "no-measure-hauls" selected fields from NA to 0
TBn[is.na(TBn$GENSPE)==TRUE,"TOTAL_WEIGHT_IN_HAUL"] <- 0
TBn[is.na(TBn$GENSPE)==TRUE,"TOTAL_NUMBER_IN_HAUL"] <- 0
TBn$GENSPE <- sp

TCn[is.na(TCn$GENSPE)==TRUE,"WEIGHT_OF_THE_FRACTION"] <- 0
TCn[is.na(TCn$GENSPE)==TRUE,"WEIGHT_OF_THE_SAMPLE_MEASURED"] <- 0
TCn[is.na(TCn$GENSPE)==TRUE,"NBSEX"] <- 0
TCn[is.na(TCn$GENSPE)==TRUE,"NBLEN"] <- 0
TCn$GENSPE <- sp


##########################################
# uniform the length classes to the same unit (mm or cm)
# WARNING: do not run until the field LENGTH_CLASSES_CODE has not been checked
if(len.unit=="cm"){

      TCn$length <- ifelse(TCn$LENGTH_CLASSES_CODE=="m", TCn$LENGTH_CLASS/10, TCn$LENGTH_CLASS)
    } else {
      
      TCn$length <- ifelse(TCn$LENGTH_CLASSES_CODE!="m", TCn$LENGTH_CLASS*10, TCn$LENGTH_CLASS)
    }
###########################################


# drop no standard coding - several no standard codes were found (i.e. @, 0, 2, 4, 5, D, f, i, m, etc), NOT needed anymore since these values not allowed in database 
#TCn <- TCn[TCn$SEX=="M" | TCn$SEX=="F" | TCn$SEX=="I" | TCn$SEX=="N",]


# assign indetermined to M and F (assumed to be 0.5 M and 0.5 F)
TCn[TCn$SEX=="I" | TCn$SEX=="N","NBLEN"] <- TCn[TCn$SEX=="I" | TCn$SEX=="N","NBLEN"]/2
TCn[TCn$SEX=="I" | TCn$SEX=="N","NBSEX"] <- TCn[TCn$SEX=="I" | TCn$SEX=="N","NBSEX"]/2

tmp <- TCn[TCn$SEX=="I" | TCn$SEX=="N",]

TCn[TCn$SEX=="I" | TCn$SEX=="N","SEX"] <- "M"
TCn <- rbind(TCn, tmp)
TCn[TCn$SEX=="I" | TCn$SEX=="N","SEX"] <- "F"


# raising the measured sample the whole catch
TCn$NBLEN.raised <- TCn$NBLEN*TCn$WEIGHT_OF_THE_FRACTION/TCn$WEIGHT_OF_THE_SAMPLE_MEASURED


# aggregate the new assigned M and F (previously I) to the original M and F
TCn2 <- summaryBy(NBLEN.raised~LENGTH_CLASS+SEX+HAUL_NUMBER+YEAR, data=TCn, id=~COUNTRY+VESSEL+MONTH.x+DAY.x+SHOOTING_TIME+HAULING_DURATION+VALIDITY+DISTANCE+
VERTICAL_OPENING+WING_OPENING+AREA+GENSPE+DEPTH+length, FUN=sum, keep.names=TRUE)


# drop wing openings=0, this was a hack to remove Spanish medits data from 1994
#TBn$WING_OPENING<-subset(TBn, WING_OPENING>0)

#calculate swept area, adjust units of measure to km, wing opening is measured in decimeteres, and distance in meters
TBn$swept<-(TBn$DISTANCE *((TBn$WING_OPENING)/10000))/1000 #Swept area in Km^2
TCn$swept<-(TCn$DISTANCE *((TCn$WING_OPENING)/10000))/1000 #Swept area in Km^2


# change too long column names (why not keeping MEDITS labels?)
colnames(TBn)[colnames(TBn)=="TOTAL_WEIGHT_IN_HAUL"] <- "PTot"
colnames(TBn)[colnames(TBn)=="TOTAL_NUMBER_IN_HAUL"] <- "NBTot"

#Keep colnames that make sense and are useful
#TBn<-subset(TBn, select= -c(CODEND_CLOSING.x, UPLOAD_DATA.y, UPLOAD_DATA.x ,COUNTRY_ACCOUNT.x,COUNTRY_ACCOUNT.y,  CODEND_CLOSING.y))

# colnames(TBn)[colnames(TBn)=="COUNTRY_INFLIE"] <- "COUNTRY"
# colnames(TCn)[colnames(TCn)=="COUNTRY_INFLIE"] <- "COUNTRY"

#REMOVE INVALID HAULS
TAn<-subset(TAn, VALIDITY=="V")
TBn<-subset(TBn, VALIDITY=="V")
TCn<-subset(TCn, VALIDITY=="V")

TA <- TAn
TB <- TBn
TC <- TCn