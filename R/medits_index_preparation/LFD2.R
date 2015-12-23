###========================================================================================####
#    Title:   LFDplyr.R
#    Release: 0.1
#                          
#    Description: Script that produces standardised length-frequency distribution from Medits data.
#    Authors:  Matteo MURENU, Alessandro MANNINI, Tristan ROUYER, Chato OSIO, Finlay SCOTT
#    Date: Created December 2015, ROMA, EWG 15-06 meeting
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

# Libraries and what not

#library(reshape2)
library(ggplot2)
library(data.table)
library(dplyr) # Not really being used - could be for speed
library(plyr) # Needed for join

#---------------------------------------
# Parameters - could be inputs for function

# Set these for your computer
# Location of TA, TB, TC and Medits Strata tables
datadir <- "../../../tables/medits/"

# Where do you want the results to go
outdir <- "../../../medits_lfd_test/data"
plotdir <- "../../../medits_lfd_test/figures"

# Have these in a data frame?
# Might be easier

#wanted <- data.frame(genus = c("MERL", "MERL", "MULL", "PAPE"),
#                     species = c("MER", "MER", "BAR", "LON"),
#                     gsa = c(9, 11, 9, 17),
#                     sex_split = FALSE,
#                     len_unit = c(rep("cm",3), "mm"))


wanted <- data.frame(genus = c("MERL", "MERL", "MULL"),
                     species = c("MER", "MER", "BAR"),
                     gsa = c(9, 11, 9),
                     sex_split = c(FALSE, TRUE,FALSE),
                     len_unit = c(rep("cm",3)))


plots <- TRUE

#---------------------------------------------------
# Data loading and editing

# Load raw data (MEDITS dataset tables) from csv and subset the GSAs and species we want
# Check the filenames of the tables and the serator (here ;)
TAn <- fread(paste0(datadir, "ta", ".csv"), sep=";", header=T)
TAn <- subset(TAn, area %in% unique(wanted$gsa))
TCn <- fread(paste0(datadir, "tc", ".csv"), sep=";", header=T)
# Insufficient subsetting as it may include unwanted rows - but enough to trim TCn down a bit
TCn <- subset(TCn, area %in% unique(wanted$gsa) & genus %in% unique(wanted$genus) & species %in% unique(wanted$species))

# For consistency
names(TAn)<-toupper(names(TAn))
names(TCn)<-toupper(names(TCn))

# Remove invalid hauls from TA 
TAn<-subset(TAn, VALIDITY=="V")

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
# Join GENUS and SPECIES columns
TCn$SP <- paste(TCn$GENUS, TCn$SPECIES)


# Do these go here - can flag up which species if put in loop below
# Checks
# Check if there are duplicated rows in the files
if (dim(unique(TAn))[1] != dim(TAn)[1]){
    stop("TA table contains duplicates\n")
}
if (dim(unique(TCn))[1] != dim(TCn)[1]){
    stop("TC table contains duplicates\n")
}
# Check Wing opening
if (!all(TAn$WING_OPENING > 0)){
    stop("WING_OPENING less than 0 in TA table\n")
}

# Load MEDITS strata dataframe
medstrata<-fread(paste0(datadir,"MEDITS_Strata.csv"), sep=",", header=T)
names(medstrata)[3:5]<-c("NSTRATE", "CODEZONE","DEPTHSTRATA")

#-------------------------------------------------------------
# The MEDITS methodology
# See page 56 of MEDITS manual 2012

# Pull MEDITS strata into TA 
# Joins by NSTRATE (which is unique in TAn)
# Check before join that all of the NSTRATE values in TA also exist in medstra - else merge will fail
if(!all(TAn$NSTRATE %in% medstrata$NSTRATE)){
    stop("Not all strata numbers in TA are in the MEDITS strata table (some missing values or -1?).\n")
}
TAn <- join(TAn, medstrata[,c("NSTRATE", "AREASTRATA", "CODESTRATA"), with=FALSE], by="NSTRATE")
# Now TAn has CODESTRATA and AREASTRATA

# Calc SWEPT AREA of each haul in TAn and adjust units of measure to km 
# Wing opening is measured in decimeteres and distance in meters, then:
TAn$DISTANCE<-as.numeric(as.character(TAn$DISTANCE))
TAn$WING_OPENING<-as.numeric(as.character(TAn$WING_OPENING))
TAn$SWEPT<-(TAn$DISTANCE *((TAn$WING_OPENING)/10000))/1000 #Swept area in Km^2
# Get total swept area of hauls in each stratum: Ak
total_swept <- ddply(TAn, .(YEAR, AREA, CODESTRATA), summarise, TOTAL_AREA_SWEPT=sum(SWEPT))

# Raise the measured sample to the whole catch in the haul
TCn$NBLEN.sample.raised <- TCn$NBLEN*TCn$WEIGHT_OF_THE_FRACTION/TCn$WEIGHT_OF_THE_SAMPLE_MEASURED

# Get weight of each strata in each area: Wk
# Each GSA is split into five depth strata and several zones
# Need to combine areas across zones, keeping strata separate
# Total area of each strata, summing across zones
strata_weight <- ddply(medstrata, .(AREA, CODESTRATA), summarise, TOTAL_STRATA_AREA = sum(AREASTRATA))
# Finally get weight of each strata in each area: Wk
strata_weight <- ddply(strata_weight, .(AREA), transform, STRATA_WEIGHT = TOTAL_STRATA_AREA / sum(TOTAL_STRATA_AREA))

# Now process wanted row by row
for (wanted_count in 1:nrow(wanted)){
    cat("Processing", as.character(wanted$genus[wanted_count]), as.character(wanted$species[wanted_count]), "in gsa" , wanted$gsa[wanted_count],"\n")
    TCsub <- subset(TCn, AREA==wanted$gsa[wanted_count] & SPECIES==wanted$species[wanted_count] & GENUS==wanted$genus[wanted_count])
    # Adjust lengths based on unit
    if(wanted$len_unit[wanted_count]=="cm"){
        TCsub$length <- ifelse(TCsub$LENGTH_CLASSES_CODE=="M", TCsub$LENGTH_CLASS/10, TCsub$LENGTH_CLASS)
    } else {
        TCsub$length <- ifelse(TCsub$LENGTH_CLASSES_CODE!="M", TCsub$LENGTH_CLASS*10, TCsub$LENGTH_CLASS)
    }
    # Check all hauls in TCsub are in TAn - else join will fail
    for (yr in unique(TCsub$YEAR)){
        if (!all(unique(subset(TCsub, YEAR==yr)$HAUL_NUMBER) %in% unique(subset(TAn, AREA==wanted$gsa[wanted_count] & YEAR==yr)$HAUL_NUMBER))){
            warning(paste("Not all hauls in TC are in TA for", as.character(wanted$genus[wanted_count]), as.character(wanted$species[wanted_count]), "in gsa" , wanted$gsa[wanted_count], "in year", yr, "\n"))
        }
    }
    # Want to sum N across hauls for each strata and length
    # Join TC with TA to bring CODESTRATA into TC
    # Joins by COUNTRY, AREA, YEAR, HAUL_NUMBER
    TCsub <- join(TCsub, TAn[,c("YEAR", "COUNTRY", "AREA", "HAUL_NUMBER", "CODESTRATA"), with=FALSE])
    # TCn now has CODESTRATA
    # Calculate the mean N per area (by length and stratum)
    # First get the total numbers in each stratum by length
    mean_n <- ddply(TCsub, .(YEAR, AREA, CODESTRATA, GENUS, SPECIES, SEX, LENGTH_CLASS), summarise, TOTAL_N = sum(NBLEN.sample.raised))
    # Bring in the total swept area in each stratum - joins by YEAR, AREA and CODESTRATA
    mean_n <- join(mean_n, total_swept)
    # Calc mean N over the area
    mean_n$MEAN_N <- mean_n$TOTAL_N / mean_n$TOTAL_AREA_SWEPT
    # Multiply the mean n by the weight of the strata
    mean_n <- join(mean_n, strata_weight)
    mean_n$WEIGHTED_MEAN_N <- mean_n$MEAN_N * mean_n$STRATA_WEIGHT
    # Sum to create the index
    if (wanted$sex_split[wanted_count]==TRUE){
        species_file_name <- paste(wanted$genus[wanted_count], wanted$species[wanted_count],wanted$gsa[wanted_count], "SEX",sep="_")
        index <- ddply(mean_n, .(AREA, YEAR, GENUS, SPECIES, SEX, LENGTH_CLASS), summarise, data = sum(WEIGHTED_MEAN_N))
    }
    else {
        species_file_name <- paste(wanted$genus[wanted_count], wanted$species[wanted_count],wanted$gsa[wanted_count], sep="_")
        index <- ddply(mean_n, .(AREA, YEAR, GENUS, SPECIES, LENGTH_CLASS), summarise, data = sum(WEIGHTED_MEAN_N))
    }
    
    # Save output
    save(index,file=paste0(outdir,"/stratified_Nlen_",species_file_name,".Rdata"))

    
    # Generate plots
    if (plots == TRUE){
        if  (wanted$len_unit[wanted_count]=="cm") {
            p <- ggplot(index, aes(y=data, x=LENGTH_CLASS/10))
            p2 <- ggplot(index, aes(y=LENGTH_CLASS/10, x=factor(YEAR)))
        }
        else {
            p <- ggplot(index, aes(y=data, x=LENGTH_CLASS))
            p2 <- ggplot(index, aes(y=LENGTH_CLASS, x=factor(YEAR)))
        } 
        # set the graph title
        if  (wanted$sex_split[wanted_count]==FALSE) {
            tt <- ggtitle(paste0(wanted$genus[wanted_count], " ", wanted$species[wanted_count],": GSA",wanted$gsa[wanted_count],"\n"))
        }
        else {
            tt <- ggtitle(paste0(wanted$genus[wanted_count], " ", wanted$species[wanted_count], ", sex"," : GSA",wanted$gsa[wanted_count],"\n"))
        } 
        # plot the histogram
        p <- p + tt + geom_bar(stat= "identity") +
            labs(x=paste0("\nLength (",wanted$len_unit[wanted_count],")"),y=expression("n/km"^"2")) +
            theme(legend.position="none")
        if  (wanted$sex_split[wanted_count]==FALSE) {
            p <- p + facet_wrap(~YEAR)
        }
        else{
            p <- p + facet_grid(YEAR~SEX)
        }
        print(p)
        ggsave(p, file=paste0(plotdir,"/stratified_Nlen_",species_file_name,".png"))
        
        # plot the boxplot
        p2 <- p2 + tt + geom_boxplot() +
            labs(y=paste0("\nLength (", wanted$len_unit[wanted_count],")\n"),x=NULL) +
            guides(fill=FALSE) + theme_bw() +
            theme(plot.title = element_text(colour = "red", size=rel(1.2))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        if  (wanted$sex_split[wanted_count]==TRUE) {
            p2 <- p2 + facet_wrap(~SEX)
        }
        print(p2)
        ggsave(p2, file=paste0(plotdir,"/bxplen_",species_file_name,".png")) 
    }
}

#--------------------------------------------------------------


# Include?
# Correct Latitude and Longitude for mapping
# Estimate mid point in Haul and convert to Decimal Degrees for plotting
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














ggplot(index, aes(x=LENGTH_CLASS, y=data)) +geom_bar(stat="identity") + facet_wrap(~YEAR)






