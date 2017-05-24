# raising_MEDITS.R
# Copyright European Union, 2016
# Author: Finlay Scott (EC JRC) <finlay.scott@ec.europa.eu>
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# Description: Script that produces standardised length-frequency distribution from Medits data.

# Based on: LFDplyr.R by Matteo MURENU, Alessandro MANNINI, Tristan ROUYER, Chato OSIO and Finlay SCOTT

# Libraries and what not
library(ggplot2)
library(data.table)
library(plyr) # Needed for join

#' Get and raise the MEDITS data for a particular species and GSA
#'
#' Uses the MEDTS TA, TC and strata tables to raise the survey data following the MEDITS methodology.
#' Requires 'medits_strata.csv' file to be in the data directory.
#' If required the output is saved in a directory
#'
#' @param genus Genus in capitals, e.g. "PAPE"
#' @param species Species in captials, e.g. "LON"
#' @param gsa Number of the GSA or GSAs, e.g. c(17, 18)
#' @param sex_split So you want the data split by sex, TRUE or FALSE
#' @param len_unit The returned unit of the length class
#' @param datadir Directory of the input data tables
#' @param outdir Directory of the produced outputs
#' @param plotdir Directory of the output plots
#' @param save_output Do you want to save the output from each species / GSA, TRUE or FALSE
#' @param plots Do you want plots or not, TRUE or FALSE
#' @param ta_name The name of the TA CSV file
#' @param tc_name The name of the TC CSV file
#' @param sep The seperator of input CSV files, e.g. ","
#' @example
#' # Needs the data to run
#' test <- raise_medits("PAPE", "LON", c(17,18), sex_split=FALSE, len_unit="cm", datadir="../../../Surveys_data/Medits/", outdir="../../../medits_out", plotdir="../../../medits_out", save_output=TRUE, plots=TRUE, ta_name="medits_ta.csv", tc_name="medits_tc.csv", sep=",")

# All species in all GSAs are returned
raise_medits <- function(genus, species, gsa, sex_split=FALSE, len_unit="cm", datadir, outdir, plotdir, plots=FALSE, ta_name="medits_ta.csv", tc_name="medits_tc.csv", sep=","){
    if(!(len_unit %in% c("mm","cm"))){
        stop("len_unit must be mm or cm")
    }
    # Data loading and editing
    # Check the filenames of the tables and the seperator etc
    TAn <- fread(paste0(datadir, ta_name), sep=sep, header=T)
    names(TAn)<-toupper(names(TAn))
    TAn <- subset(TAn, AREA %in% gsa)
    # Remove invalid hauls from TA 
    TAn<-subset(TAn, VALIDITY=="V")

    TCn <- fread(paste0(datadir, tc_name), sep=sep, header=T)
    names(TCn)<-toupper(names(TCn))
    TCn <- subset(TCn, (AREA %in% gsa) & (GENUS %in% genus) & (SPECIES %in% species))
    # Rename column names so we know what is going on
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
    TCn$LENGTH_CLASSES_CODE <-as.character(TCn$LENGTH_CLASSES_CODE)
    # Join GENUS and SPECIES columns
    TCn$SP <- paste(TCn$GENUS, TCn$SPECIES, sep=" ")

    # Checks - do these go here? - can flag up which species if put in loop below
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
    medstrata<-fread(paste0(datadir,"medits_strata.csv"), sep=";", header=T)
    names(medstrata)[3:5]<-c("NSTRATE", "CODEZONE","DEPTHSTRATA")

    # The MEDITS methodology
    # See page 56 of MEDITS manual 2012
    # Pull MEDITS strata data into TA (the columns AREASTRATA and CODESTRATA)
    # Joins by NSTRATE (which is unique in TAn)
    # Check before join that all of the NSTRATE values in TA also exist in medstra - else merge will fail
    if(!all(TAn$NSTRATE %in% medstrata$NSTRATE)){
        # Farm out to other routine
        warning("Not all strata numbers in TA are in the MEDITS strata table (some missing values or -1?). Attempting to join by DEPTH and GSA\n")
        # NSTRATE is a unique identifier - but if NSTRATE is empty in TA then we could try to join by DEPTH and GSA
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
        # Join by CODESTRATA and AREA, pulling in AREASTRATA from medstrata
        pre_dim <- dim(TAn)
        TAn <- join(TAn, medstrata[,c("AREASTRATA", "AREA", "CODESTRATA"), with=FALSE], by=c("AREA", "CODESTRATA"))
        if(dim(TAn)[1] != pre_dim[1]){
            stop("Merging medstrata and TAn has gone wrong. Number of rows of TAn have changed\n")
        }
    } else {
        # Joining by NSTRATE
        TAn <- join(TAn, medstrata[,c("NSTRATE", "AREASTRATA", "CODESTRATA"), with=FALSE], by="NSTRATE")
        # Now TAn has CODESTRATA and AREASTRATA
    }
    # Check that we have AREASTRATA for each entry
    if(any(is.na(TAn$AREASTRATA))){
        stop("Some NAs in AREASTRATA column of TAn\n")
    }

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

    # Perform raising for each GSA and species combination we have
    # Process each wanted species and GSA separately
    wanted <- expand.grid(sp = paste(genus, species, sep=" "), gsa = gsa, stringsAsFactors = FALSE)
    index_out <- data.frame()
    for (wanted_count in 1:nrow(wanted)){
        cat("Processing", wanted$sp[wanted_count], "in GSA" , wanted$gsa[wanted_count],"\n")
        TCsub <- subset(TCn, AREA==wanted$gsa[wanted_count] & SP==wanted$sp[wanted_count])
        # CHECK UNITS
        # Adjust lengths based on unit - check this
        # LENGTH_CLASSES_CODE: m for mm, 0 for 0.5 cm (??), 1 for cm
        TCsub$length <- TCsub$LENGTH_CLASS
        if(len_unit=="cm"){
            # Scale down from mm to cm
            TCsub[LENGTH_CLASSES_CODE=="m", length:= length/10]
            # Scale down from 0.5 cm to cm
            TCsub[LENGTH_CLASSES_CODE=="0", length:= length/2]
        }
        if(len_unit=="mm"){
            TCsub[LENGTH_CLASSES_CODE=="1", length:= length * 10]
            TCsub[LENGTH_CLASSES_CODE=="0", length:= length * 5]
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
        # Sum to create the index and make file name
        if (sex_split==TRUE){
            species_file_name <- paste(wanted$sp[wanted_count],wanted$gsa[wanted_count], "SEX",sep="_")
            index <- ddply(mean_n, .(AREA, YEAR, GENUS, SPECIES, SEX, LENGTH_CLASS), summarise, data = sum(WEIGHTED_MEAN_N))
        } else {
            species_file_name <- paste(wanted$sp[wanted_count],wanted$gsa[wanted_count], sep="_")
            index <- ddply(mean_n, .(AREA, YEAR, GENUS, SPECIES, LENGTH_CLASS), summarise, data = sum(WEIGHTED_MEAN_N))
        }
        species_file_name <- gsub(" ","_",species_file_name)

        # Pad out index so that each year has the same length classes - not sure it's a good idea but comparable with previous method
        unique_lengths <- unique(index$LENGTH_CLASS)
        unique_years <- unique(index$YEAR)
        if(!sex_split){
            padded_index <- expand.grid(LENGTH_CLASS = unique_lengths, YEAR = unique_years)
        } else {
            unique_sex <- unique(index$SEX)
            padded_index <- expand.grid(LENGTH_CLASS = unique_lengths, YEAR = unique_years, SEX=unique_sex)
        }
        padded_index$GENUS <- index$GENUS[1]
        padded_index$SPECIES <- index$SPECIES[1]
        padded_index$AREA <- index$AREA[1]
        # Bring in index into the padded_index
        padded_index <- join(padded_index, index)
        padded_index <- padded_index[order(padded_index$YEAR, padded_index$LENGTH_CLASS),] 
        padded_index[is.na(padded_index$data),"data"] <- 0 # Probably a bad idea
        index <- padded_index
    
        # Save output
        if (save_output == TRUE){
            save(index,file=paste0(outdir,"/stratified_Nlen_",species_file_name,".Rdata"))
        }
        index_out <- rbind(index_out, index)

        # Generate plots
        if (plots == TRUE){
            # For boxplot data
            tmp <- as.data.frame(TCn)
            box_data <-tmp[rep(row.names(tmp),round(tmp$NBLEN.sample.raised)),]  # duplicate rows
            if  (len_unit=="cm") {
                p <- ggplot(index, aes(y=data, x=LENGTH_CLASS/10))
                p2 <- ggplot(box_data, aes(y=LENGTH_CLASS/10, x=factor(YEAR)))
            } else {
                p <- ggplot(index, aes(y=data, x=LENGTH_CLASS))
                p2 <- ggplot(box_data, aes(y=LENGTH_CLASS, x=factor(YEAR)))
            } 
            # set the graph title
            if  (sex_split==FALSE) {
                tt <- ggtitle(paste0(wanted$sp, ": GSA",wanted$gsa[wanted_count],"\n"))
            } else {
                tt <- ggtitle(paste0(wanted$sp[wanted_count], ", sex"," : GSA",wanted$gsa[wanted_count],"\n"))
            } 
            # plot the histogram
            p <- p + tt + geom_bar(stat= "identity") +
                labs(x=paste0("\nLength (",len_unit,")"),y=expression("n/km"^"2")) +
                theme(legend.position="none")
            if  (sex_split==FALSE) {
                p <- p + facet_wrap(~YEAR)
            } else{
                p <- p + facet_grid(YEAR~SEX)
            }
            print(p)
            ggsave(p, file=paste0(plotdir,"/stratified_Nlen_",species_file_name,".png"))
            
            # plot the boxplot
            p2 <- p2 + tt + geom_boxplot() +
                labs(y=paste0("\nLength (", len_unit,")\n"),x=NULL) +
                guides(fill=FALSE) + theme_bw() +
                theme(plot.title = element_text(colour = "red", size=rel(1.2))) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            if  (sex_split==TRUE) {
                p2 <- p2 + facet_wrap(~SEX)
            }
            print(p2)
            ggsave(p2, file=paste0(plotdir,"/bxplen_",species_file_name,".png")) 
        }
    }
    return(index_out)
}


