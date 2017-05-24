# get_landings_length_funcs.R
# Copyright European Union, 2016
# Author: Finlay Scott (EC JRC) <finlay.scott@jrc.ec.europa.eu>
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

#' Get the length based landings data for a species and GSA
#'
#' Get the length based landings data for a species and GSA
#' If the units are a mix of cm and mm, all the data is transformed to mm
#' The returned data is by species, GSA, year, quarter, mid_length, start_length, country and fishery
#' The landings data file must have the following columns: XXXX
#' Note all requested species in all requested GSAs are returned.
#' @param landings_file Name of the landings data file, including location
#' @param species The letter three species code. 
#' @param gsa The GSA number.
#' @example
#' # Don't run
#' # dps17 <- get_length_landings("../Fisheries_data/landings.csv", "DPS", 17)
get_length_landings <- function(landings_file, species, gsa){
    landat <- fread(landings_file)
    # Area column has GSA or SA. Want just numeric.
    # Make a new numeric GSA column based on area
    landat$gsa <- as.numeric(gsub("[^0-9]", "", landat$area))
    # Clean species column - remove white space and make capitals
    landat$species <- toupper(gsub("\\s", "", landat$species))
    # Select columns of interest
    lccols <- colnames(landat)[grep("lengthclass", colnames(landat))]
    cols <- c("year", "quarter", "gsa", "species", "unit","gear","fishery","country", lccols)
    landat <- landat[, cols,with=FALSE]
    # Change lengthclass columns to just the value
    colnames(landat)[colnames(landat) %in% lccols] <- gsub("[^0-9]", "", lccols)

    # Push into helpful shape
    mdat <- melt(landat, id.vars=c("year", "quarter", "gsa", "species", "unit", "country", "gear", "fishery"), variable.name="start_length", value.name="value")
    mdat$start_length <- as.numeric(mdat$start_length)
    # Set -ve values to NA 
    mdat[(mdat$value<0) & !is.na(mdat$value),"value"] <- NA

    # Subset GSA and species - fix for multiple areas
    # Stupid data table thing with names and columns names
    sp <- species
    GSA <- gsa
    mdat <- mdat[(species%in%sp) & (gsa%in%GSA)]

    #The units for some stocks are a mix of cm and mm across the different GSAs.
    #To allow us to combine GSAs we transform all records in cm to mm.
    # Except that sometimes the unit is wrong - so ignore!
    # Can only be mm for crustaceans and cephalopod
    #if(all(c("cm","mm") %in% unique(mdat$unit))){
    #    mdat[unit=="cm"]$start_length <- mdat[unit=="cm"]$start_length * 10
    #    mdat[unit=="cm"]$unit <- "mm"
    #}
    if (all(c("cm","mm") %in% unique(mdat$unit))){
        warning("Mix of mm and cm in dataset")
    }

    # We need a mean length column (the length classes are 1cm wide).
    # Add mid lengths columns
    mdat$mid_length <- mdat$start_length + 0.5
    return(mdat)
}

plot_dist <- function(len_data, title="", max_length = 30){
    nyears <- length(unique(len_data$year))
    ncolp <- floor(sqrt(nyears))
    p <- ggplot(len_data[mid_length < max_length], aes(x=mid_length, y=value)) +
        geom_bar(stat="identity") + facet_wrap(~year, ncol=ncolp) +
        ggtitle(title) + xlab("Mid length") + ylab("Count")
    return(p)
}

