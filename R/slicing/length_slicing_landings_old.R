# length_slicing_landings.R
# Copyright 2013 Finlay Scott and Chato Osio
# Maintainer: Finlay Scott, JRC, finlay.scott@jrc.ec.europa.eu
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#--------------------------------------------------------------------------

# This script demonstrates how to read in length-based abundance data from
# the landings database and convert it to age-based data.
# Two different slicing methods are presented:
#   knife-edge

# The same techniques can be for used for length-based discards data.
# We present the knife-edge method first. This is a simple method that can
# be applied to all data sets.
# These methods are also used when slicing the medits data.

#--------------------------------------------------------------------------

# Libraries and functions

rm(list=ls())
# Load your libraries etc
library(RODBC) # odbcConnectAccess
library(ggplot2)
library(FLCore)
library(mixdist)
library(reshape2)
library(plyr)

source("length_slicing_funcs.R")

#--------------------------------------------------------------------------

# Reading in data from database, correct issues and prepare for slicing

# Details of the stock - used to pull data from the database
gsa <- "SA 7"
species <- "HKE"

# Pulling the data from the database
# Unfortunately this only works with Windows at the moment due to it being
# an Access database.
# You need to set the location of the database
database_location <- "c:/Work/Dec2013/A_D_201311"
database_name <- "B Fisheries landings at length data MED 2002-2011 2013.mdb"
ch <- odbcConnectAccess(paste(database_location,database_name,sep="/"))
# Pull out species from the database
qry1 <- paste("SELECT * FROM \"all\" WHERE AREA='", gsa, "' AND SPECIES='",
              species, "'", sep="")
dat <- sqlQuery(ch, qry1)
close(ch)

# The data is in an object dat, with lengths as columns
# Find last column before the LENGTHCLASSXX columns
last_col <- min(grep("LENGTHCLASS", x=colnames(dat)))-1
# We melt it into a more helpful shape
dm <- melt(dat, id.vars=1:last_col, variable.name="len")

# We want to split the length column so we only have the numbers
dm$len <- as.character(dm$len)
dm$len <- as.numeric(substring(dm$len, first=nchar("LENGTHCLASS")+1, last = 1000000L))
# Remove empty lengths
dm <- dm[dm$value != -1,]
# Scale the value up by 1000 if necessary (what are your units?)
dm$value <- dm$value * 1000

# Check your bin widths!
# The following will produce a table of length bins by year
# They are probably by 1 cm
dcast(ddply(dm, .(YEAR), summarise, lengths = unique(len)), lengths~YEAR)

# We now sum the numbers in each length class by Quarter and Year
dm <- ddply(dm, .(YEAR, len, QUARTER), summarise, value = sum(value))

# Plot a histogram of our data by Year and Quarter
p <- ggplot(dm) +
 geom_histogram(aes(len,weight=value),colour="darkgreen",fill="white",binwidth=4) +
 scale_x_continuous(name="Length") + scale_y_continuous(name ="Frequency") +
 facet_grid(YEAR~QUARTER)
p

# The timing of the landings through the year is important
# Here we convert the QUARTER to be a proportion through the year
# We assume that the landings happened in the middle of the quarter
# e.g. quarter 1 = 0.25/2
# A QUARTER of -1 means no quarterly information.
# It is assumed to happen in the middle of the year (0.5)
dm$timing <- NA
dm[dm$QUARTER==-1,"timing"] <- 0.5
dm[dm$QUARTER==1,"timing"] <- 0.25/2
dm[dm$QUARTER==2,"timing"] <- 0.25 + (0.5 - 0.25)/2
dm[dm$QUARTER==3,"timing"] <- 0.50 + (0.75 - 0.5)/2
dm[dm$QUARTER==4,"timing"] <- 0.75 + (1.0 - 0.75)/2


#--------------------------------------------------------------------------

# Knife-edge slicing

# The knife-edge slicing method is a simple method based on the
# von Bertalanffy growth curve. The age of the sample is calculated from
# the length using the inverse vonB equation:
# age = t0 - log(1 - L / Linf) / K
# The age is then rounded down to the nearest integer (e.g. an age 4.9
# fish is rounded down to age 4).
# A minimum age can be set (the default is 1).
# If the length is equal to or larger than Linf then the length is set
# just below Linf. This is to avoid numerical problems.
# The knife-edge function calculates the age of a fish at the start of the year.
# It therefore takes into account when during the year the sample was taken.
# This is set using the 'timing' option, which has the range 0 (start of year)
# to 1 (end of the year).

# The QUARTER column in the data from the landings database is used to set
# the timing option.

# Set some biological information about the stock
vB_pars =c(Linf=87,K=0.2,t0=-0.03)
# The maximum age you want to have in your stock - this will depend on the data
# and the stock
plusGroup <- 4
minage <- 1
ages<-minage:plusGroup

# We need to get the lengths in the middle of length class
# Here, the length classes are in 1cm bins. We therefore add 0.5 to the lengths
# so that the lengths are mid bin
dm$len_mid <- dm$len+0.5

# Apply the knife_edge function to the data
caa_qy <- ddply(dm, .(YEAR, timing), function(x,vB_pars,plusGroup,minage)
        knife_edge(x[,c("len_mid","value")], timing=x$timing[1], vB=vB, plusGroup=plusGroup, minage=minage),
        vB=vB, plusGroup=plusGroup, minage=minage)

# This gives us numbers at age by year and quarter (timing)
# We want numbers by year so we sum over the quarters
caa <- ddply(caa_qy, .(YEAR, age), summarise, value = sum(x))

# Rename the columns to something more useful for FLR
names(caa)=c("year","age","data")

# Finally turn the results into an FLQuant that can be used
# to make an FLStock
caaq <- as.FLQuant(caa[,c("age","year","data")])
plot(caaq)

#--------------------------------------------------------------------------

# FAO proportional slicing

# Set the biological parameters
vB_pars =c(Linf=87,K=0.2,t0=-0.03)

# Try to break this
vB_pars =c(Linf=50,K=0.2,t0=-0.03)

# The maximum age you want to have in your stock - this will depend on the data
# and the stock
plusGroup <- 4
minage <- 1
ages<-minage:plusGroup

head(dm)
unique(dm$len)


# Get Proportion Matrix up to Linf
prop_mat <- get_prop_matrix(K=vB_pars["K"],Linf=vB_pars["Linf"],t0=vB_pars["t0"],binwidth=1,plusgroup=plusGroup)

dm1 <- dm[dm$YEAR==2007 & dm$timing==0.125,]

# Careful - all obs with lengths > Linf are chopped off
ldatx <- expand_ldat(dm1[,c("len","value")],prop_mat)

# check dims are correct first
temp_natage <- apply(sweep(prop_mat,1,ldatx,"*"),2,sum)
# put into data.frame


#--------------------------------------------------------------------------
# Stop