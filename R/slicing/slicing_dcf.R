# Applies the simple knife-edge slicing method to the length based landings 
# data. 

# Copyright Laurence Kell, Alexander Kell, Finlay Scott, Graham Pilling,
# Valerio Bartolino, Chato Osio, Max Cardinale (2015)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#--------------------------------------------------------------------------

# Simple example of slicing the length data with the knife edge data
# to prepare the indices

# Libraries and functions
library(plyr)
library(ggplot2)
library(data.table)
source("../slicing/length_slicing_funcs.R")

# Set the GSA and species
gsa <- 9
species <- "HKE"

# Set the VB parameters in this format
vB <- c(Linf = 130, K = 0.2, t0 = -0.01)

# Read in your length based data
# Check where you saved it
dat <- fread("../../../db/tables/dcf/landings_length.csv", header=TRUE, sep=";")
# Pull out your gsa and species
dat <- subset(dat, area=gsa, species=species)

# Reshape data.frame into a more helpful shape
last_col <- min(grep("lengthclass", x=colnames(dat)))-1
dm <- melt(dat, id.vars=1:last_col, variable.name="len")
# (ignore warning)

# Coerce value to numeric
dm$value <- as.numeric(dm$value)

# We want to split the length column so we only have the numbers
dm$len <- as.character(dm$len)
dm$len <- as.numeric(substring(dm$len, first=nchar("lengthclass")+1, last = 1000000L))
# Remove empty lengths
dm <- dm[dm$value != -1,]
# Scale the value (the counts) up by 1000 if necessary (what are your units?)
dm$value <- dm$value * 1000

# We need to get the lengths in the middle of the bins
# So we need to check the bin widths
# The following will produce a table of length bins by year
# They are probably by 1 cm
dcast(ddply(dm, .(year), summarise, lengths = unique(len)), lengths~year)

# We now sum the numbers in each length class by Quarter and Year
dm <- ddply(dm, .(year, len, quarter), summarise, value = sum(value))

# Take a look at the length frequencies
# Plot a histogram of our data by Year and Quarter
p <- ggplot(dm) +
 geom_histogram(aes(len,weight=value),colour="darkgreen",fill="white",binwidth=4) +
 scale_x_continuous(name="Length") + scale_y_continuous(name ="Frequency") +
 facet_grid(year~quarter)
p

# The timing of the landings through the year is important
# Here we convert the 'quarter' column to be a proportion through the year
# We assume that the landings happened in the middle of the quarter
# e.g. quarter 1 = 0.25/2
# A quarter of -1 means no quarterly information.
# It is assumed to happen in the middle of the year (0.5)
dm$timing <- NA
dm[dm$quarter==-1,"timing"] <- 0.5
dm[dm$quarter==1,"timing"] <- 0.25/2
dm[dm$quarter==2,"timing"] <- 0.25 + (0.5 - 0.25)/2
dm[dm$quarter==3,"timing"] <- 0.50 + (0.75 - 0.5)/2
dm[dm$quarter==4,"timing"] <- 0.75 + (1.0 - 0.75)/2

# The minimum and plusgroup 
# Depends on the data and the stock
plusGroup <- 10
minage <- 1
ages<-minage:plusGroup

# We need to get the lengths in the middle of length class
# Here, the length classes are in 1cm bins. We therefore add 0.5 to the lengths
# so that the lengths are mid bin
dm$len_mid <- dm$len+0.5

# Finally, apply the knife_edge function to the data
caa_qy <- ddply(dm, .(year, timing), function(x,vB_pars,plusGroup,minage)
        knife_edge(x[,c("len_mid","value")], timing=x$timing[1], vB=vB, plusGroup=plusGroup, minage=minage),
        vB=vB, plusGroup=plusGroup, minage=minage)

# This gives us numbers at age by year and quarter (timing)
# We want numbers by year so we sum over the quarters
caa <- ddply(caa_qy, .(year, age), summarise, data = sum(data))

# Finally turn the results into an FLQuant that can be used
# to make an FLStock
library(FLCore)
caaq <- as.FLQuant(caa[,c("age","year","data")])
plot(caaq)


