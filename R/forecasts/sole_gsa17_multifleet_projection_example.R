
# Approximate multifleet projections with FLR
# An example using Sole in GSA 17
# Copyright 2014 Finlay Scott and Chato Osio
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

#------------------------------------------------------------------

# This script takes the ideas in the multifleet_generic.R script
# and uses them to make a stock specific example.
# Here we work with Sole in GSA 17 and explore what happens
# if one of the three fleets (Trawl) stops operating,

# This example is specific to Sole in GSA 17, so you should think
 #of your own scenarios.

#------------------------------------------------------------------

# Ingredients
# 1. A full FLStock i.e. results of a stock assessment.
# 2. Some way of estimating the selectivity patterns of the fleets.
#    In this example we use the historical catches by age of each fleet to estimate the partial fishing mortality.
# 3. An idea of what scenarios you would like to run.


#------------------------------------------------------------------

# Install packages if you don't have them.
# Use R 3.1

# from CRAN
#install.packages(c("copula","triangle","ggplot2"))
# from CRAN
#install.packages(c("plyr","xtable","plot3D"))
# from FLR
#install.packages(c("FLCore", "FLa4a"), repos="http://flr-project.org/R")
# from FLR
#install.packages(c("FLXSA","FLAssess","FLash","FLBRP","ggplotFL"), repos="http://flr-project.org/R")


#------------------------------------------------------------------

# Libraries and whatnot

rm(list=ls())
library(FLCore)
library(FLAssess)
library(FLa4a)
library(FLash)
library(ggplotFL)
library(plyr)

#---------------------------------------------------------

# Running the assessment

# FLIndices
idxs <- readFLIndices("../../data/sole_gsa17/TUNEFF.DAT")
# FLStock - the output from the SS3 assessment
load("../../data/sole_gsa17/stk.Rdata") 

# Use FLa4a to rerun the assessment
fmodel <- ~te(age, year, k = c(4, 10)) + s(year, k = 5, by = as.numeric(as.numeric(age == 1)))
qmodel <- list(~s(age, k=3))
rmodel <- ~factor(year)
fit <- a4aSCA(stock=sole,
        indices = idxs,
        fmodel = fmodel,
        qmodel = qmodel,
        srmodel = rmodel)
sole <- sole + fit

# Check out the results
plot(sole)

#------------------------------------------------------------------

# We are really running a 'single' fleet projection.
# The fishing mortality pattern of this single fleet is the sum of the fishing mortality pattern of the three fleets (set net, trammel net, trawl).
# To project with only 1 or 2 fleets, use the sum of the fishing mortality patterns of these fleets.
# This means that we need the fishing mortality pattern of the different fleets.
# Calculate partial fishing mortalities.

# Here we calculate the partial fishing mortalities by calculating the partial catches of the fleets.

# Sole GSA 17 has three fleets
# Set net
# Trammel net
# Trawl

#------------------------------------------------------------------

# Read in catch numbers from the fleets and make FLQuants of them
set_net <- read.csv("../../data/sole_gsa17/ITA_SET_NET.csv", header=TRUE, sep=";")
set_net <- FLQuant(t(as.matrix(set_net)[,2:8]), dimnames=list(age=0:6,year=2000:2012))
# for example
set_net

trawl <- read.csv("../../data/sole_gsa17/ITA_TRAWL.csv", header=TRUE, sep=";")
trawl <- FLQuant(t(as.matrix(trawl)[,2:8]), dimnames=list(age=0:6,year=2000:2012))

tram <- read.csv("../../data/sole_gsa17/SLO_CRO_TRAMMEL.csv", header=TRUE, sep=";")
tram <- FLQuant(t(as.matrix(tram)[,2:8]), dimnames=list(age=0:6,year=2000:2012))

# Calculate the proportion of catches from each fleet
total_catch <- set_net + trawl + tram
prop_catch_set_net <-  set_net / total_catch
prop_catch_trawl <- trawl / total_catch
prop_catch_tram <- tram / total_catch
# Put all this into an FLQuants object to make operating on them easier
prop_catches <- FLQuants(set_net = prop_catch_set_net, trawl=prop_catch_trawl, tram=prop_catch_tram)

# What do these look like?
ggplot(as.data.frame(prop_catches), aes(x=age,y=data)) + geom_line(aes(colour=qname)) + facet_wrap(~year)


# We now calculate the partial fishing mortalities of each fleet
# catch proportion of each fleet * estimated harvest rate from the stock assessment
pfs <- lapply(prop_catches, function(x) sweep(harvest(sole), 1:5, x, "*"))
# Plot these
pfs_df <- as.data.frame(pfs)
ggplot(pfs_df[pfs_df$year>2003,], aes(x=age, y=data)) + geom_line(aes(colour=qname)) + facet_wrap(~year)


# For the projections we are going to assume that the selectivities don't change over time
# So we need a single selectivity pattern for each fleet
# We take the mean over the years 2006 to 2012
pfs_mean <- lapply(pfs, function(x) apply(x[,as.character(2006:2012)] ,c(1,3:6),mean))
# Giving our selection patterns as
ggplot(as.data.frame(pfs_mean), aes(x=age, y=data)) + geom_line(aes(colour=qname)) 

# Don't worry about the actual values
# Inside the projections the fishing mortalities at age are scaled to hit the desired target.
# We are only interested in the relative shape here.

#---------------------------------------------------------

# Set up the projections
# Remember that the fwd() function uses the F pattern that
# has been set up for the future.
# This F pattern is then scaled to hit the  desired target

# Short term forecast - no stock recruitment relationship
nyears <- 3

# We are interested in projecting with and without the trawl fleet.
# This means we need to set up 2 fleets for the projection:
# 1. Selectivity pattern is for all fleets 
# 2. Selelctivity pattern is only the set net and trammel fleets

# Make the extended stocks

# Stock 1 - all the fleets
sole_allfleets <- stf(sole, nyears=3)
# By default the fishing mortality pattern in the future years is set
# as the mean of the last three years. We don't want that.
# We need to overwrite the harvest slot in the future years
# (2013 to 2015) to be the F patterns we calculated above.
# Current pattern - mean of last three years from assessment
harvest(sole_allfleets)[,as.character(2013:2015)]
# Overwrite this
harvest(sole_allfleets)[,as.character(2013:2015)] <- pfs_mean[["set_net"]] + pfs_mean[["trawl"]] + pfs_mean[["tram"]]
harvest(sole_allfleets)[,as.character(2013:2015)]

# Stock 2 - just the set net and trammel net fleets
sole_notrawl <- stf(sole, nyears=3)
# Overwrite harvest with sum of set net and trammel net fleets
harvest(sole_notrawl)[,as.character(2013:2015)] <- pfs_mean[["set_net"]] + pfs_mean[["tram"]]
harvest(sole_notrawl)[,as.character(2013:2015)]

# Compare the future fishing mortality patterns
f_df <- rbind(cbind(scenario = "notrawl", as.data.frame(pfs_mean[["set_net"]] + pfs_mean[["tram"]])), cbind(scenario = "allfleets", as.data.frame(pfs_mean[["set_net"]] + pfs_mean[["trawl"]] + pfs_mean[["tram"]])))
# Scale for comparison plots
f_df <- ddply(f_df, .(scenario), transform, scaled_f = data / max(data))
ggplot(f_df, aes(x=age, y=scaled_f)) + geom_line(aes(colour=scenario))

# We now have two stocks with different future F patterns based
# on the fleets that are operating

# Stock recruitment
# We are only performing a short term forecast so do we need to fit a stock recruitment relationship
# Instead we use the geometric mean recruitment of the last 3 years
mean_rec <- exp(mean(log(c(rec(sole)[,as.character(2010:2012)]))))

#------------------------------------------------------------------

# Status quo projection with all fleets and with no trawl fleet

# The status quo Fbar is taken as the mean of the last three years
f_allfleets_sq <- mean(c(fbar(sole)[,as.character(2010:2012)]))

# The status quo Fbar without the trawl fleet is slightly more complicated.
# We are assuming that the set net and trammel nets continue fishing at the
# same level as the last 3 years, but that there is no trawl fleet.
# We can use the partial F object we created earlier and just include the
# set net and trammel net
pfs[["set_net"]] + pfs[["tram"]]
# We want to average over the last three years, and also average over the fbar range
fbar_range <- as.character(sole@range["minfbar"]: sole@range["maxfbar"])
f_notrawl_sq <- mean((pfs[["set_net"]] + pfs[["tram"]])[fbar_range, as.character(2010:2012)])

f_allfleets_sq
f_notrawl_sq


# Set the control object for all fleets
f_allfleets_sq_ctrl <- fwdControl(data.frame(year = 2013:(2013+nyears-1), quantity = "f", val = rep(f_allfleets_sq,nyears)))
# And project
sole_allfleets_sq <- fwd(sole_allfleets, f_allfleets_sq_ctrl, sr=list(model="mean", params=FLPar(a=mean_rec)))
# Have a look
plot(sole_allfleets_sq)

# Set the control object without the trawl
f_notrawl_sq_ctrl <- fwdControl(data.frame(year = 2013:(2013+nyears-1), quantity = "f", val = rep(f_notrawl_sq, nyears)))
f_notrawl_sq_ctrl
# And project
sole_notrawl_sq <- fwd(sole_notrawl, f_notrawl_sq_ctrl, 
    sr=list(model="mean", params=FLPar(a=mean_rec)))
# Compare
plot(FLStocks(status_quo_allfleets = sole_allfleets_sq, status_quo_notrawl = sole_notrawl_sq))

# We can see that without the trawl fleet the catches are lower.
#------------------------------------------------------------------

# What happens if the set net and trammel net increase their effort to
# take advantage of the trawl fleet not operating?
# The catches will be the same as the all fleet status quo scenario.
# But only 2 fleets are operating. 

# Equal catches scenario
# Get the projected catches from the all fleet status quo position
catch_sq <- c(catch(sole_allfleets_sq)[,as.character(2013:2015)])
# Create the control object
# Now catch is the target
catch_notrawl_sq_ctrl <- fwdControl(data.frame(year = 2013:(2013+nyears-1),
                                  quantity = "catch",
                                  val = catch_sq))
catch_notrawl_sq_ctrl
# Project
sole_notrawl_equalcatch <- fwd(sole_notrawl, catch_notrawl_sq_ctrl, 
    sr=list(model="mean", params=FLPar(a=mean_rec)))


plot(FLStocks(allfleets = sole_allfleets_sq, notrawl_equalcatch = sole_notrawl_equalcatch))
# We can see that the catches are the same
# But the fishing mortality has in increased to get them
# - a consequence of the trawl fleet not operating.
# SSB is slightly affected

#-------------------------------------------------------
# Reference points
# Different selectivities / fishing mortality patterns
# give different reference points

library(FLBRP)

# By default FLBRP uses the mean fishing pattern in the last 3 years
# to calculate reference points.
# We have already set the fishing pattern for the two scenarios
# (all fleets and no trawl) when we used stf() above.

allfleets_brp <- FLBRP(sole_allfleets)
allfleets_brp <- brp(allfleets_brp)
refpts(allfleets_brp)
allfleets_f01 <- refpts(allfleets_brp)["f0.1","harvest"]

notrawl_brp <- FLBRP(sole_notrawl)
notrawl_brp <- brp(notrawl_brp)
refpts(notrawl_brp)
notrawl_f01 <- refpts(notrawl_brp)["f0.1","harvest"]

allfleets_f01
notrawl_f01
# Hardly any difference

# Project with these values if you want to
# Set ctrl objects
f_allfleets_f01_ctrl <- fwdControl(data.frame(year = 2013:(2013+nyears-1), quantity = "f",val = rep(c(allfleets_f01), nyears)))

f_notrawl_f01_ctrl <- fwdControl(data.frame(year = 2013:(2013+nyears-1), quantity = "f",val = rep(c(notrawl_f01), nyears)))

# Project
sole_allfleets_f01 <- fwd(sole_allfleets, f_allfleets_f01_ctrl, sr=list(model="mean", params=FLPar(a=mean_rec)))

sole_notrawl_f01 <- fwd(sole_notrawl, f_notrawl_f01_ctrl, sr=list(model="mean", params=FLPar(a=mean_rec)))

plot(FLStocks(f01_allfleets=sole_allfleets_f01, f01_notrawl=sole_notrawl_f01))


