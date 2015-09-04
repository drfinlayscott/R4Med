
# Generic approximate multifleet projections with FLR
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

# This script demonstrates how an approximate multifleet projection
# can be run using the FLR tools (fwd()).
# In the example a typical Short Term Forecast is run (based on the
# short_term_forecast.R script) that projects the stock forward
# under a range of different fishing mortality scenarios.

# The resulting catches of the projection are then split between 
# the fleets. The split is based on either the historical partial
# catches or fishing mortalities.

# In the projections, the future fishing mortality pattern at age (the
# selectivity is always the same) implying that the relative activity
# of each of the fleets is the same. 

# Another script is in the repository.
# This looks at how it may be possible to explore the consequences
# of one or more fleets increasing or decreasing their fishing effort.
# For example, what might happen if one of the fleets stops operating.

#------------------------------------------------------------------

# Ingredients
# 1. A full FLStock i.e. results of a stock assessment.
# 2. Some way of estimating the selectivity patterns of the fleets.
#    In this example we use the historical catches by age of each fleet
#    to estimate the partial fishing mortality.

#------------------------------------------------------------------
# Install packages if you don't have them.
# Use R >= 3.0

## from CRAN
#install.packages(c("copula","triangle","ggplot2"))
## from CRAN
#install.packages(c("plyr","xtable","plot3D"))
## from FLR
#install.packages(c("FLCore", "FLa4a"), repos="http://flr-project.org/R")
## from FLR
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
library(FLBRP)

#---------------------------------------------------------
# Running the assessment
# Read in your own stock object or rerun your assessment

# FLIndices
idxs <- readFLIndices("../../data/sole_gsa17/TUNEFF.DAT")
# FLStock - the output from the SS3 assessment
load("../../data/sole_gsa17/stk.Rdata") 

# Here we use FLa4a to rerun the assessment
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

# This is the starting point of the example.

#------------------------------------------------------------------
# We are really running a 'single' fleet projection.
# The fishing mortality pattern of this single fleet is the sum of the
# fishing mortality pattern of the fleets.
# In this example we have 3 fleets (set net, trammel net, trawl).

# We are going to run Short Term Forecast with different F levels.
# These scenarios will use the same selection pattern for the future years.
# This selection pattern represents the combined selectivity of the fleets.

# In this example we use Sole in GSA 17 which has three fleets
# Set net
# Trammel net
# Trawl

# The partial catch data is stored in separate csv files.
# We read each of them in and make FLQuant objects from them.

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
catches <- FLQuants(set_net = set_net, trawl=trawl, tram=tram)
prop_catches <- lapply(catches, function(x) x / total_catch)
# What do these look like?
ggplot(as.data.frame(prop_catches), aes(x=age,y=data)) + geom_line(aes(colour=qname)) + facet_wrap(~year)

# We also calculate the mean partial catch proportions over the period
# 2006 to 2012
# This is the period we will use to calculate the future selectivity
# of each fleet.
mean_catches <- lapply(catches, function(x) apply(x[,as.character(2006:2012)],c(1,3:6),mean))
total_mean_catches <- mean_catches[["set_net"]] + mean_catches[["trawl"]] + mean_catches[["tram"]]
# Get the proprtion of mean catches taken by each fleet
# We will use this later to partition the projected catches.
prop_mean_catches <- lapply(mean_catches, function(x) x / total_mean_catches)
ggplot(as.data.frame(prop_mean_catches), aes(x=age,y=data)) + geom_line(aes(colour=qname)) 

# We now calculate the partial fishing mortalities of each fleet
# pF = catch proportion of each fleet * estimated harvest rate from the stock assessment
pfs <- lapply(prop_catches, function(x) sweep(harvest(sole), 1:5, x, "*"))
# Plot these
pfs_df <- as.data.frame(pfs)
ggplot(pfs_df[pfs_df$year>2003,], aes(x=age, y=data)) + geom_line(aes(colour=qname)) + facet_wrap(~year)

# For the projections we are going to assume that the selectivities don't change over time
# So we need a single selectivity pattern for each fleet
# We take the mean over the years 2006 to 2012
# (same as the partial catches)
pfs_mean <- lapply(pfs, function(x) apply(x[,as.character(2006:2012)] ,c(1,3:6),mean))
# Giving our selection patterns for the future as
ggplot(as.data.frame(pfs_mean), aes(x=age, y=data)) + geom_line(aes(colour=qname)) 

# We can also calculate the proportion of total F from each fleet
total_pf <- pfs_mean[["set_net"]] + pfs_mean[["trawl"]] + pfs_mean[["tram"]]
prop_pf <- lapply(pfs_mean, function(x) x / total_pf)
ggplot(as.data.frame(prop_pf), aes(x=age, y=data)) + geom_line(aes(colour=qname)) 
# Don't worry about the actual values
# Inside the projections the fishing mortalities at age are scaled to hit the desired target.
# We are only interested in the relative shape here.


#---------------------------------------------------------------
# Short term forecast
# Taken from short_term_forecast.R in R4Med repository

# STF
# Here we run the STF for 3 years, 2013, 2014, 2015
# You can change these as appropriate
# The first year of the STF should be the next one after the final year in your stock data
# For example, the final year in the dummy stk object is 2012 so the first year of the STF is 2013
stf_years <- c(2013,2014,2015)
no_stf_years <- length(stf_years)

# Set up the future stock object.
# Here we use the default assumptions about what happens to weights and maturity in the future.
# (e.g. weights are means of the last 3 years)
# NOTE: You may want to change some of these assumptions by hand
# See the help page for stf: ?stf for more details
stf_sole <- stf(sole, nyears = no_stf_years, wts.nyears = 3)

# By default stf() sets up the future F pattern to be the mean of the last
# 3 years. We don't want that. We want to use the F pattern we calculated above
harvest(stf_sole)[,as.character(stf_years)] <- total_pf

# We can now calculate f0.1 using this F pattern
sole_brp <- brp(FLBRP(stf_sole))
refpts(sole_brp)
f01 <- c(refpts(sole_brp)["f0.1","harvest"])
f01

# We also need F status quo - the geometric mean of the last X years
# Here we use 3 years
no_sole_years <- dim(rec(sole))[2]
no_fbar_years <- 3 # Or set your own as appropriate
fbars <- fbar(sole)[,(no_sole_years - no_fbar_years + 1):no_sole_years]
fbar_status_quo <- exp(mean(log(c(fbars))))
fbar_status_quo

# Set up future recruitment to be mean of last X years
# Here we set as geometric mean of the last 3 years
no_rec_years <- 3 # Change number of years as appropriate
recs <- rec(sole)[,(no_sole_years - no_rec_years + 1):no_sole_years]
mean_rec <- exp(mean(log(c(recs))))
mean_rec

# We are going to run several F scenarios for the STF
# The scenarios are based on 'F status quo', which we calculated above as the mean F of the last X years
# An STF is for three years - you could change this but if you do you will have to hack the code below
# For a three year STF the F pattern is:
# year 1: fbar_status_quo
# year 2: fbar_status_quo * fbar_multiplier
# year 3: fbar_status_quo * fbar_multiplier
# The fbar_multiplier is the same for years 2 and 3

# We are going to run several STFs with different values for the fbar_multiplier
# The fbar_multiplier ranges from 0.1 to 2 by 0.1
fbar_multiplier <- seq(from = 0, to = 2, by = 0.1)

# We are going to build a data.frame that builds these scenarios
# Each column in the dataframe is a year
# Each row is a scenario
# Set up the fbar scenarios - note that if you project for more than 3 years you will need to add more columns / years to the matrix
fbar_scenarios <- cbind(rep(fbar_status_quo,length(fbar_multiplier)),
                        fbar_multiplier*fbar_status_quo,
                        fbar_multiplier*fbar_status_quo)
# Add the F0.1 scenario as a final scenario
fbar_scenarios <- rbind(fbar_scenarios, c(fbar_status_quo,f01,f01))
# Add some dimnames
dimnames(fbar_scenarios)[[1]] <- c(as.character(fbar_multiplier),"f0.1")

# There are various results we want to extract from the STF
# Make an empty matrix in which to store the results
stf_results <- matrix(NA,nrow = nrow(fbar_scenarios),ncol = 10)
# Change the column names to reflect years
colnames(stf_results) <- c('Ffactor','Fbar',
						   paste('Catch_', c(stf_years[1]-1,stf_years), sep=""),
						   paste('SSB_', stf_years[-1], sep=""),
						   paste('Change_SSB_', stf_years[2], "_", stf_years[3], "(%)", sep=""),
						   paste('Change_Catch_', stf_years[1]-1, "_", stf_years[2], "(%)", sep=""))


# Store the FLStock each time
sole_stf <- FLStocks()
# Loop over the scenarios
for (scenario in 1:nrow(fbar_scenarios)) {
    cat("Scenario: ", scenario, "\n")
    # Make a target object with F values for that scenario
    ctrl_target <- data.frame(year = stf_years,
                              quantity = "f",
                              val = fbar_scenarios[scenario,])
    # Set the control object - year, quantity and value for the moment
    ctrl_f <- fwdControl(ctrl_target)
    # Run the forward projection. We include an additional argument, maxF.
    # By default the value of maxF is 2.0
    # Here we increase it to 10.0 so that F is not limited
    sole_stf_fwd <- fwd(stf_sole, ctrl = ctrl_f, sr = list(model="mean", params=FLPar(a = mean_rec)), maxF = 10.0)
    ## Check it has worked - uncomment out to check scenario by scenario
    #plot(sole_stf_fwd)
    # Store the result - if you want to, comment out if unnecessary
    sole_stf[[dimnames(fbar_scenarios)[[1]][scenario]]] <- sole_stf_fwd

    # Fill results table
    stf_results[scenario,1] <- fbar_scenarios[scenario,2] / fbar_scenarios[scenario,1] # fbar status quo ratio
    stf_results[scenario,2] <- fbar(sole_stf_fwd)[,ac(stf_years[3])] # final stf year
    stf_results[scenario,3] <- catch(sole_stf_fwd)[,ac(stf_years[1]-1)] # last 'true' year
    stf_results[scenario,4] <- catch(sole_stf_fwd)[,ac(stf_years[1])] # 1st stf year
    stf_results[scenario,5] <- catch(sole_stf_fwd)[,ac(stf_years[2])] # 2nd stf year
    stf_results[scenario,6] <- catch(sole_stf_fwd)[,ac(stf_years[3])] # final stf year
    stf_results[scenario,7] <- ssb(sole_stf_fwd)[,ac(stf_years[2])] # 2nd stf year
    stf_results[scenario,8] <- ssb(sole_stf_fwd)[,ac(stf_years[3])] # final stf year
    # Change in SSB
    stf_results[scenario,9] <- (ssb(sole_stf_fwd)[,ac(stf_years[3])]-ssb(sole_stf_fwd)[,ac(stf_years[2])])/ssb(sole_stf_fwd)[,ac(stf_years[2])]*100 # change in ssb in last two stf years
    stf_results[scenario,10] <- (catch(sole_stf_fwd)[,ac(stf_years[2])]-catch(sole_stf_fwd)[,ac(stf_years[1]-1)])/catch(sole_stf_fwd)[,ac(stf_years[1]-1)]*100 # change in catch from true year, to 2nd to last stf year
}

# We have a lot of results
plot(sole_stf)
plot(lapply(sole_stf, function(x) window(x, start=2004)))

stf_results
# write.csv(stf_results, file="stf_results.csv")

#---------------------------------------------------------------
# Now we can split the catches of the results into the catches
# from each of the fleets.

# We could use the mean partial F proportions or the mean catch proportions (we calculated both of these earlier).
# They are very close so the choice makes little difference
prop_data <- rbind(
    cbind(measure = "catch_prop", as.data.frame(prop_mean_catches)),
    cbind(measure = "f_prop", as.data.frame(prop_pf))
)
ggplot(prop_data, aes(x=age,y=data)) + geom_line(aes(colour=measure)) + facet_wrap(~qname)

# Here we use the partial catches 
# Make FLQuants that go up to 2015 of the partial catch proportions
# We have already calculated the historical partial catch proportions.
# The future patial catch proportions we take as the mean of the years 2006
# to 2012 (we calculated this earlier).
future_prop_catches <- lapply(prop_catches, function(x) window(x, end=stf_years[3]))
# Copy in the future proportions
future_prop_catches[["set_net"]][,as.character(stf_years)] <- prop_mean_catches[["set_net"]]
future_prop_catches[["tram"]][,as.character(stf_years)] <- prop_mean_catches[["tram"]]
future_prop_catches[["trawl"]][,as.character(stf_years)] <- prop_mean_catches[["trawl"]]

# Calculate the future catches for each scenario
# = total catch * catch proportion of fleet
future_catch_set_net <- lapply(sole_stf, function(x) apply(catch.n(x) * future_prop_catches[["set_net"]] * catch.wt(x), 2:6, sum))
future_catch_trawl <- lapply(sole_stf, function(x) apply(catch.n(x) * future_prop_catches[["trawl"]] * catch.wt(x), 2:6, sum))
future_catch_tram <- lapply(sole_stf, function(x) apply(catch.n(x) * future_prop_catches[["tram"]] * catch.wt(x), 2:6, sum))

# Put it all together for a plot
future_catches <- rbind(
    cbind(fleet="set_net",as.data.frame(future_catch_set_net)),
    cbind(fleet="trawl",as.data.frame(future_catch_trawl)),
    cbind(fleet="tram",as.data.frame(future_catch_tram))
)
future_catches$qname <- as.character(future_catches$qname)
old_names <- future_catches$qname[future_catches$qname != "f0.1"]
paste("fsq * ", old_names,sep="")
future_catches$qname[future_catches$qname != "f0.1"] <- paste("fsq * ", old_names,sep="")

ggplot(future_catches, aes(x=year, y=data))+ geom_rect(aes(xmin=stf_years[1], xmax=stf_years[3], ymin=0,ymax=Inf), fill="pink", alpha=0.01) + geom_line(aes(colour=fleet)) + facet_wrap(~qname) + scale_x_continuous(breaks=seq(range(sole_stf)['minyear'], range(sole_stf)['maxyear'], by=2)) + coord_cartesian(xlim=c(stf_years[3]-10,stf_years[3])) 

# Now make a big table of results by exporting future_catches

# Get the mean partial Fs for the F0.1 scenario too
# We need the Fs at age for F0.1 and Fstatus quo
fbar_range <- range(stf_sole)[c('minfbar')]:range(stf_sole)[c('maxfbar')]
# fbar of fishing pattern used in forecast
fbar_fp <- mean(total_pf[ac(fbar_range),])
# f at age at F01
f_f01 <- total_pf * f01 / fbar_fp
# Check - should be same F01
mean(f_f01[ac(fbar_range),])
# Same for F status quo
f_fsq <- total_pf * fbar_status_quo / fbar_fp

# Now split these according the proportional partial Fs we calced earlier
pf_f01 <- lapply(prop_pf, function(x) x * f_f01)
pf_fsq <- lapply(prop_pf, function(x) x * f_fsq)

# Now take fbar of these - assume same fbar range for each gear - bit dodgy
pfbar_f01 <- lapply(pf_f01, function(x) mean(x[ac(fbar_range)]))
pfbar_fsq <- lapply(pf_fsq, function(x) mean(x[ac(fbar_range)]))

# Table of future catches for F0.1 scenario only
# And subset out years we ran the forecast for
out <- subset(future_catches, year %in% stf_years & qname=="f0.1")[,c("fleet", "year", "data")]
names(out)[names(out) == "data"] <- "catches"

# Add partial Fs to table
out$partial_f <- NA
# Fsqs
out[out$year==stf_years[1] & out$fleet=="set_net","partial_f"] <- pfbar_fsq[["set_net"]]
out[out$year==stf_years[1] & out$fleet=="trawl","partial_f"] <- pfbar_fsq[["trawl"]]
out[out$year==stf_years[1] & out$fleet=="tram","partial_f"] <- pfbar_fsq[["tram"]]
# F01s
out[out$year %in% stf_years[2:3] & out$fleet=="set_net","partial_f"] <- pfbar_f01[["set_net"]]
out[out$year %in% stf_years[2:3] & out$fleet=="trawl","partial_f"] <- pfbar_f01[["trawl"]]
out[out$year %in% stf_years[2:3] & out$fleet=="tram","partial_f"] <- pfbar_f01[["tram"]]

# order by year
out[order(out$year),]

#catches_f01_tab <- dcast(f01_scen, year ~ fleet)



