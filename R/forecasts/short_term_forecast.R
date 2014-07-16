# short_term_forecast.R
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

# Generic script for running short-term forecasts (STF).
# This script assumes that have already run your assessment and that you have a fully specified age-structured FLStock object

#------------------------------------------------------------------
# Libraries and data

rm(list=ls())
library(FLCore)
library(FLAssess)
library(FLash)
library(ggplotFL)
library(FLBRP)
#library(plyr)
#library(reshape2)

# Example data set - use your own
# You need a full specified FLStock object
# Here I'm loading a dummy stock object
load("../../data/stk.RData")
# Load your own data, probably using the load() function

# Quick check that the stock object is correct
summary(stk)
plot(stk)

# For the STF we would like to run a F0.1 scenario
# Use FLBRP to get F0.1
stk_brp <- brp(FLBRP(stk))
refpts(stk_brp)
f01 <- c(refpts(stk_brp)["f0.1","harvest"])
f01
# Is this number sensible?

# We also need F status quo - the geometric mean of the last X years
# Here we use 3 years
no_stk_years <- dim(rec(stk))[2]
no_fbar_years <- 3 # Or set your own as appropriate
fbars <- fbar(stk)[,(no_stk_years - no_fbar_years + 1):no_stk_years]
fbar_status_quo <- exp(mean(log(c(fbars))))

#--------------------------------------------------------------------------
# STF
# Here we run the STF for 3 years, 2013, 2014, 2015
# You can change these as appropriate
# The first year of the STF should be the next one after the final year in your stock data
# For example, the final year in the dummy stk object is 2012 so the first year of the STF is 2013
stf_nyears <- 3
final_year <- max(as.numeric(dimnames(stock.n(stk))[[2]]))
stf_years <- (final_year+1):(final_year+stf_nyears)
no_stf_years <- length(stf_years)

# Set up the future stock object.
# Here we use the default assumptions about what happens to weights, maturity and selection pattern in the future
# (e.g. weights are means of the last 3 years)
# NOTE: You may want to change some of these assumptions by hand
# See the help page for stf: ?stf for more details
stf_stk <- stf(stk, nyears = no_stf_years, wts.nyears = 3)

# Set up future recruitment to be mean of last X years
# Here we set as geometric mean of the last 3 years
no_rec_years <- 3 # Change number of years as appropriate
recs <- rec(stk)[,(no_stk_years - no_rec_years + 1):no_stk_years]
mean_rec <- exp(mean(log(c(recs))))

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

# There are various results we want to extract from the STF
# Make an empty matrix in which to store the results
stf_results <- matrix(NA,nrow = nrow(fbar_scenarios),ncol = 10)

# Update column names
colnames(stf_results) <- c('Ffactor',
                           'Fbar',
                           paste('Catch',final_year,sep="_"),
                           paste('Catch',final_year+1,sep="_"), 
                           paste('Catch',final_year+2,sep="_"),
                           paste('Catch',final_year+3,sep="_"),
                           paste('SSB',final_year+2,sep="_"),
                           paste('SSB',final_year+3,sep="_"),
                           paste('Change_SSB_',final_year+2,'-',final_year+3,'(%)',sep=""),
                           paste('Change_Catch_',final_year,'-',final_year+2,'(%)',sep=""))

# Store the FLStock each time
stk_stf <- FLStocks()
# Loop over the scenarios
for (scenario in 1:nrow(fbar_scenarios)) {
    cat("Scenario: ", scenario, "\n")
    # Make a target object withe F values for that scenario
    ctrl_target <- data.frame(year = stf_years,
                              quantity = "f",
                              val = fbar_scenarios[scenario,])
    # Set the control object - year, quantity and value for the moment
    ctrl_f <- fwdControl(ctrl_target)
    # Run the forward projection. We include an additional argument, maxF.
    # By default the value of maxF is 2.0
    # Here we increase it to 10.0 so that F is not limited
    stk_stf_fwd <- fwd(stf_stk, ctrl = ctrl_f, sr = list(model="mean", params=FLPar(a = mean_rec)), maxF = 10.0)
    ## Check it has worked - uncomment out to check scenario by scenario
    #plot(stk_stf_fwd)
    # Store the result - if you want to, comment out if unnecessary
    stk_stf[[as.character(scenario)]] <- stk_stf_fwd

    # Fill results table
    stf_results[scenario,1] <- fbar_scenarios[scenario,2] / fbar_scenarios[scenario,1] # fbar status quo ratio
    stf_results[scenario,2] <- fbar(stk_stf_fwd)[,ac(stf_years[stf_nyears])] # final stf year
    stf_results[scenario,3] <- catch(stk_stf_fwd)[,ac(final_year)] # last 'true' year
    stf_results[scenario,4] <- catch(stk_stf_fwd)[,ac(final_year+1)] # 1st stf year
    stf_results[scenario,5] <- catch(stk_stf_fwd)[,ac(final_year+2)] # 2nd stf year
    stf_results[scenario,6] <- catch(stk_stf_fwd)[,ac(final_year+3)] # final stf year
    stf_results[scenario,7] <- ssb(stk_stf_fwd)[,ac(final_year+2)] # 2nd stf year
    stf_results[scenario,8] <- ssb(stk_stf_fwd)[,ac(final_year+3)] # final stf year
    # Change in SSB
    stf_results[scenario,9] <- (ssb(stk_stf_fwd)[,ac(final_year+3)]-ssb(stk_stf_fwd)[,ac(final_year+2)])/ssb(stk_stf_fwd)[,ac(final_year+2)]*100 # change in ssb in last two stf years
    stf_results[scenario,10] <- (catch(stk_stf_fwd)[,ac(final_year+2)]-catch(stk_stf_fwd)[,ac(final_year)])/catch(stk_stf_fwd)[,ac(final_year)]*100 # change in catch from true year, to 2nd to last stf year
}

# Look at the table of results
stf_results
# export this if necessary
#write.csv(stf_results, file="stf_results.csv")

# Plotting
# Plotting is not necessary for the report but here is a crude one anyway
plot(window(stk_stf, start=2001, end=final_year+3))
stf_results
