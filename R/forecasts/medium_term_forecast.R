# medium_term_forecast.R
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

# Generic script for running medium-term forecasts (MTF).
# This script assumes that have already run your assessment and that you have a fully specified age-structured FLStock object
# Running a MTF requires a stock-recruitment relationship (SRR).
# The first part of the script tries to fit one.
# If the fit is no good, you cannot run a MTF

#------------------------------------------------------------------
# Libraries and data

rm(list=ls())
library(FLCore)
library(FLAssess)
library(FLash)
library(ggplotFL)
library(FLBRP)
library(plyr)
library(reshape2)

# Example data set - use your own
# You need a full specified FLStock object
# Here I'm loading a dummy stock object
load("../../data/stk.RData")
# Load your own data, probably using the load() function

# Quick check that the stock object is correct
summary(stk)
plot(stk)

# For the MTF we would like to run a F0.1 scenario
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
# Fitting the stock-recruitment relationship (SRR)
# As we are projecting for several years, the MTF depends on
# having a reasonable stock recruitment relationship.
# Therefore if you do not have a reasonable SRR you cannot perform an MTF

# BROKEN ATM
# Attempt to fit SRR - only if you have sufficient data
# Try different models, and try fixing different parameters
# Here we fit a Beverton-Holt, parameterised with steepness, virgin biomass and spr0.
# Other models are available (Ricker etc.)
# We fix the values of steepness and spr0
# spr0 is estimated using the spr0 method in FLR
spr0 <- spr0(stk)
# steepness we fix at 0.75 - you should try other values
steepness <- 0.75
srr_sv <- fmle(as.FLSR(stk, model = "bevholtSV"),
               fixed = list(s = steepness, spr0 = spr0),
               model="Brent")

# Is the fit any good? 
plot(srr_sv)
# You should try different levels of steepness, models etc.
# If you do not have a good SRR you cannot do a MTF

# As a last resort you could try and fit a segmented regression (!)
srr_sg <- fmle(as.FLSR(stk, model = "segreg"))
plot(srr_sg)

# This looks reasonable so we are happy to carry on with the MTF
# At the moment the SRR uses Steepness and Virgin Biomass parameterisation
# To use this with fwd() and FLBRP() we need to transform it to the standard a-b parameterisation
srr_ab <- as.FLSR(stk, model="bevholt")
params(srr_ab) <- ab(params(srr_sv), model='bevholtSV')[c('a','b'),]
# Also copy across the residuals because you need them for a stochastic projection
residuals(srr_ab) <- residuals(srr_sv)

#--------------------------------------------------------------------------

# Set up the MTF
# The MTF will run from the first year after the final year in our stock object upto 2020
final_year <- 2020
# The start year is 2013 - change this if necessary
mtf_years <- 2013:final_year
no_mtf_years <- length(mtf_years)

# Set up the future stock object.
# Here we use the default assumptions (e.g. weights are means of the last 3 years)
# NOTE: You may want to change some of these assumptions by hand
# See the help page: ?stf for more details
stk_mtf <- stf(stk, nyears = no_mtf_years, wts.nyears = 3)

# The MTF performs forecasts with 6 different fishing strategies:
# F = status quo (mean of last X years)
#     F01,
#     decrease 10% pa, 
#     decrease to F01 by 2014
#     decrease to F01 by 2015
#     decrease to F01 by 2020
no_scenarios <- 6
# We make a matrix that holds the F values of each scenario
# Rows are scenarios, columns are years
fbar_scenarios <- matrix(NA,ncol=no_mtf_years, nrow=no_scenarios)
dimnames(fbar_scenarios) <- list(scenario = 1:no_scenarios, year = mtf_years)
# Fill up the scenarios
fbar_scenarios[1,] <- fbar_status_quo # status quo
fbar_scenarios[2,] <- f01 # F01
fbar_scenarios[3,] <- fbar_status_quo * 0.9^(0:(no_mtf_years-1)) # decrease10% every year
fbar_scenarios[4,] <- c(((f01 - fbar_status_quo) / (2014 - 2013)) * (0:1) + fbar_status_quo, rep(f01, no_mtf_years - 2)) # decrease to F01 by 2014
fbar_scenarios[5,] <- c(((f01 - fbar_status_quo) / (2015 - 2013)) * (0:2) + fbar_status_quo, rep(f01, no_mtf_years - 3)) # decrease to F01 by 2015
fbar_scenarios[6,] <- ((f01 - fbar_status_quo) / (2020 - 2013)) * (0:(no_mtf_years-1)) + fbar_status_quo # decrease to F01 by 2020
# Do these look sensible?
fbar_scenarios
# Have a quick look
test <- melt(fbar_scenarios)
test$scenario <- factor(test$scenario)
ggplot(test, aes(x=year, y=value, colour=scenario)) + geom_line()

#--------------------------------------------------------------------------

# Including stochasticity from the SRR residuals

# It is possible to use the residuals from the SRR to include
# stochasticity in the projection.
# This depends on how long your time series of data is.
# How many iterations do you want? Probably no less than 500
niters <- 500
# Blow up the stock object to have many iterations
stk_mtf <- propagate(stk_mtf,niters)
# Set up some stochasticity on the SRR residuals
# Make an empty FLQuant to store the residuals
multi_rec_residuals <- FLQuant(NA, dimnames = list(year=mtf_years, iter=1:niters))
# Sample with replacement from the residuals in the SRR object
sample_years <- sample(dimnames(residuals(srr_ab))$year, niters * length(mtf_years), replace = TRUE)
# The residuals are on a log scale so we take exp() to scale them up
multi_rec_residuals[] <- exp(residuals(srr_ab)[,sample_years])

#--------------------------------------------------------------------------

# Now do the projections

# Store the results from each scenario in a list
mtf_results <- list()
# Loop over the scenarios
for (scenario in 1:nrow(fbar_scenarios)) {
    cat("Scenario: ", scenario, "\n")
    # Make the target object for the projection
    ctrl_target <- data.frame(year = mtf_years,
                              quantity = "f",
                              val = fbar_scenarios[scenario,])
    # Set the control object - year, quantity and value for the moment
    ctrl_f <- fwdControl(ctrl_target)
    # Are you doing a stochastic projection, if so use this?
    stk_mtf_fwd <- fwd(stk_mtf, ctrl = ctrl_f, sr = srr_ab, sr.residuals = multi_rec_residuals, sr.residuals.mult = TRUE, maxF = 10.0)
    # If no stochasticity, use this line instead
    #stk_mtf_fwd <- fwd(stk_mtf, ctrl = ctrl_f, sr = srr_ab)
    # Check it has worked - uncomment out to check scenario by scenario
    #plot(stk_mtf_fwd)

    # Results
    mtf_results[[scenario]] <- stk_mtf_fwd
}

# Name the results list
names(mtf_results) <- c("Constant status quo",
                        "Constant F0.1",
                        "Decrease by 10% p.a.",
                        "Decrease to F0.1 by 2014",
                        "Decrease to F0.1 by 2015",
                        "Decrease to F0.1 by 2020")

# Plot a scenario just to check it out
plot(mtf_results[[4]])

#--------------------------------------------------------------------------

# Summaries and plots

# Generate summaries of these results: SSB, Fbar, Rec and Catch
# Pull out the median and some other quantile levels (here the 10% and 90%)
quantile_levels <- c(0.1,0.5,0.9)

summary_table <- ldply(mtf_results, function(x){
    ssbs <- ssb(x)
    ssb_ql <- apply(ssbs,1:5,quantile, probs=quantile_levels)
    ssb_ql <- cbind(measure = "SSB",melt(drop(ssb_ql)))
    recs <- rec(x)
    rec_ql <- apply(recs,1:5,quantile, probs=quantile_levels)
    rec_ql <- cbind(measure = "Recruitment",melt(drop(rec_ql)))
    fbars <- fbar(x)
    fbar_ql <- apply(fbars,1:5,quantile, probs=quantile_levels)
    fbar_ql <- cbind(measure = "Fbar",melt(drop(fbar_ql)))
    catchs <- catch(x)
    catch_ql <- apply(catchs,1:5,quantile, probs=quantile_levels)
    catch_ql <- cbind(measure = "Catch",melt(drop(catch_ql)))
    output <- rbind(ssb_ql, rec_ql, fbar_ql, catch_ql)
    return(output)})

# Flatten this for plotting purposes 
stc <- dcast(summary_table, .id + measure + year ~ Var1)
# Having a '%' in your column name is a problem
# So we replace them
names(stc)[names(stc) %in% paste(quantile_levels * 100,"%",sep="")] <- paste("q",quantile_levels*100,sep="")

# A selection of plots

# Summary plot
min_year <- 2005
legend_title <- "F scenario"
p <- ggplot(stc[stc$year>=min_year,])
p <- p + geom_line(aes(x=year,y=q50, colour = .id)) + facet_wrap(~measure, scales="free")
p <- p + geom_ribbon(aes(x=year,ymin = q10, ymax = q90, fill=.id), alpha = 0.3)
p <- p + labs(colour = legend_title, fill=legend_title)
p <- p + scale_y_continuous(name="")
p <- p + theme(legend.position="bottom", legend.direction = "vertical", legend.title.align=0.5)
p <- p + guides(colour = guide_legend(nrow=2))
p
# Looks a bit crowded
# Export the plot
#png(file="example_mtf_plot.png")
#print(p)
#dev.off()
# Or use ggsave
ggsave(filename="example_mtf_plot.png", plot=p)

# Plot one scenario with a few iterations
iters <- 1:5
plotting_scenario <- "Constant F0.1"
# Make the data
iter_data <- rbind(
    cbind(measure = "SSB",as.data.frame(ssb(mtf_results[[plotting_scenario]])[,ac(min_year:final_year),,,,iters])[,c("year","iter","data")]),
    cbind(measure = "Recruitment",as.data.frame(rec(mtf_results[[plotting_scenario]])[,ac(min_year:final_year),,,,iters])[,c("year","iter","data")]),
    cbind(measure = "Fbar",as.data.frame(fbar(mtf_results[[plotting_scenario]])[,ac(min_year:final_year),,,,iters])[,c("year","iter","data")]),
    cbind(measure = "Catch",as.data.frame(catch(mtf_results[[plotting_scenario]])[,ac(min_year:final_year),,,,iters])[,c("year","iter","data")]))
p <- ggplot(stc[stc$year>=min_year & stc$.id == plotting_scenario,])
p <- p + geom_line(aes(x=year,y=q50)) + facet_wrap(~measure, scales="free")
p <- p + geom_ribbon(aes(x=year,ymin = q10, ymax = q90), alpha = 0.1)
p <- p + scale_y_continuous(name="")
p <- p + theme(legend.position = "none")
# Superimpose a few iterations
p <- p + geom_line(data=iter_data, aes(x=year,y=data, colour = iter))
p

# Export the plot
ggsave(filename="example_mtf_plot2.png", plot=p)
#png(file="example_mtf_plot2.png")
#print(p)
#dev.off()

#--------------------------------------------------------------------------
# BRP
# If you have reasonable SRR you can estimate reference points
stk_brp <- brp(FLBRP(stk, sr = srr_ab))
# plot(stk_brp) # is broken... 
refpts(stk_brp)






