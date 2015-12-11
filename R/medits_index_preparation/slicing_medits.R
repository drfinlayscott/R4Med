# Very quick and dirty example of slicing the Medits data with
# the knife edge data to prepare the indices

# Assume you have already generated the length based index csvtables using LFD.R
library(plyr)
library(ggplot2)
source("../slicing/length_slicing_funcs.r")
tabdir <- "../../../medits_test/tables/"

# Set the VB parameters
vB <- c(Linf = 130, K = 0.2, t0 = -0.01)

# Read in your length based data
len_data <- read.csv(paste0(tabdir,"stratified_Nlen_MERL_MER_all_GSA9.csv"), header=TRUE)

# Timing is part way through year
# Here month is 6 so timing = 0.5
# Slice!

# By quarter? 

age_data <- ddply(len_data, .(year), function(x) {
      knife_edge(x[,c("len","value")], vB, timing = 0.5, plusGroup=10, minage=1)
})

# Plot by age
ggplot(age_data, aes(x = age, y=data)) + geom_bar(stat="identity") + facet_wrap(~year)

# Now put this into an FLIndex for use with your assessment
library(FLCore)
flq <- as.FLQuant(age_data)
idx <- FLIndex(index=flq)


