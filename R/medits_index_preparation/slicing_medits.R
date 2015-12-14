# Very quick and dirty example of slicing the Medits data with
# the knife edge data to prepare the indices

# Assume you have already generated the length based index csvtables using LFD.R
library(plyr)
library(ggplot2)
source("../slicing/length_slicing_funcs.R")
tabdir <- "../../../medits_test/tables/"

# Set the VB parameters
vB <- c(Linf = 130, K = 0.2, t0 = -0.01)

# Read in your length based data
len_data <- read.csv(paste0(tabdir,"stratified_Nlen_MERL_MER_all_GSA9.csv"), header=TRUE)

# Timing is part way through year.
# Add a new column for timing based on month.
# e.g. month = 6,  timing = 0.5
len_data$timing = len_data$month/12

# Set the minimum age and plusgroup
minage <- 1
plusGroup <- 10

# Slice!
# Slice by year and timing
age_data <- ddply(len_data, .(year, timing), function(x) {
      knife_edge(x[,c("len","value")], vB, timing = x$timing[1], plusGroup=plusGroup, minage=minage)
})

head(age_data)
# Now sum over timing (which is now redundant)
age_data <- ddply(age_data, .(year, age), summarise, data = sum(data))


# Plot by age
ggplot(age_data, aes(x = age, y=data)) + geom_bar(stat="identity") + facet_wrap(~year)

# Now put this into an FLIndex for use with your assessment
library(FLCore)
flq <- as.FLQuant(age_data)
idx <- FLIndex(index=flq)


