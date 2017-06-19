# Copyright European Union, 2016
# Author: Finlay Scott (EC JRC) <finlay.scott@jrc.ec.europa.eu>
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# Get the length based landings data for a species and GSA
# If the units are a mix of cm and mm, all the data is transformed to mm
# The returned data is by species, GSA, year, quarter, mid_length, start_length, country and fishery
# The landings data file must have the following columns: XXXX
# Note all requested species in all requested GSAs are returned.
# @param landings_file Name of the landings data file, including location
# @param species The letter three species code. 
# @param gsa The GSA number.

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

# Used for fitting normal cumulative function
fitnorm <- function(x, length_class, value) {
   muhat <- x[1];
   sdhat <- exp(x[2])
   sum((value-pnorm(length_class,muhat,sdhat))^2)
}

# Get the LC by year
# Two methods: 
#   Take the empirical 25% of cumulative distribution
#   Fit a normal cumulative function and take the 25%
# @param dat A data.table of length distribution with columns year, mid_length and value
get_lc <- function(dat, plot_fit=FALSE, main=""){
    nyears <- length(unique(dat$year))
    nrowp <- floor(sqrt(nyears))
    ncolp <- ceiling(nyears / nrowp)
    if(plot_fit==TRUE){
        par(mfrow=c(nrowp,ncolp))
    }
    lc <- ddply(dat, .(year), function(x){
        lc_emp_25 <- NA
        lc_fit_25 <- NA
        # Only do this if 3 or more length classes
        if(sum(x$value > 0) >= 3){
            # Sort by mean_length
            x <- x[order(x$mid_length),]
            # Get cummulative dist
            cumv <- cumsum(x$value)
            # Method 1: and take the 25% length
            lc_emp_25 <- x$mid_length[which(cumv >= (max(cumv) * 0.25))[1]]
            # Method 2: Try to fit norm and get the 25%
            muinit <- x$mid_length[which.max(x$value)]
            est <- nlm(fitnorm, c(muinit,0), length_class = x$mid_length, value=cumv/max(cumv))
            estmean <- est$estimate[1]
            estsd <- exp(est$estimate[2])
            lc_fit_25 <- qnorm(0.25, estmean, estsd)
            if(plot_fit==TRUE){
                plot(x=x$mid_length, y=cumv/max(cumv), ylab="Probability", xlab="Mid length", main=main)
                lines(x=x$mid_length, y=pnorm(x$mid_length, estmean, estsd))
                lines(x=c(lc_fit_25, lc_fit_25), y=c(0,1), col="red")
                lines(x=c(lc_emp_25, lc_emp_25), y=c(0,1),  col="blue")
            }
        }
        return(data.frame(lc_emp_25 = lc_emp_25,
                          lc_fit_25 = lc_fit_25))
    })
    return(lc)
}

#plot_dist_lc <- function(sp_gsa, dat, max_length=NA){
#    dat <- dat[year %in% sp_gsa$years]
#    lcs <- get_lc(sp_gsa, dat)
#    sp <- sp_gsa[["sp"]]
#    gsas <- sp_gsa[["gsas"]]
#    p <- plot_dist(sp_gsa=sp_gsa, dat=dat, max_length=max_length)
#    lcsm <- melt(lcs, id.var="year", variable.name="Lc")
#    p <- p + geom_vline(data=lcsm, aes(xintercept=value, colour=Lc))
#    return(p)
#}

get_indicators <- function(dat, linf, lc_choice = "fit"){
    if (!(lc_choice %in% c("fit","emp"))){
        stop("lc_choice must be either 'fit' or 'emp'")
    }
    lc_col <- ifelse(lc_choice=="fit", "lc_fit_25", "lc_emp_25")
    lcs <- get_lc(dat, plot_fit=FALSE)
    colnames(lcs)[colnames(lcs) == lc_col] <- "lc"
    lopt <- linf * 2/3
    # lmean
    lmean <- ddply(dat, .(year), function(x){
                   lc <- lcs[lcs$year == x$year[1], "lc"]
                   idx <- x$mid_length > lc
                   lmeany <- sum(x[idx, "mid_length"] * x[idx, "value"], na.rm=TRUE) / sum(x[idx, "value"], na.rm=TRUE)
                   return(data.frame(lmean = lmeany))
    })
    lfem <- (0.75 * lcs[,"lc"]) + (0.25 * linf)
    # Put all together
    indicators <- cbind(lcs[,c("year", "lc")], linf=linf, lopt=lopt, lfem=lfem)
    indicators <- join(indicators, lmean)
    indicators$lmean_lopt <- indicators$lmean / indicators$lopt
    indicators$lmean_lfem <- indicators$lmean / indicators$lfem
    return(indicators)
}

plot_dist_indicators <- function(dat, linf, max_length=30, lc_choice="fit", title=""){
    indicators <- get_indicators(dat=dat, linf=linf, lc_choice=lc_choice)
    indicators_sub <- indicators[,c("year", "lc", "linf", "lfem", "lmean")]
    indim <- melt(indicators_sub, id.var="year", variable.name="Indicators")
    p <- plot_dist(len_data=dat, title=title, max_length = max_length)
    p <- p + geom_vline(data=indim, aes(xintercept=value, colour=Indicators))
    return(p)
}

#indicators can be "lmean_lopt" or "lmean_lfem"
plot_ind <- function(dat, linf, lc_choice="fit", title="", indicators = c("lmean_lopt", "lmean_lfem")){
    ind <- get_indicators(dat=dat, linf=linf, lc_choice=lc_choice)
    indm <- melt(ind[,c("year", indicators)], id.vars="year")
    ymax <- max(1.25, max(indm$value)*1.1)
    ymin <- min(0.75, min(indm$value)*0.9)
    p <- ggplot(indm, aes(x=year, y=value, group=variable)) + geom_line(aes(colour=variable)) + geom_smooth(aes(colour=variable))
    p <- p + scale_y_continuous(limits=c(ymin,ymax))
    p <- p + ggtitle(title) + xlab("Year") + ylab("") + geom_hline(aes(yintercept=1))
    p <- p + scale_x_continuous(breaks = seq(min(indm$year), max(indm$year), 1))
    if (length(indicators) == 1){
        p <- p + theme(legend.position="none") 
    }
    return(p)
}


