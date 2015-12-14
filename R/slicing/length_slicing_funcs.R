# Age slicing functions
# A collection of functions initially developed by:
# Laurence & Alexander Kell (2010)
# Valerio Bartolino, Chato Osio, Graham Pilling & Finlay Scott (2010)
# With further development by Finlay Scott <finlay.scott@cefas.co.uk>, 
# Chato Osio <giacomo-chato.osio@jrc.ec.europa.eu> and 
# Max Cardinale <massimiliano.cardinale@slu.se> (2011)
# Now collected under one roof

# Copyright Laurence Kell, Alexander Kell, Finlay Scott,
# Valerio Bartolino, Chato Osio, Max Cardinale (2013)

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

# von Bertalanffy growth function
VB <- function(age,Linf,K,t0=0)
{
    L <- Linf * (1 - exp(-K * (age - t0)))
    return(L)
}

# Inverse von Bertalanffy growth function
rVB <- function(L,Linf,K,t0=0)
{
    age <- t0 - 1/K * log(1 - L/Linf)
    return(age)
}

#--------------------------------------------------------------------------

# Deterministic age slicing (knife-edge)

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

knife_edge=function(dat,vB,timing=0.5,plusGroup=30,minage=1)
{
    len <- dat[,1]
    n <- dat[,2]
    t0 <- vB["t0"]
    K <- vB["K"]
    Linf <- vB["Linf"]
    # Calculate the expected age at length adjusted to beginning of year
    # Implement the inverse vonB growth equation with a few checks:
    # i.e. we need to ensure that len is < Linf
    #      that all ages > plusGroup get put in the plus group
    #      that all ages < minage get put in the minage
    age <- pmax(pmin(floor(t0-log(1-pmin(len/Linf,.9999999))/K+timing),plusGroup),minage)
    # Calculate frequencies
    res <- aggregate(n, list(age=age), sum)
    out <- data.frame(age=minage:plusGroup, data=0)
    out[out$age %in% res$age,"data"] <- res$x
    return(out)
}

#--------------------------------------------------------------------------

# Length slicing based upon example in FAO Fish. Tech. Man., p182
# Introduction to Tropical Fish Stock Assessment - Part 1: Manual
# Per Sparre and Siebren C. Venema
# FAO FISHERIES TECHNICAL PAPER 306/1 Rev. 2
# Food and Agriculture Organization of the United Nations
# Rome, 1998
# http://www.fao.org/docrep/w5449e/w5449e00.htm

# The method is based on calculating the proportion of each length
# bin width that should go into each age class.
# It uses a series of functions to do this (and could conceivably be
# condensed).
# A key function is get_prop_matrix() which returns a matrix of the
# proportions. This matrix can then be swept over to calculate
# numbers at age.

# Calculates the amount of overlap between the length classes (L1, L2)
# in the data and the length at ages calculated from the von Bertalanffy
# growth equation (VBL1, VBL2).
# i.e. the proportion of L1-L2 that is inside VBL1-VBL2
# There are six possible options.
# (There must be a more sophisticated way of doing this)
get_prop <- function(L1, L2, Latage1, Latage2)
{
    # Checks
    if ((L1 < 0) | (L2 < 0) | (Latage1 < 0) | (Latage2 < 0))
        stop("Lengths can not be less than 0.")    
    if ((L1 > L2) | (Latage1 > Latage2)){
        stop("First length is greater than second length.")
    }
    prop_out <- NA
    if (L1 < Latage1 & L2 < Latage2) prop_out <- 0 # outside to the left
    if (L1 > Latage2 & L2 > Latage2) prop_out <- 0 # outside to the right
    if (L1 >= Latage1 & L2 <= Latage2) prop_out <- 1 # contained inside
    if (L1 <= Latage1 & L2 <= Latage2 & L2 >= Latage1) prop_out <- (L2 - Latage1) / (L2-L1) # overlaps to the left
    if (L1 >= Latage1 & L1 <= Latage2 & L2 >=Latage2) prop_out <- (Latage2 - L1) / (L2-L1) # overlaps the right
    if (L1 <= Latage1 & L2 >= Latage2) prop_out <- (Latage2-Latage1) / (L2-L1) # L1-2 contains Latage1-2
    return (prop_out)
}

# Given VBparams, generate a vector of Lbins start values, from 0 to almost Linf-binwidth
# This allows the final bin to be close to Linf
calc_Lbin_start <- function(Linf=60, binwidth=5)
{
    Lbins <- seq(from=0,to=Linf,by=binwidth)
    # If Lbins larger or equal to Linf, drop them
    if (any(Lbins>=Linf))
        Lbins <- Lbins[-which(Lbins>=Linf)]
    return(Lbins)
}

# This function does all the hard work
get_prop_matrix <- function(K=0.2,Linf=60,t0=0,binwidth=5,plusgroup=NA)
{  
    # Start of each Lbin (from 0 to Linf-binwidth)
    Lbins <- calc_Lbin_start(Linf,binwidth)
    # Tack on end of final Lbin = almost Linf * 0.99
    # This is a hack to get around problems with dealing age at Linf
    # Sure we can do something more sophistacted than this for final Lbin (Linf-binwidth to Linf)
    Lbins <- c(Lbins,Linf*0.99999999)
    maxage <- floor(rVB(max(Lbins),Linf,K,t0))
    Latage <- VB(0:(maxage+1),Linf,K,t0)

    prop_matrix <- matrix(0,ncol=length(0:maxage), nrow=length(Lbins)-1)
    dimnames(prop_matrix) <- list(start_L=Lbins[-length(Lbins)],age=0:maxage)
    # loop over all elements and get the propotion of age at length
    # i.e. how each observed length group is spread over the ages
    for (j in 1:(dim(prop_matrix)[2])){#ages
        for (i in 1:(dim(prop_matrix)[1])){   #lengths
            prop_matrix[i,j] <- get_prop(Lbins[i], Lbins[i+1], Latage[j], Latage[j+1])
    }}

    # Temporary fix to deal with issues in the first age
    # Force rows to have sums of 1
    prop_matrix <- sweep(prop_matrix,1,apply(prop_matrix,1,sum),"/")

    # Problem with low Ls
    # Chop off rows that are NaN (low Ls) - problem with dividing by 0 and minimum age
    # fix by using t0?
    nanrows <- which(is.nan(prop_matrix[,1]))
    # set to 1 in lowest age and 0 and for all others
    prop_matrix[nanrows,] <- 0
    prop_matrix[nanrows,1] <- 1

    # Sort out plusgroup
    if (!is.na(plusgroup)){
        temp <- matrix(0,ncol=length(0:plusgroup), nrow=length(Lbins)-1)
        dimnames(temp) <- list(start_L=Lbins[-length(Lbins)],age=c(0:(plusgroup-1),paste(plusgroup,"+",sep="")))
        if (plusgroup > maxage)
            temp[,1:ncol(prop_matrix)] <- prop_matrix[]
        if (plusgroup <= maxage){
            temp[,1:plusgroup] <- prop_matrix[,1:plusgroup] # up to plus group
            temp[,plusgroup+1] <- apply(prop_matrix[,(plusgroup+1):(maxage+1)],1,sum) # the plus group
        } # same size or smaller
        prop_matrix <- temp
    }
    return(prop_matrix)
}

# Takes the length data and drops it into contiguous length classes
# (may not be contiguous due to empty bins)
# Warning: all observations with lengths > Linf are discarded!
expand_ldat <- function(ldat,pmat,force_binwidth=TRUE, warn=TRUE)
{
   
    #browser()
    if (warn==TRUE & (max(ldat[,1]) >= max(as.numeric(dimnames(pmat)$start_L))))
        warning("Max length in ldat is greater than max length in prop matrix. Removing those observations.")
    binwidth <- as.numeric(dimnames(pmat)$start_L)[2] - as.numeric(dimnames(pmat)$start_L)[1]
    ldatx <- array(0,dim= c(dim(pmat)[1],1),dimnames=list(length=rownames(pmat)))
    # Some of the length classes are not right, given the LENGTH_CLASS_CODE
    # So force round to nearest length class (bit dodgy)
    if(force_binwidth)
        ldat[,1] <- round(ldat[,1] / binwidth)*binwidth
    # sum over length classes (may be multiple entries) also has handy effect of sorting them
    ldat <- tapply(ldat[,2],ldat[,1],sum)
    #ldatx[as.numeric(rownames(ldatx)) %in% names(ldat),] <- ldat
    
    ldatx[as.numeric(rownames(ldatx)) %in% names(ldat),] <- ldat[names(ldat) %in% as.numeric(rownames(ldatx))]
    
    return(ldatx)
}

#slice <- function(sub,pmat)
#{
#    #cat("Slicing: ",unique(sub$YEAR), "\n")
#	ldat <- sub[,c("LENGTH_CLASS", "NBLEN.raised")]
#	ldatx <- expand_ldat(ldat,pmat)
#	#****** Warning - bad hack******
#    if (max(ldat[,1]) >= max(as.numeric(dimnames(prop_mat)$start_L))){
#        print("Warning! Length class >= Linf. Sample removed")
#        # Chop out counts where Length_class > Linf (dodgy)
#        ldat <- ldat[(ldat[,1] < VBparams[[sex]]["Linf"]),]
#    }
	#********************************
#	temp_natage <- apply(sweep(pmat,1,ldatx,"*"),2,sum)
#	return(data.frame(YEAR=unique(sub$YEAR),
#						HAUL_NUMBER=unique(sub$HAUL_NUMBER),
#						SEX=unique(sub$SEX),
#						AGE=0:pls.grp,
#						NATAGE=temp_natage))
#}


