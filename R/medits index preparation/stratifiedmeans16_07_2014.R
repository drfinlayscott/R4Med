#################################################################################
#    Copyright (C) <2014>  <Tristan Rouyer, Chato Osio>

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/
###################################################################################

# File updated on 16/072014 during EWG 14-09 by Chato Osio

###==================================================================================###
###
###             FUNCTION THAT PRODUCES INDICES FROM MEDITS DATA
###
###==================================================================================###
### T. ROUYER, 25/01/2013
### IFREMER SETE
###=============###
### ta, tb, tc, tt: tables from the MEDITS dataset
### species: code rubbin of the species selected
### year: year selected
### sex: sex selection NA,F,M,I
### size: size threshold for data aggregation (below/above). If 0 no account for size
### area: number of the GSA
###==================================================================================###





# Select size and number of years, already in std-indexMedits file, but here also if you want to run alone

size<-seq(min(TC$LENGTH_CLASS), max(TC$LENGTH_CLASS), 5) # size threshold for data aggregation (below/above). If 0 no account for size

year<-seq(min(TC$YEAR), max(TC$YEAR)) # number of years for which stratified numbers at length should be calculated
  
  
##==========================================##
  ## STRATUM AREA
  ## MEDITS PROTOCOL DOCUMENT (2012) P30
  ##==========================================##
# read in stratification table of the survey
stratification_scheme <- read.csv("MEDITS_Strata10122013.csv")
stratum<-stratification_scheme[stratification_scheme$AREA %in% gsa,]

    stratum.surface <- stratum$AREASTRATA
    stratum.number <- stratum$NUMBER_OF_THE_STRATUM
    total.surface <- sum(stratum.surface)
    ## WEIGHT OF EACH STRATUM
    w <- stratum.surface/total.surface
  
  ##==========================================##
  ##
  ##==========================================##
  ## COMPUTE INDICE
  ##==========================================##
  ## STORAGE
  ## indices calculated on TC
  
    
  ai <- matrix(NA,(length(size)-1),length(year))
  
  ## ABUNDANCE INDICE
  for (i in 1:length(year)) {  ## loop on years

    ## ta indices for the year
    ind.year <- which(TA$YEAR==year[i])
    TA.year <- TA[ind.year,]
    unique(TA.year$YEAR)

    ## TC indices for the year
    ind.year.TC <- which( TC$YEAR==year[i] & TC$HAUL_NUMBER %in% TA.year$HAUL_NUMBER )
    if (length(ind.year.TC)>0) {
      TC.year <- TC[ind.year.TC,]
    } else {
      print(paste('There is no daTA for ',species,' in TC in ',year[i],sep=''))
    }
    
    ## stratum for TA
    TA.strata <- TA.year$NUMBER_OF_THE_STRATUM

    ## sTAndardized number per hour This option can be turned on if needed but not tested much in current code
    #d <- rep(NA,nrow(TC.year))
    #wf <- rep(NA,nrow(TC.year))
    #for (k in 1:nrow(TC.year)) {
     # d[k] <- TA.year$HAUL_DURATION[which(TA.year$HAUL_NUMBER==TC.year$HAUL_NUMBER[k])]/60 # standardized by duration of the haul
     # wf[k] <- TC.year$WEIGHT_OF_THE_FRACTION[k]/TC.year$WEIGHT_OF_THE_SAMPLE_MEASURED[k]
    #}
    #n <- TC.year$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE*wf/d

    n <- TC.year$NBLEN.raised / TC.year$swept # these where already raised by the prop between weight of fraction and weight of sample measured
    ## variables
    xbar <- matrix(0,(length(size)-1),length(stratum.number))
    ## loop on straTA
    for (j in 1:length(stratum.number)) {
      ## get indices of the hauls from TA within the stratum j
      ind.strata <- which(TA.strata==stratum.number[j])
      ## get indices of the hauls from TC within the stratum j
      ind.strata.TC <- which(TC.year$HAUL_NUMBER %in% TA.year$HAUL_NUMBER[ind.strata])
      ## hauls from TA: used to place data in x
      TA.hauls <- TA.year$HAUL_NUMBER[order(TA.year$HAUL_NUMBER)]
      ## matrix size*haul
      x <- matrix(0,(length(size)-1),length(TA.hauls))
      ## loop on sizes
      for (siz in 1:(length(size)-1)) {
        ## index of individuals within the right size class
        ind.ind <- which(TC.year$LENGTH_CLASS[ind.strata.TC]>=size[siz]
                         & TC.year$LENGTH_CLASS[ind.strata.TC]<size[(siz+1)])
        ## if there are any
        if (length(ind.ind)>0) {
          ## get hauls
          nh <- unique(TC.year$HAUL_NUMBER[ind.strata.TC[ind.ind]])
          for (ih in 1:length(nh)) {
            ## indices of hauls within the stratum
            ind.nh <- which(TC.year$HAUL_NUMBER[ind.strata.TC[ind.ind]]==nh[ih])
            ## indices of the haul: for placing on the matrix "x"
            ind.haul <- which(TA.hauls==nh[ih])
            ## sum over hauls within the stratum
            if (length(ind.nh)>0) {
              x[siz,ind.haul] <- sum(n[ind.strata.TC[ind.ind[ind.nh]]])
            }
          }
        }
        ## average numbers per stratum and age class
        xbar[siz,j] <- mean(x[siz,])*100
      }
    }
    for (siz in 1:(length(size)-1)) {
      ai[siz,i] <- sum(xbar[siz,] * w)
    }
    rownames(ai) <- c(paste(size[-length(size)],size[-1],sep='_'))
    colnames(ai) <- paste(year,sep='')
  }
  stratified_N_atlength <- list(ai)

# convert the matrix to data frame  
stratified_N_atlength<- as.data.frame(stratified_N_atlength)  
stratified_N_atlength<-cbind(lengths = dimnames(stratified_N_atlength)[[1]], stratified_N_atlength)
dimnames(stratified_N_atlength)[[1]]<- 1:(length(size)-1)
  stratified_N_atlength$lengths<-1:(length(size)-1)
# reshape to long version from wide so that it is easier to plot
stratified_N_atlengthR <- melt(stratified_N_atlength, id= "lengths", na.rm = TRUE)
  stratified_N_atlengthR<-droplevels(stratified_N_atlengthR)
  
  
  b <- ggplot(stratified_N_atlengthR, aes(y=value, x=lengths, color=variable))+ geom_bar(stat= "identity")+facet_wrap(~variable)

# Save plot and save file
  ggsave(b, file=paste("F:/EWG14_09/stratified_N_atlength", sp , gsa ,".png", sep="") ) 
  write.csv(stratified_N_atlength, file=paste("F:/EWG14_09/stratified_N_atlength",sp,gsa,".csv", sep=""))
  
  
  

