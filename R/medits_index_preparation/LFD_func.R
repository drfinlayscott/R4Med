#========================================================================================####
#    Title:   LFD_fun.R
#    Release: 0.1
# 
#    Description: FUNCTION THAT PRODUCES STANDARDIZED LFD FROM MEDITS DATA
#    Authors:  Matteo MURENU, Alessandro MANNINI, Tristan ROUYER, Chato OSIO, Finlay SCOTT
#    Date: created on June 2015, ISPRA, EWG 15-06 meeting
#    Updates:
#      on August 2015 during EWG 15-11  by Matteo Murenu
#      on December 2015 by Finlay Scott
#      on ____ during EWG 15-xx  by _____
#
#========================================================================================###
#
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
#
###==================================================================================###


lfd <- function(TA, TB, TC, medstra, sex, len.unit, tabdir, plotdir, plots=TRUE){

    year<-sort(unique(TC$YEAR)) # number of years for which stratified numbers at length should be calculated
    # Sum TC$NBLEN.raised by CODESTRATA, LENGTH_CLASS and YEAR
    lfd.y<-tapply(TC$NBLEN.raised,list(TC$CODESTRATA,TC$LENGTH_CLASS,TC$YEAR),sum) # LFD raised array
    # Sum TA$swept by Year and codeStrata - what is swept?
    # From calling script. Total area swept TAn$swept<-(TAn$DISTANCE *((TAn$WING_OPENING)/10000))/1000 #Swept area in Km^2
    swept.y<-tapply(TA$SWEPT,list(TA$CODESTRATA,TA$YEAR),sum) # swept array
    # Do we need medstra here? Info is contained in TA after merge?
    stra.sur<-tapply(medstra[medstra$AREA==gsa[j],]$AREASTRATA,list(medstra[medstra$AREA==gsa[j],]$CODESTRATA),sum) # LFD raised array
    stra.sur.t<-sum(medstra[medstra$AREA==gsa[j],]$AREASTRATA)
    # Fraction of XXX by CODESTRATA (depth)
    # Getting area and depth - as need to raise numbers over the area
    stra.fac<-stra.sur/stra.sur.t
    country.cod<-unique(TC$COUNTRY)
    month.cod<-floor(mean(TA$MONTH))

    # Calculate indices
    # Use sweep?
    lfd.y2<-lfd.y
    for (i in 1:length(year)) {
        mth<-(match(rownames(lfd.y[,,i]),rownames(swept.y)))
        # total numbers / swept area
        lfd.y2[1:dim(lfd.y)[1],,i]<-lfd.y[1:dim(lfd.y)[1],,i]/swept.y[mth,as.character(year[[i]])] # standardized LFD array
    }

    # reshape and store calculation ####
    # Why multiply by stra.fac - proportion of total area by depth?
    # Depends on what medstra is - if actual area and depth (not swept) then it's OK.
    # lfd.y2 is simply numbers by haul. If area of depth 1 is bigger than area of depth 2 then we raise the numbers in depth 1 by bigger number
    lfd.std<-lfd.y2[1:dim(lfd.y2)[1],,]*as.vector(stra.fac[1:dim(lfd.y2)[1]])
    lfd.std.t<-melt(lfd.std)
    names(lfd.std.t)<-c("CODESTRATA","len","year","value")

    # reshape the LFD array
    # Sum over year adn length (i.e. sum CODESTRATA)
    stratified.N<-tapply(lfd.std.t$value,list(lfd.std.t$year,lfd.std.t$len),sum,na.rm=T)
    stratified.N <- melt(stratified.N, id= "lengths")

    # add some extra information #####
    names(stratified.N)<-c("year","len","value")
    #stratified.N$country<-country.cod
    stratified.N$month<-month.cod
    stratified.N$GSA<-gsa[j]
    stratified.N$species<-TC$GENSPE[1]
    stratified.N$sex<-sex
    #stratified.N<-stratified.N[c("country","GSA","year","species","month","len","value")]
    stratified.N<-stratified.N[c("GSA","year","species","month","len","value")]
    stratified.N <- stratified.N[order(stratified.N$year,stratified.N$len),] 

    # duplicate row for NBLEN>1
    tmp <- TC[, c("AREA", "YEAR","HAUL_NUMBER","GENSPE","SEX","LENGTH_CLASS","NBLEN.raised")]
    dm<-tmp[rep(row.names(tmp),round(tmp$NBLEN.raised)),]  # duplicate rows
    ## check # numbers should be the same 
    # We avoid to put here an error warning. We just check numbers.
    cat(paste0(round(sum(tmp$NBLEN.raised)), " vs ", dim(dm)[1], "      numbers should (almost) be the same "), "\n")

    # Write data
    species_file_name <- gsub(" ", "_", TC$GENSPE[1])
    write.csv(stratified.N,row.names=F, file= paste0(tabdir,"/stratified_Nlen_",species_file_name,"_",sex,"_GSA",TC$AREA[1],".csv"))

    # Generate plots
    if (plots == TRUE){
        if  (len.unit=="cm") {
            p <- ggplot(stratified.N, aes(y=value, x=len/10))
            p2 <- ggplot(dm, aes(y=LENGTH_CLASS/10, x=factor(YEAR)))
        }
        else {
            p <- ggplot(stratified.N, aes(y=value, x=len))
            p2 <- ggplot(dm, aes(y=LENGTH_CLASS, x=factor(YEAR)))
        } 
        # set the graph title
        if  (sex=="all") {
            tt <- ggtitle(paste0(TC$GENSPE[1],"    GSA",gsa[j],"\n"))
        }
        else {
            tt <- ggtitle(paste0(TC$GENSPE[1]," (",sex,")   GSA",gsa[j],"\n"))
        } 
        # plot the histogram
        p <- p + tt + geom_bar(stat= "identity")+facet_wrap(~year)+
            labs(x=paste0("\nLength (",len.unit,")"),y=expression("n/km"^"2"))+
            theme(legend.position="none")
        print(p)
        # plot the boxplot
        p2 <- p2 + tt + geom_boxplot() +
            labs(y=paste0("\nLength (",len.unit,")\n"),x=NULL)+
            guides(fill=FALSE) + theme_bw() +
            theme(plot.title = element_text(colour = "red", size=rel(1.2)))+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p2)
        ggsave(p, file=paste0(plotdir,"/stratified_Nlen_",species_file_name,"_",sex,"_GSA",TC$AREA[1],".png")) 
        ggsave(p2, file=paste0(plotdir,"/bxplen_",species_file_name,"_",sex,"_GSA",TC$AREA[1],".png")) 
    }
}

