###############################################################################
# EJ(20150610-20150703) JRC, IPSC, MAU <ernesto.jardim@jrc.ec.europa.eu>
# MSE to test the options given by MARE to EWG1509 - NWMed MAP
# Based on FLR (http://flr-project.org)
# and a4a (https://fishreg.jrc.ec.europa.eu/web/a4a)
###############################################################################
# IMPORTANT
# This code is made avaibale under a creative commons license BY-SA 4.0
# (http://creativecommons.org/licenses/by-sa/4.0/)
###############################################################################
#==============================================================================
# NOTE: The first intermediate year must be the last on the assessment so that
# the OM has information for the MPs assessment/intermediate year.
#==============================================================================

# THINGS YOU MUST DO:

# 1. Load your own stock and index data and set up plusgroup,
#     fbar range etc (if not already done): lines 40 - 54
# 2. Put in your own XSA control details: line 87
# 3. Put in your own Fmsy value: line 143

#==============================================================================
# libraries and constants
#==============================================================================
library(FLa4a)
library(FLash)
library(FLBRP)
library(FLAssess)
library(FLXSA)
library(ggplotFL)
#==============================================================================
# Read data
#==============================================================================
# hake 10
load("../../data/hke_gsa10.Rdata")

# rename your objects to stk (an FLStock) and idx (an FLIndices)
stk <- spe.stk
idx <- spe.idx

# Other stock setup if necessary
units(harvest(stk)) <- "f"
range(stk)["minfbar"] <- 1
range(stk)["maxfbar"] <- 4
range(stk)["plusgroup"] <- 6
stk <- setPlusGroup(stk, 6)

#==============================================================================
# Setup
#==============================================================================

# fixed variables
# Number of iterations
it <- 250
amx <- range(stk)["max"]
y0 <- range(stk)["minyear"] # initial data year
# data year, assessment year and initial projections year (also intermediate)
dy <- ay <- iy <- range(stk)["maxyear"]
ny <- 24 # number of years to project
fy <- iy + ny -1 # final year
vy <- ac(iy:fy) # year vector for projections
nsqy <- 3 # number of years to compute status quo metrics
trgy <- 2015
#==============================================================================
# Conditioning
#==============================================================================
# Fill in zeros with small values
catch.n(stk)[catch.n(stk)==0] <- 0.01
for (i in 1:length(idx)){
	index(idx[[i]])[index(idx[[i]])==0] <- 0.01
}

# Run your assessment.
# Here we demonstrate with XSA or a4a - use what you have

# a4a - if you have it
# fit <- sca(stk, FLIndices(idx), fit="assessment", fmodel=~s(age, k=5) + s(year, k=10))
#stk <- stk + fit
# Simulate to include uncertainty from assessment
# sstk <- stk + simulate(fit, it)

# Or using XSA - set own control
rage_SETTING <- 0 
qage_SETTING <- 5
xsa_control <- FLXSA.control(x=NULL, tol=1e-09, maxit=30, min.nse=0.3, fse=2,
							 rage=rage_SETTING, qage=qage_SETTING, shk.n=TRUE, shk.f=TRUE, shk.yrs=3, shk.ages=2,
							 window=100, tsrange=20, tspower=3, vpa=FALSE)
fit <- FLXSA(stk, idx, xsa_control)
stk <- stk + fit
# Cannot simulate from XSA, initial stock has no uncertainty from stock assessment
# All iterations are the same - underestimate uncertainty
sstk <- propagate(stk, it)

# Set up stock for projection
pstk <- stf(sstk, ny, 5, 5)
landings.n(pstk) <- propagate(landings.n(pstk), it)
discards.n(pstk) <- propagate(discards.n(pstk), it)
#------------------------------------------------------------------------------
# S/R 
#------------------------------------------------------------------------------
# Just use a fixed mean recruitment in the MSE - or use something else
# if you believe it
# S/R residuals
sr <- fmle(as.FLSR(stk, model="geomean")) # bevholt, ricker
sr.res <- window(rec(pstk), iy, fy)
sr.res[] <- sample(c(residuals(sr)), ny*it, replace=TRUE)
#------------------------------------------------------------------------------
# BRP
#------------------------------------------------------------------------------
# Get some reference points
rp <- brp(FLBRP(stk, sr))
#------------------------------------------------------------------------------
# Magic up some index variance for simulation in the future (pulled to 1st of January)
#------------------------------------------------------------------------------

for (i in 1:length(idx)){
	lst <- mcf(list(idx[[i]]@index, stock.n(stk)))
	idx.lq <- log(lst[[1]]/lst[[2]])
	idx.lq[is.infinite(idx.lq)] <- NA # fix zeros
	idx.qmu <- idx.qsig <- stock.n(iter(pstk,1))
	idx.qmu[] <- yearMeans(idx.lq)
	idx.qsig[] <- log((sqrt(yearVars(idx.lq))/yearMeans(idx.lq))^2 + 1)
	idx.q <- idx_temp <- FLQuant(NA, dimnames=dimnames(stock.n(pstk)))
	idx.q[,ac(y0:dy)] <- propagate(exp(idx.lq[,ac(y0:dy)]), it)
	idx.q[!is.na(idx.qmu)] <- rlnorm(it, idx.qmu[!is.na(idx.qmu)], idx.qsig[!is.na(idx.qmu)])
	plot(idx.q)
	idx_temp <- idx.q * stock.n(pstk)
	idx[[i]] <- FLIndex(index=idx_temp, index.q=idx.q)
	range(idx[[i]])[c("startf", "endf")] <- c(0, 0)
	plot(index(idx[[i]]))
}

#==============================================================================
# Management scenarios
#==============================================================================
#fmsy <- refpts(rp)["msy","harvest"]
#fmsy <- refpts(rp)["f0.1","harvest"]
# Insert your own Fmsy (or Fmsy proxy)
fmsy <- 0.198
# How do we come up with low B reference?
blim <- min(ssb(stk))
bpa <- blim*1.4
fupp <- 0.007801555 + 1.349401721*fmsy
flow <- 0.00296635 + 0.66021447*fmsy

#==============================================================================
# Run the MSE
#==============================================================================

# Pick your target
ftrg <- fupp
perfstats <- harvest(pstk)
perfstats[] <- NA
dt <- date()
sa <- list()
idx0 <- idx
# go fish
for(i in vy[-length(vy)]){
	gc()
	ay <- an(i)
	cat(i, " > ")
	vy0 <- 1:(ay-y0) # data years (positions vector)
	sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector)
	stk0 <- pstk[,vy0]
	# change M before assessing
	catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
	for (index_counter in 1:length(idx)){
		idx0[[index_counter]] <- idx[[index_counter]][,vy0]
		index(idx[[index_counter]])[,i] <- stock.n(pstk)[,i]*index.q(idx[[index_counter]])[,i]
	}
	
	# Run assessment using the simplest a4a model
    # Or use another if you think it's going to work
	# a4a
	 fit0 <- sca(stk0, FLIndices(idx0), fmodel=~factor(age) + factor(year))
	 stk0 <- stk0 + fit0
	# XSA
	#fit0 <- FLXSA(stk0, idx0, xsa_control)
	#stock.n(stk0) <- stock.n(fit0)
	#harvest(stk0) <- harvest(fit0)
	
	# fwd control
	fsq0 <- yearMeans(fbar(stk0)[,sqy])
	dnms <- list(iter=1:it, year=c(ay, ay+1, ay+1), c("min", "val", "max"))
	arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
	ftrg0 <- fsq0 - (fsq0-ftrg)/ifelse(trgy - ay < 1, 1, trgy - ay)
	arr0[,,"val"] <- c(fsq0, ftrg0, rep(NA, it))
	arr0[,,"min"] <- c(rep(NA, it), rep(NA, it), rep(bpa, it))
	arr0 <- aperm(arr0, c(2,3,1))
	ctrl <- fwdControl(data.frame(year=c(ay, ay+1, ay+1), quantity=c("f", "f", "ssb"),
								  val=NA))
	ctrl@trgtArray <- arr0
	stkTmp <- stf(stk0, 3)
	stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]),
				  sr.residuals.mult = TRUE)
	# OM proj
	ctrl@target <- ctrl@target[2,]
	ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
	ctrl@target["rel.year"] <- ay-1
	perfstats[1,ac(ay-1)] <- fbar(pstk)[,ac(ay-1)]
	perfstats[2,ac(ay-1)] <- fbar(stk0)[,ac(ay-1)]
	perfstats[3,ac(ay-1)] <- fsq0
	perfstats[4,ac(ay+1)] <- ftrg0
	# until 2015 keep fsq
	if(ay<2015) ctrl@trgtArray[,"val",] <- 1
	else ctrl@trgtArray[,"val",] <- c(fbar(stkTmp)[,ac(ay+1)])/c(fsq0)
	pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay+1)]),
				sr.residuals.mult = TRUE)
}

# pstk is the perceived stock
plot(window(pstk, end=2037))
# Ignore final year - not used
# Look at distribution of SSB in final year - 1 (assume projection has stabilised)
ssb(pstk)[,"2037"]
hist(ssb(pstk)[,"2037"])

# Proportion below Blim - are you less than 5%
sum(ssb(pstk)[,"2037"] < blim) / it

