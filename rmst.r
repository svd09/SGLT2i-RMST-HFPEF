# R script for estimating IPD from digitized KM figures and from 
library (survRM2)
library(IPDfromKM)
library (metaRMST)

# open the active arm KM curve
db=file.choose()
active <- read.table(db,
                     sep="\t", header=FALSE)

# open the placebo arm KM curve
db=file.choose()
comparator <- read.table(db,
                         sep="\t", header=FALSE)

# open the file containing time stamped number at risk for both arms 
db=file.choose()
no_at_risk<- read.table(db,header=FALSE)

# invert the event rate to survival free from event
active$V2 = 100-active$V2
comparator$V2 = 100-comparator$V2

### Get data from the sample dataset=======================
trisk1 <- as.integer(no_at_risk$V1)
trisk <- na.omit(trisk1)
nrisk_comparator1 <- as.integer(no_at_risk$V2)
nrisk_comparator <- na.omit(nrisk_comparator1)
nrisk_active1 <- as.integer(no_at_risk$V3)
nrisk_active <- na.omit(nrisk_active1)

### Estimate the IPD for the comparator therapy treatment group ====================
pre_comparator <- preprocess(dat=comparator, trisk=trisk,nrisk=nrisk_comparator,maxy=100)
est_comparator <- getIPD(prep=pre_comparator,armID=0,tot.events=NULL)

### Estimate the IPD for the active treatment group ====================
pre_active <- preprocess(dat=active, trisk=trisk,nrisk=nrisk_active,maxy=100)
est_active <- getIPD(prep=pre_active,armID=1,tot.events=NULL)

#combine both arms into one dataset
###for deliver
total_ipd1 <- rbind(est_comparator$IPD,est_active$IPD)
###for emperor
total_ipd2 <- rbind(est_comparator$IPD,est_active$IPD) 


### COMBINE TRIALS
deliver <- total_ipd1 # then run load trial line 5->
emperor <- total_ipd2 # then run load trial line 5->

## LABEL TRIALS
deliver$trialID = 1
emperor$trialID= 2

# combine trial ipd
overall <- rbind (deliver, emperor)

## RENAME variables for the function
names(overall) [1] <- "Time"
names(overall) [2] <- "Event"
names(overall) [3] <- "Arm"
names(overall) [4] <- "trialID"

### metaRMST from saved ipd file
### metaRMST for Primary Outcomes
obj1 <- RMSTcurves(overall, time_horizons=c(1,12,24,35), tmax=36, nboot=500,  MA_mvma = FALSE, MA_mvma_boot = FALSE,
                   MA_uni = T, MA_uni_flex = F)
pl1 <- RMSTplot(obj1, xlim=c(0,36), ylim=c(-0.1,2), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Primary Outcome", col = c("red", "blue"), trial_legend=F, MA_legend = FALSE)

metaRMSTD(overall, time_horizons=c(1,12,24,35), MA_method = "uni", nboot = 500)

### metaRMST for HHF
obj2 <- RMSTcurves(overall, time_horizons=c(1,12,24,35), tmax=36, nboot=500,  MA_mvma = FALSE, MA_mvma_boot = FALSE,
                   MA_uni = T, MA_uni_flex = F)
pl2 <- RMSTplot(obj2, xlim=c(0,36), ylim=c(-0.1,2), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Heart Failure Hospitalization", col = c("red", "blue"), trial_legend=F, MA_legend = FALSE)

metaRMSTD(overall, time_horizons=c(1,12,24,35), MA_method = "uni", nboot = 500)

###CV Mortality
obj3 <- RMSTcurves(overall, time_horizons=c(1,12,24,35), tmax=36, nboot=500,  MA_mvma = FALSE, MA_mvma_boot = FALSE,
                  MA_uni = T, MA_uni_flex = F)
pl3 <- RMSTplot(obj3, xlim=c(0,36), ylim=c(-0.1,1), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Cardiovascular Death", col = c("red", "blue"), trial_legend= F, MA_legend = FALSE)

metaRMSTD(overall, time_horizons=c(1,12,24,35), MA_method = "uni", nboot = 500)

### metaRMST for All Cause Death
obj4 <- RMSTcurves(overall, time_horizons=c(1,12,24,35), tmax=36, nboot=500,  MA_mvma = FALSE, MA_mvma_boot = FALSE,
                   MA_uni = T, MA_uni_flex = F)
pl4 <- RMSTplot(obj4, xlim=c(0,36), ylim=c(-0.5,1), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="All Cause Death", col = c("red", "blue"), trial_legend=F, MA_legend = FALSE)

metaRMSTD(overall, time_horizons=c(1,12,24,35), MA_method = "uni", nboot = 500)

###Graph
par(mfrow=c(2,2))
plot(RMSTplot(obj1, xlim=c(0,36), ylim=c(-0.1,2), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Primary Outcome", col = c("red", "blue"), trial_legend=F, MA_legend = FALSE))
plot(RMSTplot(obj2, xlim=c(0,36), ylim=c(-0.1,2), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Heart Failure Hospitalization", col = c("red", "blue"), trial_legend=F, MA_legend = FALSE))
plot(RMSTplot(obj3, xlim=c(0,36), ylim=c(-0.1,1), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Cardiovascular Death", col = c("red", "blue"), trial_legend= F, MA_legend = FALSE))
plot(RMSTplot(obj4, xlim=c(0,36), ylim=c(-0.5,1), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="All Cause Death", col = c("red", "blue"), trial_legend=F, MA_legend = FALSE))
