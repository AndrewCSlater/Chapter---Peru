library(dplyr)
library(spOccupancy)
library(coda)

plogis(.79)

######################################################
## Probability of Detection ~ Probabilty of Occurrence
load("models/DETECTION_ONLY_50k10k1k2chReps.14.Rdata")
# INTERCEPT_Reps.14
# SpLF3_DetOcc_Reps.14_buffer1Km

out = out.det
rm(out.sfMsPGOcc, out.intercept, out.det); gc()

y = out$y # 294 Species in 414 Stations
y2 = as.data.frame(rowSums(y, dims = 2, na.rm = T))
y2 = y2[colSums(y2)>0]
y3 = y2; y3[y3>0] = 1
summary(colSums(y3))

#### Species Accumulation Curves - How many sites needed to collect given number of Species ####
acum = vegan::specaccum(comm = y3, method = "random", permutations = 1000)
summary(acum)
which(acum$richness>0.95*max(acum$richness))[1] # 95% of species were found within 124 sites
which(acum$richness>0.75*max(acum$richness))[1] # 75% of species were found within 43 sites
plot(acum$sites, acum$richness,type = "l", lwd = 2,
     xlab="Number of Sites",
     ylab="Species Richness")
################################

#### Rarefacton Curve - How many samples (individual sightings) needed to collect given number of species #### 
y4 = as.data.frame(rowSums(y2))
rr = vegan::rrarefy(y4, sample = 40)
dr = vegan::drarefy(y4, sample = 40)
r = vegan::rarefy(t(y4), sample = 400)
vegan::rarecurve(t(y4), sample = c(627, 2130, 5070))
rc = vegan::rarecurve(t(y4))
rcc = rc[[1]]
which(rcc>0.95*max(rcc))[1] # 5070 of 6308 samples needed to gather 95% of all species
which(rcc>0.75*max(rcc))[1] # 2130 of 6308 samples needed to gather 75% of all species
which(rcc>0.5*max(rcc))[1] # 627 of 6308 samples needed to gather 50% of all species
#####################################################################################


summary(out)
occ = summary(out$beta.samples) ## Mean Probability of Occurrence
meanOcc = occ[[1]][,1]
det = summary(out$alpha.samples) ## Mean probability of Detection
meanDet = det[[1]][,1]

mod = lm(meanDet~meanOcc)
summary(mod)
plot(meanDet~meanOcc, xlab = "Mean Probability of Occurrence", ylab = "Mean Probability of Detection", main = "Bird Species - Peruvian Amazon")
abline(mod)

logOcc = plogis(meanOcc)
logDet = plogis(meanDet)
logmod = lm(logDet ~ logOcc)
summary(logmod)
plot(logDet~logOcc, xlab = "Mean Probability of Occurrence", ylab = "Mean Probability of Detection", main = "Bird Species - Peruvian Amazon")
abline(logmod)

#### Rarefacton Curve based on Probability of Detection using POINT and NET ####
ldet = as.data.frame(meanDet)
detInt = ldet[1:(nrow(ldet)/4),]
detPoint = ldet[((nrow(ldet)/4)+1):(nrow(ldet)/4*2),]
ldet2 = rbind(detInt, detPoint)
ldet2 = rbind(ldet2, colSums(ldet2))
# ldet2 = ldet2 * meanOcc ## Multiply Probability of Detection by Probability of Occupancy
ldet3 = plogis(ldet2)
ldet3 = ldet3 * logOcc
rownames(ldet3)[3] = "Point"
ldet3 = round(ldet3*1000, digits = 0)
ldet3 = ldet3[-2,]

rc = vegan::rarecurve(ldet3)
rcc = rc[[1]]
which(rcc>0.95*max(rcc))[1] # 1715 samples needed to observe 95% of all species
which(rcc>0.75*max(rcc))[1] # 682 samples needed to observe 75% of all species
which(rcc>0.5*max(rcc))[1] # 285 samples needed to observe 50% of all species
vegan::rarecurve(t(ldet3[2,]), sample = 1000)
vegan::rarecurve(t(ldet3[1,]), sample = 1000)
vegan::rarecurve(ldet3, sample = 1000)
#####################################################################################
ldet4 = round(scales::rescale(ldet3, to = c(1,max(ldet3))), digits = 0)


par(mfrow = c(1,2))
par(mar = c(4,4,1,1))
hist(logOcc, breaks = 15, xlab = "Probability of Occurrence", main = "")
hist(logDet, breaks = 15, xlab = "Probability of Detection", main = "")

# load library tidyverse, gridExtra and ggExtra
library(tidyverse)
library(ggExtra)
library(gridExtra)
# set theme
theme_set(theme_bw(12))

# create sample data frame
sample_data <- data.frame(logOcc, logDet)

slope = lm(logDet~logOcc)
intercept = slope$coefficients[1]
slope = slope$coefficients[2]


# create scatter plot using ggplot() function
plot <- ggplot(sample_data, aes(x=logOcc, y=logDet))+
  geom_point()+
  labs(y = "Probability of Detection", x = "Probability of Occurrence") +
  geom_abline(intercept = intercept, slope = slope, color="red", 
              linetype="dashed") +
  theme(legend.position="none")

# use ggMarginal function to create
# marginal histogram, boxplot and density plot

plot1 <- ggMarginal(plot, type="histogram")
plot2 <- ggMarginal(plot, type="boxplot")
plot3 <- ggMarginal(plot, type="density")

# combine plots in a grid
grid.arrange( plot1, plot2, plot3, ncol=3)

plot1; plot2; plot3
#######################################################

sig_sp_band = c() ## To store the species names that have significant relationships with band reflectance

load("models/OCCURRENCE_StationYear_allReps_TasseledCapBGW.Rdata")
load("models/OCCURRENCE_StationYear_12Reps_NBR_mean_500.Rdata")
out = out.occ
rm(out.occ); gc()

summary(out, level = "community")
summary(out)
waicOcc(out)

### Get Influence of DETECTION (Alpha) variables 
alp = out$alpha.samples
write.csv(bet, "model_outputs/ _alpha_detection.csv", row.names = F)

### How many species have a SIGNIFICANT relationship with the variable: ie how many have a 95%CI that excludes zero?
ncol(alp)/4
sigYN = c()                         ## Create YN significance vector
sigVAL = c()                        ## Create Influence Value vector (for if variable significantly influences a species)
aa1 = as.data.frame(alp[,285:568])  ## Subset df into all species with 1 variable
aa1 = as.data.frame(alp[,569:854])
for (j in 1:ncol(aa1)) {            ## For each species (column)
  i = quantile(aa1[,j], c(0.025,0.975))  ## Get the 2.5 & 97.5 quants
  if (i[1]<0 & i[2]>0) {ii = 0           ## If the quants surround zero Give the Value of 0
  } else {ii = 1}                       ## Else (all 95% is either above or below zero) give a value of 1 
sigYN = append(sigYN, ii)               ## Add that 1 or 0 value to the yes no vector
if (ii==1) {sigVAL = append(sigVAL, colMeans(aa1[j]))}
}  
sum(sigYN)                         ## The number of species with Significant relationships
colnames(aa1)[(sigYN)==1]          ## Find species names with significant relationships
length(which(sigVAL<0)); length(which(sigVAL>0)) # How many sp with significant relationships have Negative & Positive relationships 
hist(sigVAL)                       ## Histogram of Significant Mean Values
hist(colMeans(aa1))                ## Histogram of ALL mean values


### Get Influence of OCCURRENCE (Beta) variables 
bet = out$beta.samples
# write.csv(bet, "model_outputs/Deforest_5km_beta_occurrence.csv", row.names = F)

### How many species have a SIGNIFICANT relationship with the variable: ie how many have a 95%CI that excludes zero?
sigYN = c()                         ## Create YN significance vector
sigVAL = c()                        ## Create Influence Value vector (for if variable significantly influences a species)
# ncol(bet)/4
bb1 = as.data.frame(bet[,359:716])  ## Subset df into all species with 1 variable
#bb1 = as.data.frame(bet[,717:1074])
#bb1 = as.data.frame(bet[,1075:1432])
for (j in 1:ncol(bb1)) {            ## For each species (column)
  i = quantile(bb1[,j], c(0.025,0.975))  ## Get the 2.5 & 97.5 quants
  if (i[1]<0 & i[2]>0) {ii = 0           ## If the quants surround zero Give the Value of 0
  } else {ii = 1}                       ## Else (all 95% is either above or below zero) give a value of 1 
  sigYN = append(sigYN, ii)               ## Add that 1 or 0 value to the yes no vector
  if (ii==1) {sigVAL = append(sigVAL, colMeans(bb1[j]))}
}  
sum(sigYN)                         ## The number of species with Significant relationships
sig_sp = colnames(bb1)[(sigYN)==1]          ## Find species names with significant relationships
length(which(sigVAL<0)); length(which(sigVAL>0)) # How many sp with significant relationships have Negative & Positive relationships 
hist(sigVAL)                       ## Histogram of Significant Mean Values
hist(colMeans(bb1))                ## Histogram of ALL mean values

sig_sp = sub(".*-", replacement = "", x = sig_sp)
sig_sp = as.data.frame(cbind(sig_sp, sigVAL))
sig_sp$sigVAL[sig_sp$sigVAL>0] = "pos"
sig_sp$sigVAL[sig_sp$sigVAL<0] = "neg"
sig_sp$band = "NBR"

sig_sp_band = bind_rows(sig_sp_band, sig_sp)
write.csv(sig_sp_band, "data/Peru_6_Significant_Band-Species_SingleModels.csv", row.names = F)

#################################################################
#######  Does Band Significance vary with Species Trait?  #######
trait = read.csv("data/BIRDS_3.5_Bird_TRAITS.csv")
length(unique(sig_sp_band$sig_sp)) ## 152
sig_trait = match(x = sig_sp_band$sig_sp, table = trait$species)
sig_sp_band_trait = cbind(sig_sp_band, trait[sig_trait, c(4:9)])
library(ggplot2)
# bar plot, with each bar representing 100%
ggplot(sig_sp_band_trait, aes(x = Primary.Lifestyle, fill = band)) + 
  geom_bar(position = "stack") +
  labs(y = "Proportion")
################################################################
################################################################

plot(out$beta.samples[, 1:6], density = FALSE) # Plot Convergence of 1st 6 species
a=geweke.diag(out$beta.samples)#[1:358, ])
hist(a$z)
summary(a$z)
sum(a$z < 2 & a$z > -2, na.rm = T) / length(a$z) # 84% of species have converged well using Geweke diagnostic
summary(out$lambda.samples) # Lambda is the spatial influence per species
lam = out$lambda.samples
write.csv(lam, "model_outputs/Deforest_5km_lambda_3spatial.csv", row.names = F)
sp.prev = as.matrix(rowSums(out$y, na.rm = T))
rm(a, bet, lam); gc()

################################################
### AUC and TSS values
################################################
out = 
rm(rep14); gc()
y.rep.samples <- fitted(out)$y.rep.samples
save(y.rep.samples, file = "model_outputs/Deforest_5km_Fitted_400samples.Rdata")
str(y.rep.samples)
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out$n.post * out$n.chains
N <- nrow(out$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  # print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- pROC::auc(response = c(out$y[i, , ]), predictor = c(y.rep.samples[j, i, , ]))
  } # i (species)
} # j (iteration)
colnames(auc.vals) = rownames(out$y)
write.csv(auc.vals, "model_outputs/Deforest_5km_AUC.csv", row.names = F)
# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975)))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
mean(auc.quants[,2])
sd(auc.quants[,2])
range(auc.quants[,2])
length(which(auc.quants[,2]>0.7))
# Calculate species-specific confusion matrices ---------------------------
# Number of non-missing y values for each species
n.obs <- sum(!is.na(c(out$y[1, , ])))
# True positives
tp.samples <- matrix(NA, n.samples, N)
# True negatives
tn.samples <- matrix(NA, n.samples, N)
# False positives
fp.samples <- matrix(NA, n.samples, N)
# False negatives
fn.samples <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    tp.samples[j, i] <- sum((c(out$y[i, , ]) == 1) & (c(y.rep.samples[j, i, , ]) == 1), na.rm = TRUE)
    tn.samples[j, i] <- sum((c(out$y[i, , ]) == 0) & (c(y.rep.samples[j, i, , ]) == 0), na.rm = TRUE)
    fp.samples[j, i] <- sum((c(out$y[i, , ]) == 0) & (c(y.rep.samples[j, i, , ]) == 1), na.rm = TRUE)
    fn.samples[j, i] <- sum((c(out$y[i, , ]) == 1) & (c(y.rep.samples[j, i, , ]) == 0), na.rm = TRUE)
  }
}
# Calculate Sensitivity ---------- TRUE POSITIVE RATES
sensitivity.samples <- tp.samples / (tp.samples + fn.samples)
# Calculate Specificity ---------- TRUE NEGATIVE RATES
specificity.samples <- tn.samples / (tn.samples + fp.samples)
# Calculate quantiles of sensitivity and specificity values
sensitivity.quants <- t(apply(sensitivity.samples, 2, quantile, c(0.025, 0.5, 0.975)))
specificity.quants <- t(apply(specificity.samples, 2, quantile, c(0.025, 0.5, 0.975)))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the sensitivity and specificity
sensitivity.quants
specificity.quants
sens = mean(sensitivity.samples)
spec = mean(specificity.samples)
tss = sens+spec-1





