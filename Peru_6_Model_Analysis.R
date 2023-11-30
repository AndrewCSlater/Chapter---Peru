library(dplyr)
library(spOccupancy)
library(coda)

plogis(0.8368)

######################################################
## Probability of Detection ~ Probabilty of Occurrence
load("models/SpLF3_DetOcc_PCA1_12reps_buffer500m_MODEL.Rdata")
# SpLF3_DETonly_12reps_buffer500m_MODEL
out = out.sfMsPGOcc
rm(out.sfMsPGOcc); gc()
load("models/SpLF3_DetOcc_PCA1_12reps_buffer500m_Bayesian_P.Rdata")
# SpLF3_DETonly_12reps_buffer500m_Bayesian_P
load("models/SpLF3_DetOcc_PCA1_12reps_buffer500m_FITTED_SAMPLES.Rdata")
# SpLF3_DETonly_12reps_buffer500m_FITTED_SAMPLES

a=as.data.frame(out$beta.samples)
quantile(a[1:358,1], probs = c(0.025,0.975))

summary(out)#, level = "community")
waicOcc(out)
summary(ppc) # Bayesian P Value : Community
# Get the actual Bayesian p-values for individual species
species.bpvs <- apply(ppc$fit.y.rep > ppc$fit.y, 2, mean)
summary(species.bpvs)
hist(species.bpvs)

ppc.df <- data.frame(fit = rowMeans(ppc$fit.y),          # Take rowmeans of multispecies model to get sample mean
                     fit.rep = rowMeans(ppc$fit.y.rep), 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')

diff.fit <- ppc$fit.y.rep.group.quants[3, 1,] - ppc$fit.y.group.quants[3, 1,]
bb=order(diff.fit, decreasing = F)
which(bb==1135)
plot(diff.fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')
which(diff.fit==min(diff.fit))
(out$y[1,,1][c(764,800,820)])

### CONVERGENACE
################################################################
plot(out$beta.samples[, 1:6], density = FALSE) # Plot Convergence of 1st 6 species
a=geweke.diag(out$beta.samples)#[1:358, ])
hist(a$z)
summary(a$z)
sum(a$z < 2 & a$z > -2, na.rm = T) / length(a$z) # 75%:84% of species have converged well using Geweke diagnostic (Det Only)
# summary(out$lambda.samples) # Lambda is the spatial influence per species
lam = out$lambda.samples
# write.csv(lam, "model_outputs/Deforest_5km_lambda_3spatial.csv", row.names = F)
sp.prev = as.matrix(rowSums(out$y, na.rm = T))

################################################
### AUC and TSS values
################################################
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
write.csv(auc.vals, "model_outputs/DetOcc_PCA1_12reps_buffer500m_AUC.csv", row.names = F)
auc.vals = read.csv("model_outputs/DetOcc_PCA1_12reps_buffer500m_AUC.csv")
auc.means = colMeans(auc.vals)
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
tss.quants = sensitivity.quants + specificity.quants -1

sig_sp_band = c() ## To store the species names that have significant relationships with band reflectance
### Get Influence of DETECTION (Alpha) variables 
alp = out$alpha.samples
# write.csv(bet, "model_outputs/ _alpha_detection.csv", row.names = F)

### How many species have a SIGNIFICANT relationship with the variable: ie how many have a 95%CI that excludes zero?
ncol(alp)/4
sigYN = c()                         ## Create YN significance vector
sigVAL = c()                        ## Create Influence Value vector (for if variable significantly influences a species)
aa1 = as.data.frame(alp[,359:716])  ## Subset df into all species with 1 variable
aa1 = as.data.frame(alp[,717:1076])
for (j in 1:ncol(aa1)) {            ## For each species (column)
  i = quantile(aa1[,j], c(0.025,0.975))  ## Get the 2.5 & 97.5 quants
  if (i[1]<0 & i[2]>0) {ii = 0           ## If the quants surround zero Give the Value of 0
  } else {ii = 1}                       ## Else (all 95% is either above or below zero) give a value of 1 
  sigYN = append(sigYN, ii)               ## Add that 1 or 0 value to the yes no vector
  if (ii==1) {sigVAL = append(sigVAL, colMeans(aa1[j]))}
}  
sum(sigYN)                         ## The number of species with Significant relationships
hist(sigVAL)                       ## Histogram of Significant Mean Values
hist(colMeans(aa1))                ## Histogram of ALL mean values
sig_sp = colnames(aa1)[(sigYN)==1]        ## Find species names with significant relationships
length(which(sigVAL<0)); length(which(sigVAL>0)) # How many sp with significant relationships have Negative & Positive relationships 

sig_sp = sub(".*-", replacement = "", x = sig_sp)
sig_sp = as.data.frame(cbind(sig_sp, sigVAL))
sig_sp$PosNeg[sig_sp$sigVAL>0] = "pos"
sig_sp$PosNeg[sig_sp$sigVAL<0] = "neg"
sig_sp$band = "Afternoon"

sig_sp_band = bind_rows(sig_sp_band, sig_sp)
write.csv(sig_sp_band, "data/Peru_6_Significant_Detection-Species_12reps_buffer500m.csv", row.names = F)


#############################
sig_sp_band = c() ## To store the species names that have significant relationships with band reflectance
### Get Influence of OCCURRENCE (Beta) variables 
bet = as.data.frame(out$beta.samples)
alp = as.data.frame(out$alpha.samples)
# write.csv(bet, "model_outputs/Deforest_5km_beta_occurrence.csv", row.names = F)
head(bet)[,1]
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























### Sumarise Species by Site
y = out$y # 358 Species in 1135 Stations
y2 = as.data.frame(rowSums(y, dims = 2, na.rm = T))
y2 = y2[colSums(y2)>0]
y3 = y2; y3[y3>0] = 1
summary(colSums(y3)) # Station Richness (Species per station)
summary(rowSums(y3)) # Species Prevalence (Number of stations species found at)
sum(y2)              # Number of individual observations in total (10,005)

#### Species Accumulation Curves - How many sites needed to collect given number of Species ####
acum = vegan::specaccum(comm = t(y3), method = "random", permutations = 1000)
which(acum$richness>0.95*max(acum$richness))[1] # 95% of species were found within 908 sites
which(acum$richness>0.75*max(acum$richness))[1] # 75% of species were found within 360 sites
plot(acum$sites, acum$richness,type = "l", lwd = 2,
     xlab="Number of Sites",
     ylab="Species Richness")
################################

#### Rarefacton Curve - How many samples (individual sightings) needed to collect given number of species #### 
y4 = as.data.frame(rowSums(y2))
rr = vegan::rrarefy(y4, sample = 40)
dr = vegan::drarefy(y4, sample = 40)
r = vegan::rarefy(t(y4), sample = 400)
rc = vegan::rarecurve(t(y4))
rcc = rc[[1]]
which(rcc>0.95*max(rcc))[1] # 7882 of 10005 samples needed to gather 95% of all species
which(rcc>0.75*max(rcc))[1] # 3025 of 10005 samples needed to gather 75% of all species
which(rcc>0.5*max(rcc))[1] # 844 of 10005 samples needed to gather 50% of all species
vegan::rarecurve(t(y4), sample = c(844, 3025, 7882))
#####################################################################################

## Per Species Occ & Det probabilites
## !! NB: Length will be number of sp x number of variables (categorical variables count 1 each) + 1 (for Intercept)
occ = summary(out$beta.samples) ## Mean Probability of Occurrence
meanOcc = occ[[1]][,1]
det = summary(out$alpha.samples) ## Mean probability of Detection
meanDet = det[[1]][,1]

logOcc = plogis(meanOcc)
logDet = plogis(meanDet)


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
sample_data <- data.frame(logOcc, logDet[1:(length(logDet)/4)])
colnames(sample_data)[2] = "logDet"

slope = lm(logDet[1:(length(logDet)/4)]~logOcc)
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
plot1







