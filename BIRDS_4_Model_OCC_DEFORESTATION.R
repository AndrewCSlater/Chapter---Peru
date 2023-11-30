#####
### OCCUPANCY MODEL BIRDS USING spOCC
library(dplyr)
library(spOccupancy)
library(coda)        # for MCMC diagnostics

load("data/BIRDS_3_COUNTS_ALL_Site-X-REPLICATE-per-Species_ARRAY.Rdata")

# Turn counts to Presence/Absence
ssPnt[ssPnt>0]=1
table(ssPnt) # Quick summary 10946 observations - 1621534 Non-Observations
str(ssPnt)
sp.names = rownames(ssPnt)

### Remove Qullomayo
rem = grep(pattern = "QU", x = colnames(ssPnt))
ssPnt = ssPnt[,-rem,]
sp.names = rownames(ssPnt)

### OCCURENCE COVARIATES
### (Station Constant across Surveys)
L19 = read.csv("data/BIRDS_2_Stations_UTM19_inc_SPECIES.csv")
# L18 = read.csv("data/BIRDS_1_Stations_UTM18_Fixed_Names_Coords.csv")
# occ = full_join(L18, L19)

occ = L19
names(occ)
# occ = select(occ, 1:6, 18:20, everything(), -c("Reflectance.Band.1", "Reflectance.Band.2", "Reflectance.Band.3", "Distance.to.Quebrada", "Altitude.of.River", "Station.Notes", "Distance.of.Research.Station.to.PEM"))
occ = select(occ, 1:6, 28:30, 10, 15, 17:18, 26, 9, 19:20, 22:25)
names(occ) = c("site_N", "site_c", "station", "lon", "lat", "utm_zone", "coords", "site_station", "stat_coords", "alt", "dist_Large_river", "dist_lake_swamp", "dist_2nd_de_forest", "habitat", "prim_2nd", "dist_sml_rds", "dist_highway", "deforest_5km","deforest_2km", "deforest_1km", "human_traffic_1to10")
occ$prim_2nd[occ$prim_2nd==1] = "primary"
occ$prim_2nd[occ$prim_2nd==0] = "secondary"
occ$prim_2nd[is.na(occ$prim_2nd)] = "not_noted"

rem = which(is.na(occ$deforest_1km) | is.na(occ$deforest_2km) | is.na(occ$deforest_5km))
occ = occ[-c(rem), ]

#### Make sure Occ.cov sites (rows) are reciprocal with and in the same order as count data
occ = occ[order(occ$site_station),]

rem = which(!occ$site_station %in% colnames(ssPnt)) # Remove Occurrence vars where there are no stations in counts
occ = occ[-rem,]
rem = which(!colnames(ssPnt) %in% occ$site_station) # Remove Counts where there are no stations in Occ vars
if (length(rem)==0) {ssPnt1 = ssPnt
} else {ssPnt1 = ssPnt[,-rem,]}

rownames(occ) = occ$site_station
or = match(colnames(ssPnt1[,,1]), rownames(occ))
occ = occ[or,]

occ$habitat = gsub(pattern = " ", "", occ$habitat)
occ$habitat = gsub(pattern = "/", "_", occ$habitat)
occ$habitat = as.factor(occ$habitat)

## can't have NA in Occurence covariates
## Will change some NA to UNKNOWN

occ.covs = occ[,c(14:15, 18:20)]

### Replicates
cc=c()
for (i in 1:nrow(ssPnt1[1,,])) {    # For each row (Station)
  print(length(which(!is.na(ssPnt1)[1,i,])))
  cc = append(cc, length(which(!is.na(ssPnt1)[1,i,])))  ## get the number of Non-NA (replicates that were conducted) 
}
hist(cc, breaks = 35)
summary(cc)
quantile(cc, probs = c(0.9, 0.95, 0.975))

length((ssPnt1)[1,,]) # 14490 
length(which(is.na(ssPnt1)[1,,])) # 12152
length(which(!is.na(ssPnt1)[1,,])) # 2338

## Sites X Replicates = 14,490
## But in only 2,338 of those (16.1%) is a survey conducted

## Median number of replicates is 5, and Mean is 5.1 
### 75% of sites have 7 or fewer
### 90% have 10 or fewer
### 95% have 14 or fewer 
### 97.5% have 20 or fewer

## I am keeping only the first 14 replicates, which covers 95% of the data
ssPnt.14 = ssPnt1[,,1:14]
ssPnt.5  = ssPnt1[,,1:5]

#######
## Remove those species with NO OCCURRENCE after removing stations & replicates
rem = c()
for (i in 1:nrow(ssPnt)) {
  # i = 100
  df = ssPnt[i,,]
  site = df
  # site = rowSums(df, na.rm = T)
  # site[site>0]=1
  # if (sum(site)<2) {
  if (sum(colSums(df,na.rm = T))==0) {
    rem = append(rem, i)
  }
  # print(sum(df, na.rm = T))
}
ssPnt = ssPnt[-c(rem),,]
sp.names = rownames(ssPnt)    # 271 species remaining

rem = c()
for (i in 1:nrow(ssPnt.14)) {
  # i = 100
  df = ssPnt.14[i,,]
  site = df
  # site = rowSums(df, na.rm = T)
  # site[site>0]=1
  # if (sum(site)<2) {
  if (sum(colSums(df,na.rm = T))==0) {
    rem = append(rem, i)
  }
  # print(sum(df, na.rm = T))
}
ssPnt.14 = ssPnt.14[-c(rem),,]
sp.names.14 = rownames(ssPnt.14)    # 294 species remaining

rem = c()
for (i in 1:nrow(ssPnt.5)) {
  # i = 100
  df = ssPnt.5[i,,]
  site = df
  # site = rowSums(df, na.rm = T)
  # site[site>0]=1
  # if (sum(site)<2) {
  if (sum(colSums(df,na.rm = T))==0) {
    rem = append(rem, i)
  }
  # print(sum(df, na.rm = T))
}
ssPnt.5 = ssPnt.5[-c(rem),,]
sp.names.5 = rownames(ssPnt.5)    # 271 species remaining

## Y Counts
# Set most common and Uncommon species as first 2 in y DF to improve use of latent variables in Species Correlation Occupancy Model
o = order(apply(ssPnt.14, 1, sum, na.rm = TRUE), decreasing = T)
o = o[c(1, length(o), 3:length(o)-1)]
sp.names.14 = sp.names.14[o]

o = order(apply(ssPnt.5, 1, sum, na.rm = TRUE), decreasing = T)
o = o[c(1, length(o), 3:length(o)-1)]
sp.names.5 = sp.names.5[o]

ssPnt.14 <- ssPnt.14[sp.names.14, , ]
ssPnt.5 <- ssPnt.5[sp.names.5, , ]


surveys = list(ssPnt.14, ssPnt.5)
reps = c("Reps.14", "Reps.5")
for (j in 1:length(surveys)) {
  # j=1
y.new = surveys[[j]]

### DETECTION COVARIATES
### (Vary With Each Survey)
am.pm = read.csv("data/BIRDS_3_Det-Covs_AM-PM_per_site_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(am.pm) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(am.pm))
am.pm = am.pm[or,]
am.pm = as.matrix(am.pm)
am.pm.use = am.pm[,1:ncol(y.new[1,,])]

year = read.csv("data/BIRDS_3_Det-Covs_YEAR_per_site_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(year) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(year))
year = year[or,]
year = as.matrix(year)
year.use = year[,1:ncol(y.new[1,,])]

method = read.csv("data/BIRDS_3_Det-Covs_POINT-NET_per_site_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(method) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(method))
method = method[or,]
method = as.matrix(method)
method.use = method[,1:ncol(y.new[1,,])]

# effort = read.csv("data/BIRDS_3_Det-Covs_EFFORT_per_site_Replicate.csv", row.names = 1)
# #### Make sure Det.covs have same column names as y
# colnames(effort) = colnames(y.new[1,,])
# or = match(colnames(y.new[,,1]), rownames(effort))
# effort = effort[or,]
# effort = as.matrix(effort)
# effort = effort[,1:14]

det.covs = list(method = method.use, am.pm = am.pm.use, year = year.use)


### COORDINATES
coords = select(occ, "lon", "lat", "site_station")
XY = coords[,1:2]                                                # Select only the Projeted Coordinates
rownames(XY) = rownames(occ)

# ###############################################################
# #### Calculate & Plot Species Accumulation
# y = read.csv("data/BIRDS_USE_Station_X_Species.csv")
# y = y[y$site_stn %in% coords$site_station,]
# y = y[,-c(1:8)]
# empty = which(colSums(y)==0)
# y = y[,-empty]
# y[y>0] = 1
# accALL <- vegan::specaccum(y, method="random", permutations=100)
# plot(accALL$sites, accALL$richness,type = "l", lwd = 2,
#      xlab="Number of Sites",
#      ylab="Species Richness")
# 

#####
##### - PUT ALL VARIABLES IN A LIST - #####
data = list(y = y.new, occ.covs = occ.covs, det.covs = det.covs, coords = XY)
# det.covs = det.covs, occ.covs = occ.covs, 
str(data)
is.na(det.covs)
range(is.na(occ.covs))
# 
# ## we will supply initial values for the following parameters: alpha.comm (community-level detection coefficients), beta.comm (community-level occurrence coefficients), alpha (species-level detection coefficients), beta (species-level occurrence coefficients), tau.sq.beta (community-level occurrence variance parameters), tau.sq.alpha (community-level detection variance parameters, z (latent occurrence variables for all species).
# ## The initial values for the latent occurrence matrix are specified as a matrix with N rows corresponding to the number of species and J columns corresponding to the number of sites.
# 
# N <- dim(data$y)[1]  ## Number of Species
# 
# ms.inits <- list(alpha.comm = 0, 
#                  beta.comm = 0, 
#                  beta = 0, 
#                  alpha = 0,
#                  tau.sq.beta = 1, 
#                  tau.sq.alpha = 1, 
#                  z = apply(data$y, c(1, 2), max, na.rm = TRUE))
# 
# ## By default, the prior hyperparameter values for the community-level means have a mean of 0 and a variance of 2.72. For the community-level variance parameters, by default we set the scale and shape parameters to 0.1. Below we specify these priors explicitly.
# ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
#                   alpha.comm.normal = list(mean = 0, var = 2.72), 
#                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
#                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))


### 1 -- INTERCEPT ONLY
out.intercept <- msPGOcc(occ.formula = ~ 1, 
                         det.formula = ~ 1, 
                         data = data, 
                         # inits = ms.inits, 
                         n.samples = 6000, 
                         # priors = ms.priors, 
                         n.omp.threads = 1, 
                         verbose = TRUE, 
                         n.report = 100, 
                         n.burn = 1000,
                         n.thin = 20, 
                         n.chains = 4) #,
# k.fold = 4, k.fold.threads = 4, k.fold.only = F)
summary(out.intercept)

save(out.intercept, file = "models/1-INTERCEPT_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
load("models/1-INTERCEPT_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
##### Create & Save data
## Posterior Predictive Check
ppc.intercept = ppcOcc(out.intercept, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.intercept)
# Extract the fitted values
y.rep.intercept <- fitted(out.intercept)$y.rep.samples

# Calculate AUC using the detection-nondetection values -------------------
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out.intercept$n.post * out.intercept$n.chains
N <- nrow(data$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- pROC::auc(response = c(y.new[i, ,]), predictor = c(y.rep.intercept[j, i, , ]))
  } # i (species)
} # j (iteration)
intercept_points = list(ppc.intercept, y.rep.intercept, auc.vals)
save(intercept_points, file = "model_outputs/1-INTERCEPT_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")

load("model_outputs/1-INTERCEPT_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
ppc.intercept = intercept_points[[1]]
y.rep.intercept = intercept_points[[2]]
auc.vals = intercept_points[[3]]

# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
rownames(auc.quants) = sp.ordered
write.csv(auc.quants, "model_outputs/AUC_Intercept_POINTS.csv")
auc.quants = read.csv("model_outputs/AUC_Intercept_POINTS.csv", row.names = 1)

####### Analyse Data ########
summary(out.intercept)
# mean(out.intercept$k.fold.deviance)
waicOcc(out.intercept)
## Posterior Predictive Check
# Get the actual Bayesian p-values for individual species and community
# Species-level
species.bpvs.intercept <- apply(ppc.intercept$fit.y.rep > ppc.intercept$fit.y, 2, mean)
species.bpvs.intercept
# Community-level NB: THIS IS LITERALLY THE MEAN OF ALL SPECIES
comm.bpv <- mean(ppc.intercept$fit.y.rep > ppc.intercept$fit.y)

str(y.rep.intercept)
auc.quants
mean(auc.quants[,2])
sd(auc.quants[,2])
range(auc.quants[,2])
hist(auc.quants[,2])
length(which(auc.quants[,2]>0.7))

# Calculate species-specific confusion matrices ---------------------------
# Number of non-missing y values for each species
n.obs <- sum(!is.na(c(y.new[1, , ])))
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
    tp.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.intercept[j, i, , ]) == 1), na.rm = TRUE)
    tn.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.intercept[j, i, , ]) == 0), na.rm = TRUE)
    fp.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.intercept[j, i, , ]) == 1), na.rm = TRUE)
    fn.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.intercept[j, i, , ]) == 0), na.rm = TRUE)
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
hist(sensitivity.quants[,2])
mean(specificity.quants)
sens = mean(sensitivity.samples)
spec = mean(specificity.samples)
tss = sens+spec-1

cbind(sens, spec, tss)

##################################
### 2 -- DETECTION (AM/PM) ONLY
## specify the number of threads to use (n.omp.threads), the number of MCMC samples (n.samples), the amount of samples to discard as burn-in (n.burn), the thinning rate (n.thin), and arguments to control the display of sampler progress (verbose, n.report).
s = Sys.time()
out.det2 <- msPGOcc(occ.formula = ~ 1, 
                    det.formula = ~ effort * method,
                    data = data, 
                    inits = ms.inits, 
                    n.samples = 5000, 
                    priors = ms.priors, 
                    n.omp.threads = 1, 
                    verbose = TRUE, 
                    n.report = 100, 
                    n.burn = 2000,
                    n.thin = 15, 
                    n.chains = 4)# ,
# k.fold = 4, k.fold.threads = 4, k.fold.only = F)
ss = Sys.time()-s
summary(out.det2)
save(out.det, file = "models/1-det_Min1Sp_14Repls_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
print(ss)
load("models/1-detAMPM_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
##### Create & Save data
## Posterior Predictive Check
ppc.detAMPM = ppcOcc(out.detAMPM, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.detAMPM)
# Extract the fitted values
y.rep.detAMPM <- fitted(out.detAMPM)$y.rep.samples

# Calculate AUC using the detection-nondetection values -------------------
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out.detAMPM$n.post * out.detAMPM$n.chains
N <- nrow(data$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- pROC::auc(response = c(y.new[i, ,]), predictor = c(y.rep.detAMPM[j, i, , ]))
  } # i (species)
} # j (iteration)
detAMPM_points = list(ppc.detAMPM, y.rep.detAMPM, auc.vals)
save(detAMPM_points, file = "model_outputs/1-detAMPM_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
load("model_outputs/1-detAMPM_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")

# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
rownames(auc.quants) = sp.ordered
write.csv(auc.quants, "model_outputs/AUC_detAMPM_POINTS.csv")
auc.q.burn = read.csv("model_outputs/AUC_detAMPM_POINTS.csv", row.names = 1)

####### Analyse Data ########
summary(out.detAMPM)
# mean(out.detAMPM$k.fold.deviance)
waicOcc(out.detAMPM)
## Posterior Predictive Check
# Get the actual Bayesian p-values for individual species and community
# Species-level
species.bpvs.detAMPM <- apply(ppc.detAMPM$fit.y.rep > ppc.detAMPM$fit.y, 2, mean)
species.bpvs.detAMPM
# Community-level NB: THIS IS LITERALLY THE MEAN OF ALL SPECIES
comm.bpv <- mean(ppc.detAMPM$fit.y.rep > ppc.detAMPM$fit.y)

str(y.rep.detAMPM)
auc.quants
mean(auc.quants[,2])
sd(auc.quants[,2])
range(auc.quants[,2])
hist(auc.quants[,2])
length(which(auc.quants[,2]>0.7))

# Calculate species-specific confusion matrices ---------------------------
# Number of non-missing y values for each species
n.obs <- sum(!is.na(c(y.new[1, , ])))
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
    tp.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.detAMPM[j, i, , ]) == 1), na.rm = TRUE)
    tn.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.detAMPM[j, i, , ]) == 0), na.rm = TRUE)
    fp.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.detAMPM[j, i, , ]) == 1), na.rm = TRUE)
    fn.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.detAMPM[j, i, , ]) == 0), na.rm = TRUE)
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
hist(sensitivity.quants[,2])
mean(specificity.quants)
sens = mean(sensitivity.samples)
spec = mean(specificity.samples)
tss = sens+spec-1
###################

##################################
### 3 -- OCCURRENCE ONLY
## specify the number of threads to use (n.omp.threads), the number of MCMC samples (n.samples), the amount of samples to discard as burn-in (n.burn), the thinning rate (n.thin), and arguments to control the display of sampler progress (verbose, n.report).
s = Sys.time()

out.occ.habPRIM <- msPGOcc(occ.formula = ~ as.matrix(occ.covs$habitat, occ.covs$prim_2nd), 
                           det.formula = ~ 1, 
                           data = data, 
                           # inits = ms.inits, 
                           n.samples = 6000, 
                           # priors = ms.priors, 
                           n.omp.threads = 1, 
                           verbose = TRUE, 
                           n.report = 100, 
                           n.burn = 1000,
                           n.thin = 20, 
                           n.chains = 4)# ,
# k.fold = 4, k.fold.threads = 4, k.fold.only = F)
ss = Sys.time()-s
# summary(out.occ.habPRIM)
# save(out.occ.habPRIM, file = "1-occHABITAT_PRIMARY_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
load("models/1-occHABITAT_PRIMARY_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
waicOcc(out.occ.habPRIM)
##### Create & Save data
## Posterior Predictive Check
# ppc.habPRIM = ppcOcc(out.occ.habPRIM, fit.stat = 'freeman-tukey', group = 1)
# save(ppc.habPRIM, file = "2-ppc.habPRIM.Rdata")
load("model_outputs/2-ppc.habPRIM.Rdata")
summary(ppc.habPRIM)
# Extract the fitted values
# y.rep.habPRIM <- fitted(out.occ.habPRIM)$y.rep.samples
# save(y.rep.habPRIM, file = "2-y.rep.habPRIM.Rdata")
load("model_outputs/2-y.rep.habPRIM.Rdata")


# Calculate AUC using the detection-nondetection values -------------------
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out.occ.habPRIM$n.post * out.occ.habPRIM$n.chains
N <- nrow(data$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- pROC::auc(response = c(y.new[i, ,]), predictor = c(y.rep.habPRIM[j, i, , ]))
  } # i (species)
} # j (iteration)
detAMPM_points = list(ppc.detAMPM, y.rep.detAMPM, auc.vals)
save(detAMPM_points, file = "model_outputs/1-detAMPM_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
load("model_outputs/1-detAMPM_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")

# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
rownames(auc.quants) = sp.ordered
write.csv(auc.quants, "model_outputs/AUC_detAMPM_POINTS.csv")
auc.q.burn = read.csv("model_outputs/AUC_detAMPM_POINTS.csv", row.names = 1)

####### Analyse Data ########
summary(out.detAMPM)
# mean(out.detAMPM$k.fold.deviance)
waicOcc(out.detAMPM)
## Posterior Predictive Check
# Get the actual Bayesian p-values for individual species and community
# Species-level
species.bpvs.detAMPM <- apply(ppc.detAMPM$fit.y.rep > ppc.detAMPM$fit.y, 2, mean)
species.bpvs.detAMPM
# Community-level NB: THIS IS LITERALLY THE MEAN OF ALL SPECIES
comm.bpv <- mean(ppc.detAMPM$fit.y.rep > ppc.detAMPM$fit.y)

str(y.rep.detAMPM)
auc.quants
mean(auc.quants[,2])
sd(auc.quants[,2])
range(auc.quants[,2])
hist(auc.quants[,2])
length(which(auc.quants[,2]>0.7))

# Calculate species-specific confusion matrices ---------------------------
# Number of non-missing y values for each species
n.obs <- sum(!is.na(c(y.new[1, , ])))
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
    tp.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.detAMPM[j, i, , ]) == 1), na.rm = TRUE)
    tn.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.detAMPM[j, i, , ]) == 0), na.rm = TRUE)
    fp.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.detAMPM[j, i, , ]) == 1), na.rm = TRUE)
    fn.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.detAMPM[j, i, , ]) == 0), na.rm = TRUE)
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
hist(sensitivity.quants[,2])
mean(specificity.quants)
sens = mean(sensitivity.samples)
spec = mean(specificity.samples)
tss = sens+spec-1




# ###############################
# ####### 3 - Spatial Latent factor Model Using Forest Gradient AND Detection TIME
# s=Sys.time()
# out.POINTS_LatFact_ampm <- sfMsPGOcc(occ.formula = ~ 1,
#                             det.formula = ~ am.pm,
#                             data = data,
#                             # inits = ms.inits,
#                             n.batch = 300,
#                             n.factors = 4, 
#                             batch.length = 25,
#                             accept.rate = 0.43,
#                             # priors = ms.priors,
#                             cov.model = "exponential",
#                             # tuning = tuning.list,
#                             n.omp.threads = 1,
#                             verbose = TRUE,
#                             NNGP = TRUE,
#                             n.neighbors = 5,
#                             search.type = 'cb',
#                             n.report = 10,
#                             n.burn = 1300,
#                             n.thin = 20,
#                             n.chains = 4)
# 
# print(Sys.time()-s)
# save(out.POINTS_LatFact_ampm, file = "models/2-SpatLatFact_datAMPM_POINTS_7500its_4ch_20thn_1300brn_1240smpl.Rdata")
# load("models/models/2-SpatLatFact_datAMPM_POINTS_7500its_4ch_20thn_1300brn_1240smpl.Rdata")
# 
# rm(out.gradTime)
# summary(out.gradTimeLF)
# # mean(out.gradTimeLF$k.fold.deviance)
# waicOcc(out.gradTimeLF)
# ## Posterior Predictive Check
# ppc.sfMsP.gradTime <- ppcOcc(out.gradTimeLF, fit.stat = 'freeman-tukey', group = 1)
# summary(ppc.sfMsP.gradTime)
# # Get the actual Bayesian p-values for individual species and community
# # Species-level
# species.bpvs.gradTime <- apply(ppc.sfMsP.gradTime$fit.y.rep > ppc.sfMsP.gradTime$fit.y, 2, mean)
# species.bpvs.gradTime
# # Community-level NB: THIS IS LITERALLY THE MEAN OF ALL SPECIES
# comm.bpv <- mean(ppc.sfMsP.gradTime$fit.y.rep > ppc.sfMsP.gradTime$fit.y)
# 
# # Calculate AUC using the detection-nondetection values -------------------
# # Extract the fitted values
# y.rep.samples.gradTime <- fitted(out.gradTimeLF)$y.rep.samples
# #  Alternative using Jeff's Function
# y.rep.samples.gradTime <- fittedNew(out.gradTimeLF)$y.rep.samples
# str(y.rep.samples.gradTime)
# # Will calculate an AUC value for each species and each iteration of the 
# # MCMC in order to get an AUC estimate with uncertainty.
# n.samples <- out.gradTimeLF$n.post * out.gradTimeLF$n.chains
# N <- nrow(data$y)
# auc.vals <- matrix(NA, n.samples, N)
# for (j in 1:n.samples) {
#   print(j)
#   for (i in 1:N) {
#     auc.vals[j, i] <- pROC::auc(response = c(y.new[i, , ]), predictor = c(y.rep.samples.gradTime[j, i, , ]))
#   } # i (species)
# } # j (iteration)
# 
# # Calculate quantiles of AUC values for each species
# auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T))
# # The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# # form our credible interval for the AUC for each species.
# auc.quants
# mean(auc.quants[,2])
# sd(auc.quants[,2])
# range(auc.quants[,2])
# hist(auc.quants[,2])
# length(which(auc.quants[,2]>0.7))


### 6 - Spatial Latent factor Model Using SENTINEL1 Radar AND Detection TIME

out.SpLaFa_habPRIM <- sfMsPGOcc(occ.formula = ~ occ.covs$habitat + occ.covs$prim_2nd,
                                det.formula = ~ 1,
                                data = data,
                                # inits = ms.inits,
                                n.batch = 300,
                                n.factors = 4, 
                                batch.length = 25,
                                accept.rate = 0.43,
                                # priors = ms.priors,
                                cov.model = "exponential",
                                # tuning = tuning.list,
                                n.omp.threads = 1,
                                verbose = TRUE,
                                NNGP = TRUE,
                                n.neighbors = 5,
                                search.type = 'cb',
                                n.report = 100,
                                n.burn = 1300,
                                n.thin = 20,
                                n.chains = 4)
save(out.SpLaFa_habPRIM, file = "models/1-occHABITAT_PRIMARY_POINTS_7000its_4ch_20thn_1300brn_1240smpl.Rdata")
summary(out.SpLaFa_habPRIM)
waicOcc(out.SpLaFa_habPRIM)

## Posterior Predictive Check
ppc.SpLaFa_habPRIM <- ppcOcc(out.SpLaFa_habPRIM, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.SpLaFa_habPRIM)
# Get the actual Bayesian p-values for individual species and community
# Species-level
bayesP.SpLaFa_habPRIM <- apply(ppc.SpLaFa_habPRIM$fit.y.rep > ppc.SpLaFa_habPRIM$fit.y, 2, mean)
bayesP.SpLaFa_habPRIM
# Community-level NB: THIS IS LITERALLY THE MEAN OF ALL SPECIES
comm.bpv <- mean(ppc.SpLaFa_habPRIM$fit.y.rep > ppc.SpLaFa_habPRIM$fit.y)


# Calculate AUC using the detection-nondetection values -------------------
# Extract the fitted values
y.rep.SpLaFa_habPRIMe <- fitted(out.SpLaFa_habPRIM)$y.rep.samples
save(y.rep.SpLaFa_habPRIMe, file = "model_outputs/2-y.rep.SpLaFa_habPRIMe.Rdata")
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out.SpLaFa_habPRIM$n.post * out.SpLaFa_habPRIM$n.chains
N <- nrow(data$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- pROC::auc(response = c(y.new[i, , ]), predictor = c(y.rep.SpLaFa_habPRIMe[j, i, , ]))
  } # i (species)
} # j (iteration)

# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
auc.quants
mean(auc.quants[,2])
sd(auc.quants[,2])
range(auc.quants[,2])
hist(auc.quants[,2])
length(which(auc.quants[,2]>0.7))
# Calculate species-specific confusion matrices ---------------------------
# Number of non-missing y values for each species
n.obs <- sum(!is.na(c(y.new[1, , ])))
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
    tp.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.SpLaFa_habPRIMe[j, i, , ]) == 1), na.rm = TRUE)
    tn.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.SpLaFa_habPRIMe[j, i, , ]) == 0), na.rm = TRUE)
    fp.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.SpLaFa_habPRIMe[j, i, , ]) == 1), na.rm = TRUE)
    fn.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.SpLaFa_habPRIMe[j, i, , ]) == 0), na.rm = TRUE)
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
hist(sensitivity.quants[,2])
mean(specificity.quants)




### 7 - Spatial Latent factor Model Using SENTINEL1 Radar AND FOREST CLASS plus Detection TIME
s=Sys.time()
out.Rad_Forest_Time <- sfMsPGOcc(occ.formula = ~ as.matrix(occ.covs[,1], scale(occ.covs[,c(1,7:10)])), ## 4 Radar values VV VH mean & sd at 1 buffer value
                                 det.formula = ~ surv.time,
                                 data = data,
                                 # inits = ms.inits,
                                 n.batch = 400,
                                 n.factors = 4, 
                                 batch.length = 25,
                                 accept.rate = 0.43,
                                 # priors = ms.priors,
                                 cov.model = "exponential",
                                 # tuning = tuning.list,
                                 n.omp.threads = 3,
                                 verbose = TRUE,
                                 NNGP = TRUE,
                                 n.neighbors = 5,
                                 search.type = 'cb',
                                 n.report = 50,
                                 n.burn = 2000,
                                 n.thin = 20,
                                 n.chains = 3,
                                 k.fold = 4)

print(Sys.time()-s)
# save(out.Rad_Forest_Time, file = "models/lf4_OCCrad30m-ForestClass_DETtime_240reps_5000its_4ch_15thn_500brn_1200smpl.Rdata")
summary(out.Rad_Forest_Time)
waicOcc(out.Rad_Forest_Time)
lfload = summary(out.Rad_Forest_Time$lambda.samples)
summary(lfload[[1]][675:896,1])
## Posterior Predictive Check
ppc.sfMsP.RadTime <- ppcOcc(out.Rad_Forest_Time, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.sfMsP.RadTime)
# Get the actual Bayesian p-values for individual species and community
# Species-level
species.bpvs.gradTime <- apply(ppc.sfMsP.RadTime$fit.y.rep > ppc.sfMsP.RadTime$fit.y, 2, mean)
species.bpvs.gradTime
# Community-level NB: THIS IS LITERALLY THE MEAN OF ALL SPECIES
comm.bpv <- mean(ppc.sfMsP.RadTime$fit.y.rep > ppc.sfMsP.RadTime$fit.y)

# Calculate AUC using the detection-nondetection values -------------------
# Extract the fitted values
y.rep.samples.RadTime <- fitted(out.Rad_Forest_Time)$y.rep.samples
str(y.rep.samples.RadTime)
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out.Rad_Forest_Time$n.post * out.Rad_Forest_Time$n.chains
N <- nrow(data$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- pROC::auc(response = c(y.new[i, , ]), predictor = c(y.rep.samples.RadTime[j, i, , ]))
  } # i (species)
} # j (iteration)

# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
auc.quants
mean(auc.quants[,2])
sd(auc.quants[,2])
range(auc.quants[,2])
hist(auc.quants[,2])
length(which(auc.quants[,2]>0.7))





### 8 - Spatial Latent factor Model Using 5 PCA plus Detection TIME
### TESTING PCA - THIS USING PCA MADE FROM MCARI1 PVI VARI
s=Sys.time()
out.PCA_Time <- sfMsPGOcc(occ.formula = ~ as.matrix(scale(occ.covs[,c(2:6)])), 
                          det.formula = ~ surv.time,
                          data = data,
                          # inits = ms.inits,
                          n.batch = 1000,
                          n.factors = 4, 
                          batch.length = 25,
                          accept.rate = 0.43,
                          # priors = ms.priors,
                          cov.model = "exponential",
                          # tuning = tuning.list,
                          n.omp.threads = 3,
                          verbose = TRUE,
                          NNGP = TRUE,
                          n.neighbors = 5,
                          search.type = 'cb',
                          n.report = 100,
                          n.burn = 1000,
                          n.thin = 50,
                          n.chains = 3,
                          k.fold = 4)

print(Sys.time()-s)
# save(out.Rad_Forest_Time, file = "models/lf4_OCCrad30m-ForestClass_DETtime_240reps_5000its_4ch_15thn_500brn_1200smpl.Rdata")
summary(out.Rad_Forest_Time)
waicOcc(out.Rad_Forest_Time)
lfload = summary(out.Rad_Forest_Time$lambda.samples)
summary(lfload[[1]][675:896,1])
## Posterior Predictive Check
ppc.sfMsP.RadTime <- ppcOcc(out.Rad_Forest_Time, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.sfMsP.RadTime)
# Get the actual Bayesian p-values for individual species and community
# Species-level
species.bpvs.gradTime <- apply(ppc.sfMsP.RadTime$fit.y.rep > ppc.sfMsP.RadTime$fit.y, 2, mean)
species.bpvs.gradTime
# Community-level NB: THIS IS LITERALLY THE MEAN OF ALL SPECIES
comm.bpv <- mean(ppc.sfMsP.RadTime$fit.y.rep > ppc.sfMsP.RadTime$fit.y)

# Calculate AUC using the detection-nondetection values -------------------
# Extract the fitted values
y.rep.samples.RadTime <- fitted(out.Rad_Forest_Time)$y.rep.samples
str(y.rep.samples.RadTime)
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out.Rad_Forest_Time$n.post * out.Rad_Forest_Time$n.chains
N <- nrow(data$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- pROC::auc(response = c(y.new[i, , ]), predictor = c(y.rep.samples.RadTime[j, i, , ]))
  } # i (species)
} # j (iteration)

# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
auc.quants
mean(auc.quants[,2])
sd(auc.quants[,2])
range(auc.quants[,2])
hist(auc.quants[,2])
length(which(auc.quants[,2]>0.7))



out.ms <- msPGOcc(occ.formula = ~ as.matrix(scale(occ.covs[,1])), 
                  det.formula = ~1, 
                  data = data, 
                  inits = ms.inits, 
                  n.samples = 5000, 
                  priors = ms.priors, 
                  n.omp.threads = 1, 
                  verbose = TRUE, 
                  n.report = 100, 
                  n.burn = 500,
                  n.thin = 5, 
                  n.chains = 3)

waicOcc(out.ms) # 187127

## Spatial & Latent Factor Model - Shown to be best for Prediction by Doser 2022
dist = dist(coords)
min.dist <- min(dist)
max.dist <- max(dist)
priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
               alpha.comm.normal = list(mean = 0, var = 2.72),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
               phi.unif = list(3 / max.dist, 3 / min.dist))


## Start with 1st var & then add 1 for 5 total
s = Sys.time()
results = c()
for (i in 2:5) {
  # i = 1
  print(i)
  names(occ.covs)
  occ.ms.formula <- ~ as.matrix(scale(occ.covs[,c(2:i)]))
  
  out.sfMsPGOcc <- sfMsPGOcc(occ.formula = occ.ms.formula, 
                             det.formula = ~1, 
                             data = data, 
                             # inits = inits, 
                             n.batch = 25,           ### Runs multiple batches of given batch length which multiplied give total samples per chain 
                             batch.length = 200,     ### eg 25 * 200 = 5000 samples per chain
                             accept.rate = 0.43, 
                             priors = priors, 
                             n.factors = 2,
                             cov.model = "exponential", 
                             tuning = list(phi = 0.5), 
                             n.omp.threads = 4, 
                             verbose = TRUE, 
                             NNGP = TRUE, 
                             n.neighbors = 5, 
                             n.report = 5,         ## Batches not Samples
                             n.burn = 2000, 
                             n.thin = 10, 
                             n.chains = 3)
  
  results[i] = waicOcc(out.sfMsPGOcc)[3] # 186361
  write.csv(results, "data/WAIC_best_noncorrVars_1-5.csv", row.names = F)
  rm(out.sfMsPGOcc)
  gc()
  print(Sys.time()-s)
}



