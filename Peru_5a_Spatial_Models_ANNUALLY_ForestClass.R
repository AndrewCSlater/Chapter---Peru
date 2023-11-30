#####
### OCCUPANCY MODEL BIRDS USING spOCC
library(dplyr)
library(spOccupancy)
library(coda)        # for MCMC diagnostics

load("data/BIRDS_3_COUNTS_ALL_SiteANNUALLY-X-REPLICATE-per-Species_ARRAY.Rdata")

# Turn counts to Presence/Absence
ssPnt[ssPnt>0]=1
table(ssPnt) # Quick summary 10940 observations - 1615156 Non-Observations
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

## Try and reciprocate what will be modelled with reflectance data.
# The OCC variables do not change year on year, but split the data as if they did.

# 1. Get the ssPnt site names without years
surv_stns = colnames(ssPnt)
surv_stns = gsub("_[1-9].*$", "", surv_stns)

rem = which(!occ$site_station %in% surv_stns) # Remove Occurrence vars where there are no stations in counts
occ = occ[-rem,]
rem = which(!surv_stns %in% occ$site_station) # Remove Counts where there are no stations in Occ vars

if (length(rem)==0) {ssPnt1 = ssPnt
} else {ssPnt1 = ssPnt[,-rem,]}

surv_stns = surv_stns[-rem]

rownames(occ) = occ$site_station
or = match(surv_stns, rownames(occ))
occ = occ[or,]

occ$habitat = gsub(pattern = " ", "", occ$habitat)
occ$habitat = gsub(pattern = "/", "_", occ$habitat)
occ$habitat = as.factor(occ$habitat)

#### Can't Have NA values in Occurrence Variables...but need to keep same data as used in reflectance...
#### Possibly Use mean or median values in deforestation where NA is present ???
mean(occ$deforest_1km, na.rm = T); median(occ$deforest_1km, na.rm = T); length(which(is.na(occ$deforest_1km))) # 7.37, 4 & 291/1123
mean(occ$deforest_2km, na.rm = T); median(occ$deforest_2km, na.rm = T); length(which(is.na(occ$deforest_2km))) # 63.4 & 66   "  "
mean(occ$deforest_5km, na.rm = T); median(occ$deforest_5km, na.rm = T); length(which(is.na(occ$deforest_5km))) # 514.7 & 429 "  "
# 
# rem = which(is.na(occ$deforest_1km) | is.na(occ$deforest_2km) | is.na(occ$deforest_5km))
# occ = occ[-c(rem), ]

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
quantile(cc, probs = c(0.75, 0.9, 0.95, 0.975, 0.99))

length((ssPnt1)[1,,]) # 23583 
length(which(is.na(ssPnt1)[1,,])) # 
length(which(!is.na(ssPnt1)[1,,])) # 3265

## Sites X Replicates = 23583
## But in only 3265 of those (13.8%) is a survey conducted

## Median number of replicates is 2, and Mean is 2.9 
### 75% of sites have 4 or fewer
### 90% have 6 or fewer
### 95% have 8 or fewer 
### 97.5% have 10 or fewer
### 99% have 12 or fewer

## I am keeping only the first 14 replicates, which covers 95% of the data
ssPnt.12 = ssPnt1[,,1:12]

#######
## Remove those species with NO OCCURRENCE after removing stations & replicates
rem = c()
for (i in 1:nrow(ssPnt.12)) {
  # i = 100
  df = ssPnt.12[i,,]
  site = df
  # site = rowSums(df, na.rm = T)
  # site[site>0]=1
  # if (sum(site)<2) {
  if (sum(colSums(df,na.rm = T))==0) {
    rem = append(rem, i)
  }
  # print(sum(df, na.rm = T))
}
ssPnt.12 = ssPnt.12[-c(rem),,]
sp.names.12 = rownames(ssPnt.12)    # 358 species remaining

## Y Counts
# Set most common and Uncommon species as first 2 in y DF to improve use of latent variables in Species Correlation Occupancy Model
o = order(apply(ssPnt.12, 1, sum, na.rm = TRUE), decreasing = T)
o = o[c(1, length(o), 3:length(o)-1)]
sp.names.12 = sp.names.12[o]

ssPnt.12 <- ssPnt.12[sp.names.12, , ]

surveys = list(ssPnt.12)
reps = c("Reps.12")
# for (j in 1:length(surveys)) {
  j=1
  y.new = surveys[[j]]
  
  ### DETECTION COVARIATES
  ### (Vary With Each Survey)
  am.pm = read.csv("data/BIRDS_3_Det-Covs_AM-PM_per_siteANNUALLY_Replicate.csv", row.names = 1)
  #### Make sure Det.covs have same column names as y
  colnames(am.pm) = colnames(y.new[1,,])
  or = match(colnames(y.new[,,1]), rownames(am.pm))
  am.pm = am.pm[or,]
  am.pm = as.matrix(am.pm)
  am.pm.use = am.pm[,1:ncol(y.new[1,,])]
  
  # year = read.csv("data/BIRDS_3_Det-Covs_YEAR_per_site_Replicate.csv", row.names = 1)
  # #### Make sure Det.covs have same column names as y
  # colnames(year) = colnames(y.new[1,,])
  # or = match(colnames(y.new[,,1]), rownames(year))
  # year = year[or,]
  # year = as.matrix(year)
  # year.use = year[,1:ncol(y.new[1,,])]
  
  method = read.csv("data/BIRDS_3_Det-Covs_POINT-NET_per_siteANNUALLY_Replicate.csv", row.names = 1)
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
  
  det.covs = list(method = method.use, am.pm = am.pm.use)
  
  
  ### COORDINATES
  coords = select(occ, "lon", "lat", "site_station")
  XY = coords[,1:2]                                                # Select only the Projeted Coordinates
  rownames(XY) = rownames(occ)
  
  ### COORDINATES
  coords = select(occ, "lon", "lat", "site_c", "station", "site_station")
  coords = unique(coords)
  coords$count = 1
  df = coords
  df$csum <- ave(df$count, paste0(df$site_code, df$station),  FUN=cumsum)
  
  # ord = match(ref2$site_station, table = df$site_station)
  # df = df[ord,]
  df$lat = df$lat+df$csum
  XY = df[,1:2]                     # !! NB: Select only the Projeted Coordinates !!
  rownames(XY) = df$site_station
  
###############################################################
#### Calculate & Plot Species Accumulation
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

  occ.covs$deforest_1km[is.na(occ.covs$deforest_1km)] = median(occ.covs$deforest_1km, na.rm = T)
  occ.covs$deforest_2km[is.na(occ.covs$deforest_2km)] = median(occ.covs$deforest_2km, na.rm = T)
  occ.covs$deforest_5km[is.na(occ.covs$deforest_5km)] = median(occ.covs$deforest_5km, na.rm = T)
  occ.covs$habitat = as.character(occ.covs$habitat)
  occ.covs$habitat[is.na(occ.covs$habitat)] = "habitat-Unknown"
  occ.covs$prim_2nd[is.na(occ.covs$prim_2nd)] = "not_noted"
  
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
 
N <- dim(data$y)[1]  ## Number of Species
  # 
  # ms.inits <- list(alpha.comm = 0, 
  #                  beta.comm = 0, 
  #                  beta = 0, 
  #                  alpha = 0,
  #                  tau.sq.beta = 1, 
  #                  tau.sq.alpha = 1, 
  #                  z = apply(data$y, c(1, 2), max, na.rm = TRUE))
  # 
  dist = dist(XY)
  min.dist <- 2000 # Standard min dist is the smallest distance between 2 survey points, but where survey points are close together relative to total scale of all sites, it is better to set the min dist as transect length, or in this case most site points 100m apart all sites are 175km apart, so I will set min dist as 2000m which is in keeping with the size of a Site (which is surveyed by many stations)
  max.dist <- max(dist)
  priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                 alpha.comm.normal = list(mean = 0, var = 2.72),
                 tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                 tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                 phi.unif = list(3 / max.dist, 3 / min.dist))
  
  
### 6 - Spatial Latent factor Model
out.SpLaFa_habPRIM <- sfMsPGOcc(occ.formula = ~ 1, #occ.covs$habitat + occ.covs$prim_2nd + scale(occ.covs$deforest_1km) + scale(occ.covs$deforest_2km) + scale(occ.covs$deforest_5km)  ,
                                  det.formula = ~ det.covs$method + det.covs$am.pm,
                                  data = data,
                                  # inits = ms.inits,
                                  n.batch = 2000,
                                  n.factors = 3, 
                                  batch.length = 25,
                                  accept.rate = 0.43,
                                  priors = priors,
                                  cov.model = "exponential",
                                  # tuning = tuning.list,
                                  n.omp.threads = 1,
                                  verbose = TRUE,
                                  NNGP = TRUE,
                                  n.neighbors = 15,
                                  search.type = 'cb',
                                  n.report = 100,
                                  n.burn = 25000,
                                  n.thin = 50,
                                  n.chains = 1)
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
  
  
  
