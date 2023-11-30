### Occupancy Model

# Deforestation sites only

#####
### OCCUPANCY MODEL BIRDS USING spOCC
library(dplyr)
library(spOccupancy)
library(coda)        # for MCMC diagnostics

load("data/BIRDS_3_COUNTS_ALL_SiteANNUALLY-X-REPLICATE-per-Species_ARRAY.Rdata")

### Remove Qullomayo
rem = grep(pattern = "QU", x = colnames(ssPnt))
ssPnt = ssPnt[,-rem,]

# Turn counts to Presence/Absence
ssPnt[ssPnt>0]=1
# table(ssPnt) # Quick summary 10946 observations - 1621534 Non-Observations
str(ssPnt)
sp.names = rownames(ssPnt)
sites = colnames(ssPnt)
sites = gsub(pattern = "_[0-9].*$", replacement = "",  sites)

### OCCURENCE COVARIATES
### (Station Constant across Surveys)
L19 = read.csv("data/BIRDS_2_Stations_UTM19_inc_SPECIES.csv")
occ = L19
occ = select(occ, 1:6, 28:30, 10, 15, 17:18, 26, 9, 19:20, 22:25)
names(occ) = c("site_N", "site_c", "station", "lon", "lat", "utm_zone", "coords", "site_station", "stat_coords", "alt", "dist_Large_river", "dist_lake_swamp", "dist_2nd_de_forest", "habitat", "prim_2nd", "dist_sml_rds", "dist_highway", "deforest_5km","deforest_2km", "deforest_1km", "human_traffic_1to10")
occ$prim_2nd[occ$prim_2nd==1] = "primary"
occ$prim_2nd[occ$prim_2nd==0] = "secondary"
occ$prim_2nd[is.na(occ$prim_2nd)] = "not_noted"

rem = which(is.na(occ$deforest_1km) | is.na(occ$deforest_2km) | is.na(occ$deforest_5km))
occ = occ[-c(rem), ]

#### Make sure Occ.cov sites (rows) are reciprocal with and in the same order as count data
occ = occ[order(occ$site_station),]

rem = which(!occ$site_station %in% sites) # Remove Occurrence vars where there are no stations in counts
occ = occ[-rem,]
rem = which(!sites %in% occ$site_station) # Remove Counts where there are no stations in Occ vars
if (length(rem)==0) {ssPnt1 = ssPnt
} else {ssPnt1 = ssPnt[,-rem,]}

rownames(occ) = occ$site_station
or = match(occ$site_station, table = sites)
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
  # print(length(which(!is.na(ssPnt1)[1,i,])))
  cc = append(cc, length(which(!is.na(ssPnt1)[1,i,])))  ## get the number of Non-NA (replicates that were conducted) 
}
# hist(cc, breaks = 35)
# summary(cc)
# quantile(cc, probs = c(0.9, 0.95, 0.975))
# 
# length((ssPnt1)[1,,]) # 14490 
# length(which(is.na(ssPnt1)[1,,])) # 12152
# length(which(!is.na(ssPnt1)[1,,])) # 2338

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
  
  
  ##################################
  ### 2 -- OCCURRENCE ONLY
  ## specify the number of threads to use (n.omp.threads), the number of MCMC samples (n.samples), the amount of samples to discard as burn-in (n.burn), the thinning rate (n.thin), and arguments to control the display of sampler progress (verbose, n.report).
  out.detOcc <- msPGOcc(occ.formula = ~ deforest_5km + deforest_2km + deforest_1km + habitat + as.factor(prim_2nd), 
                        det.formula = ~ method + am.pm + year,
                        data = data, 
                        # inits = ms.inits, 
                        n.samples = 10000, 
                        # priors = ms.priors, 
                        n.omp.threads = 1, 
                        verbose = TRUE, 
                        n.report = 100, 
                        n.burn = 2000,
                        n.thin = 20, 
                        n.chains = 4)
  
  filename = paste0("/storage/hpc/44/slatera6/both_DETandOCC_",reps[j],".Rdata")
  save(out.detOcc, file = filename)
}



Min stations = 12 (1%)
# 