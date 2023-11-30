### OCCUPANCY MODEL BIRDS USING spOCC
library(dplyr)
library(spOccupancy)
library(coda)        # for MCMC diagnostics


load("data/BIRDS_3_COUNTS_ALL_SiteANNUALLY-X-REPLICATE-per-Species_ARRAY.Rdata")
# Turn counts to Presence/Absence
ssPnt[ssPnt>0]=1
### Remove Qullomayo
rem = grep(pattern = "QU", x = colnames(ssPnt))
ssPnt = ssPnt[,-rem,]
str(ssPnt)
sp.names = rownames(ssPnt)


### OCCURENCE COVARIATES
### (Station Constant across Surveys)
ref = read.csv("data/Birds_4_Survey_Reflectance_Values_SITESbyYEAR.csv")
ref$site_station = stringr::str_extract(ref$rep_ID, "[^_]*_[^_]*")
ref$site_station = paste0(ref$site_station, "_", ref$year)

#### Make sure Occ.cov sites (rows) are reciprocal with and in the same order as count data
ref = ref[order(ref$site_station),]
rem = which(!ref$site_station %in% colnames(ssPnt)) # Remove Occurrence vars where there are no stations in counts
if (length(rem>0)) {ref = ref[-rem,]}

rem = which(!colnames(ssPnt) %in% ref$site_station) # Remove Counts where there are no stations in Occ vars
if (length(rem)==0) {ssPnt1 = ssPnt
} else {ssPnt1 = ssPnt[,-rem,]}

# Reduce the reflectance data to one row for each site_station_year
ref2 = ref
ref2 = select(ref2, -rep_ID)
ref2 = ref2 %>% group_by(site_station) %>%
  summarise_all(mean) %>%
  ungroup()

rownames(ref2) = ref2$site_station
# length(unique(ref$site_station))

ord = match(colnames(ssPnt1[,,1]), rownames(ref2))
ref2 = ref2[ord,]

## can't have NA in Occurence covariates

### Replicates
cc=c()
for (i in 1:nrow(ssPnt1[1,,])) {    # For each row (Station)
  # print(length(which(!is.na(ssPnt1)[1,i,])))
  cc = append(cc, length(which(!is.na(ssPnt1)[1,i,])))  ## get the number of Non-NA (replicates that were conducted) 
}
# hist(cc, breaks = 35)
# summary(cc)          # Median=2; Mean=3
quantile(cc, probs = c(0.9, 0.95, 0.975, 0.99)) # Quants=6;8;10;12

## Median number of replicates is 5, and Mean is 5.1 
### 75% of sites have 7 or fewer
### 90% have 6 or fewer
### 95% have 8 or fewer 
### 97.5% have 10 or fewer
### 99% have 12 or fewer

## I am keeping only the first 12 replicates, which covers 99% of the data
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

# rem = c()
# for (i in 1:nrow(ssPnt1)) {
#   # i = 100
#   df = ssPnt1[i,,]
#   site = df
#   # site = rowSums(df, na.rm = T)
#   # site[site>0]=1
#   # if (sum(site)<2) {
#   if (sum(colSums(df,na.rm = T))==0) {
#     rem = append(rem, i)
#   }
#   # print(sum(df, na.rm = T))
# }
# ssPnt1 = ssPnt1[-c(rem),,]
# sp.names = rownames(ssPnt1)    # 358 species remaining

## Y Counts
# Set most common and Uncommon species as first 2 in y DF to improve use of latent variables in Species Correlation Occupancy Model
o = order(apply(ssPnt.12, 1, sum, na.rm = TRUE), decreasing = T)
o = o[c(1, length(o), 3:length(o)-1)]
sp.names.12 = sp.names.12[o]

# o = order(apply(ssPnt.10, 1, sum, na.rm = TRUE), decreasing = T)
# o = o[c(1, length(o), 3:length(o)-1)]
# sp.names.10 = sp.names.10[o]

# ssPnt1 <- ssPnt1[sp.names, , ]
ssPnt.12 <- ssPnt.12[sp.names.12, , ]

surveys = list(ssPnt1, ssPnt.12)
reps = c("Reps.all", "Reps.12")
# for (j in 1:length(surveys)) {
j=2
y.new = surveys[[j]]

# ### DETECTION COVARIATES
# ### (Vary With Each Survey)
# am.pm = read.csv("data/BIRDS_3_Det-Covs_AM-PM_per_site_Replicate.csv", row.names = 1)
# #### Make sure Det.covs have same column names as y
# colnames(am.pm) = colnames(y.new[1,,])
# or = match(colnames(y.new[,,1]), rownames(am.pm))
# am.pm = am.pm[or,]
# am.pm = as.matrix(am.pm)
# am.pm.use = am.pm[,1:ncol(y.new[1,,])]
# 
# # year = read.csv("data/BIRDS_3_Det-Covs_YEAR_per_site_Replicate.csv", row.names = 1)
# # #### Make sure Det.covs have same column names as y
# # colnames(year) = colnames(y.new[1,,])
# # or = match(colnames(y.new[,,1]), rownames(year))
# # year = year[or,]
# # year = as.matrix(year)
# # year.use = year[,1:ncol(y.new[1,,])]
# 
# method = read.csv("data/BIRDS_3_Det-Covs_POINT-NET_per_site_Replicate.csv", row.names = 1)
# #### Make sure Det.covs have same column names as y
# colnames(method) = colnames(y.new[1,,])
# or = match(colnames(y.new[,,1]), rownames(method))
# method = method[or,]
# method = as.matrix(method)
# method.use = method[,1:ncol(y.new[1,,])]
# 
# # effort = read.csv("data/BIRDS_3_Det-Covs_EFFORT_per_site_Replicate.csv", row.names = 1)
# # #### Make sure Det.covs have same column names as y
# # colnames(effort) = colnames(y.new[1,,])
# # or = match(colnames(y.new[,,1]), rownames(effort))
# # effort = effort[or,]
# # effort = as.matrix(effort)
# # effort = effort[,1:14]
# 
# det.covs = list(method = method.use, am.pm = am.pm.use) # , year = year.use
# 
# 
### COORDINATES
coords = read.csv("data/Birds_3_Surveys_SITES_ANNUAL.csv")
coords = select(coords, "lon", "lat", "site_station")
coords = unique(coords)
coords = coords[coords$site_station %in% colnames(ssPnt.12),]
m = match(coords$site_station, table = colnames(ssPnt.12))
coords = coords[m,]
XY = coords[,1:2]                                                # Select only the Projeted Coordinates
rownames(XY) = coords$site_station
# 
# #####
names(ref2)
occ_nam = colnames(ref2)
for (ii in 2:(ncol(ref2)-1)) {
  # ii=2
  print(ii)
  print(occ_nam[ii])
  occ.covs = scale(ref2[,ii])
  
# ##### - PUT ALL VARIABLES IN A LIST - #####
# data = list(y = y.new, occ.covs = occ.covs, det.covs = det.covs, coords = XY)
data = list(y = y.new, occ.covs = occ.covs, coords = XY)
# str(data)
# is.na(det.covs)
# range(is.na(data$occ.covs))

N <- dim(data$y)[1]  ## Number of Species

##################################
### 1 -- OCCURRENCE ONLY
## specify the number of threads to use (n.omp.threads), the number of MCMC samples (n.samples), the amount of samples to discard as burn-in (n.burn), the thinning rate (n.thin), and arguments to control the display of sampler progress (verbose, n.report).



s = Sys.time()
out.occ <- msPGOcc(occ.formula = ~ occ.covs, 
                   det.formula = ~ 1,
                   data = data, 
                   # inits = ms.inits, 
                   n.samples = 1000, 
                   # priors = ms.priors, 
                   n.omp.threads = 1, 
                   verbose = TRUE, 
                   n.report = 100, 
                   n.burn = 2500,
                   n.thin = 15, 
                   n.chains = 3)
print(Sys.time()-s)
# filename = paste0("models/OCCURRENCE_StationYear_12Reps_", occ_nam[[ii]] ,".Rdata")
# save(out.occ, file = filename)


filename = paste0("/storage/hpc/44/slatera6/OCCURRENCE_StationYear_12Reps_", occ_nam[[ii]] ,".Rdata")
save(out.occ, file = filename)
gc()


#######################################################
###   For - SPATIAL LATENT FACTOR models

## Spatial & Latent Factor Model - Shown to be best for Prediction by Doser 2022
dist = dist(XY)
min.dist <- 2000 # Standard min dist is the smallest distance between 2 survey points, but where survey points are close together relative to total scale of all sites, it is better to set the min dist as transect length, or in this case most site points 100m apart all sites are 175km apart, so I will set min dist as 2000m which is in keeping with the size of a Site (which is surveyed by many stations)
max.dist <- max(dist)
priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
               alpha.comm.normal = list(mean = 0, var = 2.72),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
               phi.unif = list(3 / max.dist, 3 / min.dist)) 


out.sfMsPGOcc <- sfMsPGOcc(occ.formula = ~ occ.covs,
                           det.formula = ~ 1, # method + am.pm, + year
                           data = data, 
                           n.batch = 250, # Runs multiple batches of given batch length which multiplied give total samples per chain
                           batch.length = 200,     ### eg 25 * 200 = 5000 samples per chain
                           accept.rate = 0.43, 
                           priors = priors, 
                           n.factors = 3,
                           cov.model = "exponential", 
                           tuning = list(phi = 0.5), 
                           n.omp.threads = 4, 
                           verbose = TRUE, 
                           NNGP = TRUE, 
                           n.neighbors = 15,     ## Reduce to 5 & Compare WAIC 51590 @ 15; 51954 @ 5
                           n.report = 15,         ## Batches not Samples
                           n.burn = 10000, 
                           n.thin = 100, 
                           n.chains = 1)

filename = paste0("/storage/hpc/44/slatera6/SpatLatFac_OCCURRENCE_StationYear_12Reps_", occ_nam[[ii]] ,".Rdata")
save(out.sfMsPGOcc, file = filename)
rm(out.sfMsPGOcc); gc()
} #  END OCC COVS (Reflectance Variables) loop
