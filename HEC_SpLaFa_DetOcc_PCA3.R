#####
### OCCUPANCY MODEL BIRDS USING spOCC
library(dplyr)
library(spOccupancy)
library(coda)        # for MCMC diagnostics
# library(tidyverse)

### SURVEY DATA
load("data/BIRDS_3_COUNTS_ALL_SiteANNUALLY-X-REPLICATE-per-Species_ARRAY.Rdata")

# Turn counts to Presence/Absence
ssPnt[ssPnt>0]=1
# table(ssPnt) # Quick summary 10946 observations - 1621534 Non-Observations
str(ssPnt)
sp.names = rownames(ssPnt)

### Remove Qullomayo
rem = grep(pattern = "QU", x = colnames(ssPnt))
ssPnt = ssPnt[,-rem,]
sp.names = rownames(ssPnt)

### OCCURENCE COVARIATES
ref = read.csv("data/Birds_4_Survey_Reflectance_Values_SITESbyYEAR.csv")
ref$site_station = stringr::str_extract(ref$rep_ID, "[^_]*_[^_]*")
ref$site_station = paste0(ref$site_station, "_", ref$year)

#### Make sure Occ.cov sites (rows) are reciprocal with and in the same order as count data
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

m = match(x = ref2$site_station, table = colnames(ssPnt1))
ref2 = ref2[m,]
rownames(ref2) = ref2$site_station
# write.csv(ref2, "data/Birds_4.5_Survey_Reflectance_Values_perSITEperYEAR.csv", row.names = F)

#####################################
#####################################
### Reduce and match data
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

## Y Counts
# Set most common and Uncommon species as first 2 in y DF to improve use of latent variables in Species Correlation Occupancy Model
o = order(apply(ssPnt.12, 1, sum, na.rm = TRUE), decreasing = T)
o = o[c(1, length(o), 3:length(o)-1)]
sp.names.12 = sp.names.12[o]

ssPnt.12 = ssPnt.12[sp.names.12, , ]

y.new = ssPnt.12


##############################
### DETECTION COVARIATES
### (Vary With Each Survey)
# Load - Ensure replicate names and reduce quantity to the same as used in survey data
am.pm = read.csv("data/BIRDS_3_Det-Covs_AM-PM_per_siteANNUALLY_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(am.pm) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(am.pm))
am.pm = am.pm[or,]
am.pm = as.matrix(am.pm)
am.pm.use = am.pm[,1:ncol(y.new[1,,])]

method = read.csv("data/BIRDS_3_Det-Covs_POINT-NET_per_siteANNUALLY_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(method) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(method))
method = method[or,]
method = as.matrix(method)
method.use = method[,1:ncol(y.new[1,,])]

det.covs = list(method = method.use, am.pm = am.pm.use)

########################
### Create Coordinates for each replicate
L19 = read.csv("data/Birds_3_Surveys_SITES_ANNUAL.csv")
L19 = L19[L19$site_station %in% ref2$site_station,]

###### !!  NB:  The Spatial Latent Factor Model requires UNIQUE COORDINATES for each site
## We replicate sites over years and so coordinates are replicated & this does not work
## As a work around, I add 1 to each Latitude value, cumulated for each replicate of a site
## This adds between 1 - 5 to the Latitude value & there makes each one unique
## This makes each additional site 1m distant, but the reflectance values are created from the true location & the distance correlation starts at a minimum of 2000m so te added few meters will not impact analysis
###### !! NB;     NB:     NB:     !! ######

### COORDINATES
coords = select(L19, "utm_long", "utm_lat", "site_code", "station", "site_station")
coords = unique(coords)
coords$count = 1
df = coords
df$csum <- ave(df$count, paste0(df$site_code, df$station),  FUN=cumsum)

ord = match(ref2$site_station, table = df$site_station)
df = df[ord,]
df$utm_lat = df$utm_lat+df$csum
XY = df[,1:2]                     # !! NB: Select only the Projeted Coordinates !!
rownames(XY) = df$site_station



##################################
##################################
### Reduce Dimensions of Reflectance Data
### Create a PCA vectors

df = select(ref2, -"site_station", -"year")
#calculate principal components
pca_ref2 = prcomp(x = df, scale. = T)
#reverse the signs
pca_ref2$rotation = -1 * pca_ref2$rotation
# #display principal components
# pca_ref2$rotation

#reverse the signs of the scores
pca_ref2$x = -1*pca_ref2$x

# summary(pca_ref2)
#
# #display the first six scores
# head(pca_ref2$x)
#
# biplot(pca_ref2, choices = c(1,2),  scale = 0)
#
# #devtools::install_github("vqv/ggbiplot")
# library(ggbiplot)
# g <- ggbiplot(pca_ref2,
#               choices = c(1,2),
#               obs.scale = 1,
#               var.scale = 1,
#               # groups = training$Species,
#               ellipse = TRUE,
#               circle = TRUE,
#               ellipse.prob = 0.68)
# g <- g + scale_color_discrete(name = '')
# g <- g + theme(legend.direction = 'horizontal',
#                legend.position = 'top')
# print(g)
#
# #calculate total variance explained by each principal component
# round(pca_ref2$sdev^2 / sum(pca_ref2$sdev^2) ,digits = 2)
#
# #calculate total variance explained by each principal component
# var_explained = pca_ref2$sdev^2 / sum(pca_ref2$sdev^2)
#
# #create scree plot
# qplot(c(1:26), var_explained) +
#   geom_line() +
#   xlab("Principal Component") +
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot") +
#   ylim(0, 1)

## Obtain the PCA components for each site_station
pca_occ = pca_ref2$x
# write.csv(pca_occ, "/storage/hpc/44/slatera6/SpLF3_DetOcc_PCA3_12reps_buffer500m_PCA_Components.csv", row.names = T)



#####
##### - PUT ALL VARIABLES IN A LIST - #####
data = list(y = y.new, occ.covs = pca_occ[,1:3], det.covs = det.covs, coords = XY)
# str(data)
# is.na(det.covs)
# range(is.na(pca_occ[,1:3]))

  ## we will supply initial values for the following parameters: alpha.comm (community-level detection coefficients), beta.comm (community-level occurrence coefficients), alpha (species-level detection coefficients), beta (species-level occurrence coefficients), tau.sq.beta (community-level occurrence variance parameters), tau.sq.alpha (community-level detection variance parameters, z (latent occurrence variables for all species).
  ## The initial values for the latent occurrence matrix are specified as a matrix with N rows corresponding to the number of species and J columns corresponding to the number of sites.

N <- dim(data$y)[1]  ## Number of Species

  ms.inits <- list(alpha.comm = 0,
                   beta.comm = 0,
                   beta = 0,
                   alpha = 0,
                   tau.sq.beta = 1,
                   tau.sq.alpha = 1,
                   z = apply(data$y, c(1, 2), max, na.rm = TRUE))

  ## By default, the prior hyperparameter values for the community-level means have a mean of 0 and a variance of 2.72. For the community-level variance parameters, by default we set the scale and shape parameters to 0.1. Below we specify these priors explicitly.
  ## Spatial & Latent Factor Model - Shown to be best for Prediction by Doser 2022
  dist = dist(XY)
  min.dist <- 2000 # Standard min dist is the smallest distance between 2 survey points, but where survey points are close together relative to total scale of all sites, it is better to set the min dist as transect length, or in this case most site points 100m apart all sites are 175km apart, so I will set min dist as 2000m which is in keeping with the size of a Site (which is surveyed by many stations)
  max.dist <- max(dist)
  priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                 alpha.comm.normal = list(mean = 0, var = 2.72),
                 tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                 tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                 phi.unif = list(3 / max.dist, 3 / min.dist))

##################################
### 1 -- SPATIAL LATENT FACTOR
## specify the number of threads to use (n.omp.threads), the number of MCMC samples (n.samples), the amount of samples to discard as burn-in (n.burn), the thinning rate (n.thin), and arguments to control the display of sampler progress (verbose, n.report).

out.sfMsPGOcc <- sfMsPGOcc(occ.formula = ~ PC1 + PC2 + PC3,
                           det.formula = ~ method + am.pm,
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
                           n.chains = 1,
                           k.fold = 3, k.fold.threads = 3, k.fold.seed = 3, k.fold.only = T)

# filename = paste0("/storage/hpc/44/slatera6/SpLF3_DetOcc_PCA3_12reps_buffer500m_MODEL.Rdata")
filename = paste0("/storage/hpc/44/slatera6/SpLF3_DetOcc_PCA3_12reps_buffer500m_3FOLD-XVAL.Rdata")
# /storage/hpc/44/slatera6/
save(out.sfMsPGOcc,file = filename)

################################################
### AUC and TSS values
################################################
# gc()
# y.rep.samples <- fitted(out.sfMsPGOcc)$y.rep.samples
# filename = paste0("/storage/hpc/44/slatera6/SpLF3_DetOcc_PCA3_12reps_buffer500m_FITTED_SAMPLES.Rdata")
# save(y.rep.samples, file = filename)
# # Will calculate an AUC value for each species and each iteration of the
# # MCMC in order to get an AUC estimate with uncertainty.
# n.samples <- out.sfMsPGOcc$n.post * out.sfMsPGOcc$n.chains
# N <- nrow(out.sfMsPGOcc$y)
# auc.vals <- matrix(NA, n.samples, N)
# for (j in 1:n.samples) {
#   # print(j)
#   for (i in 1:N) {
#     auc.vals[j, i] <- pROC::auc(response = c(out.sfMsPGOcc$y[i, , ]), predictor = c(y.rep.samples[j, i, , ]))
#   } # i (species)
# } # j (iteration)
# colnames(auc.vals) = rownames(out.sfMsPGOcc$y)
# write.csv(auc.vals, "/storage/hpc/44/slatera6/SpLF3_DetOcc_PCA3_12reps_buffer500m_AUC.csv", row.names = F)
# 
# load("/storage/hpc/44/slatera6/SpLF3_DetOcc_PCA3_12reps_buffer500m_MODEL.Rdata")
# ppc = ppcOcc(out.sfMsPGOcc, fit.stat = 'freeman-tukey', group = 1) # Group 1 = By Station
# save(ppc, file = "/storage/hpc/44/slatera6/SpLF3_DetOcc_PCA3_12reps_buffer500m_Bayesian_P.Rdata")





