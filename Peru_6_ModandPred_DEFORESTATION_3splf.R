library(dplyr)
library(spOccupancy)
set.seed(13)

load("data/BIRDS_3_COUNTS_ALL_SiteANNUALLY-X-REPLICATE-per-Species_ARRAY.Rdata")

### Remove Qullomayo
rem = grep(pattern = "QU", x = colnames(ssPnt))
ssPnt = ssPnt[,-rem,]

#### How many replicates per station?
reps = as.data.frame(rowSums(!is.na(ssPnt[1,,])))   ### Get number of replicates per station
names(reps) = "replicates"                   
reps$site = sub("_.*", "", rownames(reps))          ### Add a column defining the site of the station
reps$station = rownames(reps)                       ### Add a column of the station name
table(reps$site)

### For each site, select 2 stations that have been surveyed at least 3 times
rep.red = reps[reps$replicates>2,]                        # Reduce data to stations with more than 2 replicates
rep.red <-  rep.red %>% group_by(site) %>% filter(n()>1)  # Reduce data to include only sites with more than 1 station left
samp = rep.red %>% group_by(site) %>% sample_n(size = 2)  # Group data by site and select 2 random stations

div = which(colnames(ssPnt)%in%samp$station)

test = ssPnt[,div,]
train = ssPnt[,-div,]
# rm(ssPnt); gc()

# Turn counts to Presence/Absence
train[train>0]=1
# table(train) # Quick summary 9802 observations - 1437998 Non-Observations
str(train)
sp.names = rownames(train)


### OCCURENCE COVARIATES
biome = read.csv("data/Birds_4_Survey_Biomas_Habitat_SITESbyYEAR.csv")
defst = read.csv("data/Birds_4_Survey_Hansen_Deforestation_Values_SITESbyYEAR.csv")

ref = cbind(biome[c(4,1)], defst[c(1:4)])

ref2 = ref %>% group_by(site_station) %>%
  mutate(habitat = habitat, forest = deforestation_0_m) %>% 
  distinct(site_station, habitat, forest) 
length(unique(ref2$site_station)) 

ref3 = ref[,c(1,4:6)] %>% group_by(site_station) %>% # This is number of m2 deforestation in a buffer - Buffer = Pi*Radii^2 spm 
  summarise_all(mean) %>%                             # 1km buffer = 314.16Ha
  ungroup()                                           # 2km = 1256.64Ha
length(unique(ref3$site_station))                     # 5km = 7853.98Ha
ref3$deforestation_1000_m = ref3$deforestation_1000_m/10000/314.16   # Turns m2 to Ha then to proportion deforested
ref3$deforestation_2000_m = ref3$deforestation_2000_m/10000/1256.64
ref3$deforestation_5000_m = ref3$deforestation_5000_m/10000/7853.98


ref = full_join(ref2, ref3, by = "site_station")
length(unique(ref$site_station))


#### Make sure Occ.cov sites (rows) are reciprocal with and in the same order as count data
rem = which(!ref$site_station %in% colnames(train)) # Remove Occurrence vars where there are no stations in counts
if (length(rem>0)) {ref = ref[-rem,]}

rem = which(!colnames(train) %in% ref$site_station) # Remove Counts where there are no stations in Occ vars
if (length(rem)==0) {train1 = train
} else {train1 = train[,-rem,]}

ref2 = ref
m = match(x = ref2$site_station, table = colnames(train1))
ref2 = ref2[m,]
rownames(ref2) = ref2$site_station


## I am keeping only the first 12 replicates, which covers 99% of the data
train.12 = train1[,,1:12]


################
## Remove those species with less than 1% of (ie more than 10) STATIONS OCCURRENCE after removing stations & replicates
rem = c()
for (i in 1:nrow(train.12)) {
  # i = 1
  df = train.12[i,,]
  site = df
  site = rowSums(df, na.rm = T)
  site[site>0]=1
  if (sum(site)<11) {
    rem = append(rem, i)
  }
  # print(sum(df, na.rm = T))
}
train.12 = train.12[-c(rem),,]
sp.names.12 = rownames(train.12)    # 135 species remaining

# After removing sp, some sites are now empty of sightings 
# site.sum = rowSums(train.12, na.rm = T, dims = 2)
# which(colSums(site.sum)==0)

## Y Counts
# Set most common and Uncommon species as first 2 in y DF to improve use of latent variables in Species Correlation Occupancy Model
o = order(apply(train.12, 1, sum, na.rm = TRUE), decreasing = T)
o = o[c(1, length(o), 3:length(o)-1)]
sp.names.12 = sp.names.12[o]
train.12 = train.12[sp.names.12, , ]

y.new = train.12


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

###### !! NB;     NB:     NB:     !! ######
###### !!  NB:  The Spatial Latent Factor Model requires UNIQUE COORDINATES for each site
## We replicate sites over years and so coordinates are replicated & this does not work
## As a work around, I add 1 to each Latitude value, cumulated for each replicate of a site
## This adds between 1 - 5 to the Latitude value & there makes each one unique
## This makes each additional site 1m distant, but the reflectance values are created from the true location & the distance correlation starts at a minimum of 2000m so the added few meters will not impact analysis
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

##### - PUT ALL VARIABLES IN A LIST - #####
data = list(y = y.new, occ.covs = ref[2:6], det.covs = det.covs, coords = XY)
str(data)
is.na(det.covs)
range(is.na(ref[2:6]))

N <- dim(data$y)[1]  ## Number of Species

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
out.sfMsPGOcc <- sfMsPGOcc(occ.formula = ~ habitat + forest + deforestation_1000_m + deforestation_2000_m + deforestation_5000_m,
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
                           n.neighbors = 15,      
                           n.report = 15,         ## Batches not Samples
                           n.burn = 25000,
                           n.thin = 50,
                           n.chains = 1) # ,

# filename = paste0("/storage/hpc/44/slatera6/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_50000_25000_50_1_MODEL.Rdata")
filename = paste0("models/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_50000_25000_50_1_MODEL.Rdata")
save(out.sfMsPGOcc,file = filename)


y.rep.samples <- fitted(out.sfMsPGOcc)$y.rep.samples
# filename = paste0("/storage/hpc/44/slatera6/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_FITTED_SAMPLES.Rdata")
filename = paste0("models/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_FITTED_SAMPLES.Rdata")
save(y.rep.samples, file = filename)
gc()
ppc = ppcOcc(out.sfMsPGOcc, fit.stat = 'chi-squared', group = 1) # Group 1 = By Station
# save(ppc, file = "/storage/hpc/44/slatera6/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_Bayesian_P.Rdata")
save(ppc, file = "models/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_Bayesian_P.Rdata")
###### MODEL ANALYSIS & PREDICTION

load("models/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_50000_25000_50_1_MODEL.Rdata")
summary(out.sfMsPGOcc)
waicOcc(out.sfMsPGOcc)
load("models/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_Bayesian_P.Rdata")
summary(ppc)
summary(apply(ppc$fit.y.rep > ppc$fit.y, 2, mean))
hist(apply(ppc$fit.y.rep > ppc$fit.y, 2, mean))
ppc.df <- data.frame(fit = rowMeans(ppc$fit.y),          # Take rowmeans of multispecies model to get sample mean
                     fit.rep = rowMeans(ppc$fit.y.rep), 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')

diff.fit <- ppc$fit.y.rep.group.quants[3, 1,] - ppc$fit.y.group.quants[3, 1,]
bb=order(diff.fit, decreasing = F)
which(bb==1101)
plot(diff.fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')
diff.fit[732]
colnames(out.sfMsPGOcc$y)[732]

# 
############################################
############################################
### For PREDICTION 

### OCCURENCE COVARIATES
biome = read.csv("data/Birds_4_Survey_Biomas_Habitat_SITESbyYEAR.csv")
defst = read.csv("data/Birds_4_Survey_Hansen_Deforestation_Values_SITESbyYEAR.csv")

ref = cbind(biome[c(4,1)], defst[c(1:4)])

ref2 = ref %>% group_by(site_station) %>%
  mutate(habitat = habitat, forest = deforestation_0_m) %>% 
  distinct(site_station, habitat, forest) 

ref3 = ref[,c(1,4:6)] %>% group_by(site_station) %>% # This is number of m2 deforestation in a buffer - Buffer = Pi*Radii^2 spm 
  summarise_all(mean) %>%                             # 1km buffer = 314.16Ha
  ungroup()                                           # 2km = 1256.64Ha
length(unique(ref3$site_station))                     # 5km = 7853.98Ha
ref3$deforestation_1000_m = ref3$deforestation_1000_m/10000/314.16   # Turns m2 to Ha then to proportion deforested
ref3$deforestation_2000_m = ref3$deforestation_2000_m/10000/1256.64
ref3$deforestation_5000_m = ref3$deforestation_5000_m/10000/7853.98


ref = full_join(ref2, ref3, by = "site_station")
length(unique(ref$site_station))


#### Make sure Occ.cov sites (rows) are reciprocal with and in the same order as count data
rem = which(!ref$site_station %in% colnames(test)) # Remove Occurrence vars where there are no stations in counts
if (length(rem>0)) {ref = ref[-rem,]}

rem = which(!colnames(test) %in% ref$site_station) # Remove Counts where there are no stations in Occ vars
if (length(rem)==0) {test1 = test
} else {test1 = test[,-rem,]}

ref2 = ref
m = match(x = ref2$site_station, table = colnames(test1))
ref2 = ref2[m,]
rownames(ref2) = ref2$site_station

test.12 = test1[,,1:12]

sp.names.12 = rownames(test.12)    # all species remaining - But model only predicts those it is trained on

## Y Counts
# Set most common and Uncommon species as first 2 in y DF to improve use of latent variables in Species Correlation Occupancy Model
o = order(apply(test.12, 1, sum, na.rm = TRUE), decreasing = T)
o = o[c(1, length(o), 3:length(o)-1)]
sp.names.12 = sp.names.12[o]
test.12 = test.12[sp.names.12, , ]

y.new = test.12

### DETECTION COVARIATES
### (Vary With Each Survey)
# Load - Ensure replicate names and reduce quantity to the same as used in survey data
am.pm = read.csv("data/BIRDS_3_Det-Covs_AM-PM_per_siteANNUALLY_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(am.pm) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(am.pm))
am.pm = am.pm[or,]
am.pm = as.matrix(am.pm)

method = read.csv("data/BIRDS_3_Det-Covs_POINT-NET_per_siteANNUALLY_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(method) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(method))
method = method[or,]
method = as.matrix(method)

########################
### Create Coordinates for each replicate
L19 = read.csv("data/Birds_3_Surveys_SITES_ANNUAL.csv")
L19 = L19[L19$site_station %in% ref2$site_station,]
###### !! NB;     NB:     NB:     !! ######
###### !!  NB:  The Spatial Latent Factor Model requires UNIQUE COORDINATES for each site
## We replicate sites over years and so coordinates are replicated & this does not work
## As a work around, I add 1 to each Latitude value, cumulated for each replicate of a site
## This adds between 1 - 5 to the Latitude value & there makes each one unique
## This makes each additional site 1m distant, but the reflectance values are created from the true location & the distance correlation starts at a minimum of 2000m so the added few meters will not impact analysis
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
XY = df[,1:2]                     # !! NB: Select only the Projected Coordinates !!
rownames(XY) = df$site_station


#############
##### Predict Occupancy
out.sfMsPGOcc$x.names
habitatflooded = ref2$habitat
habitatflooded[habitatflooded=="flooded"]=1
habitatflooded[habitatflooded!=1]=0
habitatflooded = as.numeric(habitatflooded)
habitatterraFirma = ref2$habitat
habitatterraFirma[habitatterraFirma=="terraFirma"]=1
habitatterraFirma[habitatterraFirma!=1]=0
habitatterraFirma = as.numeric(habitatterraFirma)
forestsecondary = ref2$forest
forestsecondary[forestsecondary=="secondary"]=1
forestsecondary[forestsecondary!=1]=0
forestsecondary = as.numeric(forestsecondary)

pred.dat = data.frame(cbind(1, habitatflooded, habitatterraFirma,forestsecondary,ref2[,4:6]))
predicted.deforest = predict(object = out.sfMsPGOcc, type = "occupancy", X.0 = pred.dat, coords.0 = XY)
dimnames(predicted.deforest$psi.0.samples)[[3]] = colnames(test)
dimnames(predicted.deforest$psi.0.samples)[[2]] = out.sfMsPGOcc$sp.names
save(predicted.deforest, file = "model_outputs/Predicted_Testdata_Occupancy_DEFORESTATION.Rdata")
load("model_outputs/Predicted_Testdata_Occupancy_DEFORESTATION.Rdata")


# Turn counts to Presence/Absence
test[test>0]=1
str(test)
test1 = test[rownames(test) %in% out.sfMsPGOcc$sp.names,,]
test1 = test1[,,1:12]
test1 = rowSums(test1, dims = 2, na.rm = T)
test1[test1>0] = 1
# Put the test count species data in the same order as the model output
ord = match(rownames(predicted.deforest$psi.0.samples[1,,]), table = rownames(test1))
test1 = test1[ord,]
y.test = test1

### Get species probability of detection
det.pred = predicted.pca5.DETECTION$p.0.samples
det.pred = colMeans(det.pred, dims = 1, na.rm = T)
det.pred = read.csv("model_outputs/Det_Prob_Predicted.csv", row.names = 1)


##### AUC CALC
### AUC EXPLANATORY
train1 = rowSums(train.12, dims = 2, na.rm = T) # Complete Count Data
train1[train1>0]=1
# train1 = train1[rownames(train1)%in%rownames(out.sfMsPGOcc$y[,1,]), ] # Reduce to only species in model
# train1 = train1[, colnames(train1)%in%rownames(out.sfMsPGOcc$y[1,,])] # Reduce to only stations in model
ord = match(rownames(out.sfMsPGOcc$y[,1,]), table = rownames(train1)) # Put species into the order found in model
train1 = train1[ord,]

alp = as.data.frame(out.sfMsPGOcc$alpha.samples) # Get modeled probability of detection
alp = alp[1:nrow(train1)]                        # Select the intercept values
alp = plogis(colMeans(alp))                      # Convert to probability from logit
alp = read.csv("model_outputs/Det_Prob_Explanatory.csv", row.names = 1) # Load site by species det prob that includes number of replicates


auc.train = matrix(data = NA, nrow = out.sfMsPGOcc$n.post, ncol = nrow(train1), dimnames = list(c(1:out.sfMsPGOcc$n.post),rownames(train1))) # Create matrix to hold AUC values per posterior sample

for (i in 1:nrow(out.sfMsPGOcc$y[,1,])) { ## For each species
  # i=42
  tru = train1[i,]
  for (ii in 1:out.sfMsPGOcc$n.post) {    ## For each posterior sample
    # ii=1
    pred = out.sfMsPGOcc$z.samples[ii,i,] * alp[i,]   # Multiply occurrence by detection probabilities
    # pred = scales::rescale(pred,   to = c(0,1))          ## Try rescaling probabilities from 0-1
    # pred[pred>0.5]=1                                     ## Make prediction 0/1 based on over or under 0.5
    # pred[pred!=1]=0
    a = Metrics::auc(actual = tru, predicted = pred)  # Calculate the AUC for that species
    # auc[ii,i] = a
    auc.train[ii,i] = a                               # Add the Species AUC for that posterior sample to the df 
  }
}
hist(colMeans(auc.train, na.rm = T))
summary(colMeans(auc.train, na.rm = T))
sd(colMeans(auc.train, na.rm = T))
length(which(colMeans(auc.train, na.rm = T)>0.7))
View(auc.train)
# xx = colMeans(auc.train, na.rm = T)
# write.csv(xx, "data/Peru_99_Species_AUC_Explanatory_DEFORESTATION.csv", row.names = T)
auc.train = read.csv("data/Peru_99_Species_AUC_Explanatory_DEFORESTATION.csv")
##### AUC PREDICTIVE
# det.pred = read.csv("model_outputs/Det_Prob_Predicted.csv", row.names = 1) # Load site by species det prob that includes number of replicates
auc = matrix(data = NA, nrow = out.sfMsPGOcc$n.post, ncol = nrow(y.test), dimnames = list(c(1:out.sfMsPGOcc$n.post),rownames(y.test)))

for (i in 1:nrow(y.test)) {             ## For each species
  # i=135
  tru = y.test[i,]
  for (ii in 1:out.sfMsPGOcc$n.post) {  ## For each posterior sample
    # ii=1
    pred = predicted.deforest$psi.0.samples[ii,i,] * det.pred[i,]
    # pred = scales::rescale(pred,   to = c(0,1))           ## Try rescaling probabilities from 0-1
    # pred[pred>0.5]=1                                      ## Make prediction 0/1 based on over or under 0.5
    # pred[pred!=1]=0
    a = Metrics::auc(actual = tru, predicted = pred) #  pred.ext[ii,i,])
    auc[ii,i] = a
  }
}
hist(colMeans(auc, na.rm = T))
summary(colMeans(auc, na.rm = T))
length(which(colMeans(auc, na.rm = T)>0.7)) # 49
auc.sp = colMeans(auc, na.rm = T)
sd(auc.sp, na.rm = T)
# write.csv(auc.sp, "data/Peru_99_Species_AUC_Predicted_DEFORESTATION.csv", row.names = T)
# auc.sp = read.csv("data/Peru_99_Species_AUC_Predicted_5PCA.csv")
plot(sort(auc[,1]), type = "l")

auc.red = auc[ , colSums(is.na(auc))==0]
for (jj in 1:500) {
  sort(auc.red[jj,], na.last = NA)
  #   if (!is.na(auc.red[jj,])) {auc.red[jj,] = sort(auc.red[jj,], na.last = NA)}
}
for (jj in 1:ncol(auc)) {
  auc[,jj] = sort(auc[,jj], na.last = T)
  # if (!is.na(auc[,jj])) {auc[,jj] = sort(auc[,jj])}
}

# matplot(auc, type = "l", title(main = "Each line is a species - Range of AUC over posterior samples"))
# lines(rowMeans(auc, na.rm = T),  type = "l", lwd = 3)

# 
# abline(h=0.7, col = "red", lwd=2)
# matplot(t(auc.red), type = "l")
# lines(rowMeans(t(auc.red), na.rm = T),  type = "l", lwd = 3)
# abline(h=0.7, col = "red", lwd=2)


###### Community Similarity
## Explanatory
bray = matrix(data = NA, nrow = out.sfMsPGOcc$n.post, ncol = ncol(train1)) # Matrix to store Dissim values
tru.y = rowSums(out.sfMsPGOcc$y, na.rm = T, dims = 2)                      # Get observed communities
for (ii in 1:ncol(train1)) {            ## For each community
  # ii=1
  bc = c()                              # Empty vector to store dissim for current community
  tru = tru.y[,ii]                      # Get observed current community
  tru[tru>0]=1
  for (i in 1:500) {                      ## For each Posterior Sample
    # i=1
    d = rbinom(length(alp[,ii]), size = 1, prob=alp[,ii])   # Get the LATENT DETECTION
    pred = out.sfMsPGOcc$z.samples[i,,ii] * d               # Latents multiplied
    # pred = out.sfMsPGOcc$psi.samples[i,,ii] * alp[,ii]    # Probabilities multiplied
    # pred = scales::rescale(pred, to = c(0,1))
    # pred[pred>0.5]=1
    # pred[pred!=1]=0
    comp = rbind(tru, pred)                          # Combine Observed & Predcited communities
    dis = vegan::vegdist(comp, method = "bray", binary = F)  ## Get Dissimilarity value between the 2
    bc = c(bc, dis)                                  # Add Dissimilarity to the current community vector
  }
  bray[,ii] = bc                                     # Add community dissimilarity vector to the over matrix
}

summary(colMeans(bray))
sd(colMeans(bray), na.rm = T)
hist(colMeans(bray))

## Predictive
bray.pred = matrix(data = NA, nrow = out.sfMsPGOcc$n.post, ncol = ncol(y.test), dimnames = list(c(1:out.sfMsPGOcc$n.post),colnames(y.test))) # Matrix to store Dissim values
tru.y = y.test                          ## Get observed communities
for (ii in 1:ncol(y.test)) {            ## For each community
  # ii=1
  bc = c()                              # Empty vector to store dissim for current community
  tru = tru.y[,ii]                      # Get observed current community
  tru[tru>0]=1
  for (i in 1:500) {                    ## For each Posterior Sample
    # i=1
    d = rbinom(length(det.pred[,ii]), size = 1, prob=det.pred[,ii])   # Get the LATENT DETECTION
    pred = predicted.deforest$z.0.samples[i,,ii] * d   # Get predicted Latent or Probability of Occurrence for current site
    # pred = predicted.pca5$psi.0.samples[i,,ii]
    # pred = pred * det.pred[,ii]                    # Multiply by Probability of Detection
    # pred = scales::rescale(pred, to = c(0,1))
    # pred[pred>0.5]=1
    # pred[pred!=1]=0
    comp = rbind(tru, pred)                               # Combine Observed & Predcited communities
    dis = vegan::vegdist(comp, method = "bray", binary = F)  ## Get Dissimilarity value between the 2
    bc = c(bc, dis)                     # Add Dissimilarity to the current community vector
  }
  bray.pred[,ii] = bc                        # Add community dissimilarity vector to the over matrix
}
summary(colMeans(bray.pred))
sd(colMeans(bray.pred))
hist(colMeans(bray.pred))



##### Species RICHNESS
# Explanatory
rich = matrix(data = NA, nrow = out.sfMsPGOcc$n.post, ncol = ncol(train1), dimnames = list(c(1:out.sfMsPGOcc$n.post),colnames(train1)))         # Create matrix to store predicted richness per station per posterior sample

rich.tru = colSums(train1, na.rm = T)  # Create an observed richness per station
summary(rich.tru)
hist(rich.tru, breaks = 35)

for (ii in 1:out.sfMsPGOcc$n.post) {  ## For each posterior sample
  # ii=1
  pred = out.sfMsPGOcc$psi.samples[ii,,]   # Get the Species/Station matrix for current sample
  pred = as.matrix(pred*alp)                          # Multiply probability of occupancy by probability of detection
  pred <- matrix(rbinom((135*1101), 1, pred), 135, 1101) # Create Binomial Matrix from Preds
  # pred = scales::rescale(x = pred, to = c(0,1))
  # pred[pred>0.5]=1
  # pred[pred!=1]=0
  rich.pred = colSums(pred, na.rm = T)                  # get the site richness-es
  # comp = rbind(rich.tru, rich.pred)                     # Combine Observed & Predicted Richness
  # dis = vegan::vegdist(x = comp, method = "bray", binary = F) # Calculate Bray dissimilarity between richness-es
  rich[ii,] = rich.pred
  # rich.dif[ii,] = rich.pred - rich.tru
  # rich.dissim = c(rich.dissim, dis)
}
rich.pred = colMeans(rich)
summary(rich.pred)
sd(rich.pred)
hist(rich.pred, breaks = 35)
summary(rich.tru)
sd(rich.tru)
hist(rich.tru, breaks = 35)

## Predictive
rich.pred = matrix(data = NA, nrow = out.sfMsPGOcc$n.post, ncol = ncol(test1), dimnames = list(c(1:out.sfMsPGOcc$n.post),colnames(test1)))   # Create matrix to store predicted richness per station per posterior sample

rich.tru = colSums(test1, na.rm = T)  # Create an observed richness per station
summary(rich.tru)
hist(rich.tru, breaks = 35)

for (ii in 1:out.sfMsPGOcc$n.post) {  ## For each posterior sample
  # ii=1
  pred = predicted.deforest$psi.0.samples[ii,,]     # Get the Species/Station matrix for current sample
  pred = as.matrix(pred * det.pred)             # Multiply occupancy by mean predicted probability of detection
  pred <- matrix(rbinom((135*34), 1, pred), 135, 34) # Create Binomial Matrix from Preds
  # pred = scales::rescale(x = pred, to = c(0,1))
  # pred[pred>0.5]=1
  # pred[pred!=1]=0
  rich.pred[ii,] = colSums(pred, na.rm = T)               # get the site richness-es
  # comp = rbind(rich.tru, rich.pred)                     # Combine Observed & Predicted Richness
  # dis = vegan::vegdist(x = comp, method = "bray", binary = F) # Calculate Bray dissimilarity between richness-es
  # rich[ii,] = rich.pred
  # rich.dif[ii,] = rich.pred - rich.tru
  # rich.pred.dissim = c(rich.pred.dissim, dis)
}
rich = colMeans(rich.pred)
summary(rich)
hist(rich)
sd(rich)
summary(rich.tru)
hist(rich.tru)
sd(rich.tru)

t.test(x = rich.tru, y = rich)

