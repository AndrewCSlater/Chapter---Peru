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
ref = read.csv("data/Birds_4_Survey_Reflectance_Values_SITESbyYEAR.csv")
ref$site_station = stringr::str_extract(ref$rep_ID, "[^_]*_[^_]*")
ref$site_station = paste0(ref$site_station, "_", ref$year)

#### Make sure Occ.cov sites (rows) are reciprocal with and in the same order as count data
rem = which(!ref$site_station %in% colnames(train)) # Remove Occurrence vars where there are no stations in counts
if (length(rem>0)) {ref = ref[-rem,]}

rem = which(!colnames(train) %in% ref$site_station) # Remove Counts where there are no stations in Occ vars
if (length(rem)==0) {train1 = train
} else {train1 = train[,-rem,]}

# Reduce the reflectance data to one row for each site_station_year
ref2 = ref
ref2 = select(ref2, -rep_ID)
ref2 = ref2 %>% group_by(site_station) %>%
  summarise_all(mean) %>%
  ungroup()

m = match(x = ref2$site_station, table = colnames(train1))
ref2 = ref2[m,]
rownames(ref2) = ref2$site_station

##############
### Replicates - How many hold 99% of surveys
# cc=c()
# for (i in 1:nrow(train1[1,,])) {    # For each row (Station)
#   print(length(which(!is.na(train1)[1,i,])))
#   cc = append(cc, length(which(!is.na(train1)[1,i,])))  ## get the number of Non-NA (replicates that were conducted)
# }
# hist(cc, breaks = 35)
# summary(cc)
# quantile(cc, probs = c(0.75, 0.9, 0.95, 0.975, 0.99))

## Median number of replicates is 2, and Mean is 2.9 
### 75% of sites have 3 or fewer
### 90% have 6 or fewer
### 95% have 8 or fewer 
### 97.5% have 9 or fewer
### 99% have 12 or fewer

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
sp.names.12 = rownames(train.12)    # 134 species remaining

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


##########################################
##########################################
### Reduce Dimensions of Reflectance Data
### Create a PCA vectors
df = select(ref2, -"site_station", -"year")
#calculate principal components
pca_ref2 = prcomp(x = df, scale. = T)
#reverse the signs
pca_ref2$rotation = -1 * pca_ref2$rotation
# display principal components
# pca_ref2$rotation

#reverse the signs of the scores
pca_ref2$x = -1*pca_ref2$x

# summary(pca_ref2)

# display the first six scores
# head(pca_ref2$x)

# biplot(pca_ref2, choices = c(1,2),  scale = 0)

# devtools::install_github("vqv/ggbiplot")
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

# calculate total variance explained by each principal component
# round(pca_ref2$sdev^2 / sum(pca_ref2$sdev^2) ,digits = 2)

# calculate total variance explained by each principal component
# var_explained = pca_ref2$sdev^2 / sum(pca_ref2$sdev^2)

# create scree plot
# qplot(c(1:26), var_explained) +
#   geom_line() +
#   xlab("Principal Component") +
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot") +
#   ylim(0, 1)

## Obtain the PCA components for each site_station
pca_occ = pca_ref2$x

#####
##### - PUT ALL VARIABLES IN A LIST - #####
data = list(y = y.new, occ.covs = pca_occ[,1:5], det.covs = det.covs, coords = XY)
str(data)
is.na(det.covs)
range(is.na(pca_occ[,1:5]))

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
names(df)

out.sfMsPGOcc <- sfMsPGOcc(occ.formula = ~ 1,
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
# k.fold = 2, k.fold.threads = 1, k.fold.seed = 3)

save(out.sfMsPGOcc, file = "1.B-Training_Min11sp_INTERCEPT_50000_25000_50_1.Rdata")


