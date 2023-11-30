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

occ = select(occ, 1:5, 8, 14:15, 18:20)
occ = occ[!is.na(occ$deforest_1km),]

ref = read.csv("data/Birds_4_Survey_Reflectance_Values_SITESbyYEAR.csv")
ref$site_station = stringr::str_extract(ref$rep_ID, "[^_]*_[^_]*")

occ = full_join(ref, occ, by = "site_station")
occ = occ[!is.na(occ$blue_mean_500),]
occ = occ[!is.na(occ$habitat),]

occ$site_station = paste0(occ$site_station, "_", occ$year)
occ = select(occ,-"rep_ID")


#### Make sure Occ.cov sites (rows) are reciprocal with and in the same order as count data
rem = which(!occ$site_station %in% colnames(ssPnt)) # Remove Occurrence vars where there are no stations in counts
if (length(rem>0)) {occ = occ[-rem,]}

rem = which(!colnames(ssPnt) %in% occ$site_station) # Remove Counts where there are no stations in Occ vars
if (length(rem)==0) {ssPnt1 = ssPnt
} else {ssPnt1 = ssPnt[,-rem,]}


## Split Data into Reflectance and In-Situ
ref2 = select(occ, 1:26,28)
frst = select(occ, 27:28, 34:38)

ref3 = ref2 %>% group_by(site_station) %>%
  summarise_all(mean) %>%
  ungroup()
length(unique(ref3$site_station))
##########################################
##########################################
### Reduce Dimensions of Reflectance Data
### Create a PCA vectors
df = select(ref3, -"site_station")
#calculate principal components
pca_ref2 = prcomp(x = df, scale. = T)
#reverse the signs
pca_ref2$rotation = -1 * pca_ref2$rotation
# display principal components
# pca_ref2$rotation
#reverse the signs of the scores
pca_ref2$x = -1*pca_ref2$x
## Obtain the PCA components for each site_station
pca_occ = pca_ref2$x



frst2 = frst %>% group_by(site_station) %>%
  mutate(Habitat = habitat, Forest = prim_2nd) %>% 
  distinct(site_station, Habitat, Forest) 
length(unique(frst2$site_station)) 

frst3 = frst[,c(2,5:7)] %>% group_by(site_station) %>%
  summarise_all(mean) %>%
  ungroup()
length(unique(frst3$site_station))

frst4 = full_join(frst2, frst3, by = "site_station") 
length(unique(frst4$site_station))
 

rem = which(!colnames(ssPnt) %in% ref3$site_station) # Remove Counts where there are no stations in Occ vars
if (length(rem)==0) {ssPnt1 = ssPnt
} else {ssPnt1 = ssPnt[,-rem,]}

rem = which(!ref3$site_station %in% colnames(ssPnt1)) # Remove Occurrence vars where there are no stations in counts
if (length(rem)==0) {ref3 = ref3
} else {ref3 = ref3[-rem,]}


## Try and reciprocate what will be modelled with reflectance data.
# The OCC variables do not change year on year, but split the data as if they did.
m = match(x = colnames(ssPnt1), table = ref3$site_station)
ref3 = ref3[m,]
m = match(x = colnames(ssPnt1), table = frst4$site_station)
frst4 = frst4[m,]

## I am keeping only the first 12 replicates, which covers 95% of the data
ssPnt.12 = ssPnt1[,,1:12]

#######
## Remove those species with OCCURRENCE in fewer than 11 stations after removing stations & replicates
## This is in keeping with later models
rem = c()
for (i in 1:nrow(ssPnt.12)) {
  # i = 100
  df = ssPnt.12[i,,]
  site = df
  site = rowSums(df, na.rm = T)
  site[site>0]=1
  if (sum(site)<8) {
    # if (sum(colSums(df,na.rm = T))==0) {
    rem = append(rem, i)
  }
  # print(sum(df, na.rm = T))
}
ssPnt.12 = ssPnt.12[-c(rem),,]
sp.names.12 = rownames(ssPnt.12)    # 116 species remaining

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


method = read.csv("data/BIRDS_3_Det-Covs_POINT-NET_per_siteANNUALLY_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(method) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(method))
method = method[or,]
method = as.matrix(method)
method.use = method[,1:ncol(y.new[1,,])]


det.covs = list(method = method.use, am.pm = am.pm.use)


### COORDINATES
coords = select(occ, "lon", "lat", "site_c", "station", "site_station")
coords = unique(coords)
coords$count = 1
df = coords
df$csum <- ave(df$count, paste0(df$site_code, df$station),  FUN=cumsum)

df$lat = df$lat+df$csum
XY = df[,1:2]                     # !! NB: Select only the Projeted Coordinates !!
rownames(XY) = df$site_station


#####
# ## we will supply initial values for the following parameters: alpha.comm (community-level detection coefficients), beta.comm (community-level occurrence coefficients), alpha (species-level detection coefficients), beta (species-level occurrence coefficients), tau.sq.beta (community-level occurrence variance parameters), tau.sq.alpha (community-level detection variance parameters, z (latent occurrence variables for all species).
# ## The initial values for the latent occurrence matrix are specified as a matrix with N rows corresponding to the number of species and J columns corresponding to the number of sites.

N <- dim(y.new)[1]  ## Number of Species
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


rownames(frst4) = frst4$site_station
rownames(ref3) = ref3$site_station

##### - PUT ALL VARIABLES IN A LIST - #####
data = list(y = y.new, occ.covs = pca_occ[,1:5], det.covs = det.covs, coords = XY)

### 6 - Spatial Latent factor Model
out.SpLaFa_habPRIM <- sfMsPGOcc(occ.formula = ~ PC1 + PC2 + PC3 + PC4 + PC5,
                                  # 
                                  # occ.covs$Habitat + occ.covs$Forest + scale(occ.covs$deforest_1km) + scale(occ.covs$deforest_2km) + scale(occ.covs$deforest_5km),
                                det.formula = ~ method + am.pm,
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
save(out.SpLaFa_habPRIM, file = "1.B-WAIC_Comparison_Min11sp_5PCA_50000_25000_50_1.Rdata")
