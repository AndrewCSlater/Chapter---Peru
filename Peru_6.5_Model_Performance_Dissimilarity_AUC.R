library(spOccupancy)
library(dplyr)

## Load Model
load("models/SpLF3_DetOcc_PCA3_12reps_buffer500m_MODEL.Rdata")
# SpLF3_DETonly_12reps_buffer500m_MODEL
n.samples <- out.sfMsPGOcc$n.post * out.sfMsPGOcc$n.chains
n.reps = ncol(out.sfMsPGOcc$y[1,,])
species = rownames(out.sfMsPGOcc$y[,,1])
sites = rownames(out.sfMsPGOcc$y[1,,])

## load predicted values
load("models/SpLF3_DetOcc_PCA3_12reps_buffer500m_FITTED_SAMPLES.Rdata") ## y.rep.samples
# SpLF3_DETonly_12reps_buffer500m_FITTED_SAMPLES
dimnames(y.rep.samples) = list(c(1:n.samples), species, sites, c(1:n.reps))
str(y.rep.samples)

###########
### Create probability of occurrence = Mean of All Samples & Replicates per species per site
pred.df = matrix(NA, nrow = length(sites), ncol = length(species), dimnames = list(sites, species))
for (i in 1:length(species)) {
  #  i = 1
  sp = y.rep.samples[,i,,]
  pred.y = as.data.frame(apply(sp, c(2), mean, na.rm=T))
  pred.df[,i] = pred.y$`apply(sp, c(2), mean, na.rm = T)`
}

## THOUGHT ##
## By scaling all probabilities betweem 0-1, those above 0.5 = 1, those below = 0
library(scales)
pred.df.scaled = sapply(1:ncol(pred.df), function(i) rescale(pred.df[,i]))
colnames(pred.df.scaled) = colnames(pred.df)
pred.df.scaled[pred.df.scaled>=0.5] = 1
pred.df.scaled[pred.df.scaled!=1] = 0


## Load observed values
count = read.csv("data/BIRDS_2_Station-ANNUALLY_X_Species.csv")
count = count[count$site_stn %in% sites, ]    # Remove sites not in model
ord = match(sites, count$site_stn)            # Put rows in same order as model sites
count = count[ord,]

tru.y = count[,10:ncol(count)]
tru.y = tru.y[,colnames(tru.y) %in% species]  # Remove species not in model
ord = match(species, colnames(tru.y))         # Put columns in same order as model species
tru.y = tru.y[,ord]
rownames(tru.y) = sites


########################################
## Calculate Dissimialrity Per Species : True presence only How well does model predict Presence where it did occur

# create matix to store dissimilarities
dissim = matrix(NA,nrow = length(species), ncol = 1)

for (i in 1:length(species)) {
#  i = 1
  
### Get the sites at which species[i] was observed 
tru.count = tru.y[,i]
# tru.count[tru.count==0] = NA
rem = which(tru.count==0) ## Which sites were they not observed at - to remove those sites from the predictions
tru.count = tru.count[-rem]

dist = c() # create vector to store dissimilarities for each sample
### Using raw Probabilities
for (ii in 1:n.samples) {
#  ii = 1
sample.count = matrix(data = rowSums(y.rep.samples[ii, i,,], na.rm = T), ncol = 1, dimnames = list(sites, species[i]))
sample.count[sample.count>0]=1  
# sample.count[sample.count==0] = NA
sample.count = sample.count[-rem,]

# ### Or using Scaled Categorical Presence/Absence Probabilities
# sample.count = pred.df.scaled[,i]
# sample.count = sample.count[-rem]

comp = rbind(tru.count, sample.count)
d = vegan::vegdist(comp, method = "bray")
dist = append(dist, d)
# } # Sample Loop
# hist(dist)
# summary(dist)
dissim[i,1] = median(dist)
mean(dist)
} # Species Loop
hist(dissim)
summary(dissim)
df = as.data.frame(dissim)
df$sitesSeenAt = colSums(tru.y)
plot(df$sitesSeenAt, df$V1)
rownames(df) = colnames(tru.y)
write.csv(df, "data/Peru_6.5_Dissimilarity_per_Species_ONLYWhereItwasSeen_NotAllSites.csv")
#####################################################################

# species_df = df
# species_df$species = rownames(species_df)
# colnames(species_df)[1] = "Dissim.across.sites"
# species_df$Dissim.ScaledPA.across.sites = dissim

########################################
## Calculate Dissimialrity Per Site : True presence only for Species that were Observed

# create matix to store dissimilarities
dissim = matrix(NA,nrow = length(sites), ncol = 1)

for (i in 1:length(sites)) {
  #  i = 1
  
  ### Get the sites at which species[i] was observed 
  tru.count = tru.y[i,]
  # tru.count[tru.count==0] = NA
  rem = which(tru.count==0) ## Which sites were they not observed at - to remove those sites from the predictions
  tru.count = tru.count[,-rem]
  
  # dist = c() # create vector to store dissimilarities for each sample
  # for (ii in 1:n.samples) {
  #   #  ii = 1
  #   sample.count = matrix(data = rowSums(y.rep.samples[ii,,i,], na.rm = T), nrow = 1, dimnames = list(sites[i], species))
  #   sample.count[sample.count>0]=1  
  #   # sample.count[sample.count==0] = NA
  #   sample.count = sample.count[,-rem]
    sample.count = pred.df.scaled[i,-rem]
    
    comp = rbind(tru.count, sample.count)
    d = vegan::vegdist(comp, method = "bray")
    # dist = append(dist, d)
  # } # Sample Loop
  # hist(dist)
  # summary(dist)
  # dissim[i,1] = mean(dist)
    dissim[i,1] = d
} # Site Loop
hist(dissim)
summary(dissim)
df = as.data.frame(dissim)
df$siteRichness = rowSums(tru.y)
plot(df$siteRichness, df$V1, ylab = "Dissimilarity", xlab = "Number of Species Observed at a Site")
rownames(df) = rownames(tru.y)
write.csv(df, "data/Peru_6.5_Dissimilarity_per_Site_ONLYSpeciesObserved_NotAllSpecies.csv")
######################################

######################################


## Calculate the AUC for each species
## Based on Probability of Occurrence
## For random selections, is probability higher where species is found than is not found? 
met.df = sapply(1:ncol(tru.y), function(i) Metrics::auc(tru.y[,i], pred.df.scaled[,i]))
# met.df = sapply(1:ncol(tru.y), function(i) pROC::auc(tru.y[,i], pred.df[,i]))
auc.mean = mean(met.df)
length(which(met.df>0.7))
hist(pred.df[,1])
hist(tru.y[,1])

# species_df$AUC.scaled = met.df
write.csv(species_df, "data/Peru_99_Species_Site_and_Performance.csv", row.names = F)

auc.vals <- matrix(NA, length(species), 1)
for (j in 1:length(species)) {
  # print(j)
#  for (i in 1:N) {
    auc.vals[j,1] <- pROC::auc(response = c(tru.y[,j]), predictor = pred.df.scaled[,j])
  } # i (species)
mean(auc.vals)
summary(auc.vals)
length(which(auc.vals>0.7))
#####################

###############################################
### Compare Dissimilarity with Site Variables
env = read.csv("data/BIRDS_2_Stations_UTM19_inc_SPECIES.csv")
env = select(env, 1:6, 28:30, 10, 15, 17:18, 26, 9, 19:20, 22:25)
names(env) = c("site_N", "site_c", "station", "lon", "lat", "utm_zone", "coords", "site_station", "stat_coords", "alt", "dist_Large_river", "dist_lake_swamp", "dist_2nd_de_forest", "habitat", "prim_2nd", "dist_sml_rds", "dist_highway", "deforest_5km","deforest_2km", "deforest_1km", "human_traffic_1to10")
env$prim_2nd[env$prim_2nd==1] = "primary"
env$prim_2nd[env$prim_2nd==0] = "secondary"
env$prim_2nd[is.na(env$prim_2nd)] = "not_noted"

env = select(env, site_station, habitat, prim_2nd)

df2 = as.data.frame(cbind(sites,dissim))
df2$site_station = gsub(pattern = "_20.*$",replacement = "",x = df2$sites)

df3 = left_join(df2, env, by="site_station")
df3 = df3[!is.na(df3$habitat),]
df3$V2 = as.numeric(df3$V2)
df3$env = paste0(df3$habitat,"_",df3$prim_2nd)
df4 = df3[df3$prim_2nd!="not_noted",]
boxplot(df3$V2 ~ df3$prim_2nd)
boxplot(df4$V2 ~ df4$env)
### Dissimilarity Does Not seem to Vary with Habitat/Forest type ###
####################################################

####################################################
### Compare AUC values with Bird Traits ###
trait = read.csv("data/BIRDS_3.5_Bird_TRAITS.csv")
auc = as.data.frame(cbind(species, auc.vals))
trait = trait[trait$species %in% auc$species,]
plot_df = left_join(trait,auc, by = "species")
plot_df$V2 = as.numeric(plot_df$V2)

boxplot(plot_df$V2~plot_df$Primary.Lifestyle)
boxplot(plot_df$V2~plot_df$Habitat)


