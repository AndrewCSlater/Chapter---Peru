library(spOccupancy)
library(dplyr)
library(vegan)
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
table(reps$site)                                    ### How many Stations per Site

### For each site, select 2 stations that have been surveyed at least 3 times
rep.red = reps[reps$replicates>2,]                        # Reduce data to stations with more than 2 replicates
rep.red <-  rep.red %>% group_by(site) %>% filter(n()>1)  # Reduce data to include only sites with more than 1 station left
samp = rep.red %>% group_by(site) %>% sample_n(size = 2)  # Group data by site and select 2 random stations

div = which(colnames(ssPnt)%in%samp$station)              # What are the count column numbers of the stations that are in the randomly selected stations

test = ssPnt[,div,]    # Create Test data by selecting those stations from count data
test.12 = test[,,1:12] # Select the first 12 replicates from Test data
test.y = as.data.frame(rowSums(test.12, dims = 2, na.rm = T)) # Create Species * Station Test Data

# ## Remove those species with less than 1% of (ie more than 10) STATIONS OCCURRENCE after removing stations & replicates
# rem = which(!rownames(test.12) %in% rownames(train.12))
# test.12 = test.12[-rem,,]

## Fitted Model
load("models/1.B-Training_Min11sp_5PCA_50000_25000_50_1.Rdata")
df = out.sfMsPGOcc$y # 135 Species in 1101 Stations over 12 Replicates
# df[1,,] # Station * Replicate ie. Per Species
# df[,1,] # Species * Replicate ie. Per Station
# df[,,1] # Species * Station   ie. Per Replicate

### Create Training Species * Station Data from Model Inputs
train.y = as.data.frame(rowSums(df, dims = 2, na.rm = T))
train.y = train.y[order(rownames(train.y)),]

### Remove Species from Test dt NOt in Train data
test.y = test.y[rownames(test.y)%in%rownames(train.y),]
test.y = test.y[order(rownames(test.y)),]


### Get Effort per station in terms of Number of Replicates
effort.train = rowSums(!is.na(df[1,,]))      
effort.test = rowSums(!is.na(test.12[1,,]))

### Get Mean detection for training data
alp = as.data.frame(out.sfMsPGOcc$alpha.samples) # Get modeled probability of detection
det = alp[1:nrow(train.y)]                       # Select the intercept values
det = plogis(colMeans(det))                      # Convert to probability from logit

### Create a probability of Detection Per Species Per Station, based on Number of Replicates
eff = matrix(data = NA, nrow = length(det), ncol = length(effort.train))
for (i in 1:length(effort.train)) {    # For Each Station
  # i=1
  power = effort.train[i]              # Get the Number of replicates
  v = 1 - (1 - det)^power              # Calculate probability of Detection given number of replicates
  eff[,i] = v
}  
# write.csv(eff, "model_outputs/Det_Prob_Explanatory.csv", row.names = T)
#################################
## Species Accumulation Curves on Training & Test Data

acu.train = specaccum(comm = t(train.y), method = "random", w = effort.train, permutations = 1000)
which(acu.train$richness>0.95*max(acu.train$richness))[1] # 95% of species were found within 132 stations
which(acu.train$richness>0.75*max(acu.train$richness))[1] # 75% of species were found within 47 sites

acu.test = specaccum(comm = t(test.y), method = "random", w = effort.test, permutations = 1000)
which(acu.test$richness>0.95*max(acu.test$richness))[1] # 95% of species were found within 29 stations
which(acu.test$richness>0.75*max(acu.test$richness))[1] # 75% of species were found within 16 sites

###########################
#### Get Modelled Species Accumulations
# This is on Training Data

# create layout

filename = paste0("images/Spec_Accum_Train-Explain_Reps.png")
png(filename=filename, res=600, width = 7, height = 5, units = "in")
par(mar = c(2,2,1.5,1))
layout(matrix(c(1, 2, 3, 4), nrow = 2, 
              ncol = 2, byrow = TRUE))


latent = out.sfMsPGOcc$z.samples

pred=c()
# png(filename="images/Spec_Accum_Train-Explain.png", res=600, width = 11.69, height = 8.27, units = "in")
plot(acu.train, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlim = c(0,1000), ylim = c(0,140), 
     xlab="Number of Stations",
     ylab="Accumulated Species") # xlim = c(0,1000), xvar = "effort",

for (i in 1:nrow(latent)) { ## For each posterior sample
# i=1
  d = matrix(data = rbinom(length(eff), size = 1, prob=eff), nrow = nrow(eff), ncol = ncol(eff))
  z = latent[i,,]
  z = z * d
  sp1 = specaccum(comm = t(z), method = "random", permutations = 200) # w = effort, Use effort per station to improve
  pred = cbind(pred,sp1$richness)
}
mx = apply(X = pred, MARGIN = 1, FUN = max)
mn = apply(X = pred, MARGIN = 1, FUN = min)
md = apply(X = pred, MARGIN = 1, FUN = median)
xx = c(1:length(mx), length(mn):1)
yy = c(mx, rev(mn))
polygon(x = xx, y = yy, density = 20)
lines(md, type = "l", lwd=2, add=T)
# dev.off()

#####################
#####################

## Do the same for test data
load("model_outputs/Predicted_Testdata_Occupancy.Rdata")
load("model_outputs/Predicted_Testdata_Detection.Rdata")
latent = predicted.pca5$z.0.samples
det.pred = predicted.pca5.DETECTION$p.0.samples
det.pred = colMeans(det.pred, dims = 1, na.rm = T)

### Create a probability of Detection Per Species Per Station, based on Number of Replicates
for (i in 1:ncol(det.pred)) {    # For Each Station
  # i=1
  power = effort.test[i]              # Get the Number of replicates
  v = 1 - (1 - det.pred[,i])^power    # Calculate probability of Detection given number of replicates
  det.pred[,i] = v
} 

# write.csv(det.pred, "model_outputs/Det_Prob_Predicted.csv", row.names = T)

# png(filename="images/Spec_Accum_Test-Predicted.png", res=600, width = 11.69, height = 8.27, units = "in")
plot(acu.test, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab="Number of Stations",
     ylab="Accumulated Species", ylim = c(0,140))

pred2=c()
for (i in 1:nrow(latent)) { ## length(latent) For each posterior sample
  # i=1
  d = matrix(data = rbinom(length(det.pred), size = 1, prob=det.pred), nrow = nrow(det.pred), ncol = ncol(det.pred))
  z = latent[i,,]
  z = z * d
  sp1 = specaccum(comm = t(z), method = "random", permutations = 100) # w = effort, Use effort per station to improve
pred2 = cbind(pred2, sp1$richness)  
}
mx = apply(X = pred2, MARGIN = 1, FUN = max)
mn = apply(X = pred2, MARGIN = 1, FUN = min)
md = apply(X = pred2, MARGIN = 1, FUN = mean)
xx = c(1:length(mx), length(mn):1)
yy = c(mx, rev(mn))
polygon(x = xx, y = yy, density = 20)
lines(md, type = "l", lwd = 2)
# dev.off()

###############################
###############################

### Create Accumulations Curve for a Single Station that has multiple replicates
stn = names(which.max(effort.test)) # Get Station Name with most replicates in test data
stn.occ = ssPnt[,stn,1:12]              # Get the first 12 replicate surveys of that station from full count data
stn.occ = stn.occ[which(rownames(stn.occ)%in%rownames(train.y)),] # Remove those species not included in the model

acu.station = specaccum(comm = t(stn.occ), method = "random", permutations = 1000)
# plot(acu.station, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", ylim = c(0,50),  xlab="Number of Replicates",
#      ylab="Accumulated Species")
# 
# ### get matrix of Predicted latent occurrences for 1 Station 
# ii = which(colnames(predicted.pca5$psi.0.samples[1,,])==stn) ## Get position of the station name in Pred data
# stn.pred = predicted.pca5$z.0.samples[,,ii] ## Select the latent predictions for that station
# ii = order(colnames(predicted.pca5$psi.0.samples))
# stn.pred = stn.pred[,ii]
# 
# ### Get predicted probability of Detection for each species at that station
# det.pred = predicted.pca5.DETECTION$p.0.samples
# det.pred = colMeans(det.pred, dims = 1, na.rm = T)
# ii = which(colnames(det.pred)==stn) ## Get position of the station name in Det.Pred data
# stn.det = det.pred[,ii]             ## Select the Detection predictions for that station
# ii = order(names(stn.det))          ## Put species in Alphabetical order
# stn.det = stn.det[ii]



# png(filename="images/Spec_Accum_by_Replicates.png", res=600, width = 11.69, height = 8.27, units = "in")
plot(1, type = "n", xlab = "Replicate Surveys", 
     ylab = "Number of Species", xlim = c(0, 100),  
     ylim = c(0, 60)) 

### Perform SpecAccum for that station
stn.df = matrix(data = NA, nrow = nrow(stn.occ), ncol = ncol(stn.occ)*10) ## Create a df to store Predicted Pres/Abs per REPLICATE for EACH POSTERIOR SAMPLE
pred3=c()
for (i in 1:500) { ## nrow(stn.pred) For each posterior sample
  # i=1
  for (j in 1:120) {     ## For each REPLICATE
  # j=1
  d = rbinom(length(stn.det), size = 1, prob=stn.det) # Get the LATENT DETECTION
  z = stn.pred[i,]                                    # Get LATENT OCCUPANCY
  z = z * d                                           # MULTIPLY TOGETHER
  stn.df[,j] = z
  }
  sp1 = specaccum(comm = t(stn.df), method = "random", permutations = 100) # w = effort, Use effort per station to improve
  pred3 = cbind(pred3, sp1$richness)
# plot(sp1, ci=0, ci.type="poly", ci.col="lightgrey", lwd=1, add = T) # xvar="effort",
}
mx = apply(X = pred3, MARGIN = 1, FUN = max)
mn = apply(X = pred3, MARGIN = 1, FUN = min)
md = apply(X = pred3, MARGIN = 1, FUN = mean)
xx = c(1:length(mx), length(mn):1)
yy = c(mx, rev(mn))
polygon(x = xx, y = yy, density = 20)
lines(md, type = "l", lwd=2)
plot(acu.station, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab="Number of Replicates",
     ylab="Accumulated Species", add = T)
dev.off()















