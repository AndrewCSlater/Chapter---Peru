library(spOccupancy)

load("models/SpLF3_DETonly_12reps_buffer500m_MODEL.Rdata")
# SpLF3_DetOcc_PCA3_12reps_buffer500m_MODEL
load("models/SpLF3_DETonly_12reps_buffer500m_FITTED_SAMPLES.Rdata")
# SpLF3_DetOcc_PCA3_12reps_buffer500m_FITTED_SAMPLES

out = out.sfMsPGOcc
rm(out.sfMsPGOcc); gc()

# species = out$sp.names


### Calculate the latent Z occupancy for each species at each site
### This is the mean across all posterior samples
latZ = colMeans(out$z.samples, dims = 1) # Creates a 358 * 1135 matrix = Species by Stations
## Rows are species

auc = matrix(data = NA, nrow = 358, ncol = 400, dimnames = list(out$sp.names))

for (zz in 1:400) { # For each posterior Latent Occurence Z
#zz=1
latZ = out$z.samples[zz,,]  
# str(latZ)

### Compare the latent occupancy with the Fitted Model samples
# str(y.rep.samples)
for (smp in 1:400) { # for each posterior sample
# smp=1
f = y.rep.samples[smp,,,1] # select the first replicate of that sample 
for (s in 1:358) {    # For each species in that sample
# s=1
  a = Metrics::auc(actual = latZ[s,], predicted = f[s,])
  auc[s,smp] = a
} # End Species loop
} # End Posterior sample loop 
aucMean = as.data.frame(aucMean)
species = cbind(species, aucMean)
} # End Posterior Z Latent Occupancy loop
species = species[,-1]

aucMean = rowMeans(species)
matplot(t(species), type = "l")


aucSD = apply(auc, MARGIN = 1, FUN = sd)

mean(aucMean)
length(which(aucMean>0.6))

