#####
### OCCUPANCY MODEL BIRDS USING spOCC

library(spOccupancy)
library(coda)        # for MCMC diagnostics
# library(stars)       # for Plotting 
# library(ggplot2)

####################################
####################################      MULTI - SPECIES - OCCUPANCY - MODELS        #########################################
####################################                                                  #########################################

## spOccupancy uses nearly identical syntax for fitting multi-species occupancy models as it does for single-species models and provides the same functionality for posterior predictive checks, model assessment and selection using WAIC and k-fold cross-validation, and out of sample prediction.
# load("data/COUNTS_Site-X-REPLICATE-per-Species_ARRAY.Rdata")
load("data/BIRDS_3_COUNTS_NETS_Site-X-REPLICATE-per-Species_ARRAY.Rdata")

# Turn counts to Presence/Absence
ssNet[ssNet>0]=1
table(ssNet) # Quick summary 1009132 Non-observations - 7928 Observations
str(ssNet)
sp.names = rownames(ssNet)
colnames(ssNet[,,1])
# #######
# ## Remove those species occurring in fewer than 2 sites
# rem = c()
# for (i in 1:335) {   ## For each species
# #   # i = 100
# df = ssNet[i,,]
# site = df
# site = rowSums(df, na.rm = T)
# site[site>0]=1
# if (sum(site)<3) {
# if (sum(colSums(df,na.rm = T))<2) {
# rem = append(rem, i)
# }
# print(sum(df, na.rm = T))
# }
# }
# ssNet = ssNet[-c(rem),,]
# sp.names = sp.names[-c(rem)]

# Set most common and Uncommon species as first 2 in y DF to improve use of latent variables in Species Correlation Occupancy Model
o = order(apply(ssNet, 1, sum, na.rm = TRUE), decreasing = T)
o = o[c(1, length(o), 3:length(o)-1)]
sp.ordered = sp.names[o]
y.new <- ssNet[sp.ordered, , ]

# ### DETECTION COVARIATES
# am.pm = read.csv("data/Det-Covs_AM-PM_per_site_Replicate.csv", row.names = 1)
# 
# #### Make sure Det.covs have same column names as y
# colnames(am.pm) = colnames(y.new[1,,])
# or = match(colnames(y.new[,,1]), rownames(am.pm))
# am.pm = am.pm[or,]
# am.pm = as.matrix(am.pm)
# det.covs = list(am.pm = am.pm)


### COORDINATES
coords = read.csv("data/BIRDS_2_COORDINATES.csv")                   # Read in coordinate data
coords$Site_Station = paste0(coords$site_code, "_", coords$station) # Create a unique Site_Station name variable
keep = which(coords$Site_Station %in% colnames(ssNet[,,1]))           # Find which Station Names occur in the Count Array
coords = coords[keep,]                                              # Select only coordinates where the names are in the count array
rownames(coords) = NULL                                             # Reset rownames starting from 1
or = match(x = colnames(ssNet[,,1]), table = coords$Site_Station)     # Get the order in which the Count names match the coord names 
coords = coords[or,]                                                # Order the coords to be the same order as counts 
XY = coords[,4:5] 


#####
##### - PUT ALL VARIABLES IN A LIST - #####
data = list(y = y.new, coords = XY)
# det.covs = det.covs, occ.covs = occ.covs, 
str(data)


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
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

#########################
### 1 -- INTERCEPT ONLY
## specify the number of threads to use (n.omp.threads), the number of MCMC samples (n.samples), the amount of samples to discard as burn-in (n.burn), the thinning rate (n.thin), and arguments to control the display of sampler progress (verbose, n.report).
s = Sys.time()
out.intercept <- msPGOcc(occ.formula = ~ 1, 
                         det.formula = ~ 1, 
                         data = data, 
                         inits = ms.inits, 
                         n.samples = 6000, 
                         priors = ms.priors, 
                         n.omp.threads = 1, 
                         verbose = TRUE, 
                         n.report = 100, 
                         n.burn = 1000,
                         n.thin = 20, 
                         n.chains = 4) #,
# k.fold = 4, k.fold.threads = 4, k.fold.only = F)
ss = Sys.time()-s
print(ss)
save(out.intercept, file = "models/1-INTERCEPT_NETS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
load("models/1-INTERCEPT_NETS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
waicOcc(out.intercept)
# summary(out.intercept)
##### Create & Save data
## Posterior Predictive Check
ppc.intercept = ppcOcc(out.intercept, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.intercept)
# Extract the fitted values
y.rep.intercept <- fitted(out.intercept)$y.rep.samples

# Calculate AUC using the detection-nondetection values -------------------
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out.intercept$n.post * out.intercept$n.chains
N <- nrow(data$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- pROC::auc(response = c(y.new[i, ,]), predictor = c(y.rep.intercept[j, i, , ]))
  } # i (species)
} # j (iteration)
intercept_points = list(ppc.intercept, y.rep.intercept, auc.vals)
# save(intercept_points, file = "1-INTERCEPT_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")
# load("model_outputs/1-INTERCEPT_POINTS_6000its_4ch_20thn_1000brn_1000smpl.Rdata")

# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
rownames(auc.quants) = sp.ordered
write.csv(auc.quants, "AUC_Intercept_NETS.csv")
# auc.q.burn = read.csv("model_outputs/AUC_Intercept_POINTS.csv", row.names = 1)

####### Analyse Data ########
# summary(out.intercept)
# mean(out.intercept$k.fold.deviance)
# waicOcc(out.intercept)
## Posterior Predictive Check
# Get the actual Bayesian p-values for individual species and community
# Species-level
# species.bpvs.intercept <- apply(ppc.intercept$fit.y.rep > ppc.intercept$fit.y, 2, mean)
# species.bpvs.intercept
# Community-level NB: THIS IS LITERALLY THE MEAN OF ALL SPECIES
# comm.bpv <- mean(ppc.intercept$fit.y.rep > ppc.intercept$fit.y)
# 
# str(y.rep.intercept)
# auc.quants
# mean(auc.quants[,2])
# sd(auc.quants[,2])
# range(auc.quants[,2])
# hist(auc.quants[,2])
# length(which(auc.quants[,2]>0.7))

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
    tp.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.intercept[j, i, , ]) == 1), na.rm = TRUE)
    tn.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.intercept[j, i, , ]) == 0), na.rm = TRUE)
    fp.samples[j, i] <- sum((c(y.new[i, , ]) == 0) & (c(y.rep.intercept[j, i, , ]) == 1), na.rm = TRUE)
    fn.samples[j, i] <- sum((c(y.new[i, , ]) == 1) & (c(y.rep.intercept[j, i, , ]) == 0), na.rm = TRUE)
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
sens = mean(sensitivity.samples)
spec = mean(specificity.samples)
tss = sens+spec-1





