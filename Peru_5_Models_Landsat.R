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
# write.csv(ref2, "data/Birds_4.5_Survey_Reflectance_Values_perSITEperYEAR.csv", row.names = F)
# length(unique(ref$site_station))

ord = match(colnames(ssPnt1[,,1]), rownames(ref2))
ref2 = ref2[ord,]


######################################################################
### Find Which Variables have lowest Correlations with all others  ###
ref3 = ref2[,-c(1,28)]   ## Remove year & station
ref_cor = cor(x = ref3)  ## Create correlation matrix of all vars
ref_cor_rem = ref_cor
ref_cor_rem[upper.tri(ref_cor)] = 0 ## Remove duplicates by turning upper triangle and diganol into 0
diag(ref_cor_rem) = 0

## Remove highly correlated variables
data_new <- ref3[ , !apply(ref_cor_rem,                # From original dataframe !apply = do not select
                           2,                          # 2 = By column  
                           function(x) any(x > 0.8))]  # If any value in the correlation column > 0.8 

data_new = cbind(ref2[,c(1,28)], data_new)             # Add year & station back
head(data_new)    
write.csv(data_new, "data/Birds_4_Survey_Reflectance_Values_SITESbyYEAR_UNCORRELATED.csv", row.names = F)
#########################################################################################################

## can't have NA in Occurence covariates

### Replicates
cc=c()
for (i in 1:nrow(ssPnt1[1,,])) {    # For each row (Station)
  # print(length(which(!is.na(ssPnt1)[1,i,])))
  cc = append(cc, length(which(!is.na(ssPnt1)[1,i,])))  ## get the number of Non-NA (replicates that were conducted) 
}
hist(cc, breaks = 35)
summary(cc)          # Median=2; Mean=3
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
  
### DETECTION COVARIATES
### (Vary With Each Survey)
am.pm = read.csv("data/BIRDS_3_Det-Covs_AM-PM_per_site_Replicate.csv", row.names = 1)
#### Make sure Det.covs have same column names as y
colnames(am.pm) = colnames(y.new[1,,])
or = match(colnames(y.new[,,1]), rownames(am.pm))
am.pm = am.pm[or,]
am.pm = as.matrix(am.pm)
am.pm.use = am.pm[,1:ncol(y.new[1,,])]
  
  # year = read.csv("data/BIRDS_3_Det-Covs_YEAR_per_site_Replicate.csv", row.names = 1)
  # #### Make sure Det.covs have same column names as y
  # colnames(year) = colnames(y.new[1,,])
  # or = match(colnames(y.new[,,1]), rownames(year))
  # year = year[or,]
  # year = as.matrix(year)
  # year.use = year[,1:ncol(y.new[1,,])]
  
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
  
  det.covs = list(method = method.use, am.pm = am.pm.use) # , year = year.use
  
  
### COORDINATES
coords = select(occ, "lon", "lat", "site_station")
XY = coords[,1:2]                                                # Select only the Projeted Coordinates
rownames(XY) = rownames(occ)
  
#####
##### - PUT ALL VARIABLES IN A LIST - #####
data = list(y = y.new, occ.covs = occ.covs, det.covs = det.covs, coords = XY)
data = list(y = y.new, occ.covs = ref2)
str(data)
is.na(det.covs)
range(is.na(data$occ.covs))

N <- dim(data$y)[1]  ## Number of Species

##################################
### 1 -- OCCURRENCE ONLY
## specify the number of threads to use (n.omp.threads), the number of MCMC samples (n.samples), the amount of samples to discard as burn-in (n.burn), the thinning rate (n.thin), and arguments to control the display of sampler progress (verbose, n.report).

names(ref2)
data$occ.covs = select(data$occ.covs, TCB_SD_500,TCG_SD_500,TCW_SD_500)

s = Sys.time()
out.occ <- msPGOcc(occ.formula = ~ scale(TCB_SD_500) + scale(TCG_SD_500) + scale(TCW_SD_500), 
                   det.formula = ~ 1,
                   data = data, 
                   # inits = ms.inits, 
                   n.samples = 7500, 
                   # priors = ms.priors, 
                   n.omp.threads = 1, 
                   verbose = TRUE, 
                   n.report = 100, 
                   n.burn = 2500,
                   n.thin = 5, 
                   n.chains = 2)
print(Sys.time()-s)
filename = paste0("models/OCCURRENCE_StationYear_allReps_TasseledCapBGW.Rdata")
save(out.occ, file = filename)










