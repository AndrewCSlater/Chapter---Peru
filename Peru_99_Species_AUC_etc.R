library(dplyr)

prev = read.csv("data/Peru_99_Species_Prevealence_AUC_Influences.csv")
# prev = read.csv("data/Peru_99_Species_Prev_AUC_VarInf_Trait.csv")
x1 = read.csv("data/Peru_99_Species_AUC_Explanatory_INTERCEPT.csv")
x2 = read.csv("data/Peru_99_Species_AUC_Predicted_INTERCEPT.csv")
x3 = read.csv("data/Peru_99_Species_AUC_Explanatory_5PCA.csv")
x4 = read.csv("data/Peru_99_Species_AUC_Predicted_5PCA.csv")
x5 = read.csv("data/Peru_99_Species_AUC_Explanatory_DEFORESTATION.csv")
x6 = read.csv("data/Peru_99_Species_AUC_Predicted_DEFORESTATION.csv")


df = cbind(x1,x2,x3,x4)
d# f = cbind(x5, x6); df = df[,c(1,2,4)]; names(df) = c("species", "auc.exp.deforest", "auc.pred.deforest")
df = df[,c(1,2,4,6,8)]
names(df) = c("species", "auc.exp.int", "auc.pred.int", "auc.exp.5pca", "auc.pred.5pca")
prev = full_join(prev, df, by="species")
rm(x1,x2,x3,x4)

### Get significant relationships of OCCURRENCE (Beta) variables per species 
## ie how many have a 95%CI that excludes zero?
load("models/1.B-Training_Min11sp_5PCA_50000_25000_50_1.Rdata")
# load("models/1.B-Training_Min11sp_BiomasHabitat_HansenDeforest_50000_25000_50_1_MODEL.Rdata")

bet = (out.sfMsPGOcc$beta.samples)

sig_sp_band = matrix(data = NA, nrow = length(out.sfMsPGOcc$sp.names), ncol = length(out.sfMsPGOcc$x.names)-1, dimnames = list(out.sfMsPGOcc$sp.names, out.sfMsPGOcc$x.names[-1])) ## To store the species names that have significant relationships with band reflectance

sps = length(out.sfMsPGOcc$sp.names) ## Number of Species
grp2 = c(1:sps)
for (jj in 1:ncol(sig_sp_band)) {    ## For each variable
  #jj=1
  sps2 = sps*jj + grp2    
  
  sigYN = c()              ## Create YN significance vector
  sigVAL = c()             ## Create Influence Value vector (for if variable significantly influences a species)
  
  bb1 = as.data.frame(bet[,sps2])     ## Subset df into all species with 1 variable
  for (j in 1:ncol(bb1)) {            ## For each species (column)
    #j=1
    i = quantile(bb1[,j], c(0.025,0.975))  ## Get the 2.5 & 97.5 quants
    if (i[1]<0 & i[2]>0) {ii = 0           ## If the quants surround zero Give the Value of 0
    } else {ii = 1}                        ## Else (all 95% is either above or below zero) give a value of 1 
    sigYN = append(sigYN, ii)              ## Add that 1 or 0 value to the yes no vector: get 1 or 0 per Species
    sigVAL = append(sigVAL, colMeans(bb1[j])) ## Get mean value per species
    sigVAL = sigVAL * sigYN                ## Multiply mean by 1 or 0 to keep only significant mean values
  }  
  sig_sp_band[,jj] = sigVAL
}
sig_sp_band = as.data.frame(sig_sp_band)
sig_sp_band$species = rownames(sig_sp_band)
############ Repeat the Process for Detection Variables
alp = (out.sfMsPGOcc$alpha.samples)

sps = length(out.sfMsPGOcc$sp.names) ## Number of Species
sig_sp_det = matrix(data = NA, nrow = length(out.sfMsPGOcc$sp.names), ncol = (ncol(out.sfMsPGOcc$alpha.samples)/sps)-1, dimnames = list(out.sfMsPGOcc$sp.names, out.sfMsPGOcc$x.p.names[-1])) ## To store the species names that have significant relationships with band reflectance

grp2 = c(1:sps)
for (jj in 1:ncol(sig_sp_det)) {    ## For each variable
  #jj=1
  sps2 = sps*jj + grp2    
  
  sigYN = c()              ## Create YN significance vector
  sigVAL = c()             ## Create Influence Value vector (for if variable significantly influences a species)
  
  bb1 = as.data.frame(alp[,sps2])     ## Subset df into all species with 1 variable
  for (j in 1:ncol(bb1)) {            ## For each species (column)
    #j=1
    i = quantile(bb1[,j], c(0.025,0.975))  ## Get the 2.5 & 97.5 quants
    if (i[1]<0 & i[2]>0) {ii = 0           ## If the quants surround zero Give the Value of 0
    } else {ii = 1}                        ## Else (all 95% is either above or below zero) give a value of 1 
    sigYN = append(sigYN, ii)              ## Add that 1 or 0 value to the yes no vector: get 1 or 0 per Species
    sigVAL = append(sigVAL, colMeans(bb1[j])) ## Get mean value per species
    sigVAL = sigVAL * sigYN                ## Multiply mean by 1 or 0 to keep only significant mean values
  }  
  sig_sp_det[,jj] = sigVAL
}
sig_sp_det = as.data.frame(sig_sp_det)
sig_sp_det$species = rownames(sig_sp_det)


# prev = full_join(prev, sig_sp_det, by="species")
prev = full_join(prev, sig_sp_band, by="species")

###########  Add TRAIT data
trait = read.csv("data/BIRDS_3.5_Bird_TRAITS.csv")
trait = trait[trait$species%in%prev$species,]
prev = full_join(prev, trait, by="species")

write.csv(prev, "data/Peru_99_Species_Prev_AUC_VarInf_Trait.csv", row.names = F)

# prev = full_join(sig_sp_band, prev, by="species")
plot(prev$auc.pred.5pca, prev$auc.pred.deforest)

########################
# Look for correlations between variables
prev = read.csv("data/Peru_99_Species_Prev_AUC_VarInf_Trait.csv")
names(prev)

### AUC & Prevalence
summary(lm(prev$auc.pred.5pca ~ prev$sitesSeenAt, na.action = na.exclude))

plot(prev[,2:22])
c = cor(prev[,c(2:22,26,27,35)], use = "complete.obs")






