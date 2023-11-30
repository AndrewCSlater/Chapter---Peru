library(dplyr)

counts = read.csv("data/BIRDS_2_Station-ANNUALLY_X_Species.csv")
counts = counts[,-c(1:8)]

names(counts)[1] = "site_station"
counts = counts %>% group_by(site_station) %>% summarise_all(sum)

biome = read.csv("data/Birds_4_Survey_Biomas_Habitat_SITESbyYEAR.csv")
defst = read.csv("data/Birds_4_Survey_Hansen_Deforestation_Values_SITESbyYEAR.csv")
frst = cbind(biome[c(4,1)], defst[c(1:4)])

frst2 = frst %>% group_by(site_station) %>%
  mutate(habitat = habitat, forest = deforestation_0_m) %>% 
  distinct(site_station, habitat, forest) 

frst2 = frst2 %>% group_by(site_station) %>% filter(row_number()==1) %>% ungroup()

env = right_join(frst2, counts, by = "site_station")
env = env[!is.na(env$habitat),]



primary = env[env$forest=="primary" ,c(8:ncol(env))] # 1109
secondary = env[env$forest=="secondary" ,c(8:ncol(env))] # 26
primary = as.matrix(colSums(primary, na.rm = T))
secondary = as.matrix(colSums(secondary, na.rm = T))


flood = env[env$habitat=="flooded" ,c(8:ncol(env))] # 603
dry = env[env$habitat=="terraFirma" ,c(8:ncol(env))] # 489
agri = env[env$habitat=="agriculture" ,c(8:ncol(env))] # 43

flood = as.matrix(colSums(flood, na.rm = T))
dry = as.matrix(colSums(dry, na.rm = T))
agri = as.matrix(colSums(agri, na.rm = T))

pridry = env[env$forest=="primary" &  env$habitat=="terraFirma" ,c(8:ncol(env))] # 488
priwet = env[env$forest=="primary" &  env$habitat=="flooded" ,c(8:ncol(env))] # 601
priagri = env[env$forest=="primary" &  env$habitat=="agriculture" ,c(8:ncol(env))] # 20
twodry = env[env$forest=="secondary" &  env$habitat=="terraFirma" ,c(8:ncol(env))] # 1
twowet = env[env$forest=="secondary" &  env$habitat=="flooded" ,c(8:ncol(env))] # 2
twoagri = env[env$forest=="secondary" &  env$habitat=="flooded" ,c(8:ncol(env))] # 2


pridry = as.matrix(colSums(pridry, na.rm = T))
priwet = as.matrix(colSums(priwet, na.rm = T))
priagri = as.matrix(colSums(priagri, na.rm = T))

twodry = as.matrix(colSums(twodry, na.rm = T))
twowet = as.matrix(colSums(twowet, na.rm = T))
twoagri = as.matrix(colSums(twoagri, na.rm = T))

sp.comp = cbind(primary, secondary, flood, dry, agri)
colnames(sp.comp) = c("primary", "secondary", "flood", "dry", "agri")
# sp.comp[sp.comp>0]=1
sp.comp = t(sp.comp)

sp.comp2 = cbind(pridry, priwet, priagri, twodry, twowet, twoagri)
colnames(sp.comp2) = c("pridry", "priwet", "priagri",  "twodry", "twowet", "twoagri")
# sp.comp2[sp.comp2>0]=1
sp.comp2 = t(sp.comp2)

### Number of species per Habitat
rowSums(sp.comp)
##  Primary = 344, 2ndary = 132; Flood = 290, Dry = 264, Agri = 155
rowSums(sp.comp2)
## pridry = 258, priwet = 290, priagri = 71, twodry = 10, twowet = 7, twoagri = 7

### Calculate community similarities between Primary/Secondary/Floodplain/TerraFirma
1-vegan::vegdist(x = sp.comp, method = "jaccard")  ### NB: do 1- to get Similarity , not DISsimilarity
1-vegan::vegdist(x = sp.comp2, method = "jaccard")

#############
## Shannon Diversity Indices
## NB: Revert back to Origianl Full matrices, not the summaries used above
dryland = dry
dry[is.na(dry)] = 0
wet = flood
wet[is.na(wet)] = 0
one = primary
one[is.na(one)] = 0
two = secondary
two[is.na(two)] = 0
field = agri
field[is.na(field)] = 0

dryland = vegan::diversity(dryland)
mean(dryland) ## 4.72
wet = vegan::diversity(wet)
mean(wet) ## 4.75
one = vegan::diversity(one)
mean(one) ## 4.82
two = vegan::diversity(two)
mean(two) ## 4.68
field = vegan::diversity(field)
mean(field) ## 4.78

t.test(dryland,wet)
t.test(one,two)

## Calculate Permanova to find R2 of variation in community structure between Primary/Secondary, Floodplain/TerraFirma 
df = utm19[,c(9,26,31:ncol(utm19))]
df$Primary_1_Secondary_0[df$Primary_1_Secondary_0==0] = "secondary"
df$Primary_1_Secondary_0[df$Primary_1_Secondary_0==1] = "primary"
df$Primary_1_Secondary_0[is.na(df$Primary_1_Secondary_0)] = "unknown"
df[is.na(df)] = 0
rem = rowSums(df[,-c(1:2)])
rem = which(rem==0)
if(length(rem)>0) {df = df[-rem,]}
vegan::adonis2(df[,3:ncol(df)] ~ df$Habitat, permutations = 1000) # Permanova Test
df2 = df[df$Primary_1_Secondary_0!="unknown",]
vegan::adonis2(df2[,3:ncol(df2)] ~ df2$Primary_1_Secondary_0 + df2$Habitat, permutations = 1000) # Permanova Test


#############
## Species common between landscapes
length(which(sp.comp[1,]==1 & sp.comp[2,]==1)) # 121 Species are common to both Primary and Secondary forests 
length(which(sp.comp[3,]==1 & sp.comp[4,]==1)) # 216 Species are common to both Wet and Dry habitats
length(which(sp.comp[3,]==1 & sp.comp[5,]==1)) # 126 Species are common to both Wet and Agri habitats
length(which(sp.comp[4,]==1 & sp.comp[5,]==1)) # 120 Species are common to both Dry and Agri habitats
length(which(sp.comp[3,]==1 & sp.comp[4,]==1 & sp.comp[5,]==1)) # 107 Species are common to Wet & Dry & Agriculture habitats



length(which(sp.comp2[1,]==1 & sp.comp2[2,]==1)) # 215 Species are common to both Primary Dry & Primary Wet forests
length(which(sp.comp2[1,]==1 & sp.comp2[3,]==1)) # 47 Species are common to both Primary Dry & Primary Agri forests
length(which(sp.comp2[2,]==1 & sp.comp2[3,]==1)) # 60 Species are common to both Primary Wet & Primary Agri forests
length(which(sp.comp2[1,]==1 & sp.comp2[2,]==1 & sp.comp2[3,]==1)) # 46 Species are common to All Primary forest groups

length(which(sp.comp2[1,]==1 & sp.comp2[4,]==1)) # 4 Species are common to both Primary Dry & Secondary Dry forests






length(which(sp.comp2[1,]==1 & (sp.comp2[4,]==1 | sp.comp2[3,]==1 | sp.comp2[2,]==1))) # 81 Species are common to both Primary Dry & One of the other habitats
length(which(sp.comp2[1,]==1 & (sp.comp2[4,]==0 & sp.comp2[3,]==0 & sp.comp2[2,]==0))) # 81 Species are common to both Primary Dry & One of the other habitats
length(which(sp.comp2[3,]==1 & sp.comp2[4,]==1)) # 153 Species are common to both Secondary WET & DRY forests 







