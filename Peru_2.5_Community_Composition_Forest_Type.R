utm19 = read.csv("data/BIRDS_2_Stations_UTM19_inc_SPECIES.csv")
keep = which(complete.cases(utm19[31:ncol(utm19)]))
utm19 = utm19[keep,]

primary = utm19[utm19$Primary_1_Secondary_0==1 & !is.na(utm19$Primary_1_Secondary_0) ,c(31:ncol(utm19))] # 181
primary = as.matrix(colSums(primary, na.rm = T))
secondary = utm19[utm19$Primary_1_Secondary_0==0 & !is.na(utm19$Primary_1_Secondary_0) ,c(31:ncol(utm19))] # 196
secondary = as.matrix(colSums(secondary, na.rm = T))

flood = utm19[utm19$Habitat=="Floodplain" & !is.na(utm19$Habitat) ,c(31:ncol(utm19))] # 440
flood = as.matrix(colSums(flood, na.rm = T))
dryland = utm19[utm19$Habitat=="Terra Firme" & !is.na(utm19$Habitat) , c(31:ncol(utm19))] # 194
dryland = as.matrix(colSums(dryland, na.rm = T))

#############################################################
#### Number of 637 stations with different variables
length(which(!is.na(utm19$Habitat))) # All 637 have Habitat (Flood/Dry)
length(which(!is.na(utm19$Deforestation.within.5.km.of.site.centroid..ha.))) # 414 have Deforestation level
length(which(!is.na(utm19$Primary_1_Secondary_0))) # 377 have Forest Category - All of which also include deforestation level - calculated on next line
length(which(!is.na(utm19$Deforestation.within.5.km.of.site.centroid..ha.) & !is.na(utm19$Primary_1_Secondary_0)))
#############################################################


pridry = utm19[utm19$Primary_1_Secondary_0==1 & !is.na(utm19$Primary_1_Secondary_0) & utm19$Habitat=="Terra Firme" ,c(31:ncol(utm19))] # 17
pridry = as.matrix(colSums(pridry, na.rm = T))
priwet = utm19[utm19$Primary_1_Secondary_0==1 & !is.na(utm19$Primary_1_Secondary_0) & utm19$Habitat=="Floodplain" ,c(31:ncol(utm19))] # 164
priwet = as.matrix(colSums(priwet, na.rm = T))
twodry = utm19[utm19$Primary_1_Secondary_0==0 & !is.na(utm19$Primary_1_Secondary_0) & utm19$Habitat=="Terra Firme" ,c(31:ncol(utm19))] # 107
twodry = as.matrix(colSums(twodry, na.rm = T))
twowet = utm19[utm19$Primary_1_Secondary_0==0 & !is.na(utm19$Primary_1_Secondary_0) & utm19$Habitat=="Floodplain" ,c(31:ncol(utm19))] # 89
twowet = as.matrix(colSums(twowet, na.rm = T))

sp.comp = cbind(primary, secondary, flood, dryland)
colnames(sp.comp) = c("primary", "secondary", "flood", "dry")
# sp.comp[sp.comp>0]=1
sp.comp = t(sp.comp)

sp.comp2 = cbind(pridry, priwet, twodry, twowet)
colnames(sp.comp2) = c("pridry", "priwet", "twodry", "twowet")
# sp.comp2[sp.comp2>0]=1
sp.comp2 = t(sp.comp2)

### Number of species per Habitat
rowSums(sp.comp)
##  Primary = 225, 2ndary = 246; Flood = 319, Dry = 252
rowSums(sp.comp2)
## pridry = 85, priwet = 218, twodry = 196, twowet = 203

### Calculate community similarities between Primary/Secondary/Floodplain/TerraFirma
1-vegan::vegdist(x = sp.comp, method = "jaccard")  ### NB: do 1- to get Similarity , not DISsimilarity
1-vegan::vegdist(x = sp.comp2, method = "jaccard")

#############
## Shannon Diversity Indices
## NB: Revert back to Origianl Full matrices, not the summaries used above
dry = dryland
dry[is.na(dry)] = 0
wet = flood
wet[is.na(wet)] = 0
one = primary
one[is.na(one)] = 0
two = secondary
two[is.na(two)] = 0

dry = vegan::diversity(dry)
mean(dry) ## 2.28 2.09
wet = vegan::diversity(wet)
mean(wet) ## 2.17 1.86
one = vegan::diversity(one)
mean(one) ## 2.46 2.19
two = vegan::diversity(two)
mean(two) ## 2.15 1.96

t.test(dry,wet)
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
length(which(sp.comp[1,]==1 & sp.comp[2,]==1)) # 177 Species are common to both Primary and Secondary forests 
length(which(sp.comp[3,]==1 & sp.comp[4,]==1)) # 219 Species are common to both Primary and Secondary forests 

length(which(sp.comp2[1,]==1 & sp.comp2[2,]==1)) # 78 Species are common to both Primary Dry & Primary Wet forests
length(which(sp.comp2[1,]==1 & sp.comp2[3,]==1)) # 70 Species are common to both Primary Dry & Secondary Dry forests
length(which(sp.comp2[1,]==1 & sp.comp2[4,]==1)) # 75 Species are common to both Primary Dry & Secondary Wet forests
length(which(sp.comp2[1,]==1 & (sp.comp2[4,]==1 | sp.comp2[3,]==1 | sp.comp2[2,]==1))) # 81 Species are common to both Primary Dry & One of the other habitats
length(which(sp.comp2[1,]==1 & (sp.comp2[4,]==0 & sp.comp2[3,]==0 & sp.comp2[2,]==0))) # 81 Species are common to both Primary Dry & One of the other habitats
length(which(sp.comp2[3,]==1 & sp.comp2[4,]==1)) # 153 Species are common to both Secondary WET & DRY forests 







