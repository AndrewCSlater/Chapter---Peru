library(dplyr)
library(spOccupancy)


### REFLECTANCE COVARIATES
ref = read.csv("data/Birds_4_Survey_Reflectance_Values_SITESbyYEAR.csv")
ref$site_station = stringr::str_extract(ref$rep_ID, "[^_]*_[^_]*")   # Get the station name from the rep ID
ref$site_station = paste0(ref$site_station, "_", ref$year)           # Create station name by adding the year

# Reduce the reflectance data to one row for each site_station_year
ref = select(ref, -rep_ID, -year)
ref = ref %>% group_by(site_station) %>%
  summarise_all(mean) %>%
  ungroup()
rownames(ref) = ref$site_station

##########################################
### Reduce Dimensions of Reflectance Data
### Create a PCA vectors
df = select(ref, -"site_station")
#calculate principal components
pca_ref = prcomp(x = df, scale. = T)
#reverse the signs
pca_ref$rotation = -1 * pca_ref$rotation
# display principal components
# pca_ref$rotation

#reverse the signs of the scores
pca_ref$x = -1*pca_ref$x

# summary(pca_ref)

# display the first six scores
# head(pca_ref$x)

# biplot(pca_ref, choices = c(3,5),  scale = 0)


### FOREST COVARIATES
biome = read.csv("data/Birds_4_Survey_Biomas_Habitat_SITESbyYEAR.csv")
defst = read.csv("data/Birds_4_Survey_Hansen_Deforestation_Values_SITESbyYEAR.csv")
frst = cbind(biome[c(4,1)], defst[c(1:4)])

frst2 = frst %>% group_by(site_station) %>%
  mutate(habitat = habitat, forest = deforestation_0_m) %>% 
  distinct(site_station, habitat, forest) 

frst3 = frst[,c(1,4:6)] %>% group_by(site_station) %>% # This is number of m2 deforestation in a buffer - Buffer = Pi*Radii^2 spm 
  summarise_all(mean) %>%                             # 1km buffer = 314.16Ha
  ungroup()                                           # 2km = 1256.64Ha
                                                      # 5km = 7853.98Ha
frst3$deforestation_1000_m = frst3$deforestation_1000_m/10000/314.16   # Turns m2 to Ha then to proportion deforested
frst3$deforestation_2000_m = frst3$deforestation_2000_m/10000/1256.64
frst3$deforestation_5000_m = frst3$deforestation_5000_m/10000/7853.98
env = full_join(frst2, frst3, by = "site_station")

m = match(x = ref$site_station, table = env$site_station) # Put Stations in same order as reflectance
env = env[m,]

env$landscape = paste0(env$habitat,"_",env$forest)

env = select(env, site_station, habitat, forest)

#########################################################
### Plot PCAs grouped by Env Vars
# df2 = ref  #### Create new DF for Plot data
# 
# # env = read.csv("data/BIRDS_2_Stations_UTM19_inc_SPECIES.csv") # Load the Environmental data for the stations
# # env = select(env, 1:6, 28:30, 10, 15, 17:18, 26, 9, 19:20, 22:25)
# # names(env) = c("site_N", "site_c", "station", "lon", "lat", "utm_zone", "coords", "site_station", "stat_coords", "alt", "dist_Large_river", "dist_lake_swamp", "dist_2nd_de_forest", "habitat", "prim_2nd", "dist_sml_rds", "dist_highway", "deforest_5km","deforest_2km", "deforest_1km", "human_traffic_1to10")
# # env$prim_2nd[env$prim_2nd==1] = "primary"
# # env$prim_2nd[env$prim_2nd==0] = "secondary"
# # env$prim_2nd[is.na(env$prim_2nd)] = "not_noted"
# 
# #-----------------------
# 
# #-----------------------------
# 
# # env = select(env, site_station, habitat, prim_2nd)
# 
# df3 = left_join(df2, env, by = "site_station")
# # env2 = select(df3, site_station, year, habitat, prim_2nd)
# env2 = select(df3, site_station, year, habitat, forest)
# # env2$landscape = paste0(env2$habitat,"_",env2$prim_2nd)
# env2$landscape = paste0(env2$habitat,"_",env2$forest)
# # df3 = select(df3, -c(site_station, year, habitat, prim_2nd))
# df3 = select(df3, -c(site_station, year, habitat, forest))
# df3 = df
# 
# pca_df2 = prcomp(x = df3, scale. = T)
# #reverse the signs
# pca_df2$rotation = -1 * pca_df2$rotation
# #display principal components
# pca_df2$rotation
# 
# #reverse the signs of the scores
# pca_df2$x = -1*pca_df2$x

####### PLOT - Points Colour Coded by Group
g <- ggbiplot::ggbiplot(pca_ref,
                        choices = c(3,5),
                        obs.scale = 1,
                        var.scale = 1,
                        groups = env$landscape,
                        ellipse = TRUE,
                        circle = TRUE,
                        ellipse.prob = 0.68)
g <- g + scale_color_discrete(name = '')
g = g + theme_bw()
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')


png(filename="images/PCA1-2_HansenForest.png", res=600, width = 11.69, height = 8.27, units = "in")
print(g)
dev.off()

####### PLOT - Arrow length varies with loading
ft = factoextra::fviz_pca_var(pca_ref,axes = c(3,5),
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE     # Avoid text overlapping
)
ft = ft + theme_bw()

png(filename="images/PCA3-5_Band_Proportion.png", res=600, width = 11.69, height = 8.27, units = "in")
print(ft)
dev.off()

# figure <- ggpubr::ggarrange(g, ft,
#                     labels = c("A", "B"),
#                     ncol = 1, nrow = 2)
# figure

########################################################################################

# calculate total variance explained by each principal component
# round(pca_ref2$sdev^2 / sum(pca_ref2$sdev^2) ,digits = 2)

# create scree plot
# qplot(c(1:26), var_explained) +
#   geom_line() +
#   xlab("Principal Component") +
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot") +
#   ylim(0, 1)

#############################
### PCA Constituents

pca_ref$scale
pca_ref$center
pca_ref$sdev

cont = pca_ref$rotation^2 # Get influence of each band per the PCA component...square it to get rid of negative values
cont = round(prop.table(cont, margin = 2), digits = 2)  ## Get Proportion of Each Band per PCA Component

colSums(cont)
max.bands = unique(rownames(cont)[apply(cont,2,which.max)])  ## Get Maximum Contributor per component


pn = pca_ref$rotation  ## Turn proportional values pos or neg
pn[pn<0] = -1
pn[pn>0] = 1
cont = cont*pn

write.csv(cont, "data/Peru_6_PCA_Band_Contribution.csv", row.names = T)


pca1 = as.data.frame(sort(cont[,1], decreasing = T))
pca2 = as.data.frame(sort(cont[,2], decreasing = T))
pca3 = as.data.frame(sort(cont[,3], decreasing = T))
pca4 = as.data.frame(sort(cont[,4], decreasing = T))
pca5 = as.data.frame(sort(cont[,1], decreasing = T))




