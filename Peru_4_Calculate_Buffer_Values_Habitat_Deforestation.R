library(dplyr)
library(raster)
library(sf)
# library(rgdal)
# library(terra)

sat = raster::raster("Satellite_data/Hansen/Hansen_Yr_Forest_Loss.tif") 
sat[sat>0]=1
plot(sat)

####################################
# --- Load the Trap Point Data
site = read.csv("Data/Birds_3_Surveys.csv")
traps = site[site$utm_zone=="19L",]

yrs = sort(unique(traps$year))   # If calcualting for surveyed years
# yrs = c(2004:2020)                 # If calculating for every year

deforest = c() ## To store all years reflectance values in

for (i in 1:length(yrs)) {
  # i=1

srvy = traps[traps$year==yrs[i],] ## To get Station Reflectance only for the year of survey
# srvy = traps                ## To get station reflectance for every year (Needed for spOccupancy Multi-Season models)

xy = srvy[, 9:10]
spdf <- SpatialPointsDataFrame(coords = xy, data = srvy,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))

def = c()
b = c(0, 1000, 2000, 5000) # Buffer radii
for (ii in b) {
  # ii=1000
# Landsat Satellite Imagery
  p = spTransform(spdf, crs(sat))
  hansen = extract(x = sat,
                   y = p,
                   buffer=ii,
                   fun=sum,   ## Number of cells that are deforested
                   df=T)
  hansen = dplyr::select(hansen, -"ID")
  hansen = hansen *30 *30 ## 30*30 = Square metres per cell
  colnames(hansen) = paste0("deforestation_",ii,"_m")
if (length(def)==0) {def=hansen}
  else {def = cbind(def,hansen)}
  
  } # End buffer loop for current year
def$year = yrs[i]
def$rep_ID = srvy$rep_ID
deforest = rbind(deforest, def)
} # End years loop

deforest$deforestation_0_m[deforest$deforestation_0_m==0]="primary"
deforest$deforestation_0_m[deforest$deforestation_0_m!="primary"]="secondary"

deforest$site_station = stringr::str_extract(deforest$rep_ID, "[^_]*_[^_]*")
deforest$site_station = paste0(deforest$site_station, "_", deforest$year)

length(unique(deforest$site_station))

write.csv(deforest, "data/Birds_4_Survey_Hansen_Deforestation_Values_SITESbyYEAR.csv", row.names = F)

###############################
## Biomas Land Class
sat = raster::raster("Satellite_data/MapBiomas_Peru_Landclass/mapbiomas-peru-collection-10-2013.tif") 
plot(sat)
unique(values(sat))

site = read.csv("Data/Birds_3_Surveys.csv")
traps = site[site$utm_zone=="19L",]

yrs = sort(unique(traps$year))   # If calcualting for surveyed years
# yrs = c(2004:2020)                 # If calculating for every year

flood = c() ## To store all years reflectance values in

for (i in 1:length(yrs)) {
  # i=1
  
  srvy = traps[traps$year==yrs[i],] ## To get Station Reflectance only for the year of survey
  # srvy = traps                ## To get station reflectance for every year (Needed for spOccupancy Multi-Season models)
  
  xy = srvy[, 9:10]
  spdf <- SpatialPointsDataFrame(coords = xy, data = srvy,
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
  
# Biomas Imagery
    p = spTransform(spdf, crs(sat))
    habitat = extract(x = sat,
                   y = p,
                   # buffer=ii,
                   # fun=sum,   ## Number of cells that are deforested
                   df=T)
    habitat = dplyr::select(habitat, -"ID")
    colnames(habitat) = "habitat"

habitat$year = yrs[i]
habitat$rep_ID = srvy$rep_ID
flood = rbind(flood, habitat)
} # End years loop

flood[flood==3]="terraFirma"
flood[flood==6]="flooded"
flood[flood==33]="flooded"
flood[flood==21]="agriculture"
flood[flood==11]="flooded"
flood[flood==18]="agriculture"

flood$site_station = stringr::str_extract(flood$rep_ID, "[^_]*_[^_]*")
flood$site_station = paste0(flood$site_station, "_", flood$year)

write.csv(flood, "data/Birds_4_Survey_Biomas_Habitat_SITESbyYEAR.csv", row.names = F)

