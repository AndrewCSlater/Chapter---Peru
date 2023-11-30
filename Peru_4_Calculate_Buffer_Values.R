library(dplyr)
library(raster)
library(sf)
# library(rgdal)
# library(terra)


tifs = list.files(path = "Satellite_data/Landsat")

####################################
# --- Load the Trap Point Data
site = read.csv("Data/Birds_3_Surveys.csv")
traps = site[site$utm_zone=="19L",]

yrs = sort(unique(traps$year))   # If calcualting for surveyed years
# yrs = c(2004:2020)                 # If calculating for every year

ref = c() ## To store all years reflectance values in

for (i in 1:length(yrs)) {
  # i=1

srvy = traps[traps$year==yrs[i],] ## To get Station Reflectance only for the year of survey
# srvy = traps                ## To get station reflectance for every year (Needed for spOccupancy Multi-Season models)

sat = grep(pattern = paste0(yrs[i],".tif$"), x = tifs) # This selects the file names that end in .tif & excludes the other extensions
filename = paste0("Satellite_data/Landsat/",tifs[sat])
sat = raster::brick(filename)  

xy = srvy[, 9:10]
spdf <- SpatialPointsDataFrame(coords = xy, data = srvy,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))

b = c(500) # Buffer radii

for (ii in b) {
  # ii=500
  # Landsat Satellite Imagery
  p = spTransform(spdf, crs(sat))
  lsat_mean = extract(x = sat,
                      y = p,
                      buffer=ii,
                      fun=mean,
                      df=T)
  lsat_mean = dplyr::select(lsat_mean, -"ID")
  colnames(lsat_mean) = paste0(colnames(lsat_mean),"_mean_",ii)
  
  lsat_SD = extract(x = sat,
                    y = p,
                    buffer=ii,
                    fun=sd,
                    df=T)
  lsat_SD = dplyr::select(lsat_SD, -"ID")
  colnames(lsat_SD) = paste0(colnames(lsat_SD),"_SD_",ii)
  

lst = cbind(lsat_mean, lsat_SD)
} # End buffer Mean & SD loop for current year
lst$year = yrs[i]
lst$rep_ID = srvy$rep_ID
ref = rbind(ref, lst)
} # End years loop

write.csv(ref, "data/Birds_4_Survey_Reflectance_Values.csv", row.names = F)
write.csv(ref, "data/Birds_4_Survey_Reflectance_Values_ALL_YEARS.csv", row.names = F)
write.csv(ref, "data/Birds_4_Survey_Reflectance_Values_SITESbyYEAR.csv", row.names = F)
##########################################################
### CREATE GLCM VALUES FOR ALL BANDS & INICES & BUFFER ###
library(glcm)

### Creating GLCM over a huge raster may be excessively resource heavy, so crop rasters to points bounding envelope.
glcm_df=c()
names(sat)
# 1. Create a single layer raster from of each band of the Lsat image in turn 
for (band in c(1:length(names(sat)))) { # Chosen bands in the Raster Brick
  band = "TCW"
  ST = Sys.time()
  img = sat[[band]]
  n = names(img)
  
  # 2.Create a Raster stack of the 3 GLCM statistic values based on the current raster band
  glcm = glcm(img,
              window = c(3,3), 
              shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
              n_grey = 16, # 16 = 4-Bit, 64 = 6-Bit, 256 = 8-Bit
              statistics = c("mean", "contrast", "entropy"))
  
  plot(glcm[[1]])
  # 3.Calculate Buffer Values for the 3 GLCM statistics   
  p = spTransform(spdf, crs(glcm))
  
  for (j in b) {
    lsat_mean = extract(x = glcm,
                        y = p,
                        buffer=j,
                        fun=mean,
                        df=F)
    colnames(lsat_mean) = paste0("lsat_", n,"_mean_",j,"m_", colnames(lsat_mean))
    
    lsat_SD = extract(x = glcm,
                      y = p,
                      buffer=j,
                      fun=sd,
                      df=F)
    colnames(lsat_SD) = paste0("lsat_", n,"_SD_",j,"m_", colnames(lsat_SD))
    
    glcm_df = cbind(glcm_df, lsat_mean, lsat_SD)
  } # End Buffer extraction of GLCM
  
  filename2 = paste0("Data/golaTRAPpoints_GLCM_Landsat_2022_Jan_30-500m.csv")
  write.csv(glcm_df, filename2, row.names = F)  
  
} # End creation of all GLCM

####################################################################
### CREATE CANONICAL COMPONENTS FOR ALL GROUPS and BUFFER LEVELS ###
library(sgdm)
library(gdm)
library(dplyr)

# Load Survey Data
survey = read.csv("Data/XY_No_Satellite_Gola.csv", row.names = 1)
rownames(survey) = 1:nrow(survey)
survey$SITE_ID = rownames(survey)

# Create Count Data
y = survey[,c(20:ncol(survey))]
y[y>0]=1

# As per Pichler 2021 - Remove OTU's with fewer than 3 occurrences
# He Removes Sites with fewer than 4 OTU's - But I keep all sites that have any insects
df1 = survey
r = which(rowSums(y)>0) # Rows to keep
df1 = df1[c(r),]
c = which(colSums(y)>2) # Columns to keep
c = c+19
df1 = df1[,c(1:19,c)]
y = df1[,c(20:ncol(df1))]
y[y>0]=1

# Load Satellite data
lsat = read.csv("Data/golaTRAPpoints_Landsat_2022_March_30_500m.csv")
lsat = lsat[rownames(lsat)%in%rownames(df1),]

sat_glcm = read.csv("Data/golaTRAPpoints_GLCM_Landsat_2022_Jan_30-500m.csv")
sat_glcm = sat_glcm[rownames(sat_glcm)%in%rownames(df1),]

# hansP = read.csv("Data/golaTRAPpoints_HansenYear_maxMin_30_500m.csv")
# hansP = hansP[rownames(hansP)%in%rownames(df1),]
# # Set the Hansen data to the Max value of 30m buffer level
# hans = hansP[c(r),] # Remove the rows that are not in the count data
# hans=hans[,1]
# 
radar = read.csv("Data/golaTRAPpoints_Sentinel1_2022_q1_30_500m.csv")
radar = radar[rownames(radar)%in%rownames(df1),]

# distance = df1$distance_f

# Create Predictor (bands) and Response (count) Data Sets
df2 = df1[, c("SITE_ID", "Longitude", "Latitude")]
df2 = rename(df2, Plot_ID = SITE_ID, X=Longitude, Y=Latitude)

df2$Plot_ID = as.numeric(df2$Plot_ID) # NB: - The plot_id needs to be numeric and unique

insect = cbind(df2$Plot_ID, y) # Include plot_id in 1st column
colnames(insect)[1] = "Plot_ID"

# Split data into Bands & Indices
band_names = c("blue", "green", "red", "NIR", "SWIR", "MIR", "RE1", "RE2", "RE3")

b = c(30, 60, 90, 120, 150, 250, 500) # Buffer radii

for (buffer in b) {
  #buffer=30
  # Filter by buffer size
  raw_buffer = lsat[,grep(paste(buffer,collapse="|"),colnames(lsat))]
  # raw_buffer = radar[,grep(paste(buffer,collapse="|"),colnames(radar))]
  glcm_buffer  = sat_glcm[,grep(paste(buffer,collapse="|"),colnames(sat_glcm))]
  
  # Divide Raw & GLCM by Band & Index
  raw_bands = raw_buffer[,grep(paste(band_names,collapse="|"),colnames(raw_buffer))]
  raw_indices = raw_buffer[,-grep(paste(band_names,collapse="|"),colnames(raw_buffer))]
  glcm_bands = glcm_buffer[,grep(paste(band_names,collapse="|"),colnames(glcm_buffer))]
  glcm_indices = glcm_buffer[,-grep(paste(band_names,collapse="|"),colnames(glcm_buffer))]
  
  satellite_used = list(raw_bands, raw_indices, glcm_bands, glcm_indices)
  satellite_used_names = c("lsat_bands", "lsat_indices", "lsat_glcm-bands", "lsat_glcm-indices")
  
  # Create all combinations of groups - I already have individual Bands & Indices & combined glcm 0f the 15 possible combinations, so am creating the other 12
  a = list(satellite_used[[1]], satellite_used[[2]], satellite_used[[3]], satellite_used[[4]], cbind(satellite_used[[1]], satellite_used[[2]]), cbind(satellite_used[[1]], satellite_used[[3]]), cbind(satellite_used[[1]], satellite_used[[4]]), cbind(satellite_used[[2]], satellite_used[[3]]), cbind(satellite_used[[2]], satellite_used[[4]]), cbind(satellite_used[[3]], satellite_used[[4]]), cbind(satellite_used[[1]], satellite_used[[2]], satellite_used[[3]]), cbind(satellite_used[[1]], satellite_used[[2]], satellite_used[[4]]), cbind(satellite_used[[1]], satellite_used[[3]], satellite_used[[4]]), cbind(satellite_used[[2]], satellite_used[[3]], satellite_used[[4]]), cbind(satellite_used[[1]], satellite_used[[2]], satellite_used[[3]], satellite_used[[4]]))
  
  bb = c(satellite_used_names[[1]], satellite_used_names[[2]], satellite_used_names[[3]], satellite_used_names[[4]], paste0(satellite_used_names[[1]], satellite_used_names[[2]]), paste0(satellite_used_names[[1]], satellite_used_names[[3]]), paste0(satellite_used_names[[1]], satellite_used_names[[4]]), paste0(satellite_used_names[[2]], satellite_used_names[[3]]), paste0(satellite_used_names[[2]], satellite_used_names[[4]]), paste0(satellite_used_names[[3]], satellite_used_names[[4]]), paste0(satellite_used_names[[1]], satellite_used_names[[2]], satellite_used_names[[3]]), paste0(satellite_used_names[[1]], satellite_used_names[[2]], satellite_used_names[[4]]), paste0(satellite_used_names[[1]], satellite_used_names[[3]], satellite_used_names[[4]]), paste0(satellite_used_names[[2]], satellite_used_names[[3]], satellite_used_names[[4]]), paste0(satellite_used_names[[1]], satellite_used_names[[2]], satellite_used_names[[3]], satellite_used_names[[4]]))
  
  
  for (s in 1:15) {  ## For each combination of satellite data
    #  s=15
    sat = a[[s]]
    
    # There must not be any NA in either matrix
    # This replaces any NA in the satellite data with the column (Band) mean value
    for(ii in 1:ncol(sat)){
      sat[is.na(sat[,ii]), ii] <- mean(sat[,ii], na.rm = TRUE)
    }
    
    bbb = bb[[s]] # Current Group Name
    
    
    ########################################
    # Create 1 Canon from the 4 Radar values
    # sat = raw_buffer
    
    # The band data must have variation and cannot have a zero Standard Deviation
    # Find columns with SD of 0 (There is no variability in the column)
    r = which(apply(sat,2,sd)==0)
    
    use = if(is.na(r[2])) {
      sat
    } else {
      sat[,-c(r)]
    }
    
    bands = cbind(df2, use) # Include plot_id and X Y in first 3 columns
    
    ### RUN THE SGDM PROCESS
    # 1 - Test models using a range of penalisation values from 0.6 - 1.0 on both 
    # predictor and response variables - the output is a 5x5 vector of RMSE values
    mod = sgdm.param(predData = bands, bioData = insect, k = 10)
    gc()
    file = paste0("Data/sgdm/ModParams_gola_k10_",bbb,"_",buffer,".csv")
    # file = paste0("Data/sgdm/ModParams_gola_k1_RADAR_",buffer,".csv")
    write.csv(mod, file = file)
    
    # 2 - Retrieve the best performing model based on the highest vector value -
    # (using the 'sgdm.best' function, with  output option 'm')
    #  This provides the appropriate penalisation values to be used in the model
    best_mod = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "m", k = 10)
    
    # 3 - Extract Sparse Canonical Components relating to the best model
    # (using the 'sgdm.best' function, with  output option 'c')
    best_mod_cpnts = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "c", k = 10)
    file = paste0("Data/sgdm/components/Can_Cpnts_k10_",bbb,"_",buffer,".csv")
    # file = paste0("Data/sgdm/components/Can_Cpnts_k1_RADAR_",buffer,".csv")
    write.csv(best_mod_cpnts, file = file, row.names = F)
    
    ## Check variable contributions
    # 8 - Extract the canonical vectors of the best model -
    # (using the 'sgdm.best' function, with output option 'v')
    best_mod_vects = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "v", k = 10)
    
    # Add the band names to the vector
    best_mod_vects = as.data.frame(best_mod_vects)
    rownames(best_mod_vects) = colnames(bands[,-c(1:3)])
    best_mod_vects$band = colnames(bands[,-c(1:3)])
    
    # Save this value to use to calculate Canonical Components across all pixels of an image
    file = paste0("Data/sgdm/vectors_to_create_Can_Cpnts_k10_",bbb,"_",buffer,".csv")
    # file = paste0("Data/sgdm/vectors_to_create_Can_Cpnts_k1_RADAR_",buffer,".csv")
    write.csv(best_mod_vects, file = file, row.names = T)
    gc()
  } # End of Each Satellite Group Loop
} # End of Each Buffer Size 



