library(readxl)
library(reshape2)
library(tidyr)
library(dplyr)
library(rgdal)
library(openxlsx)
library(lubridate)

# ##### READ DATA IN - PERU BIRDS #####
# #####################################################################################
points_in = read_excel("data/RAW New2 - FF Cleaned Master Bird Database - COPY.xlsx", "Point Counts")
nets_in = as.data.frame(read_excel("data/RAW New2 - FF Cleaned Master Bird Database - COPY.xlsx", "Banding"))
station_data = as.data.frame(read_excel("data/RAW New2 - FF Cleaned Master Bird Database - COPY.xlsx", "19L Station Data", range = cell_rows(1:748))) # Specify rows as they are all blank below but read by R
net_effort = read_excel("data/RAW New2 - FF Cleaned Master Bird Database - COPY.xlsx", "Daily Effort")

# ### FILTER OUT SIGHTINGS WITH NO LOCATION and from Quellomayo ###
# ### ie - no Coordinate Value AND no Station Value
# ##       and Remove un required columns
pi = points_in
pi = filter(pi, !is.na(Station) | !is.na(Easting))
pi = filter(pi, Site_Code  != "QUE" & Site_Code  != "Que")
names(pi)
pi = select(pi, 1:6, 9, 11:16, 18, 21)



ni = nets_in
ni$Easting[is.na(ni$Northing)]=NA
ni$Northing[is.na(ni$Easting)]=NA
ni = filter(ni, !is.na(Station) & Station!="NA" & Station!="-" & Station!="-9" |
              !is.na(Easting) & Easting!="NA" & !is.na(Northing) & Northing!="NA")
ni = filter(ni, Site_Code  != "QUE" & Site_Code  != "Que")
names(nets_in)
ni = select(ni, 7, 9, 28, 31:34, 41:42)

st = filter(station_data, !is.na(Easting) & Easting!="NA" & !is.na(Northing) & Northing!="NA" & Northing != "no data")
names(st)
st = select(st, 1:5, 8:9, 14:25)


# ##### There are issues with Station Names being different across the Count & Station
# ### Data sheets. I will create unique identifiers for each Station in Count & Station DFs
# ### and Merge the 2 data sets by the unique identifier


# ## Add Location IDs from Joined Site & Station
# ## Add Coordinate IDs from Joined Easting & Northing
ni$coord_ID = paste0(ni$Easting,"_",ni$Northing)
ni$location = paste0(ni$Site_Code,"_",ni$Station)
length(unique(ni$coord_ID)) # 594 Unique Coordinates in Nets
length(unique(ni$location)) # 834 Unique Locations in Nets

pi$coord_ID = paste0(pi$Easting,"_",pi$Northing)
pi$location = paste0(pi$Site_Code,"_",pi$Station)
length(unique(pi$coord_ID)) # 201 Unique Coordinates in Points
length(unique(pi$location)) # 215 Unique Locations in Points

st$location = paste0(st$Site_Code,"_",st$Station)
st$coord_ID = paste0(st$Easting,"_",st$Northing)
length(unique(st$coord_ID)) # 730 Unique Coordinates in Stations (Although 742 Station Lines)
# This is because there are a few Stations that share coordinate's, I think because they changed names, or have multiple names for whatever reason?
length(unique(st$location)) # 741 Unique Locations in Station - FPE1 in IBAA is replicated but with 2 different coordinates


n = ni
p = pi

# Add a Net_Point Column to n & p for future differentiation
n$net_point = as.character("net")
p$net_point = as.character("point")

# Net counts dont have a "count" column, so I will add one to include a numeric count value
# All values will be one because each sighting as a single capture
# Some Point Count values are NA, so I will give them all as a count of 1
n$Count = as.numeric(1)
p$Count[is.na(p$Count)]=1 # NA's in count need to be 1

# Add an am_pm column to Both as sometimes there were surveys done in the morning and again in afternoon, and this will be needed to define which effort data information gets included
str(n)
n$am_pm = format(n$Capture_Time, "%p") # Create column by formatting the time column %p value
unique(n$am_pm) # New column has AM, PM and NA
n$am_pm[is.na(n$am_pm)]="UK" # Replaces NA with UK (unknown)

str(p)
p$am_pm = format(p$`Start Time`, "%p") # Create column by format time column %p value
unique(p$am_pm) # New column has AM, PM and NA
p$am_pm[is.na(p$am_pm)]="UK" # Replaces NA with UK (unknown)



# Add Pnt, Net or St in front of all Point, Net and Station column names
# for ease of differentiation once merged
colnames(p) = paste0("Pnt_",colnames(p))
# Add un prefixed Location & Coordinate columns for use in merging
p$location = p$Pnt_location
p$coord_ID = p$Pnt_coord_ID

colnames(n) = paste0("Net_",colnames(n))
# Add un prefixed Location & Coordinate columns for use in merging
n$location = n$Net_location
n$coord_ID = n$Net_coord_ID

colnames(st) = paste0("St_",colnames(st))
# Add un prefixed Location & Coordinate columns for use in merging
st$location = st$St_location
st$coord_ID = st$St_coord_ID



# #####     Transform all Counts to Wide     #####
pw = p
# Add a unique Site_Station_Date_am/pm column, that identifies individual surveys
pw$unique = paste0(pw$location,"_",pw$Pnt_Date,"_",pw$Pnt_am_pm)
# Calculate the count of each species for each survey
pw = group_by(pw, unique, Pnt_Species) %>%
  mutate(Count = sum(Pnt_Count)) %>%
  ungroup()
# Remove space between species names so that they can become column headers
pw$Pnt_Species = gsub(pw$Pnt_Species, pattern = " ", replacement = "_")
# Remove rows with Duplicated Species & Survey 
pw1 = distinct(pw, unique, Pnt_Species, .keep_all = T)
# Pivot the table wider to put the Species names as columns based on the unique column
pw2 = pivot_wider(pw1, names_from = Pnt_Species,values_from = Pnt_Count,values_fill = 0)
names(pw2)
length(unique(pw2$unique)) # There are 408 individual Survey rows

#### DO THE SAME FOR NET COUNTS ####
####################################
####                            ####




# Add Columns Based Solely on LOCATION and NOT INDIVIDUAL SURVEYS
# Then merge them with Station Data by the Location column
pw3 = pw2
pw3$location = paste0(pw3$Pnt_Site_Code,"_",pw3$Pnt_Station)
pul_st = left_join(pw3, st, by = "location")
# Keep those with matches to make a full DF later
pu_L_matched = filter(pul_st, !is.na(St_location)) # 1684 of 2407

# Remove those with matches to try and match the remaining by coordinate
puc = filter(pul_st, is.na(St_location)) # 723 of 2407

# Remove the Station Data and rejoin again by Coordinate ID
names(puc)
puc2 = select(puc, -266:-287)
puc2$coord_ID = puc2$Pnt_coord_ID
puc_st = left_join(puc2, st, by = "coord_ID") 
# Keep those with matches to add to the Matched Df started from Matched Locations
pu_c_matched = filter(puc_st, !is.na(St_coord_ID) & St_coord_ID != "NA_NA") # 633 of 723

names(pu_c_matched)
pu_c_matched = select(pu_c_matched, -18:-20, -266:-271, -286:-288)
names(pu_L_matched)
pu_L_matched = select(pu_L_matched, -18:-20, -266:-270, -285:-287)

# Bind the Matched Data together to get a Matched Data Set
P_matched = bind_rows(pu_L_matched,pu_c_matched)
write.csv(P_matched, file = "data/Points_with_Station.csv", row.names = F)

# ## Remove those with matches - the remainder have no reciprocal station data
# There are 90 surveys, over 17 Locations, without any attached Station Data
p_no_st = filter(puc_st, is.na(St_coord_ID) | St_coord_ID == "NA_NA") # 90 Surveys
p_no_st_un = distinct(p_no_st, Pnt_location, .keep_all = T) # 17 Locations
write.csv(p_no_st, file = "data/Points_without_reciprocal_Station.csv", row.names = F)








p_st_loc = left_join(p,stco, by="location")
a = filter(p_st_loc, !is.na(Station.y)) # 2202 Lines of 3110 have reciprocal station data
b = filter(p_st_loc, is.na(Station.y)) # 908 Lines with no reciprocal data

# Take the lines with no location reciprocal to see if there is a reciprocal from coordinates.
c = select(b, 1:13,17) %>% rename(coord_ID = "coord_ID.x")
d = left_join(c,stco,by="coord_ID") %>% filter(coord_ID != "NA_NA")
# Make sure all column names are the same between DFs to allow rbind
names(a)
a1 = select(a, 1:13, 17)
names(a1) = gsub(pattern = ".x", replacement = "", x = names(a1))
names(d)
d1 = select(d, 1:13, 18)
names(d1) = gsub(pattern = ".x", replacement = "", x = names(d1))
d1 = rename(d1, Habitat = "Habitat.y")

e = rbind(a1,d1) # I now have all elements of Point Count Data joined with Station Data where either a Unique Location or Coordinate name is reciprocated 
length(unique(e$coord_ID)) # 201 = Same as at start
length(unique(e$location)) # 204 = original was 215

# e has 3029 Observations and the Original counts (without Quell & NAs was 3110)
# Select the 81 rows that are missing to analyse them
aj = anti_join(p,e,by="location")
# Reduce the 81 lines to only 1 per site = 11 Missing Sites
aju = select(aj, 1:2, 4:6, 12:13) %>% distinct(Station, .keep_all = TRUE)
# Create 
puloc = distinct(p, location, .keep_all = T)
write.csv(puloc, file = "data/pointsIN_unique_Locations.csv", row.names = F)
# pucod = distinct(p, coord_ID, .keep_all = T)
# pu = rbind(puloc, pucod) %>% distinct()
# all_equal(target = puloc ,current = pu, ignore_col_order = T,ignore_row_order = T)

write.csv(aj, file = "data/point_no_reciprocal_station.csv", row.names = F)
write.csv(e, file = "data/point_all_with_reciprocal_station.csv", row.names = F)
write.csv(p, file = "data/point_base_input.csv", row.names = F)
write.csv(stco, file = "data/station_base_input.csv", row.names = F)



# ##### COUNT DATA #####




# ##### Create an All Count Data Frame #####
n_count = select(nets_in, unique_id, Site_Code, Station, coord_ID, Date, am_pm, net_point, Count, Species)
p_count = select(points_in, unique_id, Site_Code, Station, coord_ID, Date, am_pm, net_point, Count, Species)

# ##### Group by Unique Survey then Species to sum Count Per Species Per Survey
n_count2 = group_by(n_count, unique_id, Species, Site_Code, Station, coord_ID, Date, am_pm, net_point) %>%
  summarise(sum(Count)) %>% rename(count = "sum(Count)") %>% ungroup()
p_count2 = group_by(p_count, unique_id, Species, Site_Code, Station, coord_ID, Date, am_pm, net_point) %>%
  summarise(sum(Count)) %>% rename(count = "sum(Count)") %>% ungroup()
# ## rbind them to form 1 dataframe ##
all_counts = rbind(n_count2,p_count2)
all_counts_wide = pivot_wider(all_counts, names_from = Species, values_from = count, values_fill = 0)
dupes = all_counts_wide[which(duplicated(all_counts_wide[,c("unique_id")])==T),] # Should be zero!

# ### Visualize Some Count Histograms ###
g= all_counts_wide[,7:ncol(all_counts_wide)]
g = g[, colSums(g != 0) > 0] # Removes species with zero sightings

range(colSums(g))# Shows Range of counts per species
hist(log(colMeans(g>0))) # Proportion of site_dates each species was encountered
hist(colSums(g)) # Actual overall count of each species
hist(rowSums(g>0)) # Species Richness per site_date & the Frequency of that amount of richness

write.csv(all_counts_wide, file = "data/WrangleNEW_Bird_Counts_All_Wide.csv", row.names = F)

# ######################################
# ##### Spatial Coordinates #####

# Create Site Locations from Mean Coordinates
###########################################################################
# ### CONVERT UTM to LON(X Axis) LAT(Y Axis) ###
library(proj4)
# Rename erroneous CoOrd Data with NA and remove all NA as they can't be included in conversion
station_data$Easting[station_data$Easting=="NA" | station_data$Easting=="no data"]=NA
station_data$Northing[station_data$Northing=="NA" | station_data$Northing=="no data"]=NA
utm1 = station_data
utm1 = filter(utm1, !is.na(Easting)) # Remove NA as they vant be included in conversion

# Select just the Coordinate Data and make sure it's a numeric format
utm1 = select(utm1,"Easting", "Northing")
utm1$Easting = as.numeric(utm1$Easting)
utm1$Northing = as.numeric(utm1$Northing)

sputm = SpatialPoints(utm1, proj4string = CRS("+proj=utm + zone=19L + south + datum=WGS84 +units=m"))
spgeo = spTransform(sputm, CRS("+proj=longlat + datum=WGS84"),inverse=T) # I think Inverse=T refers to Southern Hemisphere
xy=as.data.frame(spgeo@coords)
names(xy)
xy = rename(xy, x = Easting, y = Northing)
plot(xy)
# #####################################################################################








all_counts_wide$site_station = paste0(all_counts_wide$Site_Code,"_",all_counts_wide$Station)
station_data$site_station = paste0(station_data$Site_Code,"_",station_data$Station)

# There should be same number of unique station locations in Station and Count data?
# But the are a lot more Count Data stations than in the Site Data
length(unique(all_counts_wide$site_station)) #1045
length(unique(station_data$site_station)) #746

a = full_join(station_data, all_counts_wide, by = "site_station")
names(a)
b=select(a, 26:27, 1:6, 24, 28:32)
write.csv(b, file = "data/StationMERGECount.csv", row.names = F)

