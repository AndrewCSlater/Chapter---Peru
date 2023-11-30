library(readxl)
library(dplyr)
library(lubridate)
# library(reshape2)
# library(tidyr)
# library(rgdal)
# library(openxlsx)

###########   UNIFY SITE NAMES ACROSS BIRD DATA   #################

##### READ DATA IN - PERU BIRDS #####
#####################################################################################
####################
##  STATION DATA  ##
utm18 = read_excel("data/RAW New2 - FF Cleaned Master Bird Database.xlsx", "18L Station Data", range = cell_rows(1:41))
utm19 = read_excel("data/RAW New2 - FF Cleaned Master Bird Database.xlsx", "19L Station Data", range = cell_rows(1:772))

## Get Site Names & Codes
df = rbind(utm18[,1:2], utm19[,1:2])
# Look to see if any replcates with different specllings etc
sort(unique(df$`Site Name`))
length(unique(df$`Site Name`)) # 27 unique names
length(unique(df$`Site Code`)) # 27 unique codes
nc = paste0(df$`Site Name`,"_", df$`Site Code`) # combine names & codes to see if there are any anomalies
length(unique(nc))             # 27 Name-Code combinations
sort(unique(nc))
### STATION DATA NAMES ALL OK

### Reconcile with Site Names of COUNT Data
points_in = read_excel("data/RAW New2 - FF Cleaned Master Bird Database.xlsx", "Point Counts")
nets_in = as.data.frame(read_excel("data/RAW New2 - FF Cleaned Master Bird Database.xlsx", "Banding", guess_max = 10000)) # Guess max lets R sample more rows to guess the column format - It was importing Time Open as logi and losing data at default level
pdf = rbind(points_in[,1:2], nets_in[,32:33])
length(unique(pdf$`Site Name`)) # 34 Names from Net & Point Counts
length(unique(pdf$`Site Code`)) # 33 Codes from Net & Point Counts
cnc = paste0(pdf$`Site Name`, "_", pdf$`Site Code`)
length(unique(cnc))             # 34 Unique combinations of count site names & codes
## There is a disparity between number of Codes & Names & a difference between number of sites in Station & Count data 

### Look for naming anomalies
# Get unique count site names
un_ct_nm = sort(unique(pdf$`Site Name`))

# Find which unique names are not in Station site names
un_ct_nm[which(!un_ct_nm %in% df$`Site Name`)]
unique(df$`Site Name`)[which(!unique(df$`Site Name`) %in% pdf$`Site Name`)]
unique(df$`Site Name`)[which(unique(df$`Site Name`) %in% pdf$`Site Name`)]
sort(unique(c(unique(pdf$`Site Name`), unique(df$`Site Name`))))

#### From Inspection
## 1 - Remove "-9" as these count sites have no location
## 2 - Change "EI- Farm" & "EI Farm" to "EI-Farm"
## 3 - Change "Monte Amazonico" to "Monte Amazonica"

nets_in$`Site Name` = gsub(pattern = "Monte Amazonico", replacement = "Monte Amazonica", x = nets_in$`Site Name`)
nets_in$`Site Name` = gsub(pattern = "EI- Farm", replacement = "EI-Farm", x = nets_in$`Site Name`)
nets_in$`Site Name` = gsub(pattern = "EI Farm", replacement = "EI-Farm", x = nets_in$`Site Name`)
nets_in = nets_in[nets_in$`Site Name` != "-9", ]

points_in$`Site Name` = gsub(pattern = "Monte Amazonico", replacement = "Monte Amazonica", x = points_in$`Site Name`)
points_in$`Site Name` = gsub(pattern = "EI- Farm", replacement = "EI-Farm", x = points_in$`Site Name`)
points_in$`Site Name` = gsub(pattern = "EI Farm", replacement = "EI-Farm", x = points_in$`Site Name`)
points_in = points_in[points_in$`Site Name` != "-9", ]

### Re-examine the data
pdf = rbind(points_in[,1:2], nets_in[,32:33])
length(unique(pdf$`Site Name`)) # 31 Names from Net & Point Counts
length(unique(pdf$`Site Code`)) # 32 Codes from Net & Point Counts
cnc = paste0(pdf$`Site Name`, "_", pdf$`Site Code`)
length(unique(cnc))             # 32 Unique combinations of count site names & codes
sort(unique(cnc))
sort(unique(nc))

#### From Inspection
## 4 - Change CODE "EI-FARM" to "EI-Farm"
nets_in$`Site Code` = gsub(pattern = "EI-FARM", replacement = "EI-Farm", x = nets_in$`Site Code`)
points_in$`Site Code` = gsub(pattern = "EI-FARM", replacement = "EI-Farm", x = points_in$`Site Code`)

### Re-examine the data
pdf = rbind(points_in[,1:2], nets_in[,32:33])
length(unique(pdf$`Site Name`)) # 31 Names from Net & Point Counts
length(unique(pdf$`Site Code`)) # 31 Codes from Net & Point Counts
cnc = paste0(pdf$`Site Name`, "_", pdf$`Site Code`)
length(unique(cnc))             # 31 Unique combinations of count site names & codes
sort(unique(cnc))
sort(unique(nc))

sort(unique(cnc)[!unique(cnc) %in% unique(nc)])
sort(unique(nc)[!unique(nc) %in% unique(cnc)])

#### From Inspection
## 5 - Change CODE "AZUL" to "AZU"
## 6 - Change CODE "IBA" to "IBAA"
nets_in$`Site Code` = gsub(pattern = "AZUL", replacement = "AZU", x = nets_in$`Site Code`)
points_in$`Site Code` = gsub(pattern = "AZUL", replacement = "AZU", x = points_in$`Site Code`)
nets_in$`Site Code` = gsub(pattern = "IBA", replacement = "IBAA", x = nets_in$`Site Code`)
points_in$`Site Code` = gsub(pattern = "IBA", replacement = "IBAA", x = points_in$`Site Code`)

### Re-examine the data
pdf = rbind(points_in[,1:2], nets_in[,32:33])
length(unique(pdf$`Site Name`)) # 31 Names from Net & Point Counts
length(unique(pdf$`Site Code`)) # 31 Codes from Net & Point Counts
cnc = paste0(pdf$`Site Name`, "_", pdf$`Site Code`)
length(unique(cnc))             # 31 Unique combinations of count site names & codes
sort(unique(cnc))
sort(unique(nc))

sort(unique(cnc)[!unique(cnc) %in% unique(nc)])
## There are now 4 Sites in the COUNTS df that DO NOT occur in the STATIONS df - If no coordinates occur for these they will be removed
sort(unique(nc)[!unique(nc) %in% unique(cnc)])
## There are now no Site Name - Codes in the STATION df that do not exist in the COUNTS df's

pss = paste0(points_in$`Site Name`, "_", points_in$`Site Code`)
sort(unique(pss)[!unique(pss) %in% unique(nc)])
#### ALL POINT COUNT SITES & CODES are now represented in the STATION DATA

nets_in$`Station GPS (Easting)`[nets_in$`Site Name` == "Bamboo"]
nets_in$`Station GPS (Easting)`[nets_in$`Site Name` == "Pampas Del Heath"]
nets_in$`Station GPS (Easting)`[nets_in$`Site Name` == "CBP-ICA"]
nets_in$`Station GPS (Easting)`[is.na(nets_in$`Site Name`)]
### None of these Sites are in the STATION df or have their own Coordinates, so will be removed

nets_in = nets_in[nets_in$`Site Name` != "Bamboo" & nets_in$`Site Name` != "Pampas Del Heath" & nets_in$`Site Name` != "CBP-ICA" & !is.na(nets_in$`Site Name`),]
nss = paste0(nets_in$`Site Name`, "_", nets_in$`Site Code`)
sort(unique(nss)[!unique(nss) %in% unique(nc)])
#### ALL NET COUNT SITES & CODES are now represented in the STATION DATA

length(unique(nc))  # 27
length(unique(pss)) # 12
length(unique(nss)) # 27

###### SAVE THE UPDATED DF's
write.csv(points_in, "data/BIRDS_1_Counts_Points_Fixed_SitesNames.csv", row.names = F)
write.csv(nets_in, "data/BIRDS_1_Counts_Nets_Fixed_SitesNames.csv", row.names = F)

##################################################################################
##################################################################################

#####   Find Unique Station/Coordinate Locations per Site
zone18 = utm18
zone19 = utm19

## 1. Make Coordinates numeric & get rid of decimal points
zone18$`UTM Easting` = round(as.numeric(zone18$`UTM Easting`))
zone18$`UTM Northing` = round(as.numeric(zone18$`UTM Northing`))
zone19$`UTM Easting` = round(as.numeric(zone19$`UTM Easting`))
zone19$`UTM Northing` = round(as.numeric(zone19$`UTM Northing`))

## 1. Remove Stations with No Coordinates
zone18 = zone18[!is.na(zone18$`UTM Northing`),]
zone19 = zone19[!is.na(zone19$`UTM Northing`),]

# There should be no overlap between UTM 18 & 19 Sites
length(which(zone18$`Site Name` %in% zone19$`Site Name`)) # OK
length(which(zone18$`Site Code` %in% zone19$`Site Code`)) # OK
length(which(zone18$`UTM Easting` %in% zone19$`UTM Easting`)) # OK
length(which(zone18$`UTM Northing` %in% zone19$`UTM Northing`)) #OK


zone18$coords = paste0(zone18$`UTM Easting`, "_", zone18$`UTM Northing`)
zone18$stations = paste0(zone18$`Site Code`, "_", zone18$`Station Code`)
zone18$stat_co = paste0(zone18$stations, "_", zone18$coords)
length(unique(zone18$stations)) # 39 Unique station names
length(unique(zone18$coords))   # 38 Unique coordinates
length(unique(zone18$stat_co))  # 39 Unique station coordinate combinations
## More Station Names than Coordinates
## Check for Replicated Names
ed = zone18 %>%
  group_by(stations) %>%
  filter(n()>1) %>%
  ungroup()
## MARIO3 & QMario3 have the same coords - I will remove QMario3
zone18 = zone18[zone18$`Station Code` != "QMario3",]
##### ZONE 18 NOW HAS NO REPLAICATED COORDINATES OR STATION NAMES
write.csv(zone18, "data/BIRDS_1_Stations_UTM18_Fixed_Names_Coords.csv", row.names = F)

###

zone19$coords = paste0(zone19$`UTM Easting`, "_", zone19$`UTM Northing`)
zone19$stations = paste0(zone19$`Site Code`, "_", zone19$`Station Code`)
zone19$stat_co = paste0(zone19$stations, "_", zone19$coords)
length(unique(zone19$stations)) # 740 Unique station names
length(unique(zone19$coords))   # 730 Unique coordinates
length(unique(zone19$stat_co))  # 742 Unique station coordinate combinations
## More Station Names than Coordinates
## Check for Replicated Names
ed = zone19 %>%
  group_by(stations) %>%
  filter(n()>1) %>%
  ungroup()
## FPE1 & TT1.2 have two sets of coords each
# zone19$`Station Code` = make.unique(zone19$`Station Code`)
zone19$stations = paste0(zone19$`Site Code`, "_", zone19$`Station Code`)
zone19$stations = make.unique(zone19$stations)
zone19$stat_co = paste0(zone19$stations, "_", zone19$coords)

ed = zone19 %>%
  group_by(coords) %>%
  filter(n()>1) %>%
  ungroup()
unique(ed$coords)
## Only 12 sets of Unique Coords from 24 entries
rem = sort(ed$stations)
rem = rem[c(2:5, 10, 12, 14, 16, 18, 20, 23:24)]
zone19 = zone19[!zone19$stations %in% rem,]

length(unique(zone19$stations)) # 730 Unique station names
length(unique(zone19$coords))   # 730 Unique coordinates
length(unique(zone19$stat_co))  # 730 Unique station coordinate combinations
##### ZONE 19 NOW HAS NO REPLAICATED COORDINATES OR STATION NAMES
write.csv(zone19, "data/BIRDS_1_Stations_UTM19_Fixed_Names_Coords.csv", row.names = F)


#######################################################################
#######################################################################

zone19 = read.csv("data/BIRDS_1_Stations_UTM19_Fixed_Names_Coords.csv")
zone18 = read.csv("data/BIRDS_1_Stations_UTM18_Fixed_Names_Coords.csv")

## Stations & Coords of POINT Count
points_in = read.csv("data/BIRDS_1_Counts_Points_Fixed_SitesNames.csv")
names_point = c("site_name", "site_code", "date", "station", "utm_long", "utm_lat", "utm_zone", "observer", "time_start", "time_end", "duration_mins", "temp_start", "temp_end", "cloud_pct", "humidity", "name_Family",  "name_common", "name_Species", "name_code", "seen_heard_both", "number",  "distance_to_animal", "notes")
names(points_in) = names_point
points_in$utm_long = round(as.numeric(points_in$utm_long), digits = 0)
points_in$utm_lat = round(as.numeric(points_in$utm_lat), digits = 0)

nets_in = read.csv("data/BIRDS_1_Counts_Nets_Fixed_SitesNames.csv")
nn = names(nets_in)
nn = c("Family", "Species..English.", "Species..Scientific.", "Species.Code", "Date..dd.mm.yyyy.", "Year", "Month", "Capture.Time", "Site.Name", "Site.Code", "Station", "Net", "Station.GPS..Easting.", "Station.GPS..Northing.", "Station.GPS.Zone", "Station.Altitude", "Temp.Open..C.", "Temp.Mid..C.", "Temp.Close..C.", "Cloud.Cover.Open", "Cloud.Cover.Mid", "Cloud.Cover.Close", "Precip.Open", "Precip.Mid", "Precip.Close", "Wind.Open", "Wind.Mid", "Wind.Close", "First.Net.Time.Open", "Last.Net.Time.Open", "Time.open", "First.Net.Time.Close", "Last.Net.Time.Close", "Time.Closed", "Sample.Time..hrs.", "Number.of.Nets..12m.")
nets_in = select(nets_in, all_of(nn))
nets_in = select(nets_in, "Site.Name", "Site.Code", "Station", "Station.GPS..Easting.", "Station.GPS..Northing.", "Station.GPS.Zone", "Family", "Species..English.", "Species..Scientific.", "Species.Code", everything())
names(nets_in)[1:10] = c("site_name", "site_code", "station", "utm_long", "utm_lat", "utm_zone", "name_Family",  "name_common", "name_Species", "name_code")
nets_in$utm_long = round(as.numeric(nets_in$utm_long), digits = 0)
nets_in$utm_lat = round(as.numeric(nets_in$utm_lat), digits = 0)

# pts = points_in[,c(2,4,5,6)]
# nts = nets_in[,c(2:5)]

points_in$Site_Station = paste0(points_in$site_code, "_", points_in$station)
points_in$coords = paste0(points_in$utm_long, "_", points_in$utm_lat)
points_in$Site_coords = paste0(points_in$Site_Station, "_", points_in$coords)

nets_in$Site_Station = paste0(nets_in$site_code, "_", nets_in$station)
nets_in$coords = paste0(nets_in$utm_long, "_", nets_in$utm_lat)
nets_in$Site_coords = paste0(nets_in$Site_Station, "_", nets_in$coords)

#######
# Find Stations that have No Coords AND Name is not in SITES df
length(which((is.na(nets_in$utm_long) | is.na(nets_in$utm_lat)) & (!nets_in$Site_Station %in% zone18$stations | !nets_in$Site_Station %in% zone19$stations))) ## There are 1132 Net count lines without coords AND matching station name - These will be removed
rem = which((is.na(nets_in$utm_long) | is.na(nets_in$utm_lat)) & (!nets_in$Site_Station %in% zone18$stations | !nets_in$Site_Station %in% zone19$stations))
nets_in = nets_in[-rem,]

length(which((is.na(points_in$utm_long) | is.na(points_in$utm_lat)) & (!points_in$Site_Station %in% zone18$stations | !points_in$Site_Station %in% zone19$stations))) ## There are 44 Point count lines without coords AND matching station name - These will be removed
rem = which((is.na(points_in$utm_long) | is.na(points_in$utm_lat)) & (!points_in$Site_Station %in% zone18$stations | !points_in$Site_Station %in% zone19$stations))
points_in = points_in[-rem,]


# ALL COUNTS SHOULD NOW HAVE COORDINATES AND/OR RECIPROCAL SITE DATA WITH COORDINATES
write.csv(nets_in, "data/BIRDS_1.1_Counts_Nets_removed_NoCOORDS.csv", row.names = F)
write.csv(points_in, "data/BIRDS_1.1_Counts_Points_removed_NoCOORDS.csv", row.names = F)
#################

nets_in = read.csv("data/BIRDS_1.1_Counts_Nets_removed_NoCOORDS.csv")
points_in = read.csv("data/BIRDS_1.1_Counts_Points_removed_NoCOORDS.csv")
# pts = points_in[,c(1,2,4,5,6)]
# nts = nets_in[,c(1:5)]
# 
# nets_in$Site_Station = paste0(nets_in$site_code, "_", nets_in$station)
# nets_in$coords = paste0(nets_in$utm_long, "_", nets_in$utm_lat)
# 
# points_in$Site_Station = paste0(points_in$site_code, "_", points_in$station)
# points_in$coords = paste0(points_in$utm_long, "_", points_in$utm_lat)


### Out Of Interest ###
# Look at number of Coordinates 7 Stations from SITE data that are in COUNT data

length(unique(zone18$coords)) # 38
length(which(unique(zone18$coords) %in% points_in$coords | unique(zone18$coords) %in% nets_in$coords)) # 34 --- ie 4 Coords don't have counts
not = which(!unique(zone18$coords) %in% points_in$coords & !unique(zone18$coords) %in% nets_in$coords)
zone18$stations[not] # MARIO1,2,3 QMario1

length(unique(zone19$coords)) # 730
length(which(unique(zone19$coords) %in% points_in$coords | unique(zone19$coords) %in% nets_in$coords)) # 625 --- ie 105 Coords don't have counts
not = which(!unique(zone19$coords) %in% points_in$coords & !unique(zone19$coords) %in% nets_in$coords)
zone19$stations[not] # 


length(unique(zone18$stations)) # 38
length(which(unique(zone18$stations) %in% points_in$Site_Station | unique(zone18$stations) %in% nets_in$Site_Station)) # 33 --- ie 5 Stations don't have counts
not = which(!unique(zone18$stations) %in% points_in$Site_Station & !unique(zone18$stations) %in% nets_in$Site_Station)
zone18$stations[not] # MARIO1,2,3 QMario1, QP7

length(unique(zone19$stations)) # 730
length(which(unique(zone19$stations) %in% points_in$Site_Station | unique(zone19$stations) %in% nets_in$Site_Station)) # 408 --- ie 322 Stations don't have counts
not = which(!unique(zone19$stations) %in% points_in$Site_Station & !unique(zone19$stations) %in% nets_in$Site_Station)
zone19$stations[not] # 


##########
##########   Find where Count Coordinates match Site Coordinates, but Names do not   ###########
##########

wrong_name = which(!is.na(match(x = points_in$coords, table = zone18$coords)) & is.na(match(x = points_in$Site_Station, table = zone18$stations)))
### Find what the right names of those coordinates should be
right_name = match(x = points_in$coords[c(wrong_name)], table = zone18$coords)
### Replace the wrong names with the right names
points_in$site_name[wrong_name] = zone18$Site.Name[right_name]
points_in$site_code[wrong_name] = zone18$Site.Code[right_name]
points_in$station[wrong_name] = zone18$Station.Code[right_name]

### Repeat for zone19

### Find where Count Coordinates match Site Coordinates, but Names do not
wrong_name = which(!is.na(match(x = points_in$coords, table = zone19$coords)) & is.na(match(x = points_in$Site_Station, table = zone19$stations)))
### Find what the right names of those coordinates should be
right_name = match(x = points_in$coords[c(wrong_name)], table = zone19$coords)
### Replace the wrong names with the right names
points_in$site_name[wrong_name] = zone19$Site.Name[right_name]
points_in$site_code[wrong_name] = zone19$Site.Code[right_name]
points_in$station[wrong_name] = zone19$Station.Code[right_name]

points_in$Site_Station = paste0(points_in$site_code, "_", points_in$station)
## All Point Counts with Coordinates in Site data now have names matching with Site Data

### Repeat both for Net Counts ###

### Find where Count Coordinates match Site Coordinates, but Names do not
wrong_name = which(!is.na(match(x = nets_in$coords, table = zone18$coords)) & is.na(match(x = nets_in$Site_Station, table = zone18$stations)))
### Find what the right names of those coordinates should be
right_name = match(x = nets_in$coords[c(wrong_name)], table = zone18$coords)
### Replace the wrong names with the right names
nets_in$site_name[wrong_name] = zone18$Site.Name[right_name]
nets_in$site_code[wrong_name] = zone18$Site.Code[right_name]
nets_in$station[wrong_name] = zone18$Station.Code[right_name]

### Repeat for zone19

### Find where Count Coordinates match Site Coordinates, but Names do not
wrong_name = which(!is.na(match(x = nets_in$coords, table = zone19$coords)) & is.na(match(x = nets_in$Site_Station, table = zone19$stations)))
### Find what the right names of those coordinates should be
right_name = match(x = nets_in$coords[c(wrong_name)], table = zone19$coords)
### Replace the wrong names with the right names
nets_in$site_name[wrong_name] = zone19$Site.Name[right_name]
nets_in$site_code[wrong_name] = zone19$Site.Code[right_name]
nets_in$station[wrong_name] = zone19$Station.Code[right_name]

nets_in$Site_Station = paste0(nets_in$site_code, "_", nets_in$station)

wrong_name = which(!is.na(match(x = nets_in$coords, table = zone19$coords)) & is.na(match(x = nets_in$Site_Station, table = zone19$stations)))
### Find what the right names of those coordinates should be
right_name = match(x = nets_in$coords[c(wrong_name)], table = zone19$coords)
#### All counts where their coordinates match with Sites coordinates now have common names.


##########
##########   Find where Count Names match Site Names, but Coordinates do not   ###########
##########
wrong_coord = which(is.na(match(x = points_in$coords, table = zone18$coords)) & !is.na(match(x = points_in$Site_Station, table = zone18$stations)))
### ZERO Point Counts in ZONE 18
wrong_coord = which(is.na(match(x = points_in$coords, table = zone19$coords)) & !is.na(match(x = points_in$Site_Station, table = zone19$stations)))
### 13 Point Counts in ZONE 19
right_coord = match(x = points_in$Site_Station[c(wrong_coord)], table = zone19$stations)
### Replace the wrong names with the right names
points_in$utm_long[wrong_coord] = zone19$UTM.Easting[right_coord]
points_in$utm_lat[wrong_coord] = zone19$UTM.Northing[right_coord]

points_in$coords = paste0(points_in$utm_long, "_", points_in$utm_lat)

### Repeat both for Net Counts ###
wrong_coord = which(is.na(match(x = nets_in$coords, table = zone18$coords)) & !is.na(match(x = nets_in$Site_Station, table = zone18$stations)))
### ZERO Net Counts in ZONE 18
wrong_coord = which(is.na(match(x = nets_in$coords, table = zone19$coords)) & !is.na(match(x = nets_in$Site_Station, table = zone19$stations)))
### 101 Net Counts in ZONE 19
right_coord = match(x = nets_in$Site_Station[c(wrong_coord)], table = zone19$stations)
### Replace the wrong names with the right names
nets_in$utm_long[wrong_coord] = zone19$UTM.Easting[right_coord]
nets_in$utm_lat[wrong_coord] = zone19$UTM.Northing[right_coord]

nets_in$coords = paste0(nets_in$utm_long, "_", nets_in$utm_lat)
#### All counts where their NAMES match with Sites NAMES now have common COORDINATES

####
# ALL COUNTS THAT HAVE RECIPROCAL NAMES/COORDINATES IN SITE DATA ARE NOW THE SAME ACROSS ALL DATA

## NOW CHECK FOR ANOMOLIES IN THOSE COUNTS AT LOCATIONS THAT ARE NOT IN SITE DATA

# Point counts that don't have reciprocal coords or names in sites
pt_out = points_in[!points_in$coords %in% zone18$coords & !points_in$coords %in% zone19$coords & !points_in$Site_Station %in% zone18$stations & !points_in$Site_Station %in% zone19$stations,] # 323 Lines excluded
length(unique(pt_out$coords)) # 13 Coordinates
length(unique(pt_out$Site_Station)) # 13 Site_Station
length(unique(paste0(pt_out$Site_Station,"_",pt_out$coords))) # 13 Site_Coord Combos
# Point Counts are OK - ie they have the same number of Coords & Site names & they don't cross match to make a greater number

# Net counts that don't have reciprocal coords or names in sites
net_out = nets_in[!nets_in$coords %in% zone18$coords & !nets_in$coords %in% zone19$coords & !nets_in$Site_Station %in% zone18$stations & !nets_in$Site_Station %in% zone19$stations,] # 49 Lines Excluded
length(unique(net_out$coords)) # 6 Coordinates
length(unique(net_out$Site_Station)) # 9 Site_Station
length(unique(paste0(net_out$Site_Station,"_",net_out$coords))) # 9 Site_Coord Combos

unique(net_out$Site_Station) %in% zone18$stations | unique(net_out$Site_Station) %in% zone19$stations
unique(net_out$coords) %in% zone18$coords | unique(net_out$coords) %in% zone19$coords

reps = net_out %>%
  group_by(Site_Station) %>%
  filter(n()>1) %>%
  ungroup()

# LIM_Fiery-3, LIM_Fiery-2, LIM_Fiery-1 all have the same coords
# LIM_BA-2, LIM_BA-2 both have the same coords

### Change LIM_Fiery to LIM_Fiery-1-3
### Change LIM_BA to LIM_BA-2-3

nets_in$station[nets_in$Site_Station=="LIM_Fiery-3" | nets_in$Site_Station=="LIM_Fiery-2" | nets_in$Site_Station=="LIM_Fiery-1"] = "Fiery-1-3"
nets_in$station[nets_in$Site_Station=="LIM_BA-2" | nets_in$Site_Station=="LIM_BA-3"] = "BA-2-3"

nets_in$Site_Station = paste0(nets_in$site_code, "_", nets_in$station)

# Repeat Test
net_out = nets_in[!nets_in$coords %in% zone18$coords & !nets_in$coords %in% zone19$coords & !nets_in$Site_Station %in% zone18$stations & !nets_in$Site_Station %in% zone19$stations,] # 49 Lines Excluded
length(unique(net_out$coords)) # 6 Coordinates
length(unique(net_out$Site_Station)) # 6 Site_Station
length(unique(paste0(net_out$Site_Station,"_",net_out$coords))) # 6 Site_Coord Combos
# Nets Outside of Site Data are now OK

####
### Out of Interest ###
# Check if Names/Coords are reciprocated between Nets & Sites

unique(net_out$coords) %in% pt_out$coords
unique(net_out$Site_Station) %in% pt_out$Site_Station

unique(pt_out$coords) %in% net_out$coords
unique(pt_out$Site_Station) %in% net_out$Site_Station
## Nothing Reciprocal between Nets/Points that are not in Sites

## Check all Unique Coords & Stations in Nets & Points
length(unique(points_in$coords)) # 242
length(unique(points_in$Site_Station)) # 242
length(unique(paste0(points_in$coords, points_in$Site_Station))) # 242
# There are the same number of unique coordinates, station names & station-coord combos

length(unique(nets_in$coords)) # 640
length(unique(nets_in$Site_Station)) # 649
length(unique(paste0(nets_in$coords, nets_in$Site_Station))) # 660
### There are 640 Coordinates but 649 Sites & 660 Combos

## Get all Sites with repeated names
reps = nets_in %>%
  group_by(Site_Station) %>%
  filter(n()>1) %>%
  ungroup()

## Select Unique Site/Station/Coordinate sets 
un = unique(reps[,c(1:5, 37,38)])

### Reduce the Unique data to find replicated Site_stations
rep_un = un %>%
  group_by(Site_Station) %>%
  filter(n()>1) %>%
  ungroup()

length(unique(rep_un$Site_Station)) ## There are 11 Stations that have More than 1 set of coords
rep_stat = unique(rep_un$Site_Station) # Get those repeat station names

### For each of those stations with multple coords
for (i in 1:length(rep_stat)) {
df = nets_in[nets_in$Site_Station==rep_stat[i],]
keep = table(unlist(df$coords)) # Get a count of how many times each coordinate occurs
use = names(keep[which.max(keep)]) # Keep the one that occurs most
# Replace the Given Long & Lat with the most used Long & Lat
nets_in$utm_long[nets_in$Site_Station==rep_stat[i]] = as.numeric(substr(use, 1, 6))
nets_in$utm_lat[nets_in$Site_Station==rep_stat[i]] = as.numeric(substr(use, 8, nchar(use)))
# Recreate the coords value using the new long lat data
nets_in$coords = paste0(nets_in$utm_long,"_",nets_in$utm_lat)
}

length(unique(nets_in$coords)) # 640
length(unique(nets_in$Site_Station)) # 649
length(unique(paste0(nets_in$coords, nets_in$Site_Station))) # 649
### There are still 640 Coordinates and 649 Sites
### But there are now also only 649 Site-Coord combos.
## This means that 1 set of coordinates may have more than 1 Station Name attached to it, which is fine because each station only has 1 set of coords



##########
########## Merge the Net & Point count data to check for disparities between names & coords
names(nets_in)
nt = nets_in[,c(1:6)]
names(points_in)
pt = points_in[,c(1,2,4:7)]

df = rbind(nt,pt)
df$coords = paste0(df$utm_long,"_",df$utm_lat)
df$site_stn = paste0(df$site_code, "_", df$station)

# Check for site & coordinate disparity now the point & net data are combined
rep = unique(df)  # Get all unique sites (Names & Coords) = 700

# I don't mind coordinates having multiple names
# But Names must only have 1 set of coords
# Group the unique Locations by Name & pull out those that are repeated
rep2 = rep %>%
  group_by(site_stn) %>%
  filter(n()>1) %>%
  ungroup()
rep_n = unique(rep2$site_stn) # There are 6 Names with two coordinates
rep_n %in% zone19$stations     # All names occur in the SITE df

for (i in 1:length(rep_n)) {
  nets_in$utm_long[nets_in$Site_Station==rep_n[i]] = zone19$UTM.Easting[zone19$stations==rep_n[i]]
  nets_in$utm_lat[nets_in$Site_Station==rep_n[i]] = zone19$UTM.Northing[zone19$stations==rep_n[i]]
  points_in$utm_long[points_in$Site_Station==rep_n[i]] = zone19$UTM.Easting[zone19$stations==rep_n[i]]
  points_in$utm_lat[points_in$Site_Station==rep_n[i]] = zone19$UTM.Northing[zone19$stations==rep_n[i]]
}
nets_in$coords = paste0(nets_in$utm_long, "_", nets_in$utm_lat)
points_in$coords = paste0(points_in$utm_long, "_", points_in$utm_lat)

###########################
# Remove SPACES from Species Names
nets_in$name_Species = gsub(pattern = " ", replacement = "_", x = nets_in$name_Species)
points_in$name_Species = gsub(pattern = " ", replacement = "_", x = points_in$name_Species)

# ALL COUNTS SHOULD NOW HAVE COORDINATES AND/OR RECIPROCAL SITE DATA WITH COORDINATES
# And Site Names that do not have multiple coordinates
write.csv(nets_in, "data/BIRDS_1.2_Counts_Nets_NamesTOcoords_unique.csv", row.names = F)
write.csv(points_in, "data/BIRDS_1.2_Counts_Points_NamesTOcoords_unique.csv", row.names = F)

###############################################################################################

nets_in = read.csv("data/BIRDS_1.2_Counts_Nets_NamesTOcoords_unique.csv")
points_in = read.csv("data/BIRDS_1.2_Counts_Points_NamesTOcoords_unique.csv")

####
#### ADD THE MISSING STATIONS & COORDINATES to the SITES df

# Point Count Sites not in SITE data
nocopo = which(!points_in$coords %in% zone18$coords & !points_in$coords %in% zone19$coords)
miss_point = points_in[nocopo,]
# Net Count Sites not in SITE data
nocone = which(!nets_in$coords %in% zone18$coords & !nets_in$coords %in% zone19$coords)
miss_net = nets_in[nocone,]

## Check if any missing Point & Net Count Sites are reciprocated between both counts
miss_net$coords %in% miss_point$coords # No missing net coords are in missing site coords
miss_point$coords %in% miss_net$coords

miss_net$Site_Station %in% miss_point$Site_Station # No missing net Site_Station are in missing site Site_Station
miss_point$Site_Station %in% miss_net$Site_Station

### Add the missing Data to the SITE df
miss_net = unique(miss_net[,1:6])
colnames(miss_net) = colnames(zone19)[1:6]
# Missing Nets are all Zone19
zone19 = bind_rows(zone19, miss_net)
miss_point = unique(miss_point[,c(1,2,4:7)])
colnames(miss_point) = colnames(zone19)[1:6]
# Missing Points are all Zone18
zone18 = bind_rows(zone18, miss_point)

zone18$coords = paste0(zone18$UTM.Easting, "_", zone18$UTM.Northing)
zone19$coords = paste0(zone19$UTM.Easting, "_", zone19$UTM.Northing)

zone18$stations = paste0(zone18$Site.Code, "_", zone18$Station.Code)
zone19$stations = paste0(zone19$Site.Code, "_", zone19$Station.Code)

colnames(zone18)[1:6] = c("site_name", "site_code", "station", "utm_long", "utm_lat", "utm_zone")
colnames(zone19)[1:6] = c("site_name", "site_code", "station", "utm_long", "utm_lat", "utm_zone")

write.csv(zone18, "data/BIRDS_1.2_Stations_UTM18_ADDED_Missing_Count_Sites.csv", row.names = F)
write.csv(zone19, "data/BIRDS_1.2_Stations_UTM19_ADDED_Missing_Count_Sites.csv", row.names = F)

###############################################################################################
###############################################################################################

#### Add Lat Lon instead of UTM coordinates
# Create Spatial Coordinates from UTM data
xy18 = zone18
xy19 = zone19
library(proj4)
# Find the projection proj4 description of the coordinate reference system
# from https://spatialreference.org/ref/?search=UTM+zone+19S
# UTM zone 19 South = (espg 32719)
# Get the proj4 Description from https://spatialreference.org/ref/epsg/32719/proj4/
proj19 <- "+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj18 <- "+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Transform the data
lon_lat18 <- project(xy18[,4:5], proj18, inverse = TRUE)
lon_lat19 <- project(xy19[,4:5], proj19, inverse = TRUE)

# convert to a data frame
x_y18 = data.frame(site_station = xy18$stations, lon = lon_lat18$x, lat = lon_lat18$y)
x_y19 = data.frame(site_station = xy19$stations, lon = lon_lat19$x, lat = lon_lat19$y)
## Add new coords to Sites df
zone18 = tibble::add_column(zone18, x_y18[,2:3], .after = "utm_zone")
zone19 = tibble::add_column(zone19, x_y19[,2:3], .after = "utm_zone")
coords = rbind(zone18[,c(21,20,7,8)], zone19[,c(29,28,7,8)])


write.csv(coords, "data/BIRDS_USE_1_Coordinates_Peru_Stations_LatLon.csv", row.names = F)
write.csv(zone18, "data/BIRDS_1.3_Stations_UTM18_add_LAtLonCoords.csv", row.names = F)
write.csv(zone19, "data/BIRDS_1.3_Stations_UTM19_add_LAtLonCoords.csv", row.names = F)
