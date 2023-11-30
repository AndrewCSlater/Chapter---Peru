library(dplyr)


utm19 = read.csv("data/BIRDS_1_Stations_UTM19_Fixed_Names_Coords.csv")
u19 = utm19[,c(1:6,26:28)]
utm18 = read.csv("data/BIRDS_1_Stations_UTM18_Fixed_Names_Coords.csv")
u18 = utm18[,c(1:6,18:20)]
n = c("site_name", "site_code", "station", "utm_long", "utm_lat", "utm_zone", "coords", "stations", "stn_coords")
names(u18) = n
names(u19) = n
site = rbind(u19,u18)
length(unique(site$stations)) # 768 Station Names
length(unique(site$coords)) # 768 Coordinates
length(unique(site$stn_coords)) # 768 Station_Coord combos


counts = read.csv("data/BIRDS_1.5_Combined_NetPoint_Counts.csv")
counts$stations = paste0(counts$site_code, "_", counts$station)
counts$coords = paste0(counts$utm_long, "_", counts$utm_lat)
counts$stn_coords = paste0(counts$stations, "_", counts$coords)
length(unique(counts$stations)) # 773 Station Names
length(unique(counts$coords)) # 687 Coordinates
length(unique(counts$stn_coords)) # 803 Station_Coord combos!!

length(which(!unique(counts$stations) %in% site$stations))     # 138 of 773 station names are not in Site data
length(which(!unique(counts$coords) %in% site$coords))         # 28 of 687 coordinates are not in Site data
length(which(!unique(counts$stn_coords) %in% site$stn_coords)) # 193 of 803 station coord combos are not in Site data


### 1 - Where coordinates Match, I will Copy Site Station Name to Count Station Name
swap = match(counts$coords, site$coords)
for (i in 1:length(swap)) {
  # i=1
  
  if (!is.na(swap[i])) { counts$station[i] = site$station[swap][i];
  counts$site_name[i] = site$site_name[swap][i];
  counts$site_code[i] = site$site_code[swap][i];
  counts$site_name[i] = site$site_name[swap][i]}
}

counts$stations = paste0(counts$site_code, "_", counts$station)
counts$stn_coords = paste0(counts$stations, "_", counts$coords)

### 2 - Possible wrong length of UTM coords
swap = which(nchar(counts$utm_long)!=6)
swap = which(nchar(counts$utm_lat)!=7)    ###  Coord is too short
unique(counts$stations[swap]) %in% site$stations 
err = unique(counts$stations[swap])       ### Names of stations with short coords

for (i in 1:nrow(counts)) {       ## For each row of count data; If the station is one with short coords    
  #i=1                                      
  if (counts$stations[i] %in% err) {
    a = match(counts$stations[i], site$stations);  ### Match the coord with the site data
    counts$utm_lat[i] = site$utm_lat[a];           ### and Copy it across to count data
    }
}
counts$coords = paste0(counts$utm_long, "_", counts$utm_lat)
counts$stn_coords = paste0(counts$stations, "_", counts$coords)

length(which(!unique(counts$stations) %in% site$stations))     # 21 of 692 station names are not in Site data
length(which(!unique(counts$coords) %in% site$coords))         # 25 of 685 coordinates are not in Site data
length(which(!unique(counts$stn_coords) %in% site$stn_coords)) # 36 of 694 station coord combos are not in Site data

length(unique(counts$stations)) # 692 Station Names
length(unique(counts$coords)) # 685 Coordinates
length(unique(counts$stn_coords)) # 694 Station_Coord combos
### This implies that Some Station Names are associated with more than 1 coordinate.
### Coords can have multiple names but names must have only 1 coord

stn = counts[,21:23]
stn = distinct(stn)
stn = stn %>% group_by(stations) %>%
  filter(n()>1) %>%
  ungroup
multi = unique(stn$stations) # Get the unique stations that have multiple coords

# table(counts$coords[counts$stations == multi[[1]]])
# site$coords[site$stations == multi[[1]]]

for (j in 1:length(multi)) { # For each Station with mutliple coords

for (i in 1:nrow(counts)) {  # Run through each row of counts
 # i=1
if (counts$stations[i] == multi[[j]] & counts$coords[i] != site$coords[site$stations == multi[[j]]]) { # Where the count station = the multi coord station & the coords differ from the site coords for that station
  counts$station[i] = paste0(counts$station[i], "_Alt_Coord")}    # Add an Alt_coord text to the site name
} # End Counts
} # End Multi coord stations

counts$stations = paste0(counts$site_code, "_", counts$station)
counts$stn_coords = paste0(counts$stations, "_", counts$coords)

length(unique(counts$stations)) # 694 Station Names
length(unique(counts$coords)) # 685 Coordinates
length(unique(counts$stn_coords)) # 694 Station_Coord combos
# There are now the same number of Station-Coord combos as Station Names - IE Each station now has only 1 set of coords
#### 

a = which(!counts$stations %in% site$stations)
an = unique(counts$stations[a])  # These 23 count station names do not exist in site data 
### As matching Coordinates have already had Station Names unified, those Stations that do not have names in Site Data must also not have coords in site data & are therefore new & unique

b = which(!counts$coords %in% site$coords)
bn = unique(counts$stations[b])  # These are the 34 count station names where the coordinates do not exist in site data
# In Some, Neither the Names Or Coordinates exist & are entirely new data (Those in the variable 'an' above)
# In others the Station Name exists but has Different coordinates - This needs to be rectified by changing Count Station Name & Adding to Site Data

### Those observations with Non-Matching Coordinates but Matching Names have been named out of character or have the wrong coordinates...As I do not know this I will keep the sightings but make the Station names unique & they will then also be unique Sites

## Single out the station names that have wrong/alternative coordinates
c = which(!bn %in% an)
c = bn[c] ## These are the station names which occur in both Site & Count but where the coordinates in Count do not exist in site

######## Examine the Stations & their Site & Count Coordinates
for (i in 1:length(c)) {
# print(table(counts$coords[counts$stations == c[[i]]])) # Each Station with wrong coords ONLY has 1 set of coords
print(site$coords[site$stations == c[[i]]])
print(unique(counts$coords[counts$stations == c[[i]]]))
}

## Visual Inspection shows that "SFO_TT4", "TNS_TNS8", "EI_TT14" have a Single Digit out & on Google Earth shows that the Site Data is correct
## LIM Hermit Counts have 3 coordinates for 9 Stations and again I will take the Site coordinates to be correct

for (j in 1:length(c)) { # For each Station with mutliple coords
  # j=1
  for (i in 1:nrow(counts)) {  # Run through each row of counts
    # i=1
    if (counts$stations[i] == c[[j]]) { # Where the count station = the wrong coord station
      counts$utm_long[i] = site$utm_long[site$stations==c[[j]]]; # the count long & lat gets replaced by the site long & lat
      counts$utm_lat[i] = site$utm_lat[site$stations==c[[j]]]}    
  } # End Counts
} # End wrong coord stations

counts$coords = paste0(counts$utm_long, "_", counts$utm_lat)
counts$stn_coords = paste0(counts$stations, "_", counts$coords)


####################
### SOME SITE STATIONS HAVE DUAL NAMES THAT DO NOT SHOW UP IN COUNTS
## eg LIM_Fiery-1 or Nun 1 in Site is not in Count - BUT - LIM_Fiery-1 IS in Count
## Their Count Coords are replicated amongst sites, so I will swap them to the Sites coords

for (i in 1:length(an)) { # For each station that has no name in sites
  
ii = grep(pattern = an[i], x = site$stations) # Check to see if the count station name is included as part of a site station
if (length(ii>0)) {                           # If it is
 counts$utm_long[counts$stations==an[i]] = site$utm_long[ii];
counts$utm_lat[counts$stations==an[i]] = site$utm_lat[ii]
  }
}

counts$coords = paste0(counts$utm_long, "_", counts$utm_lat)
counts$stn_coords = paste0(counts$stations, "_", counts$coords)

a = which(!counts$stations %in% site$stations)
an = unique(counts$stations[a])  # These 23 count station names do not exist in site data 
### As matching Coordinates have already had Station Names unified, those Stations that do not have names in Site Data must also not have coords in site data & are therefore new & unique

b = which(!counts$coords %in% site$coords)
bn = unique(counts$stations[b])  # These are the 15 count station names where the coordinates do not exist in site data

c = an[!an %in% bn] ### Count Stations whose coordinates exist but names don't in Sites data - This is OK

write.csv(counts, "data/BIRDS_1.75_Combined_Counts_Unified_Station_Coords.csv", row.names = F)





