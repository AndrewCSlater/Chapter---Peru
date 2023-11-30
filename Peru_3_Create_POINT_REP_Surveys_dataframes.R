library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)

##### For POINT Counts
points = read.csv("data/BIRDS_1.3_Counts_Points_Corrected_Species_Spelling.csv")
names(points)
points = select(points,  "Site_Station", "date", "observer", "time_start", "duration_mins", "temp_start", "temp_end", "cloud_pct", "humidity", "name_Species")

#### Correct data structure
points$date = as.Date(points$date)
points$time_start = as.POSIXct(points$time_start)                      
points = tibble::add_column(points, AM.PM = format(points$time_start, "%p"), .after = "time_start")
points = tibble::add_column(points, rep_ID = paste0(points$Site_Station, "_", points$date, "_", points$AM.PM), .before = "Site_Station")
points$temp_start = as.numeric(points$temp_start)
points[points== -9] = NA

### To use AM.PM counts as a Detection Covariate, all surveys need a value as NA is not accepted
### So Turn NA into UK (unknown)
points$AM.PM[is.na(points$AM.PM)] = "uk"

points$count = 1

# Create a unique rows df to pivot wider
dist_pnt = distinct(points)

## Turn long data into wide data...1 column per species...1 row per survey session
wide = pivot_wider(dist_pnt, names_from = name_Species, values_from = count, values_fill = 0)
names(wide)
length(unique(wide$rep_ID)) # Only 512 Unique ID's from 520 rows. Some Stations were counted more than once in the AM or PM on a given day
wide %>%
  group_by(rep_ID) %>%
  filter(n()>1) %>%
  ungroup()                 # Shows 15 Surveys in total from 7 Site Days

# Make the replicate name Unique
wide$rep_ID = make.unique(wide$rep_ID)
length(unique(wide$rep_ID)) # Now all 520 surveys have unique replicate numbers

range(table(wide$Site_Station)) # Stations surveyed between 1 & 11 times
median(table(wide$Site_Station)) # 2
hist(table(wide$Site_Station), breaks = 20)

### Group by Transects and sum all sightings
### NB: These are not abundance counts, but the number of individual surveys in which the given species was recorded
tt = wide %>% group_by(Site_Station) %>%
  summarise(across(c(12:ncol(wide)-1),sum))


### Group all surveys by Unique Replicate
tt_day = wide %>% group_by(rep_ID, Site_Station, date, AM.PM) %>%
  summarise(across(c(12:ncol(wide)-4),sum))
length(unique(tt_day$Site_Station)) # There are 520 Surveys spread over 242 Unique Locations
tt_day = tt_day[order(tt_day$date),]

## Create DETECTION COVARIATES as a Per Site Per Survey matrix
# AM or PM  
morn = tt_day[,c(1:2,4)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 11 so rename the Unique surveys as 1
morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-11 for each station
  group_by(Site_Station) %>%
  mutate(rep_ID = cumsum(rep_ID)) %>%
  ungroup

## Pivot the data so that Each Station has Values of AM PM or NA for each of the 11 possible replicates
morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = AM.PM))
rownames(morn2) = morn2$Site_Station
morn2 = morn2[,-1]

# morn2 <- morn2 %>%
#   group_by(Site_Station) %>%
#   fill(everything(), .direction = "downup") %>%
#   slice(1)

write.csv(morn2, "data/Det-Covs_AM-PM_per_site_Replicate.csv", row.names = T)

##############
##### Create Matrix per Species for Analysis
df = tt_day
df = ungroup(df)

length(unique(paste0(points$date,points$AM.PM))) # There are 111 individual survey date/am-pm combinations
range(table(df$Site_Station))                    # There are a maximum of 11 surveys conducted at any Station


points = as.character(unique(df$Site_Station)) # Variable of Transect Names
df1 = df[FALSE,] # Create a blank DF based on the columns of the original

for (i in 1:length(points)) {   # For each transect
  #i=1
  sb = c(1:nrow(df[df$Site_Station==points[i],]))  # Create a vector to name the survey replicates (1-number of replicates)
  tr = df[df$Site_Station==points[i],]                 # Subset data to current Transect
  tr$rep_ID = sb                                # Replace Transect Names with the survey number just created
  
  # NB: Currently there are many points with fewer than 11 surveys
  # table(df$Site_Station)
  # So the shortfall must be made up of NA values so that a Site*Survey per Species Array can be created (Dimensions need to be common across all) 
  if (length(sb)<11) { 
    m = as.data.frame(matrix(data = NA, nrow = 11-length(sb), ncol = ncol(tr)))
    m[,2] = points[i]                ## Station
    m[,1] = c(length(sb)+1:nrow(m))  ## Survey ID
    colnames(m) = colnames(tr)
    tr = rbind(tr,m)
  }
  ## Add the newly named and lengthened set of Transect data to the new DF 
  df1 = rbind(df1,tr)
}

# Create empty list or array to hold all species matrices
site.surv.spec = list()

## For each individual Species
for (i in 5:ncol(df1)) {
  # i = 5
  df2 = df1[,c(1:2,i)]
  # bird = colnames(df[i])
  colnames(df2)[3] = "count"
  ## Pivot the data wider so that each survey is a column (there should be 240 surveys per site), this reduces the rows to 1 per site & whether a bird was seen or not at that site & survey date. If there was not a survey at a given site & date, NA is inserted.  
  # df2 = pivot_wider(data = df1, names_from = Date, values_from = count)
  df3 = df2 %>% group_by(Site_Station) %>%
    pivot_wider( names_from = rep_ID, values_from = count)
  df3 = df3 %>% ungroup() %>% filter(!if_all(2:ncol(df3), is.na)) # Get rid of the all NA row
  transects = df3[,1]
  df3 = as.matrix(df3[,-1])

# assign(x = bird, value = df2)
site.surv.spec = append(site.surv.spec, list(df3))
}
names(site.surv.spec) = colnames(df[5:ncol(df)])

ssPnt <- array (
  data = do.call(rbind, lapply(site.surv.spec, as.vector)), 
  dim = c(length(site.surv.spec), dim(site.surv.spec[[1]])),
  dimnames = list(colnames(df[5:ncol(df)]), transects$Site_Station, colnames(df3))
)
str(ssPnt)
save(ssPnt, file = "data/BIRDS_3_COUNTS_POINTS_Site-X-REPLICATE-per-Species_ARRAY.Rdata")


#######################################################################################
#### REPEAT FOR NET COUNTS
#######################################################################################

nets = read.csv("data/BIRDS_1.3_Counts_Nets_Corrected_Species_Spelling.csv")
names(nets)
nets = select(nets,  "Site_Station", "Date..dd.mm.yyyy.", "Capture.Time", "name_Species")
names(nets) = c("Site_Station", "date", "time", "name_Species")

#### Correct data structure
nets$date = as.Date(nets$date)
nets$time =  as.POSIXct(
  as.numeric((nets$time)) * 24 * 60 * 60,
  origin = "1904-01-01 00:00:00", tz="GMT"
)
nets = tibble::add_column(nets, AM.PM = format(nets$time, "%p"), .after = "time")
nets$AM.PM[is.na(nets$AM.PM)] = "uk"
nets = tibble::add_column(nets, rep_ID = paste0(nets$Site_Station, "_", nets$date, "_", nets$AM.PM), .before = "Site_Station")
nets = select(nets,-time)
nets$count = 1

# Create a unique rows df to pivot wider
dist_net = distinct(nets)

## Turn long data into wide data...1 column per species...1 row per survey session
wide = pivot_wider(dist_net, names_from = name_Species, values_from = count, values_fill = 0)
names(wide)
length(unique(wide$rep_ID)) # 3036 Unique ID's from 3036 rows.

range(table(wide$Site_Station)) # Stations surveyed between 1 & 27 times
median(table(wide$Site_Station)) # 4
hist(table(wide$Site_Station), breaks = 20)

### Group by Transects and sum all sightings
### NB: These are not abundance counts, but the number of individual surveys in which the given species was recorded
tt = wide %>% group_by(Site_Station) %>%
  summarise(across(c(5:ncol(wide)-1),sum))

### Group all surveys by Unique Replicate
tt_day = wide %>% group_by(rep_ID, Site_Station, date, AM.PM) %>%
  summarise(across(c(5:ncol(wide)-4),sum))
length(unique(tt_day$Site_Station)) # There are 3036 Surveys spread over 650 Unique Locations

##############
##### Create Matrix per Species for Analysis
df = tt_day
df = ungroup(df)

length(unique(paste0(nets$date,nets$AM.PM))) # There are 1146 individual survey date/am-pm combinations
range(table(df$Site_Station))                # There are a maximum of 27 surveys conducted at any Station


nts = as.character(unique(df$Site_Station)) # Variable of Transect Names
df1 = df[FALSE,] # Create a blank DF based on the columns of the original

for (i in 1:length(nts)) {   # For each transect
  #i=1
  sb = c(1:nrow(df[df$Site_Station==nts[i],]))  # Create a vector to name the survey replicates (1-number of replicates)
  tr = df[df$Site_Station==nts[i],]             # Subset data to current Transect
  tr$rep_ID = sb                                # Replace Transect Names with the survey number just created
  
  # NB: Currently there are many points with fewer than 27 surveys
  # table(df$Site_Station)
  # So the shortfall must be made up of NA values so that a Site*Survey per Species Array can be created (Dimensions need to be common across all) 
  if (length(sb)<27) { 
    m = as.data.frame(matrix(data = NA, nrow = 27-length(sb), ncol = ncol(tr)))
    m[,2] = nts[i]                ## Station
    m[,1] = c(length(sb)+1:nrow(m))  ## Survey ID
    colnames(m) = colnames(tr)
    tr = rbind(tr,m)
  }
  ## Add the newly named and lengthened set of Transect data to the new DF 
  df1 = rbind(df1,tr)
}

# Create empty list or array to hold all species matrices
site.surv.spec = list()

## For each individual Species
for (i in 5:ncol(df1)) {
  # i = 5
  df2 = df1[,c(1:2,i)]
  # bird = colnames(df[i])
  colnames(df2)[3] = "count"
  ## Pivot the data wider so that each survey is a column (there should be 27 surveys per site), this reduces the rows to 1 per site & whether a bird was seen or not at that site & survey date. If there was not a survey at a given site & date, NA is inserted.  
  # df2 = pivot_wider(data = df1, names_from = Date, values_from = count)
  df3 = df2 %>% group_by(Site_Station) %>%
    pivot_wider( names_from = rep_ID, values_from = count)
  df3 = df3 %>% ungroup() %>% filter(!if_all(2:ncol(df3), is.na)) # Get rid of the all NA row
  transects = df3[,1]
  df3 = as.matrix(df3[,-1])
  
  # assign(x = bird, value = df2)
  site.surv.spec = append(site.surv.spec, list(df3))
}
names(site.surv.spec) = colnames(df[5:ncol(df)])

ssNet <- array (
  data = do.call(rbind, lapply(site.surv.spec, as.vector)), 
  dim = c(length(site.surv.spec), dim(site.surv.spec[[1]])),
  dimnames = list(colnames(df[5:ncol(df)]), transects$Site_Station, colnames(df3))
)
str(ssNet)
save(ssNet, file = "data/BIRDS_3_COUNTS_NETS_Site-X-REPLICATE-per-Species_ARRAY.Rdata")









