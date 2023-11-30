#####
### OCCUPANCY MODEL BIRDS USING spOCC
library(dplyr)
library(spOccupancy)
library(coda)        # for MCMC diagnostics
library(tidyr)
# library(tibble)
# library(lubridate)

##### For POINT Counts
counts = read.csv("data/BIRDS_2_Counts_All_Corrected_Species_Spelling.csv")
names(counts)
counts  = counts[,-c(21:23)]
counts$count = 1
counts = tibble::add_column(counts, site_station = paste0(counts$site_code, "_", counts$station), .after = "am.pm")
counts = tibble::add_column(counts, coords = paste0(counts$utm_long, "_", counts$utm_lat), .after = "site_code")
counts = tibble::add_column(counts, rep_ID = paste0(counts$site_station, "_", counts$date, "_", counts$am.pm), .after = "coords")

# Create a unique rows df to pivot wider
counts2 = dplyr::select(counts, -c("point_mins", "net_hrs","nets_qty","effort_point","effort_net_time","effort_net_nets","effort_net_total"))
# counts2 = counts

counts2 = counts2 %>%
  group_by(rep_ID) %>%
  mutate(effort = mean(effort)) %>%
  ungroup

dist_pnt = distinct(counts2) # 10940

dist_pnt$am.pm[is.na(dist_pnt$am.pm)] = "Unknown"
# dist_pnt$net_hrs[is.na(dist_pnt$net_hrs) & dist_pnt$method=="net"] = "Unknown"
# dist_pnt$nets_qty[is.na(dist_pnt$nets_qty) & dist_pnt$method=="net"] = "Unknown"

##################################################################
### Get Communities caught in Nets & Points & those common to both
net = dist_pnt[dist_pnt$method=="net" & dist_pnt$utm_zone=="19L",]
net_sp = unique(net$species)              ### 272 NET Species
length(unique(net$rep_ID))                ### 2893 Surveys
pnt = dist_pnt[dist_pnt$method=="point" & dist_pnt$utm_zone=="19L",]
pnt_sp = unique(pnt$species)              ### 216 POINT Species
length(unique(pnt$rep_ID))                ### 404 Surveys

both = which(net_sp %in% pnt_sp)          ### 130 species common to both types
both = net_sp[both]                       ### 142 species unique to NETS; 83 species unique to POINTS
both = sort(both)                         ### 355 Species total

## Repeat for A.M. and P.M. surveys
am = dist_pnt[dist_pnt$am.pm=="AM" & dist_pnt$utm_zone=="19L",]
am_sp = unique(am$species)                   ### 349 MORNING Species (of 354)
length(unique(am$rep_ID))                    ### 2968 Surveys
pm = dist_pnt[dist_pnt$am.pm=="PM" & dist_pnt$utm_zone=="19L",]
pm_sp = unique(pm$species)                   ### 102 AFTERNOON Species (of 354)
length(unique(pm$rep_ID))                    ### 215 Surveys

both = which(am_sp %in% pm_sp)            ### 97 species common to both types
both = am_sp[both]                        ### 252 species unique to MORNING; 5 species unique to AFTERNOON
both = sort(both)                         ### 354 Species total
##################################################################

## Turn long data into wide data...1 column per species...1 row per survey session
wide = tidyr::pivot_wider(dist_pnt, names_from = species, values_from = count, values_fill = 0) # 3566
names(wide)
length(unique(wide$rep_ID)) # 3566 Unique ID's from 3566 rows.

range(table(wide$site_station)) # Stations surveyed between 1 & 35 times
median(table(wide$site_station)) # 4
hist(table(wide$site_station), breaks = 20)

## Create Unique Surveys data frame
library(terra)
crs <- "+proj=utm +zone=19 +south"
p1 <- vect(as.data.frame(wide[ ,6:7]), geom=c("utm_long", "utm_lat"), crs=crs)
p1
p2 <- project(p1, "+proj=longlat")
p3 = as.data.frame(geom(p2))
# convert to a data frame
p3 = data.frame(lon = p3$x, lat = p3$y)
## Add new coords to Sites df
wide2 = tibble::add_column(wide, p3, .after = "utm_zone")
wide2 = wide2[,1:17]
write.csv(wide2, "data/Birds_3_Surveys.csv", row.names = F)


### Group by Transects and sum all sightings
### NB: These are not abundance counts, but the number of individual surveys in which the given species was recorded
tt = wide %>% group_by(site_station) %>%
  summarise(across(c(25:ncol(wide)-1),sum))

names(wide)

### Group all surveys by Unique Replicate
tt_day = wide %>% group_by(rep_ID, site_station, year, month, am.pm, effort, method) %>%
  summarise(across(c(16:ncol(wide)-7),sum))
length(unique(tt_day$site_station)) # There are 3566 Surveys spread over 694 Unique Locations



## Create DETECTION COVARIATES as a Per Site Per Survey matrix
# AM or PM
morn = tt_day[,c(1:2,5)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 11 so rename the Unique surveys as 1
morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-34 for each station
  group_by(site_station) %>%
  mutate(rep_ID = cumsum(rep_ID)) %>%
  ungroup
## Pivot the data so that Each Station has Values of AM PM or NA for each of the 35 possible replicates
morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = am.pm))
rownames(morn2) = morn2$site_station
morn2 = morn2[,-1]
table(unlist(morn2))
write.csv(morn2, "data/BIRDS_3_Det-Covs_AM-PM_per_site_Replicate.csv", row.names = T)

# Point or Net
morn = tt_day[,c(1:2,7)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 11 so rename the Unique surveys as 1
morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-35 for each station
  group_by(site_station) %>%
  mutate(rep_ID = cumsum(rep_ID)) %>%
  ungroup
## Pivot the data so that Each Station has Values of AM PM or NA for each of the 35 possible replicates
morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = method))
rownames(morn2) = morn2$site_station
morn2 = morn2[,-1]

table(unlist(morn2))

length(morn2[morn2=="net"]) # = 23778
length(morn2[morn2=="point"]) # = 21236
write.csv(morn2, "data/BIRDS_3_Det-Covs_POINT-NET_per_site_Replicate.csv", row.names = T)

# EFFORT
morn = tt_day[,c(1:2,6)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 35 so rename the Unique surveys as 1
# zero = which(morn$nets_qty=="Unknown" | morn$net_hrs=="Unknown")
# morn$effort[zero] = "Unknown"
# morn = morn[, -c(3:4)]
morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-35 for each station
  group_by(site_station) %>%
  mutate(rep_ID = cumsum(rep_ID)) %>%
  ungroup
## Pivot the data so that Each Station has Values of AM PM or NA for each of the 35 possible replicates
morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = effort))
rownames(morn2) = morn2$site_station
morn2 = morn2[,-1]
write.csv(morn2, "data/BIRDS_3_Det-Covs_EFFORT_per_site_Replicate.csv", row.names = T)


# YEAR
morn = tt_day[,c(1:2,3)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 11 so rename the Unique surveys as 1
morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-35 for each station
  group_by(site_station) %>%
  mutate(rep_ID = cumsum(rep_ID)) %>%
  ungroup
## Pivot the data so that Each Station has Values of AM PM or NA for each of the 35 possible replicates
morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = year))
rownames(morn2) = morn2$site_station
morn2 = morn2[,-1]
write.csv(morn2, "data/BIRDS_3_Det-Covs_YEAR_per_site_Replicate.csv", row.names = T)

# MONTH
morn = tt_day[,c(1:2,4)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 11 so rename the Unique surveys as 1
morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-35 for each station
  group_by(site_station) %>%
  mutate(rep_ID = cumsum(rep_ID)) %>%
  ungroup
## Pivot the data so that Each Station has Values of AM PM or NA for each of the 35 possible replicates
morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = month))
rownames(morn2) = morn2$site_station
morn2 = morn2[,-1]
write.csv(morn2, "data/BIRDS_3_Det-Covs_MONTH_per_site_Replicate.csv", row.names = T)


##############
##### Create Matrix per Species for Analysis
df = tt_day
df = ungroup(df)

length(unique(paste0(counts$date,counts$am.pm))) # There are 1234 individual survey date/am-pm combinations
range(table(df$site_station))                    # There are a maximum of 35 surveys conducted at any Station


counts = as.character(unique(df$site_station)) # Variable of Transect Names
df1 = df[FALSE,] # Create a blank DF based on the columns of the original

for (i in 1:length(counts)) {   # For each transect
  #i=1
  sb = c(1:nrow(df[df$site_station==counts[i],]))  # Create a vector to name the survey replicates (1-number of replicates)
  tr = df[df$site_station==counts[i],]                 # Subset data to current Transect
  tr$rep_ID = sb                                # Replace Transect Names with the survey number just created
  
  # NB: Currently there are many counts with fewer than 35 surveys
  # table(df$Site_Station)
  # So the shortfall must be made up of NA values so that a Site*Survey per Species Array can be created (Dimensions need to be common across all) 
  if (length(sb)<35) { 
    m = as.data.frame(matrix(data = NA, nrow = 35-length(sb), ncol = ncol(tr)))
    m[,2] = counts[i]                ## Station
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
for (i in 10:ncol(df1)) {
  # i = 8
  df2 = df1[,c(1:2,i)]
  # bird = colnames(df[i])
  colnames(df2)[3] = "count"
  ## Pivot the data wider so that each survey is a column (there should be 240 surveys per site), this reduces the rows to 1 per site & whether a bird was seen or not at that site & survey date. If there was not a survey at a given site & date, NA is inserted.  
  # df2 = pivot_wider(data = df1, names_from = Date, values_from = count)
  df3 = df2 %>% group_by(site_station) %>%
    pivot_wider( names_from = rep_ID, values_from = count)
  df3 = df3 %>% ungroup() %>% filter(!if_all(2:ncol(df3), is.na)) # Get rid of the all NA row
  transects = df3[,1]
  df3 = as.matrix(df3[,-1])

# assign(x = bird, value = df2)
site.surv.spec = append(site.surv.spec, list(df3))
}
names(site.surv.spec) = colnames(df[10:ncol(df)])

ssPnt <- array (
  data = do.call(rbind, lapply(site.surv.spec, as.vector)), 
  dim = c(length(site.surv.spec), dim(site.surv.spec[[1]])),
  dimnames = list(colnames(df[10:ncol(df)]), transects$site_station, colnames(df3))
)
str(ssPnt)
save(ssPnt, file = "data/BIRDS_3_COUNTS_ALL_Site-X-REPLICATE-per-Species_ARRAY.Rdata")
