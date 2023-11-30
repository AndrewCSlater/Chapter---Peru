
### OCCUPANCY MODEL BIRDS USING spOCC
library(dplyr)
library(spOccupancy)
library(coda)        # for MCMC diagnostics
library(tidyr)
# library(tibble)
# library(lubridate)

plogis(.79)

##### For POINT Counts
counts = read.csv("data/BIRDS_2_Counts_All_Corrected_Species_Spelling.csv")
counts  = counts[,-c(11:17,21:23)]
counts$count = 1

counts = tibble::add_column(counts, site_station = paste0(counts$site_code, "_", counts$station, "_", counts$year), .after = "am.pm")
counts = tibble::add_column(counts, coords = paste0(counts$utm_long, "_", counts$utm_lat), .after = "site_code")
counts = tibble::add_column(counts, rep_ID = paste0(counts$site_station, "_", counts$date, "_", counts$am.pm), .after = "coords")

counts2 = counts

counts2 = counts2 %>%
  group_by(rep_ID) %>%
  mutate(effort = mean(effort)) %>%
  ungroup

dist_pnt = distinct(counts2) # 10940

dist_pnt$am.pm[is.na(dist_pnt$am.pm)] = "Unknown"

# ##################################################################
# ### Get Communities caught in Nets & Points & those common to both
# net = dist_pnt[dist_pnt$method=="net" & dist_pnt$utm_zone=="19L",]
# net_sp = unique(net$species)              ### 272 NET Species
# length(unique(net$rep_ID))                ### 2893 Surveys
# pnt = dist_pnt[dist_pnt$method=="point" & dist_pnt$utm_zone=="19L",]
# pnt_sp = unique(pnt$species)              ### 216 POINT Species
# length(unique(pnt$rep_ID))                ### 404 Surveys
# 
# both = which(net_sp %in% pnt_sp)          ### 130 species common to both types
# both = net_sp[both]                       ### 142 species unique to NETS; 83 species unique to POINTS
# both = sort(both)                         ### 355 Species total
# 
# ## Repeat for A.M. and P.M. surveys
# am = dist_pnt[dist_pnt$am.pm=="AM" & dist_pnt$utm_zone=="19L",]
# am_sp = unique(am$species)                   ### 349 MORNING Species (of 354)
# length(unique(am$rep_ID))                    ### 2968 Surveys
# pm = dist_pnt[dist_pnt$am.pm=="PM" & dist_pnt$utm_zone=="19L",]
# pm_sp = unique(pm$species)                   ### 102 AFTERNOON Species (of 354)
# length(unique(pm$rep_ID))                    ### 215 Surveys
# 
# both = which(am_sp %in% pm_sp)            ### 97 species common to both types
# both = am_sp[both]                        ### 252 species unique to MORNING; 5 species unique to AFTERNOON
# both = sort(both)                         ### 354 Species total
##################################################################

## Turn long data into wide data...1 column per species...1 row per survey session
wide = tidyr::pivot_wider(dist_pnt, names_from = species, values_from = count, values_fill = 0) # 3566
# names(wide)
length(unique(wide$rep_ID)) # 3566 Unique ID's from 3566 rows.

range(table(wide$site_station)) # Stations surveyed between 1 & 21 times
median(table(wide$site_station)) # 2
hist(table(wide$site_station), breaks = 20)
reps = as.data.frame(table(wide$site_station))
summary(reps[,2]) ## 1.000   1.000   2.000   2.921   4.000  21.000
quantile(reps[,2], probs = c(0.8, 0.9, 0.95, 0.975, 0.99)) ## 4     6     8    10    12


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
write.csv(wide2, "data/Birds_3_Surveys_SITES_ANNUAL.csv", row.names = F)

### Group by Transects and sum all sightings
### Group all surveys by Unique Replicate
tt_day = wide %>% group_by(rep_ID, site_station, year, month, am.pm, effort, method) %>%
  summarise(across(c(16:ncol(wide)-7),sum))
length(unique(tt_day$site_station)) # There are 3566 Surveys spread over 1221 Unique Locations (Each Station is Considered a Separate Place each Year)

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
write.csv(morn2, "data/BIRDS_3_Det-Covs_AM-PM_per_siteANNUALLY_Replicate.csv", row.names = T)

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
write.csv(morn2, "data/BIRDS_3_Det-Covs_POINT-NET_per_siteANNUALLY_Replicate.csv", row.names = T)

# # EFFORT
# morn = tt_day[,c(1:2,6)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
# morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 35 so rename the Unique surveys as 1
# # zero = which(morn$nets_qty=="Unknown" | morn$net_hrs=="Unknown")
# # morn$effort[zero] = "Unknown"
# # morn = morn[, -c(3:4)]
# morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-35 for each station
#   group_by(site_station) %>%
#   mutate(rep_ID = cumsum(rep_ID)) %>%
#   ungroup
# ## Pivot the data so that Each Station has Values of AM PM or NA for each of the 35 possible replicates
# morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = effort))
# rownames(morn2) = morn2$site_station
# morn2 = morn2[,-1]
# write.csv(morn2, "data/BIRDS_3_Det-Covs_EFFORT_per_site_Replicate.csv", row.names = T)
# 
# 
# # YEAR
# morn = tt_day[,c(1:2,3)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
# morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 11 so rename the Unique surveys as 1
# morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-35 for each station
#   group_by(site_station) %>%
#   mutate(rep_ID = cumsum(rep_ID)) %>%
#   ungroup
# ## Pivot the data so that Each Station has Values of AM PM or NA for each of the 35 possible replicates
# morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = year))
# rownames(morn2) = morn2$site_station
# morn2 = morn2[,-1]
# write.csv(morn2, "data/BIRDS_3_Det-Covs_YEAR_per_site_Replicate.csv", row.names = T)
# 
# # MONTH
# morn = tt_day[,c(1:2,4)]  ## Select the Unique Survey ID, the Station Name & Whether the survey was AM or PM
# morn$rep_ID = 1           ## The Site Count Replicates are named 1 - 11 so rename the Unique surveys as 1
# morn = morn %>%           ## By getting the cumulative Count of the Unique ID we then get replicate numbers between 1-35 for each station
#   group_by(site_station) %>%
#   mutate(rep_ID = cumsum(rep_ID)) %>%
#   ungroup
# ## Pivot the data so that Each Station has Values of AM PM or NA for each of the 35 possible replicates
# morn2 = as.data.frame(pivot_wider(morn, names_from = c(rep_ID), values_from = month))
# rownames(morn2) = morn2$site_station
# morn2 = morn2[,-1]
# write.csv(morn2, "data/BIRDS_3_Det-Covs_MONTH_per_site_Replicate.csv", row.names = T)


##############
##### Create Matrix per Species for Analysis
df = tt_day
df = ungroup(df)

length(unique(paste0(counts$date,counts$am.pm))) # There are 1234 individual survey date/am-pm combinations
range(table(df$site_station))                    # There are a maximum of 21 surveys conducted at any Station

counts = as.character(unique(df$site_station)) # Variable of Transect Names
df1 = df[FALSE,] # Create a blank DF based on the columns of the original

for (i in 1:length(counts)) {   # For each transect
  #i=1
  sb = c(1:nrow(df[df$site_station==counts[i],]))  # Create a vector to name the survey replicates (1-number of replicates)
  tr = df[df$site_station==counts[i],]                 # Subset data to current Transect
  tr$rep_ID = sb                                # Replace Transect Names with the survey number just created
  
  # NB: Currently there are many counts with fewer than 21 surveys
  # table(df$site_station)
  # So the shortfall must be made up of NA values so that a Site*Survey per Species Array can be created (Dimensions need to be common across all) 
  if (length(sb)<21) { 
    m = as.data.frame(matrix(data = NA, nrow = 21-length(sb), ncol = ncol(tr)))
    m[,2] = counts[i]                ## Station
    m[,1] = c(length(sb)+1:nrow(m))  ## Survey ID
    colnames(m) = colnames(tr)
    tr = rbind(tr,m)
  }
  ## Add the newly named and lengthened set of Transect data to the new DF 
  df1 = rbind(df1,tr)
}

###############################################################
##### Code Added 28/11/2023 - Find number of surveys & species found by stations with different habitats & methods
names(df1)
df2 = select(df1, 2,5,7,8:ncol(df1))
df2 = df2[!is.na(df2$am.pm) | !is.na(df2$method),]
rem = grep(pattern = "QUE", x = df2$site_station)
df2 = df2[-rem,]
df3 = df2 %>% group_by(site_station, am.pm, method) %>% summarise_all(sum) %>% ungroup()
length(which(df3$method=="net")) # 1223
length(which(df3$method=="point")) # 246

stations = colnames(train1) ### IMPORTED FROM Peru_6_ModandPred_5PCA_3splf
rem = which(!df2$site_station %in% stations)
rem = df2$site_station[rem]
df3 = df2[!df2$site_station %in% rem,]

df3Time = df3[,-3]
df3Time = df3Time %>% group_by(site_station, am.pm) %>% summarise_all(sum) %>% ungroup()
keep = which(!duplicated(df3Time$site_station))
keep2 = df3Time[keep,]
length(which(keep2$am.pm=="AM")) # 1091
length(which(keep2$am.pm=="PM")) # 6
length(which(keep2$am.pm=="Unknown")) # 38
cc = keep2[,3:ncol(keep2)]
cc[cc>0]=1
keep2 = cbind(keep2[,2],cc)
morn = (which(colSums(keep2[keep2$am.pm=="AM",-1])>0)) # 349
aft = (which(colSums(keep2[keep2$am.pm=="PM",-1])>0)) # 7
unkn = (which(colSums(keep2[keep2$am.pm=="Unknown",-1])>0)) # 58
length(which(aft %in% morn)) # 6
length(which(unkn %in% morn)) # 57
length(which(aft %in% unkn)) # 5


df3Type = df3[,-2]
df3Type = df3Type %>% group_by(site_station, method) %>% summarise_all(sum) %>% ungroup()
keep = which(!duplicated(df3Type$site_station))
keep2 = df3Type[keep,]
length(which(keep2$method=="point")) # 86
length(which(keep2$method=="net")) # 1049
cc = keep2[,3:ncol(keep2)]
cc[cc>0]=1
keep2 = cbind(keep2[,2],cc)
nn = (which(colSums(keep2[keep2$`keep2[, 2]`=="net",-1])>0)) # 272
pp = (which(colSums(keep2[keep2$`keep2[, 2]`=="point",-1])>0)) # 161
length(which(pp %in% nn))
#######################################################################


# Create empty list or array to hold all species matrices
site.surv.spec = list()

## For each individual Species
for (i in 8:ncol(df1)) {
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
names(site.surv.spec) = colnames(df[8:ncol(df)])

ssPnt <- array (
  data = do.call(rbind, lapply(site.surv.spec, as.vector)), 
  dim = c(length(site.surv.spec), dim(site.surv.spec[[1]])),
  dimnames = list(colnames(df[8:ncol(df)]), transects$site_station, colnames(df3))
)
str(ssPnt)
save(ssPnt, file = "data/BIRDS_3_COUNTS_ALL_SiteANNUALLY-X-REPLICATE-per-Species_ARRAY.Rdata")
# load("data/BIRDS_3_COUNTS_ALL_SiteANNUALLY-X-REPLICATE-per-Species_ARRAY.Rdata")
