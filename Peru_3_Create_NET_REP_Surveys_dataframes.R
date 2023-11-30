library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)

#######################################################################################
####  FOR NET COUNTS
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
length(unique(tt_day$Site_Station)) # There are 3036 Surveys spread over 648 Unique Locations
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

write.csv(morn2, "data/Det-Covs_NETS_AM-PM_per_site_Replicate.csv", row.names = T)



####################
# Read In NET SURVEY EFFORT Data
library(readxl)
net.effort = read_excel("data/RAW New2 - FF Cleaned Master Bird Database.xlsx", "Daily Effort", range = cell_rows(1:1281))
names(net.effort)
net.effort = select(net.effort, "SITE CODE", "DATE", "NUMBER OF NETS", "MORNING OR AFTERNOON", "NET EFFORT TIME...12", "NET STATION #1", "NET STATION #2", "NET STATION #3", "NET STATION #4", "NET STATION #5...17", "NET STATION #5...18", "...19")
names(net.effort) = c("site", "date", "n_nets", "AM.PM", "net_hours", "stn1", "stn2", "stn3", "stn4", "stn5", "stn6", "stn7")

net.effort$AM.PM[net.effort$AM.PM=="MORNING"] = "AM"
net.effort$AM.PM[net.effort$AM.PM=="AFTERNOON"] = "PM"
net.effort$AM.PM[net.effort$AM.PM=="ALL-DAY"] = "both"

net.effort$net_hours = as.numeric(net.effort$net_hours)
net.effort$n_nets = as.numeric(net.effort$n_nets)
### Some surveys do not record Number of Nets used or Hours the nets were open : ie survey EFFORT
length(which(is.na(net.effort$net_hours))) # 131 Surveys have no net hours recorded
hist(net.effort$net_hours)
length(which(is.na(net.effort$n_nets))) # 50 Surveys do not record the number of nets
hist(net.effort$n_nets)
### We are not overly interested in measuring the impact of net size & time, but using them may improve our ability to determine where species occur, so I will keep all data & Replace NA values with MEAN values as the data seems normally distributed
net.effort$net_hours[is.na(net.effort$net_hours)] = mean(net.effort$net_hours, na.rm = T)
net.effort$n_nets[is.na(net.effort$n_nets)] = mean(net.effort$n_nets, na.rm = T)

net.effort$date = as.Date(net.effort$date)
length(which(net.effort$date %in% dist_net$date)) # 1153 of 1280 lines have dates reciprocated in the Count Dates
dist_net$dateTime = paste0(dist_net$date,"_", dist_net$AM.PM)
net.effort$dateTime = paste0(net.effort$date,"_", net.effort$AM.PM)
length(which(net.effort$dateTime %in% dist_net$dateTime)) # 1049 of 1280 lines have dates-AM/PM reciprocated in the Count Dates

net.effort = select(net.effort, 1:5, 13, everything())
net.effort = tibble::add_column(net.effort, Site_Station = paste0(points$Site_Station, "_", points$date, "_", points$AM.PM), .before = "Site_Station")





###### Up to 7 Stations are used per row. These need to be split into 1 per row
st1 = net.effort[,c(1:6,7)]
st2 = net.effort[,c(1:6,8)]
st3 = net.effort[,c(1:6,9)]
st4 = net.effort[,c(1:6,10)]
st5 = net.effort[,c(1:6,11)]
st6 = net.effort[,c(1:6,12)]
st7 = net.effort[,c(1:6,13)]

## Create a Site Code value in Count Data 
dist_net$site = sub(pattern = "_.*", replacement = "", x = dist_net$Site_Station)
length(unique(dist_net$site)) # 25 Unique Sites
length(unique(net.effort$site)) # 28 Unique Sites ### What are different?
# Visual inspection
unique(net.effort$site)[!unique(net.effort$site) %in% dist_net$site]
unique(dist_net$site)
net.effort$site = gsub(pattern = "-WALTER'S HOUSE", replacement = "", x = net.effort$site)
net.effort$site = gsub(pattern = "AZUL", replacement = "AZU", x = net.effort$site)

st1$Site_Station = paste0(st1$site, "_", st1$stn1)
st2$Site_Station = paste0(st2$site, "_", st2$stn2)
st3$Site_Station = paste0(st3$site, "_", st3$stn3)
st4$Site_Station = paste0(st4$site, "_", st4$stn4)
st5$Site_Station = paste0(st5$site, "_", st5$stn5)
st6$Site_Station = paste0(st6$site, "_", st6$stn6)
st7$Site_Station = paste0(st7$site, "_", st7$stn7)






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









