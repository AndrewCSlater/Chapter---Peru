library(dplyr)

#################
# THERE IS AN ISSUE WITH INDIVIDUAL COUNTS APPEARING TWICE WITH DIFFERENT EFFORTS _ THIS NEEDS SORTING




## Stations & Coords of POINT Count
points_in = read.csv("data/BIRDS_1_Counts_Points_Fixed_SitesNames.csv")
names_point = c("site_name", "site_code", "date", "station", "utm_long", "utm_lat", "utm_zone", "observer", "time_start", "time_end", "duration_mins", "temp_start", "temp_end", "cloud_pct", "humidity", "name_Family",  "name_common", "name_Species", "name_code", "seen_heard_both", "number",  "distance_to_animal", "notes")
names(points_in) = names_point
points_in$utm_long = round(as.numeric(points_in$utm_long), digits = 0)
points_in$utm_lat = round(as.numeric(points_in$utm_lat), digits = 0)
points_in = select(points_in, "site_name", "site_code", "station", "utm_long", "utm_lat", "utm_zone", "date", "time_start", "duration_mins", "name_Species")
points_in = tibble::add_column(points_in, method = "point", .after = "utm_zone")

nets_in = read.csv("data/BIRDS_1_Counts_Nets_Fixed_SitesNames.csv")
nn = names(nets_in)
nn = c("Species..Scientific.", "Date..dd.mm.yyyy.", "Year", "Month", "Capture.Time", "Site.Name", "Site.Code", "Station", "Station.GPS..Easting.", "Station.GPS..Northing.", "Station.GPS.Zone", "First.Net.Time.Open", "Last.Net.Time.Open", "Time.open", "Median.Time.Opened", "First.Net.Time.Close", "Last.Net.Time.Close", "Time.Closed", "Median.Time.Closed", "Sample.Time..hrs.", "Number.of.Nets..12m.")
nets_in = select(nets_in, all_of(nn))
nets_in = select(nets_in, "Site.Name", "Site.Code", "Station", "Station.GPS..Easting.", "Station.GPS..Northing.", "Station.GPS.Zone", "Species..Scientific.", everything())
names(nets_in)[1:7] = c("site_name", "site_code", "station", "utm_long", "utm_lat", "utm_zone", "name_Species")
nets_in$utm_long = round(as.numeric(nets_in$utm_long), digits = 0)
nets_in$utm_lat = round(as.numeric(nets_in$utm_lat), digits = 0)
nets_in = tibble::add_column(nets_in, method = "net", .after = "utm_zone")

#### Combine Nets & Points
df = bind_rows(points_in, nets_in)

## Change erroneous values to NA
df[df=="NA"] = NA
df[df=="-" | df=="-9"] = NA

### Find and remove counts which have no station name or coordiantes
rem = which(is.na(df$station) & is.na(df$utm_lat) | is.na(df$utm_long))
df = df[-rem,]
rownames(df) = NULL

### Find and remove counts which have no date
rem = which(is.na(df$date) & is.na(df$Date..dd.mm.yyyy.))
df = if(length(rem>0)) {df[-rem,]
} else {df
    } 

### Correct & Unify date
swap = which(is.na(df$date))
df$date[swap] = df$Date..dd.mm.yyyy.[swap]
df = select(df, -Date..dd.mm.yyyy.)
rownames(df) = NULL
df$Year = as.integer(substr(df$date, 1, 4)) # create a numeric year for each count by taking first 4 digits of the date

### Correct & unify Time to AM.PM
df$time_start = as.POSIXct(df$time_start)                      
df = tibble::add_column(df, AM.PM = format(df$time_start, "%p"), .after = "time_start")

sort(unique(df$Capture.Time))
df$Capture.Time[df$Capture.Time == "######" | df$Capture.Time == "*" | df$Capture.Time == "0" |  df$Capture.Time == "1st" | df$Capture.Time == "2nd" | df$Capture.Time == "3rd" | df$Capture.Time == "4th"] = NA
df$AM.PM[df$Capture.Time == "AFTERNOON"] = "PM"

### There are several time formats that I need to unify
df$Capture.Time = gsub(pattern = ":", replacement = ".", x =  df$Capture.Time) # Change : to . so string can be converted to number
df$Capture.Time = as.numeric(df$Capture.Time)                                  # Convert to number
# Some numbers include whole days (integer) and fractions of day (decimal)
swap = which(!is.integer(df$Capture.Time) & df$Capture.Time>1000)              # Find which numbers are big enough to need converting     
df$Capture.Time[swap] = df$Capture.Time[swap] - as.integer(df$Capture.Time)[swap] #  Remove the number of days (integer) to leave fraction of day

## Turn Fractions into Hours of the day
swap = which(df$Capture.Time<1) 
df$Capture.Time[swap] = round(df$Capture.Time[swap]*24, digits = 2)

### Some times are 24 hr without decimal
swap = which(df$Capture.Time > 24) # Find which numbers are bigger than 24
df$Capture.Time[swap] = df$Capture.Time[swap]/100 # Divide by 100 to turn into decimal time

# sort(unique(as.numeric(df$Capture.Time[!is.na(df$Capture.Time)])))

swap =  which(df$Capture.Time < 12 & !is.na(df$Capture.Time))
df$AM.PM[swap] = "AM"
swap = which(df$Capture.Time >= 12 & !is.na(df$Capture.Time))
df$AM.PM[swap] = "PM"
 
names(df)
df2 = select(df, 1:6, 8, 13:14, 10, 7, 11, 24:25, 12)
df2$Number.of.Nets..12m. = gsub(pattern = " x .*", replacement = "", x = df2$Number.of.Nets..12m.)
df2$Number.of.Nets..12m. = gsub(pattern = " X .*", replacement = "", x = df2$Number.of.Nets..12m.)
df2$Number.of.Nets..12m. = as.integer(df2$Number.of.Nets..12m.)

df2$Sample.Time..hrs. = 24*df2$Sample.Time..hrs.

#####
# Create an effort column by scaling number of nets, net time & count time
df2$duration_mins[is.na(df2$duration_mins)] = 0
df2$point_effort = scales::rescale(df2$duration_mins)

df2$Sample.Time..hrs.[df2$Sample.Time..hrs.==0] = NA
df2$net_effort_time = df2$Sample.Time..hrs.
df2$net_effort_time[is.na(df2$Sample.Time..hrs.) & !is.na(df2$Number.of.Nets..12m.)] = mean(df2$Sample.Time..hrs., na.rm = T)

df2$net_effort_nets = df2$Number.of.Nets..12m.
df2$net_effort_nets[is.na(df2$Number.of.Nets..12m.) & !is.na(df2$Sample.Time..hrs.)] = mean(df2$Number.of.Nets..12m., na.rm = T)

df2$net_effort = df2$net_effort_time * df2$net_effort_nets
df2$net_effort = scales::rescale(df2$net_effort)

df2$effort = pmax(df2$net_effort, df2$point_effort, na.rm = T)

df2$Month = as.integer(substr(df2$date, 6, 7))

names(df2)
n = c( "site_name", "site_code", "station", "utm_long", "utm_lat", "utm_zone", "date", "year", "month", "am.pm", "method", "point_mins", "net_hrs", "nets_qty", "species", "effort_point", "effort_net_time", "effort_net_nets", "effort_net_total", "effort")
names(df2) = n
df2 = select(df2, "site_name", "site_code", "station", "utm_long", "utm_lat", "utm_zone", "date", "year", "month", "am.pm",  "point_mins", "net_hrs", "nets_qty", "effort_point", "effort_net_time", "effort_net_nets", "effort_net_total", "effort", "method", "species")

### In a few cases the effort recorded per survey varies between recordings.
### Take the mean effort per Survey (Station - Date - AM.PM)
df3 = df2 %>% group_by(site_code, station, date, am.pm) %>%
  mutate(effort = mean(effort)) %>%
           ungroup

write.csv(df3, "data/BIRDS_1.5_Combined_NetPoint_Counts.csv", row.names = F)

