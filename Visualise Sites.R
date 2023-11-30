library(openxlsx)
library(ggplot2)
library(rgdal) # Convert UTM to Lon Lat
library(geosphere) # Measure distances between Lon_Lat points
library(here)
here()

point = read.csv(file = "data/wrangle2_point_count_wide_speciesXsite_date.csv")
net = read.csv(file = "data/wrangle2_net_count_wide_speciesXsite_date.csv")
sites = read_excel("data/RAW - Fauna Forever - Bird Counts.xlsx", sheet =   4)

p = select(point, 1:7)
n = select(net, 1:5)
s = select(sites, 1:5)

s = rename(s, station = "Station Code")

p$net_point = as.factor("point")
p$count_date = as.Date(p$count_date, format = "%d/%m/%Y")
n$net_point = as.factor("net")
n$count_date = as.Date(n$count_date, format = "%Y-%m-%d")

# To bind point & net together they must have the same columns and names
p = rename(p, Site_Code = "Site.Code", Site_Name = "Site.Name")
n$UTM.Easting = as.integer(NA)
n$UTM.Northing = as.integer(NA)
all = rbind(p,n)

# Some Station Names are repeated over several Sites and so a Unique Site_Station code is needed
all$unique = as.factor(paste0(all$Site_Code,"_",all$station))
s$unique =  as.factor(paste0(s$`Site Code`,"_",s$station))

# Add Month & Year columns for future analysis
all$Year = format(all$count_date, format = "%Y")
all$Month = format(all$count_date, format = "%m")

# Filter Net & Point data out to compare their unique values vs total unique values
nn=filter(all, net_point == "net")
pp=filter(all, net_point == "point")
length(unique(nn$unique))
length(unique(pp$unique))

# Make sure all names and codes are as Factors
all$Site_Code = as.factor(all$Site_Code)
all$Site_Name = as.factor(all$Site_Name)
all$station = as.factor(all$station)
all$unique = as.factor(all$unique)

length(unique(all$unique))

length(unique(all$Site_Name))
length(unique(p$Site_Name))
length(unique(n$Site_Name))

length(unique(all$station))
length(unique(n$station))
length(unique(p$station))

length(unique(all$count_date))
length(unique(p$count_date))
length(unique(n$count_date))

# Combine Count Data with Location Coordinates
all2 = full_join(all,s, by = "unique")
names(all2)

# Create a column of North & East - If there is no coordinates from the Sites spreadsheet
# Then copy across coordinates from the Points column, because all point counts had coordinates integral in data
all2$Easting = ifelse(is.na(all2$`UTM Easting`),all2$UTM.Easting,all2$`UTM Easting`)
all2$Northing = ifelse(is.na(all2$`UTM Northing`),all2$UTM.Northing,all2$`UTM Northing`)
# Remove the extra coord columns to leave the united one
all2 = select(all2, -"UTM Easting", -"UTM Northing", -"UTM.Easting", -"UTM.Northing")

# Do the same for site name & station
all2$`Site Name`= as.character(all2$`Site Name`)
all2$Site_Name = as.character(all2$Site_Name)
all2$site_name = ifelse(is.na(all2$Site_Name),all2$`Site Name`,all2$Site_Name)

all2$station.y = as.character(all2$station.y)
all2$station.x = as.character(all2$station.x)
all2$station = ifelse(is.na(all2$station.y),all2$station.x,all2$station.y)

# Remove any lines without a Coordinate Reference
qq = all2
qq$Northing=ifelse(qq$Northing=="NA",NA,qq$Northing)
qq= filter(qq, !is.na(qq$Northing))
length(unique(qq$unique))
length(unique(qq$count_date))
length(unique(qq$site_name))

#####################################################
### Add LONG LAT Coordinates from UTM Coordinates ###
utm1 = data.frame(as.integer(qq$Easting), as.integer(qq$Northing))

sputm = SpatialPoints(utm1, proj4string = CRS("+proj=utm + zone=19L + south + datum=WGS84 +units=m"))
spgeo = spTransform(sputm, CRS("+proj=longlat + datum=WGS84"),inverse=T) # I think Inverse=T refers to Southern Hemisphere
xy=as.data.frame(spgeo@coords)
qq$Lon = xy$as.integer.qq.Easting.
qq$Lat = xy$as.integer.qq.Northing.
plot(qq$Lon, qq$Lat)

####################################
#### VISUALISATION OF LOCATIONS ####
stations = select(qq, 6:9,13:18)
stations$Northing = as.numeric(stations$Northing)
stations$Easting = as.numeric(stations$Easting)
stations$Year = as.factor(stations$Year)
stations$Month = as.factor(stations$Month)

stations = filter(stations, !is.na(stations$Year))

##### Get Mean Location of All Sites #####
st = group_by(stations, net_point, site_name) %>%
  summarise(mean(Easting, na.rm = T), mean(Northing, na.rm = T), mean(Lon, na.rm = T), mean(Lat, na.rm = T)) %>%
  rename(Easting = "mean(Easting, na.rm = T)", Northing = "mean(Northing, na.rm = T)", Lon = "mean(Lon, na.rm = T)", Lat =  "mean(Lat, na.rm = T)")
st = filter(st, ! is.na(net_point))

##### Create a list of Site Mean Coordinates #####
GE = ungroup(st) %>%  select(2,5,6)
GE = group_by(GE, site_name) %>%
  slice_head() %>% ungroup() %>% as.data.frame()
rownames(GE)=GE$site_name
GE = GE[,-1]
write.csv(GE, file = "data/Coords_Peru_Long_LAt.csv")

##### Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent kms #####
st = st[order(st$Easting),]
ww = as.data.frame(st$Easting)
ww$E = as.integer(ww$`st$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
st$Easting = ww$E

st = st[order(st$Northing),]
ww = as.data.frame(st$Northing)
ww$E = as.integer(ww$`st$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
st$Northing = ww$E

length(unique(st$site_name))
length(unique(st$Easting))

####################################################
######## Plot Mean Coordinates of All Sites ########
ggplot(st, aes(x=Easting, y=Northing)) +
  geom_point(aes(shape = net_point)) + # Show dots
  geom_text( 
    label=st$site_name, 
    size = 2,
    nudge_x = 70, nudge_y = 5, 
    check_overlap = T) +
  ggtitle("Point Counts - Mean Coordinates") + 
      xlab("East of Minimum") + ylab("North of Minimum") 

###################################################  
##### Filter Each Site to Plot their Stations #####
qu = filter(stations, site_name == "Quellomayo")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = net_point, main = "Quellomayo", data = qu)
##########
qu = filter(stations, site_name == "Boca Pariamanu Native Community" | site_name == "Boca Parimanu Native Community")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = net_point, data = qu, main = "Boca Pariamanu Native Community")
##########
qu = filter(stations, site_name == "Amazon Rainforest Conservation Centre" | site_name == "Las Piedras Amazon Centre" | site_name == "Las Piedras Biodiversity Station")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = net_point, data=qu, main = "Amazon - Rainforest CC & Las Piedras Bio & Amazon ")
##########
qu = filter(stations, site_name == "Los Amigos Biological Station" | site_name == "Reserva Amazonica")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = Year, data=qu, main = "Los Amigos & Reserva Amazonica - All Nets")
##########
qu = filter(stations, site_name == "Tambopata Research Centre")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = Year, data=qu, main = "Tambopata Research Centre, all Points")
##########
qu = filter(stations, site_name == "Wilderness International 1" | site_name == "Hacienda Tambopata" | site_name == "Neotropical Field Station" | site_name == "Collpas Tambopata Inn")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = site_name, data=qu, main = "Four Sites - Point,Point,Net,Net")
##########
qu = filter(stations, site_name == "Secret Forest")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = net_point, data=qu, main = "Secret Forest")
##########
qu = filter(stations, site_name =="Sachavacayoc Centre" | site_name == "El Gato Homestay Baltimore")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = Year, data=qu, main = "El Gato - Point_Net, Sacha - Net")
##########
qu = filter(stations, site_name == "Tambopata Lodge" | site_name == "Venado Station" | site_name == "Explorer's Inn" | site_name == "Explorer's Inn - Farm" | site_name == "Explorer's Inn Farm")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = site_name, data=qu, main = "Venado, Tambo, Ex Inn & Farm")
##########
qu = filter(stations, site_name == "Saona Lodge")
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
qu = qu[order(qu$Easting),]
ww = as.data.frame(qu$Easting)
ww$E = as.integer(ww$`qu$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Easting = ww$E

qu = qu[order(qu$Northing),]
ww = as.data.frame(qu$Northing)
ww$E = as.integer(ww$`qu$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
qu$Northing = ww$E
qplot(Easting, Northing,  colour = Year,  shape = net_point, data=qu, main = "Saona Lodge")
##########




##############################################################################
#####     From ALEX meeting 15-4-21     ######################################

names(all)
table(all$Site_Name)
table(all$Site_Name, all$Year)

zdat = as.data.frame(matrix(table(all$Site_Name, all$Year), nrow=length(unique(all$Site_Name))))
names(zdat) = sort(unique(all$Year)) # Adds Column Names (Year)
rownames(zdat) = sort(unique(all$Site_Name)) # Adds Row Names (Site_Name)
zdat 

# How many sites per year?
colSums(zdat>0)

# How many surveys (net & point) per year
colSums(zdat)

sum(zdat)

# How many times have there been surveys in a given Year/Month at a different Site?
(table(all$Site_Name, all$Year, all$Month)>0) # Site_Name = Rows, Year = Columns, Month = Different Table

######################################################

secretF.surveys = all[all$Site_Name=="Secret Forest",]

zdat = table(secretF.surveys$Year, secretF.surveys$Month)
zdat = as.data.frame(matrix( zdat, nrow=length(unique(secretF.surveys$Year))))
names(zdat) = sort(unique(secretF.surveys$Month))
rownames(zdat) = sort(unique(secretF.surveys$Year))
zdat

secretF.surveys = droplevels(secretF.surveys)
zdat = table(secretF.surveys$Year, secretF.surveys$station)
zdat = as.data.frame(matrix( zdat, nrow=length(unique(secretF.surveys$Year))))
names(zdat) = sort(unique(secretF.surveys$station))
rownames(zdat) = sort(unique(secretF.surveys$Year))

# Number of surveys per year, per station at Secret Forest
zdat
# Number of surveys per station
summary(colSums(zdat))



###################################################################
######## Import Carbon Store TIF from Asner - 16-4-21 #############
       #  Project Site Coords onto it  #

library(raster) # To import TIF
library(rgdal) # To convert and reproject CRS

tif = raster("data/peru_acd.tif") # Above Ground Carbon Density in Mg(Tonnes)/Ha
tif
plot(tif)

xy = as.data.frame(st[,5:6])
# NB: The SpatialPointDataFrame MUST be its own CRS, NOT what you want it projected to.
spdf <- SpatialPointsDataFrame(coords = xy, data = xy,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))

p <- spTransform(spdf, crs(tif)) # This is where you transorfm the CRS from what it is in SPDF into what you want it to be (in this case the CRS of tif)    

points(p$Lon, p$Lat) # This projects the converted points onto the current tif image
# NB: the Column Names are those of the original xy data.

# To Extract TIF VALUES from defined Vector Points
CHa = extract(tif, p) # Values of tif, at points p
st$C_Ha = CHa

### Plot Sites & Carbon Values ###
# Reduce the UTM numbers to run from 1 to Max - and to /1000 to represent km??
st = st[order(st$Easting),]
ww = as.data.frame(st$Easting)
ww$E = as.integer(ww$`st$Easting`-ww[1,1]+1)
ww$E = ww$E/1000
st$Easting = ww$E

st = st[order(st$Northing),]
ww = as.data.frame(st$Northing)
ww$E = as.integer(ww$`st$Northing`-ww[1,1]+1)
ww$E = ww$E/1000
st$Northing = ww$E
qu=filter(st, !is.na(C_Ha))
qplot(Easting, Northing,  colour = C_Ha,  data=qu, main = "Peru Sites - AG Carbon per Hectare",xlab="East from Minimum", ylab = "North from Minimum")

###########################################





