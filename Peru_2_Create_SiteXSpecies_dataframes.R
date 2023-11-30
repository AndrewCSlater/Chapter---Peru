library(dplyr)
library(tidyr)

# Load Count Data
count = read.csv("data/BIRDS_1.75_Combined_Counts_Unified_Station_Coords.csv")
sort(unique(count$species))

zone19 = read.csv("data/BIRDS_1_Stations_UTM19_Fixed_Names_Coords.csv")
zone18 = read.csv("data/BIRDS_1_Stations_UTM18_Fixed_Names_Coords.csv")


##### 1 - Check and Correct Species Names
### There are many species which are unidentified or misspelled

my_patterns = c('_sp$', 'sp\\.', 'SP$', '\\?', '^NA$', '\\-9', '\\-$', 'unknown', 'ident', 'nombre', 'Nota_de_voz', "Nota De Voz", 'VIDEO', 'band_lost', "Band Lost", "-18") # Vector of patterns found within species names that need to be removed

# Find all patterns within the vector
rem = grep(pattern = paste(my_patterns, collapse='|'), x = count$species, ignore.case = T) # Which lines contain the patterns for removal
count2 = count[-c(rem),]
rownames(count2) = NULL
rem = which(is.na(count2$species))
count2 = count2[-c(rem),]
sort(unique(count2$species))



# Remove brackets () from names and capitalise Genus & small letter species names
count2$species = gsub(pattern = "\\(", replacement = "",
                               gsub(pattern = "\\)", replacement = "", x = count2$species))
count2$species = gsub(pattern = "__", replacement = "_", x = count2$species)
count2$species = gsub(pattern = "cf. ", replacement = "", x = count2$species)
count2$species = gsub(pattern = "cf._", replacement = "", x = count2$species)
count2$species = gsub(pattern = "  ", replacement = " ", x = count2$species)
count2$species = gsub(pattern = " ", replacement = "_", x = count2$species)
count2$species = stringr::str_to_title(count2$species)

#####################################
# Compare similar spellings of species

u = sort(unique(count2$species)) # 516
# compare all strings against each other
d <- adist(u)
# Do not list combinations of similar words twice
d[lower.tri(d)] <- NA
# Say your threshold below which you want to consider strings similar is 4 edits:
a <- which(d > 0 & d < 4, arr.ind = TRUE)
pairs <- cbind(u[a[,1]], u[a[,2]])
pairs
## Will have to curate the results yourself to avoid accidental equalization of unequal factors.
## Do this reproducably by using a named vector as a translation dictionary. For example:
dict <- c(
  # incorrect spellings          correct spellings
  # -------------------------    ----------------------------
  "Arremon_tacilurnus"        =  "Arremon_taciturnus",
  "Automolus_rufipeleatus"    =  "Automolus_rufipileatus",
  "Chloroeryle_aenea"         =  "Chloroceryle_aenea",
  "Conophophaga_peruviana"    =  "Conopophaga_peruviana",
  "Cyanoloxia_rothschiidii"   =  "Cyanoloxia_rothschildii",
  "Cyhorhinus_arada"          =  "Cyphorhinus_arada",
  "Dendrexetastes_rufigala"   =  "Dendrexetastes_rufigula",
  "Dendrocicla_merula"        =  "Dendrocincla_merula",
  "Dendrocolaptes_picummus"   =  "Dendrocolaptes_picumnus",
  "Epinecrohylla_haematonota" =  "Epinecrophylla_haematonota",
  "Epinecrohylla_leucophthalma" = "Epinecrophylla_leucophthalma",
  "Epinecrphylla_haematonota" =  "Epinecrophylla_haematonota",
  "Eupetonema_macrura"        =  "Eupetomena_macroura",
  "Eupetonema_macroura"       =  "Eupetomena_macroura", 
  "Galbula_cyanascens"        =  "Galbula_cyanescens",
  "Geoteygon_montana"         =  "Geotrygon_montana",
  "Leptotila_rufaxila"        =  "Leptotila_rufaxilla",
  "Leptotilla_verreauxi"      =  "Leptotila_verreauxi",
  "Megarhynchus_pitangua"     =  "Megarynchus_pitangua",
  "Mitraphanes_olivaceus"     =  "Mitrephanes_olivaceus",
  "Monasa_morpheus"           =  "Monasa_morphoeus",
  "Mymotherula_axillaris"     =  "Myrmotherula_axillaris",
  "Myrmotherulamenetriesii"   =  "Myrmotherula_menetriesii",
  "Mytrephanes_olivaceus"     =  "Mitrephanes_olivaceus",
  "Onychorhynchos_coronatus"  =  "Onychorhynchus_coronatus",
  "Oryzoboros_angolensis"     =  "Oryzoborus_angolensis",
  "Phaetornis_guy"            =  "Phaethornis_guy",
  "Phaetornis_hispidus"       =  "Phaethornis_hispidus",
  "Phaetornis_stuarti"        =  "Phaethornis_stuarti",
  "Philydor_erythopterum"     =  "Philydor_erythropterum",
  "Picummus_dorbignyanus"     =  "Picumnus_dorbignyanus",
  "Piprites_choris"           =  "Piprites_chloris",
  "Platyrhinchus_platyrhynchos" = "Platyrinchus_platyrhynchos",
  "Rhamphocelus_carbo"        =  "Ramphocelus_carbo",
  "Rhegmatorhina_melanostica" =  "Rhegmatorhina_melanosticta",
  "Sclerurus_caudacatus"      =  "Sclerurus_caudacutus",
  "Sclerurus_cauducatus"      =  "Sclerurus_caudacutus",
  "Selenidera_reinwardii"     =  "Selenidera_reinwardtii",
  "Sittasomus_griseicapilus"  =  "Sittasomus_griseicapillus",
  "Thamnophilius_aethiops"    =  "Thamnophilus_aethiops",
  "Thraupis_episcopis"        =  "Thraupis_episcopus",
  "Troglodites_aedon"         =  "Troglodytes_aedon",
  "Turdus_hauxwelli"          =  "Turdus_hauxwellii",
  "Vireo_olivaceous"          =  "Vireo_olivaceus",
  "Colibri_corunscas"         =  "Colibri_coruscans",
  "Columbius_talpacoti"       =  "Columbina_talpacoti",
  "Eupetonema_macroura"       =  "Eupetomena_macroura",
  "Formicarius_rufrifroms"    =  "Formicarius_rufifrons",
  "Glaucis_hirsuta"           =  "Glaucis_hirsutus",
  "Hylophylax_naevia"         =  "Hylophylax_naevius",
  "Laniocera_hyphophyrra"     =  "Laniocera_hypopyrra",
  "Machareopterus_pyrocephalus" = "Machaeropterus_pyrocephalus",
  "Ncytidromus_albicollis"    =  "Nyctidromus_albicollis",
  "Phaethornis_phillipii"     =  "Phaethornis_philippii",
  "Phylidor_pyrrhodes"        =  "Philydor_pyrrhodes",
  "Schiffornis_turdinus"      =  "Schiffornis_turdina",
  "Sclaretia_naevia"          =  "Sclateria_naevia",
  "Troglodytes_solticialis"   =  "Troglodytes_solstitialis",
  "Coryspingus_cucullatus"    = "Coryphospingus_cucullatus",
  "Myiothlypsis_coronatus"    = "Myiothlypis_coronata",
  "Myrmotherula_menestresi"   = "Myrmotherula_menetriesii"
  # "Philydor_erythrocercum"    = "Philydor_erythropterum", 
  #"Ramphotrigon_fuscicauda"   = "Ramphotrigon_ruficauda",
  #"Thamnophilus_doliatus"     = "Thamnophilus_palliatus",
  #"Tinamus_major"             = "Tinamus_tao",
  #"Piaya_cayana"              = "Tityra_cayana"
)

# The correct levels need to be included, to
dict <- c(dict, setNames(u,u))

# Then convert factor column to character by using as.character and apply the dictionary on the original character vector df2$name_Species:
count2$species <- dict[count2$species]

u = sort(unique(count2$species)) # 456
# compare all strings against each other
d <- adist(u)
# Do not list combinations of similar words twice
d[lower.tri(d)] <- NA
# Say your threshold below which you want to consider strings similar is 3 edits:
a <- which(d > 0 & d < 2, arr.ind = TRUE)
pairs <- cbind(u[a[,1]], u[a[,2]])
pairs

write.csv(count2, "data/BIRDS_2_Counts_All_Corrected_Species_Spelling.csv", row.names = F)
sort(unique(count2$species))
count2 = read.csv("data/BIRDS_2_Counts_All_Corrected_Species_Spelling.csv")

##### 2 - Create a Site X Species Presence/Absence df
names(count2)
df = count2[,c(1:6,20)]
df$coords = paste0(df$utm_long,"_",df$utm_lat)
df$site_stn = paste0(df$site_code, "_", df$station)

# # To get station by year rather than station overall counts 
# df.yr = count2[,c(1:6,8,20)]
# df.yr$coords = paste0(df.yr$utm_long,"_",df.yr$utm_lat)
# df.yr$site_stn = paste0(df.yr$site_code, "_", df.yr$station, "_", df.yr$year)

length(unique(df$species)) # 456
length(unique(df$site_code)) # 25
length(unique(df$coords)) # 694
length(unique(df$site_stn)) # 694

##############################################################
# Add a COUNT column for each sighting
# df = df.yr
df$count = 1
df2 = distinct(df)
wide = as.data.frame(pivot_wider(data = df2, names_from = species, values_from = count, values_fill = 0))
sp = wide[,10:ncol(wide)]
sp = sp[,order(colnames(sp))]
wide = cbind(wide[,1:9], sp)


library(proj4)
# Find the projection proj4 description of the coordinate reference system
# from https://spatialreference.org/ref/?search=UTM+zone+19S
# UTM zone 19 South = (espg 32719)
# Get the proj4 Description from https://spatialreference.org/ref/epsg/32719/proj4/
proj19 <- "+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj18 <- "+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Transform the data
l8 = which(wide$utm_zone=="18L")
l9 = which(wide$utm_zone=="19L")
lon_lat18 <- project(wide[l8, 4:5], proj18, inverse = TRUE)
lon_lat19 <- project(wide[l9, 4:5], proj19, inverse = TRUE)

# convert to a data frame
x_y18 = data.frame(stations = wide$site_stn[l8], lon = lon_lat18$x, lat = lon_lat18$y)
x_y18$row = l8
x_y19 = data.frame(stations = wide$site_stn[l9], lon = lon_lat19$x, lat = lon_lat19$y)
x_y19$row = l9
xy89 = rbind(x_y18, x_y19)
xy89 = xy89[order(xy89$row),]


## Add new coords to Sites df
wide = tibble::add_column(wide, xy89[,2:3], .after = "utm_zone")

write.csv(wide, "data/BIRDS_2_Station-ANNUALLY_X_Species.csv", row.names = F)
df = wide[wide$utm_zone=="19L",c(11:ncol(wide))]
prev19 = as.data.frame(colSums(df))
prev19$species = rownames(prev19)
colnames(prev19)[1] = c("station_prevalence")
write.csv(prev19, "data/BIRDS_2_Species_Prevalence_UTM19.csv", row.names = F)

## Explore Count Data
## Species Prevalence by Site
summary(colSums(wide[9:ncol(wide)]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    1.00    4.00   18.36   19.25  302.00 
hist(colSums(wide[9:ncol(wide)]), breaks = 24)

## Species Richness per Site
summary(rowSums(wide[9:ncol(wide)]))
# Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 1.00    6.00    9.00    12.06   16.75   67.00 


### EXAMINE COUNTS BY SITE RATHER THAN STATION
site = wide[,c(2,9:ncol(wide))]
site_counts = as.data.frame(site %>%
  group_by(site_code) %>%
  summarise_all(sum))
rownames(site_counts) = site_counts$site_code
site_counts = site_counts[-1]
site_counts[site_counts>0] = 1

write.csv(site_counts, "data/BIRDS_2_!SITE!_X_Species.csv", row.names = T)

## Species Prevalence by Site
summary(colSums(site_counts))
# Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
# 1.00    1.00    2.00    4.73    7.00   23.00 
hist(colSums(site_counts))

## Species Richness per Site
summary(rowSums(site_counts))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00     41   82.00    86.24  124.00  193.00 


### Add Site Counts to SITES data
names(wide)[names(wide)=="site_stn"] = "stations"
count18 = left_join(zone18, wide[,c(8:ncol(wide))], by = "stations")
count19 = left_join(zone19, wide[,c(8:ncol(wide))], by = "stations")

library(proj4)
# Find the projection proj4 description of the coordinate reference system
# from https://spatialreference.org/ref/?search=UTM+zone+19S
# UTM zone 19 South = (espg 32719)
# Get the proj4 Description from https://spatialreference.org/ref/epsg/32719/proj4/
proj19 <- "+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj18 <- "+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Transform the data
lon_lat18 <- project(count18[, 4:5], proj18, inverse = TRUE)
lon_lat19 <- project(count19[, 4:5], proj19, inverse = TRUE)

# convert to a data frame
x_y18 = data.frame(stations = count18$stations, lon = lon_lat18$x, lat = lon_lat18$y)
x_y19 = data.frame(stations = count19$stations, lon = lon_lat19$x, lat = lon_lat19$y)

## Add new coords to Sites df
count18 = tibble::add_column(count18, x_y18[,2:3], .after = "GPS.Zone")
count19 = tibble::add_column(count19, x_y19[,2:3], .after = "GPS.Zone")

write.csv(count18, "data/BIRDS_2_Stations_UTM18_inc_SPECIES.csv", row.names = F)
write.csv(count19, "data/BIRDS_2_Stations_UTM19_inc_SPECIES.csv", row.names = F)

names(count19)
xy = rbind(count18[,c(1:6, 18:19)], count19[,c(1:6, 26:27)])

library(proj4)
# Find the projection proj4 description of the coordinate reference system
# from https://spatialreference.org/ref/?search=UTM+zone+19S
# UTM zone 19 South = (espg 32719)
# Get the proj4 Description from https://spatialreference.org/ref/epsg/32719/proj4/
proj19 <- "+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj18 <- "+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Transform the data
l8 = which(count2$utm_zone=="18L")
l9 = which(count2$utm_zone=="19L")
lon_lat18 <- project(count2[l8, 4:5], proj18, inverse = TRUE)
lon_lat19 <- project(count2[l9, 4:5], proj19, inverse = TRUE)

# convert to a data frame
x_y18 = data.frame(stations = count2$stations[l8], lon = lon_lat18$x, lat = lon_lat18$y)
x_y18$row = l8
x_y19 = data.frame(stations = count2$stations[l9], lon = lon_lat19$x, lat = lon_lat19$y)
x_y19$row = l9
xy89 = rbind(x_y18, x_y19)
xy89 = xy89[order(xy89$row),]


## Add new coords to Sites df
count2 = tibble::add_column(count2, xy89[,2:3], .after = "utm_zone")

write.csv(count2, "data/BIRDS_2_Counts_All_Corrected_Species_Spelling.csv", row.names = F)
#################################################################################



