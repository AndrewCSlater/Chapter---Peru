library(remotes)
install_github('MDecuy/AVOCADO')  
library(AVOCADO)


??AVOCADO



#Loading the raster stack
MDD <- stack(system.file("Landsat_NDMI_2010-01-01_2018-01-01.tif", package = "AVOCADO"))
#Extract the dates from the brick
lan.info <- getSceneinfo(names(MDD))
#example output should look like: LT51670651984137     TM  167  65 1984-05-16
lan.dates <- as.Date(lan.info$date)
#example output should look like: "1984-05-16"

a = stack(system.file())