require(ggplot2)
require(ggspatial)
require(sf)
library(rgdal)
library(terra)
require(rnaturalearth)
require(rnaturalearthdata)

## Bring In Survey Point Data ##
site = read.csv("data/Model_Input_Fixed_Site_Covariates.csv")


#get an object which has South America
SAm = ne_countries(scale = "medium", continent = "South America", returnclass = "sf")
plot(SAm[6])

theme_set(theme_bw()) #optional, dark-on-light theme is best for maps

chart = ggplot(data=site) +
  geom_sf(data=SAm, colour="black") +
  # set background colour and text size  
  theme(panel.background=element_rect(fill='aliceblue'), text=element_text(size=12)) +
  # Add a map scale bar
  annotation_scale(location = "br",
                   bar_cols = c("grey60", "white"),
                   width_hint = 0.25) +
  # Add North Arrow
  annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(1, "cm"), pad_y = unit(1, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20")) +
  # set the bounding coordinates of the area of the map you want to include
  coord_sf(crs = st_crs(5373), xlim = c(-70.5, -68.5), ylim = c(-11, -14), expand = T) +
  # , datum = NA) + # This Removes the Coordinate Axis Values
  # specify what your coordinate columns are called, and if applicable which column you want to use to colour-code the points
  
  geom_point(data = site, mapping=aes(x=lon, y=lat, color = Site_Code)) +
  # put labels for your axes and legend
  labs(x="", y="", fill = "Site_Code") +
  theme(legend.text=element_text(size = 8)) +
  theme(legend.title=element_text(size = 10))
chart
