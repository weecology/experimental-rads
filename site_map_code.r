# Rmaps for data
library(maps)
library(mapdata)
library(maptools)
library(scales)
library(RgoogleMaps)
library(ggplot2)
library(RColorBrewer)
library(mapproj)


wd = "C:\\Users\\sarah\\Documents\\GitHub\\experimental-rads"
#wd = "/Users/sarah/Desktop/WE_Dropbox_stuff/ExperimentalMacroecology/"
setwd(wd)

sites = read.csv("sites_analysis_data.csv", na.strings = 'NULL', colClasses = "character")

lat = as.numeric(sites$latitude)
lon = as.numeric(sites$longitude)

map("world", col = "black", fill = TRUE, bg = "white", lwd = 0.05)
points(lon, lat, col = "yellowgreen", pch = 19, cex = 0.75)
axis(side = 1)
axis(side = 2)

#read in a shapefile (saved from blm website)
world<- readShapePoly("continent.shp")

# Create a lat-long dataframe from the maps package
world <- map_data("world")
worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill="white", colour="black")

# Use cartesian coordinates
worldmap
# With default mercator projection
worldmap + coord_map()
# Other projections
worldmap + coord_map("cylindrical")
worldmap + coord_map("azequalarea",orientation=c(-36.92,174.6,0))
worldmap + coord_map("ortho")

#add points to ggplot2 map
pts <- data.frame(x = lon, y = lat)
worldMap <-ggplot(world) + geom_path(aes(x,y))
worldMap + geom_point(data = pts, aes(x, y))

