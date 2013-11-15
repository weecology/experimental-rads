# Rmaps for data
library(maps)
library(mapdata)
library(maptools)
library(scales)
library(RgoogleMaps)
library(ggplot2)
library(RColorBrewer)
library(mapproj)
library(ggmap)


wd = "C:\\Users\\sarah\\Documents\\GitHub\\experimental-rads"
setwd(wd)


comms = read.csv("community_analysis_data.csv", na.strings = 'NULL')
comps = read.csv("orderedcomparisons.csv")
  names(comps)<-c('ref', 'controID','expID')
sites = read.csv("sites_analysis_data.csv", na.strings = 'NULL', colClasses = "character")
expers = read.csv("experiments_analysis_data.csv")


#---------------------------
#Just grab loc, taxa and experiment type from data used in the study
taxa = c()
type = c()
lon = as.numeric()
lat = as.numeric()

for (iRow in 1:nrow(comps)){
  control = comps[iRow,2]  #find control in pair
  experiment = comps[iRow,3]  # find experiment in pair
  # Check that < 10% of individuals are unidentified. If meets criteria, continue
  if (percent_unidSpp(control, comms) == "OK" & percent_unidSpp(experiment, comms) == "OK"){
    a1 = sort(as.numeric(comms[which(comms[,2] == control & comms[,7] != 0), 8])) #vector of control abundances
    a2 = sort(as.numeric(comms[which(comms[,2] == experiment & comms[,7] != 0), 8])) #vector of exp abundances
    # Check that there are at least 5 species and 30 individuals in each community, If yes, proceed.
    if (length(a1) > 4 & length(a2) > 4 & sum(a1) > 29 & sum(a2) > 29){
      cloc = as.numeric(sites[which(sites[,2] == control), c(7,6)])
      eloc = as.numeric(sites[which(sites[,2] == experiment), c(7,6)])
      lon = append(lon, cloc[1])
      lon = append(lon, eloc[1])
      lat = append(lat, cloc[2])
      lat = append(lat, eloc[2])
      #add each of these twice
      taxa = append(taxa, as.character(expers[which(expers[,2]==control),7]))# find taxonomic group from experiments table
      type = append(type, as.character(expers[which(expers[,2]==control),4])) # find experiment type from experiments table
      taxa = append(taxa, as.character(expers[which(expers[,2]==experiment),7]))# find taxonomic group from experiments table
      type = append(type, as.character(expers[which(expers[,2]==experiment),4])) # find experiment type from experiments table
}}}

#collapse taxon types into broader categories so there aren't so many factors
taxa[taxa=='carabid']<-'insect'
taxa[taxa=='lepidopteran']<- 'insect'
taxa[taxa=='odonate']<- 'insect'
taxa[taxa=='orthoptera']<-'insect'
taxa[taxa=='orthoptera ']<-'insect'
taxa[taxa=='beetle']<-'insect'
taxa[taxa=='microarthropods']<-'microarthropod'
taxa[taxa=='reptile']<-'herpetofauna'


data = as.data.frame(cbind(lon, lat))

#-------------------------------------------------
##### from : http://stackoverflow.com/questions/16028659/plots-on-a-map-using-ggplot2
#Get world map info
world_map <- map_data("world")

#Creat a base plot
p <- ggplot() + coord_fixed()

#Add map to base plot
base_world <- p + geom_polygon(data=world_map,aes(x=long, y=lat,group=group), fill = "white", col = "black") + theme_bw()

base_world + geom_point(data=data, aes(x=lon, y=lat, col = taxa), alpha = 0.5, cex = 5) +
              theme(legend.position = "right") + element_blank()




