# All files should be organized as such:
# RAD: SiteID | experiment/control | year | genus | species | id2species | abundance
# Code will compare S, N and composition among treatments at each site,
# output and compare RADs, parameters for each treatment at each site.

library(ggplot2)

#---------------------------------------------------------------------------------
#          setup - select wd, import data, source code,  file to collect results
#---------------------------------------------------------------------------------

#wd = "/Users/sarah/Documents/GitHub/experimental-rads"
wd = "C:\\Users\\sarah\\Documents\\GitHub\\experimental-rads"
setwd(wd)

source("ExpRADsFunctions.R")   #Run the code containing the functions

comms = read.csv("community_analysis_data.csv", na.strings = 'NULL')
#comps = read.csv("comparison_analysis_data.csv") #this file is unordered. Looks less nice when plotted.
comps = read.csv("orderedcomparisons.csv")
  names(comps)<-c('ref', 'controID','expID')
expers = read.csv("experiments_analysis_data.csv")
sites = read.csv("sites_analysis_data.csv")

#--------------------------------------------------------------
#          generate values and comparisons between the sites 
#--------------------------------------------------------------
#open plotting pdf window
pdf("allRADs_orderedbytaxa.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(5,4), mar=c(1.5,2,3,1), oma=c(1,1,1,1))

#descriptive variables
refID = c()
cID = c()
eID = c()
Cshape = c()
Eshape = c()
CS = as.numeric()
CN = as.numeric()
Jc = as.numeric()
ES = as.numeric()
EN = as.numeric()
Je = as.numeric()
taxon = c()
etype = c()
m2 = as.numeric()
# comparison variables
BCJ = as.numeric()
BCrad = as.numeric()
BCS = as.numeric()
BCN = as.numeric()
BCcomp = as.numeric()
percS = as.numeric()
percN = as.numeric()

c = as.numeric()
e = as.numeric()
r2 = as.numeric()

compc = as.numeric()
compe = as.numeric()
compr2 = as.numeric()

for (iRow in 1:nrow(comps)){
  control = comps[iRow,2]  #find control in pair
  experiment = comps[iRow,3]  # find experiment in pair
  taxa = as.character(expers[which(expers[,2]==control),7])# find taxonomic group from experiments table
  type = as.character(expers[which(expers[,2]==control),4]) # find experiment type from experiments table
  extent = sites[which(sites[,2]==control),11]
  ref = as.character(comps[iRow,1])
  # Check that < 10% of individuals are unidentified. If meets criteria, continue
  if (percent_unidSpp(control, comms) == "OK" & percent_unidSpp(experiment, comms) == "OK"){
    a1 = sort(as.numeric(comms[which(comms[,2] == control & comms[,7] != 0), 8])) #vector of control abundances
    a2 = sort(as.numeric(comms[which(comms[,2] == experiment & comms[,7] != 0), 8])) #vector of exp abundances
    # Check that there are at least 5 species and 30 individuals in each community, If yes, proceed.
    if (length(a1) > 4 & length(a2) > 4 & sum(a1) > 29 & sum(a2) > 29){
      # record all values in a comparison matrix
      relcon = relabund(a1) #make the lengths the same!
      relexp = relabund(a2)  #make the lengths the same!
      comparison_matrix = abundMerge(relcon, relexp)
      c = append(c, comparison_matrix[,1])
      e = append(e, comparison_matrix[,2])
      r2 = append(r2, rsquare(comparison_matrix[,1], comparison_matrix[,2]))
      # record all composition-specific values in a comparison matrix
      comparison = subset(comms[which(comms$siteID == control | comms$siteID == experiment),])
      comparison = reshape_data(comparison) #table species & abundance in paired communities 
      comparison = comparison[,c(2:ncol(comparison))]
      comparison[1,] = comparison[1,]/sum(comparison[1,]) #convert to relabundance
      comparison[2,] = comparison[2,]/sum(comparison[2,]) #convert to relabundance
      compc = append(compc, as.numeric(comparison[1,]))
      compe = append(compe, as.numeric(comparison[2,]))
      compr2 = append(compr2, rsquare(as.numeric(comparison[1,]), as.numeric(comparison[2,])))
      # find categorical shapes (logseries vs. lognormal)
      if(expers[which(expers[,2]==control),10] == 1) { #is it raw abundance data?
        d = dist.test(a1, a2)
        Cshape = append(Cshape, d$con)
        Eshape = append(Eshape, d$exp)
      }
      else {      #if mean abundance, can't get the data, (needs integers)
        Cshape = append(Cshape, "ERROR")
        Eshape = append(Eshape, "ERROR")
      }
      #plot the compared data
      RAD_plot(control, experiment, a1, a2, taxa)
      # descriptors
      con_s = length(a1)
      con_n = sum(a1)
      con_j = SimpE(comms[which(comms[,2]==control),])
      exp_s = length(a2)
      exp_n = sum(a2)
      exp_j = SimpE(comms[which(comms[,2]==experiment),])
      # get summary statistics from comparisons
      BCJ = append(BCJ, BCdist(matrix(c(con_j, exp_j), nrow = 1, ncol = 2)))
      BCrad = append(BCrad, BCdist(abundMerge(relabund(a1), relabund(a2))))
      BCS = append(BCS, BCdist(matrix(c(con_s, exp_s), nrow = 1, ncol = 2))) 
      BCN = append(BCN, BCdist(matrix(c(con_n, exp_n), nrow = 1, ncol = 2)))
      BCcomp = append(BCcomp, BCdist(subset(comms[which(comms$siteID == control | comms$siteID == experiment),])))
      percS = append(percS, ((exp_s - con_s)/con_s)*100) 
      percN = append(percN, ((exp_n - con_n)/con_n)*100)
      taxon = append(taxon, taxa)
      etype = append(etype, type)
      # record summary descriptive variables
      refID = append(refID, ref)
      cID = append(cID, control)
      eID = append(eID, experiment)
      CS = append(CS, con_s)
      CN = append(CN, con_n)
      Jc = append(Jc, round(con_j,4))
      ES = append(ES, exp_s)
      EN = append(EN, exp_n)
      Je = append(Je, round(exp_j,4))
      m2 = append(m2, extent)
    }}}

dev.off()


#collapse taxon types into broader categories so there aren't so many factors
taxon[taxon=='carabid']<-'insect'
taxon[taxon=='lepidopteran']<- 'insect'
taxon[taxon=='odonate']<- 'insect'
taxon[taxon=='orthoptera']<-'insect'
taxon[taxon=='orthoptera ']<-'insect'
taxon[taxon=='beetle']<-'insect'
taxon[taxon=='microarthropods']<-'microarthropod'
taxon[taxon=='reptile']<-'herpetofauna'


#------------------------------------- 
#                 results
#-------------------------------------
#root mean squared error for the variables. Usually used as standard deviation of model prediction error, but can be
# used as an indicator of the degree of change between control (obs) and the experiment (sim)
comp_rmse = round(rmse(compe, compc),4)
n_rmse = round(rmse(EN, CN),4)
s_rmse = round(rmse(ES, CS),4)
relabun_rmse = round(rmse(e,c),4)
j_rmse = round(rmse(Je, Jc),4)

comp_r2 = round(rsquare(compc, compe),4)
n_r2 = round(rsquare(CN,EN),4)
s_r2 = round(rsquare(CS,ES),4)
relabun_r2 = round(rsquare(c,e),4)
j_r2 = round(rsquare(Jc,Je),4)

#count the communities displaying various shapes. This does take into account duplicates. Should be correct
shapes = count_RAD_shapes(cID, eID, Cshape, Eshape)


#fiddle plots to see if small-scale sites pick up more variability
par(mfrow=c(3,2))

plot(BCcomp, m2, pch=19)
plot(abs(percS), m2, pch=19, xlim=c(0,150))
plot(abs(percN), m2, pch=19, xlim=c(0,300))
plot(BCJ, m2, pch=19)
plot(BCrad, m2, pch=19)


#------null modeling results
BCnull = c()
Nnull = c()
Snull = c()
Jnull = c()
RADnull = c()

for (iRow in 1:nrow(comps)){
  control = comps[iRow,2]  #find control in pair
  experiment = comps[iRow,3]  # find experiment in pair
  # Check that < 10% of individuals are unidentified. If meets criteria, continue
  if (percent_unidSpp(control, comms) == "OK" & percent_unidSpp(experiment, comms) == "OK"){
    a1 = sort(as.numeric(comms[which(comms[,2] == control & comms[,7] != 0), 8])) #vector of control abundances
    a2 = sort(as.numeric(comms[which(comms[,2] == experiment & comms[,7] != 0), 8])) #vector of exp abundances
    # Check that there are at least 5 species and 30 individuals in each community, If yes, proceed.
    if (length(a1) > 4 & length(a2) > 4 & sum(a1) > 29 & sum(a2) > 29){
      comparison = subset(comms[which(comms$siteID == control | comms$siteID == experiment),])
      comparison$spname = paste(comparison$genus, comparison$species, sep="_")
      comparison$grp = as.numeric(comparison$siteID==control) #control site == 1, make sure to order comparisons, control[1,], experiment[2,]
      data = comparison[,c(2,8,9,10)]
      sitexsp <- cast(data, grp + siteID  ~ spname, value = "abundance", fun = mean)
        sitexsp[is.na(sitexsp)] <- 0
        sitexsp = arrange(sitexsp, desc(grp))
      BCresult = NullCommunityBC(sitexsp[,-c(1,2)], 100)
        BCnull = append(BCnull, BCresult)
      Nresult = NullCommunityN(sitexsp[,-c(1,2)], 100)
        Nnull = append(Nnull, Nresult)
      Sresult = NullCommunityS(sitexsp[,-c(1,2)], 100)
        Snull = append(Snull, Sresult)
      Jresult = NullCommunityJ(sitexsp[,-c(1,2)], 100)
        Jnull = append(Jnull, Jresult)
    }}}

Nindex = Nnull=="Sign." 
Sindex = Snull=="Sign."
BCindex = BCnull=="Sign."

#------------------------------------------------------------- 
#                 map site locations, color coded by taxa
#-------------------------------------------------------------
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
#Get world map info
world_map <- map_data("world")

#Creat a base plot
p <- ggplot() + coord_fixed()

#Add map to base plot
base_world <- p + geom_polygon(data=world_map,aes(x=long, y=lat,group=group), fill = "white", col = "black") + theme_bw()

#add points to map
site_map <- base_world + geom_point(data=data, aes(x=lon, y=lat, col = taxa), alpha = 0.5, cex = 5) +
  theme(legend.position = "right") + element_blank()
