# All files should be organized as such:
# RAD: SiteID | experiment/control | year | genus | species | id2species | abundance
# Code will compare S, N and composition among treatments at each site,
# output and compare RADs, parameters for each treatment at each site.

library(ggplot2)
library(gridExtra)

#---------------------------------------------------------------------------------
#          setup - select wd, import data, source code,  file to collect results
#---------------------------------------------------------------------------------

#"wd = "/Users/sarah/Documents/GitHub/experimental-rads/"
wd = "C:\\Users\\sarah\\Documents\\GitHub\\experimental-rads"
setwd(wd)

source("ExpRADsFunctions.R")   #Run the code containing the functions

comms = read.csv("data/community_analysis_data.csv", na.strings = 'NULL')
#comps = read.csv("comparison_analysis_data.csv") #this file is unordered. Looks less nice when plotted.
comps = read.csv("data/orderedcomparisons.csv")
  names(comps)<-c('ref', 'controID','expID')
expers = read.csv("data/experiments_analysis_data.csv")
sites = read.csv("data/sites_analysis_data.csv")

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
ranklr = as.numeric()

compc = as.numeric()
compe = as.numeric()
compr2 = as.numeric()
complr = as.numeric()

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
        tr = as.data.frame(t(comparison_matrix))
      rlr = sapply(tr, function(x) LogRatio(x[2],x[1]) )
      c = append(c, comparison_matrix[,1])
      e = append(e, comparison_matrix[,2])
      r2 = append(r2, rsquare(comparison_matrix[,1], comparison_matrix[,2]))
      ranklr = append(ranklr, mean(rlr, na.rm=TRUE))
      
      # record all composition-specific values in a comparison matrix
      comparison = subset(comms[which(comms$siteID == control | comms$siteID == experiment),])
      comparison = reshape_data(comparison) #table species & abundance in paired communities 
      comparison = comparison[,c(2:ncol(comparison))]
      comparison[1,] = comparison[1,]/sum(comparison[1,]) #convert to relabundance
      comparison[2,] = comparison[2,]/sum(comparison[2,]) #convert to relabundance
        lr = sapply(comparison, function(x) LogRatio(x[2], x[1]) )
        complr = append(complr, mean(lr, na.rm=TRUE))
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


#log ratio differences between values
s_df = as.data.frame(t(cbind(CS,ES)))
n_df = as.data.frame(t(cbind(CN,EN)))
e_df = as.data.frame(t(cbind(Jc, Je)))

Slr = sapply(s_df, function(x) LogRatio(x[2], x[1]) )
Nlr = sapply(n_df, function(x) LogRatio(x[2], x[1]) )
Elr = sapply(e_df, function(x) LogRatio(x[2], x[1]) )

#count the communities displaying various shapes. This does take into account duplicates. Should be correct
shapes = count_RAD_shapes(cID, eID, Cshape, Eshape)


#fiddle plots to see if small-scale sites pick up more variability
par(mfrow=c(3,2))

plot(BCcomp, m2, pch=19)
plot(abs(percS), m2, pch=19, xlim=c(0,150))
plot(abs(percN), m2, pch=19, xlim=c(0,300))
plot(BCJ, m2, pch=19)
plot(BCrad, m2, pch=19)

#put the results together for later plotting and comparison
diversity = data.frame(taxa, etype, CS, ES, CN, EN, Je, Jc)
composition = data.frame(compc, compe)
relabundance = data.frame(c,e)

#----------------------------------------------------------------------- 
#                 FIGURE 1. map site locations, color coded by taxa
#-----------------------------------------------------------------------
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
  data = cbind(data, taxa = as.factor(taxa))

#Get world map info
mapWorld <- borders("world", fill="gray65", col="gray65")

base_world = ggplot() + mapWorld + theme_bw()

site_map = base_world + 
  geom_point(data=data, aes(x=lon, y=lat, col = taxa, fille = taxa, shape = taxa), cex = 5) +
  theme(legend.position = "right") + element_blank() + xlab("Longitude") + ylab("Latitude") +
  scale_colour_grey(start=0, end=0.3) + scale_shape_discrete(solid=F) + 
  scale_x_continuous(breaks = seq(-180, 180, by=30)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 15)) +
  theme(text = element_text(size=20))
site_map


#----------------------------------------------------------------------- 
#                 FIGURE 2. compare control and manipulated data
#-----------------------------------------------------------------------
#plot results along 1:1 line

#Fig 2A - plots for composition
compchange = ggplot(data=composition, aes(compc, compe)) + geom_point(alpha=0.5, size=3) + 
  xlab("species relative abundance") + ylab("species relative abundance") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_bw() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("A")

#Fig 2B - plots for abundance
abunchange = ggplot(data=diversity, aes(CN, EN)) + geom_point(alpha=0.5, size=3) + 
  xlab("total abundance") + ylab("total abundance") + 
  scale_x_log10(breaks = c(30, 100, 500, 2500, 6500), limits = c(30,6500)) +
  scale_y_log10(breaks = c(30, 100, 500, 2500, 6500), limits = c(30,6500)) + theme_bw() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("B")

# Fig 2C - plots for S
schange = ggplot(data=diversity, aes(CS, ES)) + geom_point(alpha=0.5, size=3) + 
  xlab("species richness") + ylab("species richness") + 
  scale_x_log10(breaks = c(5, 10, 25, 50, 100, 200), limits = c(5, 200)) +
  scale_y_log10(breaks = c(5, 10, 25, 50, 100, 200), limits = c(5, 200)) + theme_bw() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("C")

# Fig 1D - plots for evenness
evenchange = ggplot(data=diversity, aes(Jc, Je)) + geom_point(alpha=0.5, size=3) + 
  xlab("Simpson's evenness") + ylab("Simpson's evenness") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_bw() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("D")

# Fig 1E.  compare relative abundance at each rank in all paired sites
rankabunchange = ggplot(data=relabundance, aes(c, e)) + geom_point(alpha=0.5, size=3) + 
  xlab("rank relative abundance") + ylab("rank relative abundance") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_bw() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("E")
  
grid.arrange(compchange, abunchange, schange, evenchange, rankabunchange, nrow=2)


#------------------------------------------------------------------------------------------- 
#                 FIGURE 3. histograms of log-ratio difference in treatment vs. controls
#-------------------------------------------------------------------------------------------

logratio = data.frame()



