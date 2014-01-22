# All files should be organized as such:
# RAD: SiteID | experiment/control | year | genus | species | id2species | abundance
# Code will compare S, N and composition among treatments at each site,
# output and compare RADs, parameters for each treatment at each site.

library(ggplot2)
library(gridExtra)
library(GGally)

#---------------------------------------------------------------------------------
#          setup - select wd, import data, source code,  file to collect results
#---------------------------------------------------------------------------------

wd = "/Users/sarah/Documents/GitHub/experimental-rads/"
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

#build up dataframe
desc = data.frame(cID = 1, eID = 1, CS = 1, ES = 1, CN = 1, EN = 1, Jc = 1, Je = 1, m2 = 1)
charvars = data.frame(refID = NA, Cshape = NA, Eshape = NA, taxon = NA, etype = NA)
compvals = data.frame(BCcomp = 1, BCN = 1, BCS = 1, BCJ = 1, BCrad = 1, percS = 1, percN = 1)
rankvals = data.frame(r2 = 1, ranklr = 1, sdrlr=1)
compositionvals = data.frame(compr2 = 1, complr = 1, sdlr = 1)
outcount = 1
c = NULL
e = NULL
compc = NULL
compe = NULL
complrvals = NULL
ranklrvals = NULL

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
      # record all values in a comparison matrix - to compare relative abundance at each rank
      a1 = sort(a1, decreasing = TRUE)
      a2 = sort(a2, decreasing = TRUE)
      relcon = relabund(a1) 
      relexp = relabund(a2)  
      comparison_matrix = abundMerge(relcon, relexp) #makes the lengths the same
        tr = as.data.frame(t(comparison_matrix))  #row 1 is control, row 2 is experiment
      rlr = sapply(tr, function(x) LogRatio(x[2],x[1]) )
      ranklrvals = append(ranklrvals, rlr)
      rankvals[outcount,] = c(rsquare(comparison_matrix[,1], comparison_matrix[,2]), mean(rlr), sd(rlr))
      c = append(c, comparison_matrix[,1])
      e = append(e, comparison_matrix[,2])
      
      # record all composition-specific values in a comparison matrix, reshape into siteXspecies matrix
      comparison = reshape_data(subset(comms[which(comms$siteID == control | comms$siteID == experiment),]))
      #make sure dataframe is ordered control (row 1), experiment (row 2)
      cntrl = comparison[which(comparison$siteID == control),c(2:ncol(comparison))]
        cntrl = cntrl/sum(cntrl)
      exprm = comparison[which(comparison$siteID == experiment),c(2:ncol(comparison))]
        exprm = exprm/sum(exprm)
      comparison = rbind(cntrl, exprm) #make sure control is row 1
      #take the log ratio of raw abundance
        lr = sapply(comparison, function(x) LogRatio(x[2], x[1]) )
      complrvals = append(complrvals, lr)
      compositionvals[outcount,] = c(rsquare(as.numeric(comparison[1,]), as.numeric(comparison[2,])), mean(lr), sd(lr))
      compc = append(compc, as.numeric(comparison[1,]))
      compe = append(compe, as.numeric(comparison[2,]))
     
      # find categorical shapes (logseries vs. lognormal)
      if(expers[which(expers[,2]==control),10] == 1) { #is it raw abundance data?
        d = dist.test(a1, a2)
        charvars[outcount,] = c(ref, d$con, d$exp, taxa, type)
      }
      else {      #if mean abundance, can't get the data, (needs integers)
        charvars[outcount,] = c(ref, "ERROR", "ERROR", taxa, type)
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
      
      BCJ = BCdist(matrix(c(con_j, exp_j), nrow = 1, ncol = 2))
      BCrad = BCdist(abundMerge(relabund(a1), relabund(a2)))
      BCS = BCdist(matrix(c(con_s, exp_s), nrow = 1, ncol = 2)) 
      BCN = BCdist(matrix(c(con_n, exp_n), nrow = 1, ncol = 2))
      BCcomp = BCdist(subset(comms[which(comms$siteID == control | comms$siteID == experiment),]))
      percS = ((exp_s - con_s)/con_s)*100
      percN = ((exp_n - con_n)/con_n)*100

      # record summary descriptive variables
      desc[outcount,] = c(control, experiment, con_s, exp_s, con_n, exp_n, con_j, exp_j, extent)
      # get summary statistics from comparisons
      compvals[outcount,] = c(BCcomp, BCN, BCS, BCJ, BCrad, percS, percN)
      outcount = outcount + 1
    }}}

dev.off()


#collapse taxon types into broader categories so there aren't so many factors
charvars$taxon[charvars$taxon=='carabid']<-'insect'
charvars$taxon[charvars$taxon=='lepidopteran']<- 'insect'
charvars$taxon[charvars$taxon=='odonate']<- 'insect'
charvars$taxon[charvars$taxon=='orthoptera']<-'insect'
charvars$taxon[charvars$taxon=='orthoptera ']<-'insect'
charvars$taxon[charvars$taxon=='beetle']<-'insect'
charvars$taxon[charvars$taxon=='microarthropods']<-'microarthropod'
charvars$taxon[charvars$taxon=='reptile']<-'herpetofauna'

results = cbind(charvars, desc, compvals, rankvals, compositionvals)
#------------------------------------- 
#                 results
#-------------------------------------
#root mean squared error for the variables. Usually used as standard deviation of model prediction error, but can be
# used as an indicator of the degree of change between control (obs) and the experiment (sim)
comp_rmse = round(rmse(compe, compc),4)
n_rmse = round(rmse(results$EN, results$CN),4)
s_rmse = round(rmse(results$ES, results$CS),4)
relabun_rmse = round(rmse(e,c),4)
j_rmse = round(rmse(results$Je, results$Jc),4)

comp_r2 = round(rsquare(compc, compe),4)
n_r2 = round(rsquare(results$CN,results$EN),4)
s_r2 = round(rsquare(results$CS,results$ES),4)
relabun_r2 = round(rsquare(c,e),4)
j_r2 = round(rsquare(results$Jc,results$Je),4)

#log ratio differences between values
s_df = as.data.frame(t(cbind(results$CS, results$ES)))
n_df = as.data.frame(t(cbind(results$CN, results$EN)))
e_df = as.data.frame(t(cbind(results$Jc, results$Je)))

Slr = sapply(s_df, function(x) LogRatio(x[2], x[1]) )
Nlr = sapply(n_df, function(x) LogRatio(x[2], x[1]) )
Elr = sapply(e_df, function(x) LogRatio(x[2], x[1]) )

#count the communities displaying various shapes. This does take into account duplicates. Should be correct
shapes = count_RAD_shapes(results$cID, results$eID, results$Cshape, results$Eshape)

#fiddle plots to see if small-scale sites pick up more variability
par(mfrow=c(3,2))

plot(results$BCcomp, desc$m2, pch=19)
plot(abs(results$percS), desc$m2, pch=19, xlim=c(0,150))
plot(abs(results$percN), desc$m2, pch=19, xlim=c(0,300))
plot(results$BCJ, desc$m2, pch=19)
plot(results$BCrad, desc$m2, pch=19)

#put the results together for later plotting and comparison
diversity = cbind(charvars[,c(4:5)], desc)
composition = data.frame(compc, compe)
relabundance = data.frame(c,e)

taxa = diversity$taxon
etype = diversity$etype
complr = compositionvals$complr
ranklr = rankvals$ranklr
lograt = data.frame(complr, ranklr, Nlr, Slr, Elr, taxa, etype)

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

ggsave(site_map, file = "site_map.jpeg", dpi = 300, width = 9, height = 4.5)

#----------------------------------------------------------------------- 
#                 FIGURE 2. compare control and manipulated data
#-----------------------------------------------------------------------
#plot results along 1:1 line

#Fig 2A - plots for composition
compchange = ggplot(data=composition, aes(compc, compe)) + geom_point(alpha=0.5, size=3) + 
  xlab("species relative abundance") + ylab("species relative abundance") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("A")

#Fig 2B - plots for abundance
abunchange = ggplot(data=diversity, aes(CN, EN)) + geom_point(alpha=0.5, size=3) + 
  xlab("total abundance") + ylab("total abundance") + 
  scale_x_log10(breaks = c(30, 100, 500, 2500, 6500), limits = c(30,6500)) +
  scale_y_log10(breaks = c(30, 100, 500, 2500, 6500), limits = c(30,6500)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("B")

# Fig 2C - plots for S
schange = ggplot(data=diversity, aes(CS, ES)) + geom_point(alpha=0.5, size=3) + 
  xlab("species richness") + ylab("species richness") + 
  scale_x_log10(breaks = c(5, 10, 25, 50, 100, 200), limits = c(5, 200)) +
  scale_y_log10(breaks = c(5, 10, 25, 50, 100, 200), limits = c(5, 200)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("C")

# Fig 1D - plots for evenness
evenchange = ggplot(data=diversity, aes(Jc, Je)) + geom_point(alpha=0.5, size=3) + 
  xlab("Simpson's evenness") + ylab("Simpson's evenness") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("D")

# Fig 1E.  compare relative abundance at each rank in all paired sites
rankabunchange = ggplot(data=relabundance, aes(c, e)) + geom_point(alpha=0.5, size=3) + 
  xlab("rank relative abundance") + ylab("rank relative abundance") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("E")
  
  
grid.arrange(compchange, abunchange, schange, evenchange, rankabunchange, nrow=2)


#------------------------------------------------------------------------------------------- 
#                 FIGURE 3. histograms of log-ratio difference in treatment vs. controls
#-------------------------------------------------------------------------------------------
#plot histograms of the log-ratio results

#Fig 3A - plots for composition
comphist = ggplot(data=lograt, aes(complr)) + geom_histogram() + 
  xlab("mean log-ratio of species relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-1.5,1.5, by=0.5), limits = c(-1.5,1.5)) + theme_classic() +  
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("A")

nhist = ggplot(data=lograt, aes(Nlr)) + geom_histogram() + 
  xlab("log-ratio of total abundance") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-1.5,1.5, by=0.5), limits = c(-1.5,1.5)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("B")

shist = ggplot(data=lograt, aes(Slr)) + geom_histogram() + 
  xlab("log-ratio of species richness") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-1.5,1.5, by=0.5), limits = c(-1.5,1.5)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("C")

ehist = ggplot(data=lograt, aes(Elr)) + geom_histogram() + 
  xlab("log-ratio of Simpson's evenness") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-1.5,1.5, by=0.5), limits = c(-1.5,1.5)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("D")

rankhist = ggplot(data=lograt, aes(ranklr)) + geom_histogram() + 
  xlab("mean log-ratio of rank relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-1.5,1.5, by=0.5), limits = c(-1.5,1.5)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("A")

grid.arrange(comphist, nhist, shist, ehist, rankhist, nrow=2)


#------------------------------------------------------------------------------------------- 
#                 plot the standard deviations of the log-ratios for composition and rank
#                   and the raw log-ratio values from them - mean is maybe not representative?
#-------------------------------------------------------------------------------------------
compossd = ggplot(data=results, aes(sdlr)) + geom_histogram() + 
  xlab("sd of mean log-ratio of composition relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,6, by=0.5), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,20, by=5), limits = c(0,20)) +
  theme(text = element_text(size=20)) + ggtitle("A")

ranksd = ggplot(data=results, aes(sdrlr)) + geom_histogram() + 
  xlab("sd of mean log-ratio of rank relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,6, by=0.5), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,20, by=5), limits = c(0,20)) +
  theme(text = element_text(size=20)) + ggtitle("B")

grid.arrange(compossd, ranksd, nrow = 1)

# raw vals:
composrawlr = ggplot(data=data.frame(complrvals), aes(complrvals)) + geom_histogram(binwidth=0.5) + 
  xlab("all log-ratio of composition relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-9,9, by=1), limits = c(-9,9)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,1000, by=100), limits = c(0,1000)) +
  theme(text = element_text(size=20)) + ggtitle("A")

rankrawlr = ggplot(data=data.frame(ranklrvals), aes(ranklrvals)) + geom_histogram(binwidth=0.5) + 
  xlab("all log-ratio of rank relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-9,9, by=1), limits = c(-9,9)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,1000, by=100), limits = c(0,1000)) +
  theme(text = element_text(size=20)) + ggtitle("B")

grid.arrange(composrawlr, rankrawlr, nrow = 1)


#------------------------------------------------------------------------------------------- 
#                 FIGURE 4. pairs plots of the log-ratios for the 5 variables 
#-------------------------------------------------------------------------------------------
logratios = cbind(BC=results$BCcomp, lograt[,c(1,3,4,5,2,6)])

ggpairs(logratios, colour = "taxa")

logratios = cbind(BC=results$BCcomp, lograt[,c(1,3,4,5,2)], results$sdlr, results$sdrlr)
pairs(logratios, pch = 19)
