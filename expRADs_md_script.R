# All files should be organized as such:
# RAD: SiteID | experiment/control | year | genus | species | id2species | abundance
# Code will compare S, N, evenness, rank, and composition among treatments at each site,
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

# import the data tables
comms = read.csv("data/community_analysis_data.csv", na.strings = 'NULL')
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
rankvals_p1 = data.frame(r2_p1 = 1, ranklr_p1 = 1, sdrlr_p1 = 1, absranklr_p1 = 1, abssdrlr_p1 = 1)
rankvals_p01 = data.frame(r2 = 1, ranklr = 1, sdrlr = 1, absranklr = 1, abssdrlr =1)
compositionvals_p1 = data.frame(compr2_p1 = 1, complr_p1 = 1, sdlr_p1 = 1, abscomplr_p1 = 1, abssdlr_p1 = 1, abscomplrzero_p1 = 1)
compositionvals_p01 = data.frame(compr2 = 1, complr = 1, sdlr = 1, abscomplr = 1, abssdlr = 1, abscomplrzero = 1)
outcount = 1
c_p01 = NULL
e_p01 = NULL
c_p1 = NULL
e_p1 = NULL
c = NULL
e = NULL
compc_p1 = NULL
compe_p1 = NULL
compc_p01 = NULL
compe_p01 = NULL
compc = NULL
compe = NULL
complrvals_p01 = NULL
ranklrvals_p01 = NULL
complrvals_p1 = NULL
ranklrvals_p1 = NULL

for (iRow in 1:nrow(comps)){
  control = comps[iRow,2]  #find control in pair
  experiment = comps[iRow,3]  # find experiment in pair
  taxa = as.character(expers[which(expers[,2]==control),7])# find taxonomic group from experiments table
  type = as.character(expers[which(expers[,2]==control),4]) # find experiment type from experiments table
  extent = sites[which(sites[,2]==control),11]
  ref = as.character(comps[iRow,1])
  
  # Check that < 10% of individuals are unidentified. If meets criteria, continue
  if (percent_unidSpp(control, comms) == "OK" & percent_unidSpp(experiment, comms) == "OK"){
    a1 = sort(as.numeric(comms[which(comms[,2] == control & comms[,7] != 0), 8]), decreasing = TRUE) #control abundances
    a2 = sort(as.numeric(comms[which(comms[,2] == experiment & comms[,7] != 0), 8]), decreasing = TRUE) #exp abundances
    
    # Check that there are at least 5 species and 30 individuals in each community, If yes, proceed.
    if (length(a1) > 4 & length(a2) > 4 & sum(a1) > 29 & sum(a2) > 29){
      # Sort the abundances
      #record all values in a comparison matrix - to compare relative abundance at each RANK
      comparison_matrix = abundMerge(a1, a2) #makes the lengths the same, fills with zeroes
      
      #add a constant 1 to all values before taking relative abundance at each RANK
      comparison_p1 = comparison_matrix + 1
      relcon = relabund(comparison_p1[,1]) 
      relexp = relabund(comparison_p1[,2])  
      tr = data.frame(rbind(relcon,relexp),row.names=NULL) #row 1 is control, row 2 is experiment
      rlr = sapply(tr, function(x) LogRatio(x[2],x[1]) )
      ranklrvals_p1 = append(ranklrvals_p1, rlr)
      rankvals_p1[outcount,] = c(rsquare(as.numeric(tr[1,]), as.numeric(tr[2,])), median(rlr), sd(rlr), median(abs(rlr)), sd(abs(rlr)))
      c_p1 = append(c_p1, relcon)
      e_p1 = append(e_p1, relexp)
      
      #add a constant 0.01 to all values before taking relative abundance at each RANK
      comparison_p01 = comparison_matrix + 0.01
      relcon = relabund(comparison_p01[,1]) 
      relexp = relabund(comparison_p01[,2])  
      tr = data.frame(rbind(relcon,relexp),row.names=NULL) #row 1 is control, row 2 is experiment
      rlr = sapply(tr, function(x) LogRatio(x[2],x[1]) )
      ranklrvals_p01 = append(ranklrvals_p01, as.numeric(rlr))
      rankvals_p01[outcount,] = c(rsquare(as.numeric(tr[1,]), as.numeric(tr[2,])), median(rlr), sd(rlr), median(abs(rlr)), sd(abs(rlr)))
      c_p01 = append(c_p01, relcon)
      e_p01 = append(e_p01, relexp)
      
      #store actual relative abundance values
      realcon = relabund(comparison_matrix[,1])
      realexp = relabund(comparison_matrix[,2])
      c = append(c, realcon)
      e = append(e, realexp)
      
      # record all COMPOSITION-specific values in a comparison matrix, reshape into siteXspecies matrix
      comparison = reshape_data(subset(comms[which(comms$siteID == control | comms$siteID == experiment),])) 
      #dataframe is ordered control (r1), experiment (r2)
      
      #add a constant 1, to account for zeroes and take the relative abundance
      cntrl_p1 = comparison[which(comparison$siteID == control),c(2:ncol(comparison))] + 1 
        cntrl_p1 = cntrl_p1/sum(cntrl_p1)
      exprm_p1 = comparison[which(comparison$siteID == experiment),c(2:ncol(comparison))] + 1 
        exprm_p1 = exprm_p1/sum(exprm_p1)
      
      #add a constant 0.01, to account for zeroes and take the relative abundance
      cntrl_p01 = comparison[which(comparison$siteID == control),c(2:ncol(comparison))] + 0.01 
        cntrl_p01 = cntrl_p01/sum(cntrl_p01)
      exprm_p01 = comparison[which(comparison$siteID == experiment),c(2:ncol(comparison))] + 0.01 
        exprm_p01 = exprm_p01/sum(exprm_p01)
      
      #Keep the zeroes and take the relative abundance, to compare only species in common
      cntrl_nz = comparison[which(comparison$siteID == control),c(2:ncol(comparison))] 
        cntrl_nz = cntrl_nz/sum(cntrl_nz)
      exprm_nz = comparison[which(comparison$siteID == experiment),c(2:ncol(comparison))]
        exprm_nz = exprm_nz/sum(exprm_nz)
      
      comparison_p1 = rbind(cntrl_p1, exprm_p1) 
      comparison_p01 = rbind(cntrl_p01, exprm_p01)
      comparison_nz = rbind(cntrl_nz, exprm_nz)
      
      #take the log ratio of raw abundance and record the results
      lr_p1 = sapply(comparison_p1, function(x) LogRatio(x[2], x[1]) )
      lr_p01 = sapply(comparison_p01, function(x) LogRatio(x[2], x[1]) )     
      lr_nz = sapply(comparison_nz, function(x) LogRatio_noZero(x[2], x[1]))
      
      complrvals_p1 = append(complrvals_p1, as.numeric(lr_p1))
      compositionvals_p1[outcount,] = c(rsquare(as.numeric(comparison_p1[1,]), as.numeric(comparison_p1[2,])), median(lr_p1), sd(lr_p1), median(abs(lr_p1)), sd(abs(lr_p1)), median(abs(lr_nz),na.rm=TRUE))
      compc_p1 = append(compc_p1, as.numeric(comparison_p1[1,]))
      compe_p1 = append(compe_p1, as.numeric(comparison_p1[2,]))
      
      complrvals_p01 = append(complrvals_p01, as.numeric(lr_p01))
      compositionvals_p01[outcount,] = c(rsquare(as.numeric(comparison_p01[1,]), as.numeric(comparison_p01[2,])), median(lr_p01), sd(lr_p01), median(abs(lr_p01)), sd(abs(lr_p01)), median(abs(lr_nz),na.rm=TRUE))
      compc_p01 = append(compc_p01, as.numeric(comparison_p01[1,]))
      compe_p01 = append(compe_p01, as.numeric(comparison_p01[2,]))
     
      compc = append(compc, cntrl_nz)
      compe = append(compe, exprm_nz)
      
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

results = cbind(charvars, desc, compvals, rankvals_p01, compositionvals_p01, rankvals_p1, compositionvals_p1)

# calculate the log ratio differences between values
s_df = as.data.frame(t(cbind(results$CS, results$ES)))
n_df = as.data.frame(t(cbind(results$CN, results$EN)))
e_df = as.data.frame(t(cbind(results$Jc, results$Je)))

Slr = sapply(s_df, function(x) LogRatio(x[2], x[1]) )
Nlr = sapply(n_df, function(x) LogRatio(x[2], x[1]) )
Elr = sapply(e_df, function(x) LogRatio(x[2], x[1]) )

#put the results together for later plotting and comparison
diversity = cbind(charvars[,c(4:5)], desc)
#composition relative abundance
composition = data.frame(compc_p01, compe_p01, compc_p1, compe_p1)
#rank relative abundance
relabundance = data.frame(c_p01, e_p01, c_p1, e_p1)

taxa = diversity$taxon
etype = diversity$etype
complr_p01 = compositionvals_p01$complr
abscomplr_p01 = compositionvals_p01$abscomplr
ranklr_p01 = rankvals_p01$ranklr
absranklr_p01 = rankvals_p01$absranklr
complr_p1 = compositionvals_p1$complr_p1
abscomplr_p1 = compositionvals_p1$abscomplr_p1
ranklr_p1 = rankvals_p1$ranklr_p1
absranklr_p1 = rankvals_p1$absranklr_p1
lograt = data.frame(complr_p01, ranklr_p01, abscomplr_p01, absranklr_p01, Nlr, Slr, Elr, taxa, etype, 
                    complr_p1, ranklr_p1, abscomplr_p1, absranklr_p1)

# for pairs plot
logratios = lograt[,c(3,5,6,7,4,12,13)]
logratios$Nlr = abs(logratios$Nlr)
logratios$Slr = abs(logratios$Slr)
logratios$Elr = abs(logratios$Elr)
names(logratios) = c("composition", "abundance", "richness", "evenness", "rank", "composition_p1", "rank_p1")


#------------------------------------- 
#                 results
#-------------------------------------

# Find the R2 values for the paired data
comp_r2 = round(rsquare(compc, compe),4)
comp_p01_r2 = round(rsquare(compc_p01, compe_p01),4)
comp_p1_r2 = round(rsquare(compc_p1, compe_p1),4)
n_r2 = round(rsquare(results$CN,results$EN),4)
s_r2 = round(rsquare(results$CS,results$ES),4)
rank_r2 = round(rsquare(c,e),4)
rank_p01_r2 = round(rsquare(c_p01,e_p01),4)
rank_p1_r2 = round(rsquare(c_p1,e_p1),4)
j_r2 = round(rsquare(results$Jc,results$Je),4)

# count the communities displaying various shapes. This does take into account duplicates. Should be correct
shapes = count_RAD_shapes(results$cID, results$eID, results$Cshape, results$Eshape)

# Print the range of values
for (col in 1:ncol(logratios)){
  print(paste(names(logratios[col]), ":", min(logratios[,col]), "-", max(logratios[,col]), sep = " "))
}


# For TABLE C2
# Print the upper and lower quantiles mean, median and standard deviation of the absolute log ratio
for (col in 1:ncol(logratios)){
  mn = round(mean(logratios[,col]),3)
  m = round(median(logratios[,col]),3)
  stdev = round(sd(logratios[,col]),3)
  q = round(quantile(logratios[,col], c(0.05, 0.95)),3)
  print(paste(names(logratios[col]), ": mean = ", mn, ", median = ", m, ", sd = ", stdev, sep=""))
  print(paste(names(logratios[col]), ": the 5% and 95% quantiles are:", q[[1]], "and", q[[2]]))
}

#the mean and median value when composition includes only species that were present in both communities
mn = round(mean(compositionvals_p01$abscomplrzero),3)
md = round(median(compositionvals_p01$abscomplrzero),3)
std = round(sd(compositionvals_p01$abscomplrzero),3)
print(paste("when only consider species that were present in both communities, composition log ratio is median = ",
            md, "sd =", std, "and mean =", mn))

# mean rank log ratio for communities where species richness change was zero
rankzero = logratios[which(logratios$richness == 0),]
mn = round(mean(rankzero$rank),3)
md = round(median(rankzero$rank),3)
std = round(sd(rankzero$rank),3)
print(paste("When richness doesn't change, rank log ratio (constant = 0.01) drops to median =", md,
            "sd =", std, "and mean =", mn))

# mean rank log ratio for communities where species richness change was zero
mn = round(mean(rankzero$rank_p1),3)
md = round(median(rankzero$rank_p1),3)
std = round(sd(rankzero$rank_p1),3)
print(paste("When richness doesn't change, rank log ratio (constant = 1) is: median =", md,
            "sd =", std, "and mean =", mn))

# Print Pearson's correlation coefficients
print ("Pearson's correlation coefficients")
print(paste("composition vs. richness:", round(cor.test(logratios$composition, logratios$richness)$estimate[[1]],3)))
print(paste("composition vs. evenness:", round(cor.test(logratios$composition, logratios$evenness)$estimate[[1]],3)))
print(paste("composition vs. rank:", round(cor.test(logratios$composition, logratios$rank)$estimate[[1]],3)))

print(paste("abundance vs. richness:", round(cor.test(logratios$abundance, logratios$richness)$estimate[[1]],3)))
print(paste("abundance vs. evenness:", round(cor.test(logratios$abundance, logratios$evenness)$estimate[[1]],3)))
print(paste("abundance vs. rank:", round(cor.test(logratios$abundance, logratios$rank)$estimate[[1]],3)))

# Print other summarizing information
print(paste(length(logratios$abundance[logratios$abundance==0]), "pairs did not have a change in abundance"))
print(paste(length(logratios$richness[logratios$richness==0]), "pairs did not have a change in richness"))
print(paste(length(logratios$evenness[logratios$evenness==0]), "pairs did not have a change in evenness"))

print(paste(length(logratios$abundance[logratios$abundance >= 0.6931472]), "pairs at least doubled or halved abundance"))
print(paste(length(logratios$richness[logratios$richness >= 0.6931472]), "pairs at least doubled or halved richness"))
print(paste(length(logratios$evenness[logratios$evenness >= 0.6931472]), "pairs at least doubled or halved evenness"))

lm1 = lm(Slr + Elr + Nlr + abscomplr_p01 + absranklr_p01 ~ taxa, logratios)
summary(lm1)
      
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

ggsave(site_map, file = "site_map.eps", dpi = 600, width = 9, height = 4.5)


#------------------------------------------------------------------------------------------- 
#                 FIGURE 2. histograms of log-ratio difference in treatment vs. controls
#-------------------------------------------------------------------------------------------
#plot histograms of the log-ratio results for constant=-0.01, ABSOLUTE VALUE
setEPS()
postscript(file = "Figure_2.eps", height = 6, width = 10)

abscomp = ggplot(data=logratios, aes(composition)) + geom_histogram(binwidth=0.25) + 
  xlab("median population-level |log ratio|") + ylab("frequency") +
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("A")

nhist = ggplot(data=logratios, aes(abundance)) + geom_histogram(binwidth=0.25) + 
  xlab("total abundance |log ratio|") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("B")

shist = ggplot(data=logratios, aes(richness)) + geom_histogram(binwidth=0.25) + 
  xlab("species richness |log ratio|") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("C")

ehist = ggplot(data=logratios, aes(evenness)) + geom_histogram(binwidth=0.25) + 
  xlab("Simpson's evenness |log ratio|") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("D")

absrank = ggplot(data=logratios, aes(rank)) + geom_histogram(binwidth=0.25) + 
  xlab("median rank |log ratio|") + ylab("frequency") +
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("E")

grid.arrange(abscomp, nhist, shist, ehist, absrank, nrow=2)

dev.off()

#------------------------------------------------------------------------------------------- 
#                 FIGURE 3 - pairs plots of the log-ratios and correlatin coefficients
#-------------------------------------------------------------------------------------------
#this is a hack to ggpairs to get a white background on plots
setEPS()
postscript(file = "Figure_3.eps", height = 7, width = 7)

pairs = ggpairs(logratios[,c(1:5)], upper = "blank")
text1 = ggally_text("composition", size = 6) + theme_classic()
text2 = ggally_text(paste("Cor:", round(cor.test(logratios$composition, logratios$abundance)$estimate[[1]],3)), size=6) + theme_classic()
text3 = ggally_text(paste("Cor:", round(cor.test(logratios$composition, logratios$richness)$estimate[[1]],3)), size=6) + theme_classic()
text4 = ggally_text(paste("Cor:", round(cor.test(logratios$composition, logratios$evenness)$estimate[[1]],3)), size=6) + theme_classic()
text5 = ggally_text(paste("Cor:", round(cor.test(logratios$composition, logratios$rank)$estimate[[1]],3)), size=6) + theme_classic()

pairs = putPlot(pairs, text1, 1, 1)
pairs = putPlot(pairs, text2, 1, 2)
pairs = putPlot(pairs, text3, 1, 3)
pairs = putPlot(pairs, text4, 1, 4)
pairs = putPlot(pairs, text5, 1, 5)

text6 = ggally_text("abundance", size = 6) + theme_classic()
text7 = ggally_text(paste("Cor:", round(cor.test(logratios$abundance, logratios$richness)$estimate[[1]],3)), size=6) + theme_classic()
text8 = ggally_text(paste("Cor:", round(cor.test(logratios$abundance, logratios$evenness)$estimate[[1]],3)), size=6) + theme_classic()
text9 = ggally_text(paste("Cor:", round(cor.test(logratios$abundance, logratios$rank)$estimate[[1]],3)), size=6) + theme_classic()

pairs = putPlot(pairs, text6, 2, 2)
pairs = putPlot(pairs, text7, 2, 3)
pairs = putPlot(pairs, text8, 2, 4)
pairs = putPlot(pairs, text9, 2, 5)

text10 = ggally_text("richness", size = 6) + theme_classic()
text11 = ggally_text(paste("Cor:", round(cor.test(logratios$richness, logratios$evenness)$estimate[[1]],3)), size=6) + theme_classic()
text12 = ggally_text(paste("Cor:", round(cor.test(logratios$richness, logratios$rank)$estimate[[1]],3)), size=6) + theme_classic()

pairs = putPlot(pairs, text10, 3, 3)
pairs = putPlot(pairs, text11, 3, 4)
pairs = putPlot(pairs, text12, 3, 5)

text13 = ggally_text("evenness", size = 6) + theme_classic()
text14 = ggally_text(paste("Cor:", round(cor.test(logratios$evenness, logratios$rank)$estimate[[1]],3)), size=6) + theme_classic()

pairs = putPlot(pairs, text13, 4, 4)
pairs = putPlot(pairs, text14, 4, 5)

text15 = ggally_text("rank", size = 6) + theme_classic()
pairs = putPlot(pairs, text15, 5, 5)

p1 = ggplot(logratios, aes(composition, abundance)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p2 = ggplot(logratios, aes(composition, richness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p3 = ggplot(logratios, aes(composition, evenness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p4 = ggplot(logratios, aes(composition, rank)) + geom_point() + theme_classic() + xlab("") + ylab("")
p5 = ggplot(logratios, aes(abundance, richness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p6 = ggplot(logratios, aes(abundance, evenness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p7 = ggplot(logratios, aes(abundance, rank)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p8 = ggplot(logratios, aes(richness, evenness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p9 = ggplot(logratios, aes(richness, rank)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p10 = ggplot(logratios, aes(evenness, rank)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")

pairs = putPlot(pairs, p1, 2, 1)
pairs = putPlot(pairs, p2, 3, 1)
pairs = putPlot(pairs, p3, 4, 1)
pairs = putPlot(pairs, p4, 5, 1)
pairs = putPlot(pairs, p5, 3, 2)
pairs = putPlot(pairs, p6, 4, 2)
pairs = putPlot(pairs, p7, 4, 3)
pairs = putPlot(pairs, p8, 5, 2)
pairs = putPlot(pairs, p9, 5, 3)
pairs = putPlot(pairs, p10, 5, 4)

pairs

dev.off()

#------------------------------------------------------------------------------------------- 
#                 APPENDIX FIGURE C1. Results from previous version of the ms
#-------------------------------------------------------------------------------------------
#plot histograms of the results
setEPS()
postscript(file = "Figure_C1.eps", height = 5, width = 10)

#Appendix Panel figures - plots for composition
bchist = ggplot(data=compvals, aes(BCcomp)) + geom_histogram(binwidth=0.05) + 
  xlab("species composition Bray-Curtis") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,1, by=0.25), limits = c(0,1)) + theme_classic() +  
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("A")

percnhist = ggplot(data=compvals, aes(abs(percN))) + geom_histogram(binwidth=25) + 
  xlab("|percent change| abundance") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,500, by=100), limits = c(0,500)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("B")

percshist = ggplot(data=compvals, aes(abs(percS))) + geom_histogram(binwidth=25) + 
  xlab("|percent change| richness") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,500, by=100), limits = c(0,500)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("C")

bcjhist = ggplot(data=compvals, aes(abs(BCJ))) + geom_histogram(binwidth=0.05) + 
  xlab("Simpson's evenness Bray-Curtis") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,1, by=0.25), limits = c(0,1)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("D")

bcradhist = ggplot(data=compvals, aes(abs(BCrad))) + geom_histogram(binwidth=0.05) + 
  xlab("rank Bray-Curtis") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,1, by=0.25), limits = c(0,1)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("E")

grid.arrange(bchist, percnhist, percshist, bcjhist, bcradhist, nrow=2)

dev.off()

#-----------------------------------------------------------------------------------------------
#            APPENDIX FIGURE C2. plot log ratio values before taking absolute value
#-----------------------------------------------------------------------------------------------
#plot histograms of the log-ratio results
setEPS()
postscript(file = "Figure_C2.eps", height = 7, width = 10)

#Appendix Panel figures - plots for composition
comphist = ggplot(data=lograt, aes(complr_p01)) + geom_histogram(binwidth=0.25) + 
  xlab("median population-level log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-6,6, by=1), limits = c(-6,6)) + theme_classic() +  
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("A")

comphistP1 = ggplot(data=lograt, aes(complr_p1)) + geom_histogram(binwidth=0.25) + 
  xlab("median population-level log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-6,6, by=1), limits = c(-6,6)) + theme_classic() +  
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("B")

nhist = ggplot(data=lograt, aes(Nlr)) + geom_histogram(binwidth=0.25) + 
  xlab("total abundance log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-6,6, by=1), limits = c(-6,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("C")

shist = ggplot(data=lograt, aes(Slr)) + geom_histogram(binwidth=0.25) + 
  xlab("species richness log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-6,6, by=1), limits = c(-6,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("D")

ehist = ggplot(data=lograt, aes(Elr)) + geom_histogram(binwidth=0.25) + 
  xlab("Simpson's evenness log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-6,6, by=1), limits = c(-6,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("E")

rankhist = ggplot(data=lograt, aes(ranklr_p01)) + geom_histogram(binwidth=0.25) + 
  xlab("median rank log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-6,6, by=1), limits = c(-6,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("F")

rankhistP1 = ggplot(data=lograt, aes(ranklr_p1)) + geom_histogram(binwidth=0.25) + 
  xlab("median rank log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-6,6, by=1), limits = c(-6,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=12)) + ggtitle("G")

grid.arrange(comphist, comphistP1, nhist, shist, ehist, rankhist, rankhistP1, nrow=3)

dev.off()

#-----------------------------------------------------------------------------------------------
#            APPENDIX FIGURE C3. plot all of the raw log ratio values for compostion and rank
#-----------------------------------------------------------------------------------------------
setEPS()
postscript(file = "Figure_C3.eps", height = 6, width = 10)

composrawlr = ggplot(data=data.frame(complrvals_p01), aes(complrvals_p01)) + geom_histogram(binwidth=0.25) + 
  xlab("all log-ratio of composition relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-12,12, by=2), limits = c(-12,12)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,700, by=100), limits = c(0,700)) +
  theme(text = element_text(size=12)) + ggtitle("A")

rankrawlr = ggplot(data=data.frame(ranklrvals_p01), aes(ranklrvals_p01)) + geom_histogram(binwidth=0.25) + 
  xlab("all log-ratio of rank relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-12,12, by=2), limits = c(-12,12)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,700, by=100), limits = c(0,700)) +
  theme(text = element_text(size=12)) + ggtitle("B")

composrawlrP1 = ggplot(data=data.frame(complrvals_p1), aes(complrvals_p1)) + geom_histogram(binwidth=0.25) + 
  xlab("all log-ratio of composition relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-12,12, by=2), limits = c(-12,12)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,700, by=100), limits = c(0,700)) +
  theme(text = element_text(size=12)) + ggtitle("C")

rankrawlrP1 = ggplot(data=data.frame(ranklrvals_p1), aes(ranklrvals_p1)) + geom_histogram(binwidth=0.25) + 
  xlab("all log-ratio of rank relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-12,12, by=2), limits = c(-12,12)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,700, by=100), limits = c(0,700)) +
  theme(text = element_text(size=12)) + ggtitle("D")

grid.arrange(composrawlr, rankrawlr, composrawlrP1, rankrawlrP1, nrow = 2)

dev.off()

#----------------------------------------------------------------------- 
#            APPENDIX FIGURE C-4. compare control and manipulated data
#-----------------------------------------------------------------------
#plot results along 1:1 line
df = data.frame(as.numeric(compc),as.numeric(compe))
names(df) = c("compc", "compe")

pdf(file = "Figure_C4.pdf", height = 6, width = 10)

#Fig A - plots for composition
compchange = ggplot(df, aes(compc, compe)) + geom_point(alpha=0.5, size=3) + 
  xlab("species relative abundance") + ylab("species relative abundance") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=12)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("A")

#Fig B - plots for abundance
abunchange = ggplot(data=diversity, aes(CN, EN)) + geom_point(alpha=0.5, size=3) + 
  xlab("total abundance") + ylab("total abundance") + 
  scale_x_log10(breaks = c(30, 100, 500, 2500, 6500), limits = c(30,6500)) +
  scale_y_log10(breaks = c(30, 100, 500, 2500, 6500), limits = c(30,6500)) + theme_classic() +
  theme(text = element_text(size=12)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("B")

# Fig C - plots for S
schange = ggplot(data=diversity, aes(CS, ES)) + geom_point(alpha=0.5, size=3) + 
  xlab("species richness") + ylab("species richness") + 
  scale_x_log10(breaks = c(5, 10, 25, 50, 100, 200), limits = c(5, 200)) +
  scale_y_log10(breaks = c(5, 10, 25, 50, 100, 200), limits = c(5, 200)) + theme_classic() +
  theme(text = element_text(size=12)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("C")

# Fig D - plots for evenness
evenchange = ggplot(data=diversity, aes(Jc, Je)) + geom_point(alpha=0.5, size=3) + 
  xlab("Simpson's evenness") + ylab("Simpson's evenness") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=12)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("D")

# Fig E - compare relative abundance at each rank in all paired sites
rankabunchange = ggplot(data=relabundance, aes(c, e)) + geom_point(alpha=0.5, size=3) + 
  xlab("rank relative abundance") + ylab("rank relative abundance") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=12)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("E")

grid.arrange(compchange, abunchange, schange, evenchange, rankabunchange, nrow=2)

dev.off()

#----------------------------------------------------------------------------------
#             APPENDIX FIGURE C-5. Compare species population response
#----------------------------------------------------------------------------------
setEPS()
postscript(file = "Figure_C5.eps", height = 4, width = 10)

#considering all species + 0.01
all = ggplot(compositionvals_p01, aes(abscomplr)) + geom_histogram(binwidth=0.5) + theme_classic() +
  theme(text = element_text(size=20)) + ggtitle("A") + xlab("median |log ratio|") +
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) 

#considering all species + 1
allP1 = ggplot(compositionvals_p1, aes(abscomplr_p1)) + geom_histogram(binwidth=0.5) + theme_classic() +
  theme(text = element_text(size=20)) + ggtitle("B") + xlab("median |log ratio|") +
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) 

#considering only species present in both control and experiment plots
nonzero = ggplot(compositionvals_p01, aes(abscomplrzero)) + geom_histogram(binwidth=0.5) + theme_classic() +
  theme(text = element_text(size=20)) + ggtitle("C")  + xlab("median |log ratio|") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6))

grid.arrange(all, allP1, nonzero, nrow = 1)

dev.off()


#----------------------------------------------------------------------------------
#             APPENDIX FIGURE C-6. Compare rank abundance distribution response
#----------------------------------------------------------------------------------
setEPS()
postscript(file = "Figure_C6.eps", height = 6, width = 10)

#considering all species + 0.01
all = ggplot(logratios, aes(rank)) + geom_histogram(binwidth=0.25) + theme_classic() +
  theme(text = element_text(size=12)) + ggtitle("A") + xlab("rank median |log ratio|") +
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) 

#considering only species present in both control and experiment plots + 0.01
zerochange = ggplot(rankzero, aes(rank)) + geom_histogram(binwidth=0.25) + theme_classic() +
  theme(text = element_text(size=12)) + ggtitle("B")  + xlab("rank median |log ratio|") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) 

#considering all species + 1
allP1 = ggplot(logratios, aes(rank_p1)) + geom_histogram(binwidth=0.25) + theme_classic() +
  theme(text = element_text(size=12)) + ggtitle("C") + xlab("rank median |log ratio|") +
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6))

#considering only species present in both control and experiment plots + 1
zerochangeP1 = ggplot(rankzero, aes(rank_p1)) + geom_histogram(binwidth=0.25) + theme_classic() +
  theme(text = element_text(size=12)) + ggtitle("D")  + xlab("rank median |log ratio|") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) 

grid.arrange(all, zerochange, allP1, zerochangeP1, nrow = 2)

dev.off()


#---------------------------------------------------------------------------------------------- 
#             APPENDIX FIGURE C-7. pairs plots of the log-ratios and correlation coefficients
#----------------------------------------------------------------------------------------------
setEPS()
postscript(file = "Figure_C7.eps", height = 7, width = 7)

#this is a hack to ggpairs to get a white background on plots
pairs = ggpairs(logratios[,c(6,2,3,4,7)], upper = "blank")
text1 = ggally_text("composition", size = 6) + theme_classic()
text2 = ggally_text(paste("Cor:", round(cor.test(logratios$composition_p1, logratios$abundance)$estimate[[1]],3)), size=6) + theme_classic()
text3 = ggally_text(paste("Cor:", round(cor.test(logratios$composition_p1, logratios$richness)$estimate[[1]],3)), size=6) + theme_classic()
text4 = ggally_text(paste("Cor:", round(cor.test(logratios$composition_p1, logratios$evenness)$estimate[[1]],3)), size=6) + theme_classic()
text5 = ggally_text(paste("Cor:", round(cor.test(logratios$composition_p1, logratios$rank_p1)$estimate[[1]],3)), size=6) + theme_classic()

pairs = putPlot(pairs, text1, 1, 1)
pairs = putPlot(pairs, text2, 1, 2)
pairs = putPlot(pairs, text3, 1, 3)
pairs = putPlot(pairs, text4, 1, 4)
pairs = putPlot(pairs, text5, 1, 5)

text6 = ggally_text("abundance", size = 6) + theme_classic()
text7 = ggally_text(paste("Cor:", round(cor.test(logratios$abundance, logratios$richness)$estimate[[1]],3)), size=6) + theme_classic()
text8 = ggally_text(paste("Cor:", round(cor.test(logratios$abundance, logratios$evenness)$estimate[[1]],3)), size=6) + theme_classic()
text9 = ggally_text(paste("Cor:", round(cor.test(logratios$abundance, logratios$rank_p1)$estimate[[1]],3)), size=6) + theme_classic()

pairs = putPlot(pairs, text6, 2, 2)
pairs = putPlot(pairs, text7, 2, 3)
pairs = putPlot(pairs, text8, 2, 4)
pairs = putPlot(pairs, text9, 2, 5)

text10 = ggally_text("richness", size = 6) + theme_classic()
text11 = ggally_text(paste("Cor:", round(cor.test(logratios$richness, logratios$evenness)$estimate[[1]],3)), size=6) + theme_classic()
text12 = ggally_text(paste("Cor:", round(cor.test(logratios$richness, logratios$rank_p1)$estimate[[1]],3)), size=6) + theme_classic()

pairs = putPlot(pairs, text10, 3, 3)
pairs = putPlot(pairs, text11, 3, 4)
pairs = putPlot(pairs, text12, 3, 5)

text13 = ggally_text("evenness", size = 6) + theme_classic()
text14 = ggally_text(paste("Cor:", round(cor.test(logratios$evenness, logratios$rank_p1)$estimate[[1]],3)), size=6) + theme_classic()

pairs = putPlot(pairs, text13, 4, 4)
pairs = putPlot(pairs, text14, 4, 5)

text15 = ggally_text("rank", size = 6) + theme_classic()
pairs = putPlot(pairs, text15, 5, 5)

p1 = ggplot(logratios, aes(composition_p1, abundance)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p2 = ggplot(logratios, aes(composition_p1, richness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p3 = ggplot(logratios, aes(composition_p1, evenness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p4 = ggplot(logratios, aes(composition_p1, rank_p1)) + geom_point() + theme_classic() + xlab("") + ylab("")
p5 = ggplot(logratios, aes(abundance, richness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p6 = ggplot(logratios, aes(abundance, evenness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p7 = ggplot(logratios, aes(abundance, rank_p1)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p8 = ggplot(logratios, aes(richness, evenness)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p9 = ggplot(logratios, aes(richness, rank_p1)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")
p10 = ggplot(logratios, aes(evenness, rank_p1)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) + xlab("") + ylab("")

pairs = putPlot(pairs, p1, 2, 1)
pairs = putPlot(pairs, p2, 3, 1)
pairs = putPlot(pairs, p3, 4, 1)
pairs = putPlot(pairs, p4, 5, 1)
pairs = putPlot(pairs, p5, 3, 2)
pairs = putPlot(pairs, p6, 4, 2)
pairs = putPlot(pairs, p7, 4, 3)
pairs = putPlot(pairs, p8, 5, 2)
pairs = putPlot(pairs, p9, 5, 3)
pairs = putPlot(pairs, p10, 5, 4)

pairs

dev.off()
