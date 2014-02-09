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
rankvals = data.frame(r2 = 1, ranklr = 1, sdrlr = 1, absranklr = 1, abssdrlr = 1)
compositionvals = data.frame(compr2 = 1, complr = 1, sdlr = 1, abscomplr = 1, abssdlr = 1, abscomplrzero = 1)
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
        tr = tr + 0.0001 #a a small constant to account for zeroes in the data
      rlr = sapply(tr, function(x) LogRatio(x[2],x[1]) )
      ranklrvals = append(ranklrvals, rlr)
      rankvals[outcount,] = c(rsquare(comparison_matrix[,1], comparison_matrix[,2]), median(rlr), sd(rlr), median(abs(rlr)), sd(abs(rlr)))
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
      comparison = comparison + 0.0001  #a a small constant to account for zeroes in the data
      #take the log ratio of raw abundance
        lr = sapply(comparison, function(x) LogRatio(x[2], x[1]) )
        lr_nz = sapply(comparison, function(x) LogRatio_noZero(x[2], x[1]))
      complrvals = append(complrvals, lr)
      compositionvals[outcount,] = c(rsquare(as.numeric(comparison[1,]), as.numeric(comparison[2,])), median(lr), sd(lr), median(abs(lr)), sd(abs(lr)), median(abs(lr_nz),na.rm=TRUE))
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

# calculate the log ratio differences between values
s_df = as.data.frame(t(cbind(results$CS, results$ES)))
n_df = as.data.frame(t(cbind(results$CN, results$EN)))
e_df = as.data.frame(t(cbind(results$Jc, results$Je)))

Slr = sapply(s_df, function(x) LogRatio(x[2], x[1]) )
Nlr = sapply(n_df, function(x) LogRatio(x[2], x[1]) )
Elr = sapply(e_df, function(x) LogRatio(x[2], x[1]) )

#put the results together for later plotting and comparison
diversity = cbind(charvars[,c(4:5)], desc)
composition = data.frame(compc, compe)
relabundance = data.frame(c,e)

taxa = diversity$taxon
etype = diversity$etype
complr = compositionvals$complr
abscomplr = compositionvals$abscomplr
ranklr = rankvals$ranklr
absranklr = rankvals$absranklr
lograt = data.frame(complr, ranklr, abscomplr, absranklr, Nlr, Slr, Elr, taxa, etype)

# for pairs plot
logratios = lograt[,c(3,5,6,7,4)]
logratios$Nlr = abs(logratios$Nlr)
logratios$Slr = abs(logratios$Slr)
logratios$Elr = abs(logratios$Elr)
names(logratios) = c("composition", "abundance", "richness", "evenness", "rank")


#------------------------------------- 
#                 results
#-------------------------------------

# Find the R2 values for the paired data
comp_r2 = round(rsquare(compc, compe),4)
n_r2 = round(rsquare(results$CN,results$EN),4)
s_r2 = round(rsquare(results$CS,results$ES),4)
relabun_r2 = round(rsquare(c,e),4)
j_r2 = round(rsquare(results$Jc,results$Je),4)

# count the communities displaying various shapes. This does take into account duplicates. Should be correct
shapes = count_RAD_shapes(results$cID, results$eID, results$Cshape, results$Eshape)

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
mn = round(mean(compositionvals$abscomplrzero),3)
md = round(median(compositionvals$abscomplrzero),3)
print(paste("when only consider species that were present in both communities, composition log ratio is median = ",
            md, "and mean =", mn))

# mean rank log ratio for communities where species richness change was zero
rankzero = logratios[which(logratios$richness == 0),]
mn = round(mean(rankzero$rank),3)
md = round(median(rankzero$rank),3)
print(paste("When richness doesn't change, rank log ratio drops to median =", md,
            "and mean =", mn))

# Print Pearson's correlation coefficients
print ("Pearson's correlation coefficients")
print(paste("composition vs. richness:", round(cor.test(logratios$composition, logratios$richness)$estimate[[1]],3)))
print(paste("composition vs. evenness:", round(cor.test(logratios$composition, logratios$evenness)$estimate[[1]],3)))
print(paste("composition vs. rank:", round(cor.test(logratios$composition, logratios$rank)$estimate[[1]],3)))

print(paste("abundance vs. richness:", round(cor.test(logratios$abundance, logratios$richness)$estimate[[1]],3)))
print(paste("abundance vs. evenness:", round(cor.test(logratios$abundance, logratios$evenness)$estimate[[1]],3)))
print(paste("abundance vs. rank:", round(cor.test(logratios$abundance, logratios$rank)$estimate[[1]],3)))

print(paste(length(logratios$abundance[logratios$abundance==0]), "pairs did not have a change in abundance"))
print(paste(length(logratios$richness[logratios$richness==0]), "pairs did not have a change in richness"))
print(paste(length(logratios$evenness[logratios$evenness==0]), "pairs did not have a change in evenness"))

print(paste(length(logratios$abundance[logratios$abundance >= 0.6931472]), "pairs at least doubled or halved abundance"))
print(paste(length(logratios$richness[logratios$richness >= 0.6931472]), "pairs at least doubled or halved richness"))
print(paste(length(logratios$evenness[logratios$evenness >= 0.6931472]), "pairs at least doubled or halved evenness"))

      
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


#------------------------------------------------------------------------------------------- 
#                 FIGURE 2. histograms of log-ratio difference in treatment vs. controls
#-------------------------------------------------------------------------------------------
#plot histograms of the log-ratio results, ABSOLUTE VALUE

abscomp = ggplot(data=results, aes(abscomplr)) + geom_histogram() + 
  xlab("composition median absolute log ratio") + ylab("frequency") +
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=16)) + ggtitle("A")

nhist = ggplot(data=lograt, aes(abs(Nlr))) + geom_histogram() + 
  xlab("total abundance absolute log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=16)) + ggtitle("B")

shist = ggplot(data=lograt, aes(abs(Slr))) + geom_histogram() + 
  xlab("species richness absolute log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=16)) + ggtitle("C")

ehist = ggplot(data=lograt, aes(abs(Elr))) + geom_histogram() + 
  xlab("Simpson's evenness absolute log ratio") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=16)) + ggtitle("D")

absrank = ggplot(data=results, aes(absranklr)) + geom_histogram() + 
  xlab("rank median absolute log ratio") + ylab("frequency") +
  scale_x_continuous(breaks = seq(0,6, by=1), limits = c(0,6)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=16)) + ggtitle("E")

grid.arrange(abscomp, nhist, shist, ehist, absrank, nrow=2)


#------------------------------------------------------------------------------------------- 
#                 FIGURE 3 - pairs plots of the log-ratios and correlatin coefficients
#-------------------------------------------------------------------------------------------
#this is a hack to ggpairs to get a white background on plots
pairs = ggpairs(logratios, upper = "blank")
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




#------------------------------------------------------------------------------------------- 
#                 APPENDIX FIGURE. similar to Fig 2, but with directional data, histograms of log-ratio difference in treatment vs. controls
#-------------------------------------------------------------------------------------------
#plot histograms of the log-ratio results + Bray-Curtis

#Appendix Panel figures - plots for composition
BChist = ggplot(data=results, aes(BCcomp)) + geom_histogram() + 
  xlab("Bray-Curtis dissimilarity") + ylab("frequency") +
  scale_x_continuous(breaks = seq(0,1, by=0.2), limits = c(0,1)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("A")

comphist = ggplot(data=lograt, aes(complr)) + geom_histogram() + 
  xlab("mean log-ratio of species relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-3,3, by=1), limits = c(-3,3)) + theme_classic() +  
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("B")

nhist = ggplot(data=lograt, aes(Nlr)) + geom_histogram() + 
  xlab("log-ratio of total abundance") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-3,3, by=1), limits = c(-3,3)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("C")

shist = ggplot(data=lograt, aes(Slr)) + geom_histogram() + 
  xlab("log-ratio of species richness") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-3,3, by=1), limits = c(-3,3)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("D")

ehist = ggplot(data=lograt, aes(Elr)) + geom_histogram() + 
  xlab("log-ratio of Simpson's evenness") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-3,3, by=1), limits = c(-3,3)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("E")

rankhist = ggplot(data=lograt, aes(ranklr)) + geom_histogram() + 
  xlab("mean log-ratio of rank relative abundances") + ylab("frequency") + 
  scale_x_continuous(breaks = seq(-3,3, by=1), limits = c(-3,3)) + theme_classic() + 
  scale_y_continuous(breaks = seq(0,40, by=10), limits = c(0,40)) +
  theme(text = element_text(size=20)) + ggtitle("F")

grid.arrange(BChist, comphist, nhist, shist, ehist, rankhist, nrow=2)

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



#----------------------------------------------------------------------------------
#                 RESULTS - mean, standard deviation, and correlation 
#----------------------------------------------------------------------------------







#----------------------------------------------------------------------- 
#            APPENDIX FIGURE C-4. compare control and manipulated data
#-----------------------------------------------------------------------
#plot results along 1:1 line

#Fig A - plots for composition
compchange = ggplot(data=composition, aes(compc, compe)) + geom_point(alpha=0.5, size=3) + 
  xlab("species relative abundance") + ylab("species relative abundance") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("A")

#Fig B - plots for abundance
abunchange = ggplot(data=diversity, aes(CN, EN)) + geom_point(alpha=0.5, size=3) + 
  xlab("total abundance") + ylab("total abundance") + 
  scale_x_log10(breaks = c(30, 100, 500, 2500, 6500), limits = c(30,6500)) +
  scale_y_log10(breaks = c(30, 100, 500, 2500, 6500), limits = c(30,6500)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("B")

# Fig C - plots for S
schange = ggplot(data=diversity, aes(CS, ES)) + geom_point(alpha=0.5, size=3) + 
  xlab("species richness") + ylab("species richness") + 
  scale_x_log10(breaks = c(5, 10, 25, 50, 100, 200), limits = c(5, 200)) +
  scale_y_log10(breaks = c(5, 10, 25, 50, 100, 200), limits = c(5, 200)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("C")

# Fig D - plots for evenness
evenchange = ggplot(data=diversity, aes(Jc, Je)) + geom_point(alpha=0.5, size=3) + 
  xlab("Simpson's evenness") + ylab("Simpson's evenness") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("D")

# Fig E - compare relative abundance at each rank in all paired sites
rankabunchange = ggplot(data=relabundance, aes(c, e)) + geom_point(alpha=0.5, size=3) + 
  xlab("rank relative abundance") + ylab("rank relative abundance") + 
  scale_x_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits = c(0,1)) + theme_classic() +
  theme(text = element_text(size=20)) + 
  geom_abline(intercept = 0, slope = 1) + ggtitle("E")

grid.arrange(compchange, abunchange, schange, evenchange, rankabunchange, nrow=2)
