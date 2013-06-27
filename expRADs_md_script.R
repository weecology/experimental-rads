# All files should be organized as such:
# RAD: SiteID | experiment/control | year | genus | species | id2species | abundance
# Code will compare S, N and composition among treatments at each site,
# output and compare RADs, parameters for each treatment at each site.

#---------------------------------------------------------------------------------
#          setup - select wd, import data, source code,  file to collect results
#---------------------------------------------------------------------------------

wd = "/Users/sarah/Documents/GitHub/experimental-rads"
setwd(wd)

source("ExpRADsFunctions.R")   #Run the code containing the functions

comms = read.csv("community_analysis_data.csv", na.strings = 'NULL')
comps = read.csv("comparison_analysis_data.csv")
names(comps)<-c('ref', 'controID','expID')
expers = read.csv("experiments_analysis_data.csv")

#--------------------------------------------------------------
#          generate values and comparisons between the sites 
#--------------------------------------------------------------
#open plotting pdf window
pdf("allRADs.pdf", 7, 10, paper = "letter", pointsize = 10)
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

for (iRow in 1:nrow(comps)){
  control = comps[iRow,2]  #find control in pair
  experiment = comps[iRow,3]  # find experiment in pair
  taxa = as.character(expers[which(expers[,2]==control),7])# find taxonomic group from experiments table
  type = as.character(expers[which(expers[,2]==control),4]) # find experiment type from experiments table
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
n_rmse = round(rmse(EN, CN),4)
s_rmse = round(rmse(ES, CS),4)
relabun_rmse = round(rmse(e,c),4)
j_rmse = round(rmse(Je, Jc),4)

n_r2 = round(rsquare(CN,EN),4)
s_r2 = round(rsquare(CS,ES),4)
relabun_r2 = round(rsquare(c,e),4)
j_r2 = round(rsquare(Jc,Je),4)

#count the communities displaying various shapes. This Does take into account duplicates. Should be correct
shapes = count_RAD_shapes(cID, eID, Cshape, Eshape)


#--------------------------------------------------------
#          Variance partitioning analysis and results 
#--------------------------------------------------------

#### Standardize the variables
stdz_bc_rad = standardize(BCrad)
stdz_bc_s = standardize(BCS)
stdz_bc_n = standardize(BCN)
stdz_perc_s = standardize(abs(percS))
stdz_perc_n = standardize(abs(percN))
stdz_bc_comp = standardize(BCcomp)
taxon = as.factor(taxon)
etype = as.factor(etype)

# variance partitioning of variable impact on RADs, using Bray-Curtis S and N
lm_full = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_bc_s + stdz_bc_n + stdz_bc_s:stdz_bc_n)
lm_comp = lm(stdz_bc_rad ~ stdz_bc_comp)
lm_sn = lm(stdz_bc_rad ~ stdz_bc_s + stdz_bc_n + stdz_bc_s:stdz_bc_n)

r2_full = round(summary.lm(lm_full)$r.squared,4)
r2_comp = round(summary.lm(lm_comp)$r.squared,4)
r2_sn = round(summary.lm(lm_sn)$r.squared,4)

var_comp = round(summary.lm(lm_full)$r.squared - summary.lm(lm_sn)$r.squared, 4)
var_sn = round(summary.lm(lm_full)$r.squared-summary.lm(lm_comp)$r.squared,4)

### variance partitioning of variable impact on RADs using Bray-Curtis S and N, composition and TAXONOMIC GROUP
lm_bc_full = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_bc_s + stdz_bc_n + stdz_bc_s:stdz_bc_n + taxon)
lm_bc_comm_vars = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_bc_s + stdz_bc_n + stdz_bc_s:stdz_bc_n)
lm_bc_taxa = lm(stdz_bc_rad ~ taxon)

r2_bcfull = round(summary.lm(lm_bc_full)$r.squared,4)
r2_bccomm_vars = round(summary.lm(lm_bc_comm_vars)$r.squared,4)
r2_bctaxa = round(summary.lm(lm_bc_taxa)$r.squared,4)

var_bccommvars = round(summary.lm(lm_bc_full)$r.squared - summary.lm(lm_bc_taxa)$r.squared, 4)
var_bctaxa = round(summary.lm(lm_bc_full)$r.squared-summary.lm(lm_bc_comm_vars)$r.squared,4)

### variance partitioning of variable impact on RADs using Bray-Curtis S and N, composition and EXPERIMENT TYPE 
lm_bc_full = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_bc_s + stdz_bc_n + stdz_bc_s:stdz_bc_n + etype)
lm_bc_comm_vars = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_bc_s + stdz_bc_n + stdz_bc_s:stdz_bc_n)
lm_bc_etype = lm(stdz_bc_rad ~ etype)

r2_bcfull = round(summary.lm(lm_bc_full)$r.squared,4)
r2_bccomm_vars = round(summary.lm(lm_bc_comm_vars)$r.squared,4)
r2_bcetype = round(summary.lm(lm_bc_etype)$r.squared,4)

var_bccommvars = round(summary.lm(lm_bc_full)$r.squared - summary.lm(lm_bc_etype)$r.squared, 4)
var_bcetype = round(summary.lm(lm_bc_full)$r.squared-summary.lm(lm_bc_comm_vars)$r.squared,4)

# variance partitioning of variable impact on RADs, using abs % difference S and N
lm_perc_full = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_perc_s + stdz_perc_n + stdz_perc_s:stdz_perc_n)
lm_perc_comp = lm(stdz_bc_rad ~ stdz_bc_comp)
lm_perc_sn = lm(stdz_bc_rad ~ stdz_perc_s + stdz_perc_n + stdz_perc_n:stdz_perc_s)

r2_pfull = round(summary.lm(lm_perc_full)$r.squared,4)
r2_pcomp = round(summary.lm(lm_perc_comp)$r.squared,4)
r2_psn = round(summary.lm(lm_perc_sn)$r.squared,4)

var_pcomp = round(summary.lm(lm_perc_full)$r.squared - summary.lm(lm_perc_sn)$r.squared, 4)
var_psn = round(summary.lm(lm_perc_full)$r.squared-summary.lm(lm_perc_comp)$r.squared, 4)

### variance partitioning of variable impact on RADs, using abs % difference S and N, composition and TAXONOMIC GROUP
lm_perc_full = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_perc_s + stdz_perc_n + stdz_perc_s:stdz_perc_n + taxon)
lm_perc_comm_vars = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_perc_s + stdz_perc_n + stdz_perc_s:stdz_perc_n)
lm_perc_taxa = lm(stdz_bc_rad ~ taxon)

r2_pfull = round(summary.lm(lm_perc_full)$r.squared, 4)
r2_pcomm_vars = round(summary.lm(lm_perc_comm_vars)$r.squared, 4)
r2_ptaxa = round(summary.lm(lm_perc_taxa)$r.squared, 4)

var_pcommvars = round(summary.lm(lm_perc_full)$r.squared - summary.lm(lm_perc_taxa)$r.squared, 4)
var_ptaxa = round(summary.lm(lm_perc_full)$r.squared - summary.lm(lm_perc_comm_vars)$r.squared, 4)

### variance partitioning of variable impact on RADs, using abs % difference S and N, composition and EXPERIMENT TYPE 
lm_perc_full = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_perc_s + stdz_perc_n + stdz_perc_s:stdz_perc_n + etype)
lm_perc_comm_vars = lm(stdz_bc_rad ~ stdz_bc_comp + stdz_perc_s + stdz_perc_n + stdz_perc_s:stdz_perc_n)
lm_perc_etype = lm(stdz_bc_rad ~ etype)

r2_full = round(summary.lm(lm_perc_full)$r.squared, 4)
r2_comm_vars = round(summary.lm(lm_perc_comm_vars)$r.squared, 4)
r2_etype = round(summary.lm(lm_perc_etype)$r.squared, 4)

var_pcommvars = round(summary.lm(lm_perc_full)$r.squared - summary.lm(lm_perc_etype)$r.squared, 4)
var_etype = round(summary.lm(lm_perc_full)$r.squared - summary.lm(lm_perc_comm_vars)$r.squared, 4)
