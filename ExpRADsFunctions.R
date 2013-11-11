##### This is the source code for use with Experimental_RADs.R, the Experiments and Macroecology project

require(vegan)
require(BiodiversityR)
require(plotrix)
require(graphics)
require(CCA)
require(VGAM)
require(nlme)
require(lme4)
#require(languageR)
require(poilog)
require(scatterplot3d)
require(hydroGOF)

# All files should be organized as such:
# RAD: SiteID | experiment/control | year | genus | species | id2species | abundance
# Code will compare S, N and composition among treatments at each site,
# output and compare RADs, parameters for each treatment at each site.

reshape_data = function(dat){
  # puts data in "wide" format CATEGORIZED BY SITEID for use with composition analysis
  
  dat2 = aggregate(list(abun = dat$abundance),by=list(siteID = dat$siteID, genus = dat$genus,
                                                  species=dat$species),FUN=sum)
  dat2$genus_spp = paste(dat2$genus, dat2$species, sep = "_")
  dat2 = dat2[,c(1,5,4)]
  
  dat3 = reshape(dat2, idvar = c('siteID'), timevar = 'genus_spp', direction = 'wide')
  dat3[is.na(dat3)] = 0
  
  return(dat3)
}


SimpE = function(spp_data){
  # input a matrix where species are columns and rows are sites, with two sites. Calculates inverse simpson's diversity
  # and returns Simpson's evenness, where 1 represents a very even community and 0 represents a very uneven community
  spp_matrix = reshape_data(spp_data)
  D = diversity(spp_matrix[,2:ncol(spp_matrix)], index = 'invsimpson') #leave out first column, which is siteID
  J = D / specnumber(spp_matrix[,2:ncol(spp_matrix)])
  return (J)
}


BCdist = function(two_spp_data){
  # input a matrix where species are columns and rows are sites, with two sites. Compares the two rows using bray-curtis metric
  # returns the bray-curtis dissimilarity value, where 0 is the same and 1 is completely different
  if (ncol(two_spp_data) == 2){
  data = t(two_spp_data) # transpose the data matrix
  }
  else {
    data = reshape_data(two_spp_data) 
    data = data[,c(2:ncol(data))]
    data[1,] = data[1,]/sum(data[1,])
    data[2,] = data[2,]/sum(data[2,])
  }
  BC = vegdist(data, method = 'bray')
  return(BC)
}

CCA_composition = function(species_matrix, experiment){
  # input species matrix with species as columns and rows as site variables and experiment data (in same order) as
  # a factor variable. Makes a CCA plot with centroids and returns a significance test on the pCCA axes.
  pcca <- cca(species_matrix ~ experiment)
  CCA1 = as.vector(pcca$CCA$centroids[,1])
  CCA2 = as.vector(pcca$CCA$centroids[,2])
  pl <- plot(pcca, display = c("sp",'cn'), type = "n",xlab='', ylab='', xaxt='n')
  sel <- orditorp(pcca, dis = 'sp', lab = F, pcol = 'gray40', pch =19)    
  text(CCA1, CCA2, labels = c('1','2'), col='firebrick3', cex=1.25)
  
  pccaResult = anova(pcca)
  return (pccaResult)
}


percent_unidSpp = function(siteID, dataframe) {
  # input a numeric code for siteID and the dataframe, determines if the data can be used. 
  # col 2 of the dataframe is siteID, col 6 indicates level of identification & col 8 is the abundance
  # each row of the dataframe is a species-level record.
  # Data can be used if unidentified individuals are < 10% of the total number of individuals
  NAindex = which(is.na(dataframe[which(dataframe[,2] == siteID),8]))
  speciesData = as.numeric(dataframe[which(dataframe[,2] == siteID & dataframe[,7] != 0),8])
  unidData = as.numeric(dataframe[which(dataframe[,2] == siteID & dataframe[,7] == 0),8])
  # if there are NAs in dataframe, then P/A data, data is NOT OK
  if (length(NAindex) > 0) {
    return('NA')
  }
  #if there are no unid species, data is OK
  else if (length(unidData) == 0){
    return('OK')
  }
  # if data passes both checks, calculate percent unidentified individuals
  else {
    percent_unid = sum(unidData)/(sum(speciesData) + sum(unidData))
    if (percent_unid <= 0.10) {
      return ('OK') }
    else {
      return (round(percent_unid * 100,2)) 
    }
  }
}


get_lambda_sad=function(S,N){
  ## Function: Calculate lambda_1 using (B.4) from Harte et al. 2008
  # This is the untruncated form of the log-series MLE
  if (N<=0){
    print("Error: N must be greater than 0.")
    return(NA)}
  else if (S<=0){
    print("Error: S must be greater than 0.")
    return(NA)}
  else if (S/N>=1){
    print("Error: N must be greater than S.")
    return(NA)}
  else {
    ## Solve for lambda 1 in Harte et al. 2008 based on (B.4)
    y=function(x) 1/log(1/(1-x))*x/(1-x)-N/S
    ## There should be one root between 0 and 1 since y is monotonously decreasing
    bound=10^-10  ## Define lower and upper bound for x
    p=uniroot(y,lower=bound,upper=1-bound)$root   ## Parameter for log-series
    lambda_sad = p
    return(lambda_sad)
  }
}


chi.test = function(v1, v2){
  ## Function: chi square test of two vectors     ## FIX ME, MAKE SURE I WORK CORRECTLY!
  dat_comb=c(v1, v2)
  q=unique(as.numeric(quantile(dat_comb, probs=seq(0, 1, 0.2)))) ## 5 bins
  count.v1=hist(v1, breaks=q,plot=FALSE)$counts
  count.v2=hist(v2, breaks=q,plot=FALSE)$counts
  p.chi=chisq.test(cbind(count.v1, count.v2), simulate.p.value=TRUE)$p.value
  return(p.chi)
}


sad_shape = function(abund){
  # Function: input a vector of abundances. Calculates MLE of log-series and of poisson log-normal.
  # Gives a numeric weight for which model is better supported. Outputs 'logs' or 'logn'.
  #mle for untruncated log-series
  p = get_lambda_sad(length(abund), sum(abund)) # inputs: S, N
  lik.logs = sum(dlog(abund, prob = p, log = TRUE))
  #mle for untruncated poisson log-normal
  logn_vals = as.list(poilogMLE(abund, startVals = c(mean(log(abund)), sd(log(abund))))$par)
  lik.logn = sum(log(dpoilog(abund, logn_vals$mu, logn_vals$sig)))
  n = length(abund)
  k1 = 1
  k2 = 2
  AICc.logs = 2 * k1 - 2 * lik.logs + 2 * k1 * (k1 + 1) / (n - k1 - 1)
  AICc.logn = 2 * k2 - 2 * lik.logn + 2 * k2 * (k2 + 1) / (n - k2 - 1)
  AICc.min = min(AICc.logs, AICc.logn)
  weight.logs = exp(-(AICc.logs - AICc.min) / 2)
  weight.logn = exp(-(AICc.logn - AICc.min) / 2)
  weight = weight.logs / (weight.logs + weight.logn)  ## weight of data coming from a log-series
  if (weight > 0.5){
    type = 'logs'}
  else {type = 'logn'}
  return(list(type = type,w = weight))
}


dist.test=function(v1, v2){
  ## Function to compare the predicted abundances with empirical abundances
  ## Using both chi-square test and weights b/w log-series and log-normal
  ## pred is value (list) taken from SAD_rank.r
  ## dat is one row in Sarah's winter_tot.csv
  require(poilog)
  ## chi-square test
  p.chi = chi.test(v1,v2)
  ## weights
  weight1 = sad_shape(v1)
  weight2 = sad_shape(v2)
  
  return(list(p_chi = p.chi, con = weight1$type, exp = weight2$type))
}


rsquare = function(con, exp){
  # get fit to 1:1 line
  return (1 - sum((con - exp) ** 2) / sum((con - mean(con)) ** 2))
}


RAD_plot_and_data = function(siteID1, siteID2, abundance1, abundance2, taxa){
  #Input a numeric site ID, a list of species and a list of nonzero abundance. Makes a logscale Rank-abundance
  #graph labeled with the string idenfier and outputs species richness, total abundance, and parameter p value.  
  #gets poisson log-normal MLE parameters mu and sigma
  cont_par = as.list(poilogMLE(abundance1, startVals = c(mean(log(abundance1)), sd(log(abundance1))))$par)
  expt_par = as.list(poilogMLE(abundance2, startVals = c(mean(log(abundance2)), sd(log(abundance2))))$par)
  #plot relative abundance, rank abunduance distributions
  relabun1 = relabund(abundance1)
  relabun2 = relabund(abundance2)
  plot(NA, NA, xlab = '', ylab = '', xlim = c(0, max(c(length(abundance1), length(abundance2)))), 
        ylim = c(0,max(c(max(relabun1),max(relabun2)))), cex.axis = 0.75)
  lines(c(1:length(abundance1)), relabun1, type = 'l', pch = 20, cex = 2, lwd = 2)
  lines(c(1:length(abundance2)), relabun2, type = 'l', pch = 20, cex = 2, lwd = 2, lty = 2, col = 'deeppink3')
    legend('topright', c('control', 'experiment'), bty = 'n', lty = c(1,2), lwd = 2, col = c('black', 'deeppink3'), seg.len = 3, cex = 0.75)
      #legend('topright', c('control', paste('mu, sigma = ', round(cont_par$mu,2), ', ', round(cont_par$sig,2)),
      #                     'experiment', paste('mu, sigma = ', round(expt_par$mu,2),', ', round(expt_par$sig,2))), 
      #      bty = 'n', lty = c(1,NA,2,NA), lwd = 2, col = c('black', NA, 'deeppink3', NA), cex = 0.75)
    mtext('Relative Abundance', side = 2, line = 2, cex = 0.5)
    mtext('Rank', side = 1, line = 2, cex = 0.5)
    mtext(paste(taxa, ': ', siteID1, 'vs.', siteID2, sep = ' '), side = 3, line = .5, cex = 0.75)
  
  SAD_results = c(siteID1, cont_par$mu, cont_par$sig, siteID2, expt_par$mu, expt_par$sig)
  return (SAD_results)
}
 
RAD_plot = function(siteID1, siteID2, abundance1, abundance2, taxa){
  #Input a numeric site ID, a list of species and a list of nonzero abundance. Makes a relative abundance Rank-abundance
  #graph labeled with the string idenfier

  #plot relative abundance, rank abunduance distributions
  relabun1 = relabund(abundance1)
  relabun2 = relabund(abundance2)
  plot(NA, NA, xlab = '', ylab = '', xlim = c(0, max(c(length(abundance1), length(abundance2)))), 
       ylim = c(0,max(c(max(relabun1),max(relabun2)))), cex.axis = 0.75)
  lines(c(1:length(abundance1)), relabun1, type = 'l', pch = 20, cex = 2, lwd = 2)
  lines(c(1:length(abundance2)), relabun2, type = 'l', pch = 20, cex = 2, lwd = 2, lty = 2, col = 'deeppink3')
  legend('topright', c('control', 'experiment'), bty = 'n', lty = c(1,2), lwd = 2, col = c('black', 'deeppink3'), seg.len = 3, cex = 0.75)
  #legend('topright', c('control', paste('mu, sigma = ', round(cont_par$mu,2), ', ', round(cont_par$sig,2)),
  #                     'experiment', paste('mu, sigma = ', round(expt_par$mu,2),', ', round(expt_par$sig,2))), 
  #      bty = 'n', lty = c(1,NA,2,NA), lwd = 2, col = c('black', NA, 'deeppink3', NA), cex = 0.75)
  mtext('Relative Abundance', side = 2, line = 2, cex = 0.5)
  mtext('Rank', side = 1, line = 2, cex = 0.5)
  mtext(paste(taxa, ': ', siteID1, 'vs.', siteID2, sep = ' '), side = 3, line = .5, cex = 0.75)
}


relabund = function(abundances){
  # get relative abundances from a vector of raw abundances, rounded to 3 decimal places and sorted in decreasing order
  relativeabund = sort((abundances/sum(abundances)), decreasing = TRUE) 
  return (relativeabund)
}


shortestVector = function(v1, v2, maxLen){
  # identify the shorter vector
  if (length(v1) < maxLen){
    return(list(small = v1, large = v2))
  }
  else {
    return(list(large = v1, small = v2))
  }
}


abundMerge = function(r1, r2){
  # put two relative abundances into a matrix together for euclidean distance analysis
  ranks = c(1:max(c(length(r1),length(r2)))) # makes a list of ranks from 1 to the length of the longer vector
  rSizes = shortestVector(r1, r2, max(ranks)) # returns a list of the vectors, labeling each "small" or "large"
  diff = max(ranks) - length(rSizes$small)  # the difference in the number of ranks
  rSizes$small = append(rSizes$small, rep(0,diff)) # adds the difference to the smaller vector as zeroes - so they are the same length
  # put the data (now equal lengths) into a matrix where col 1 = control data and col 2 = experiment data
  matrix = matrix(data = c(as.numeric(unlist(rSizes[1])), as.numeric(unlist(rSizes[2]))), nrow = max(ranks), ncol = 2)
  return (matrix)
}


Euclidean = function(twosite_matrix){
  # find Euclidean distance between the RADs, returns ED rounded to 4 decimal places
  if (ncol(twosite_matrix) == 2){
    ED = sqrt(sum((twosite_matrix[,1] - twosite_matrix[,2]) ^ 2))
  }
  else {
    data = reshape_data(twosite_matrix) 
    data = data[,c(2:ncol(data))]
    data[1,] = data[1,]/sum(data[1,])
    data[2,] = data[2,]/sum(data[2,])
    data2 = t(data)
  ED = sqrt(sum((data2[,1] - data2[,2]) ^ 2))
  }
  return (round(ED, 4))
}


standardize = function(vector){
  # input a vector of values. Takes the mean (meanx) and standard deviation (sdx) of the vector, 
  #    and outputs a new vector containing the standardized values (stdz_vector)
  meanx = sum(vector)/length(vector)
  sdx = sd(vector)
  stdz_vector = as.numeric()
  for (i in 1:length(vector)){
    newVal = (vector[i] - meanx)/sdx
    stdz_vector = append(stdz_vector, newVal)
  }
return (stdz_vector)
}



analyzeSN = function(matrix){
  # proportion difference in S and N, with sign
  s1 = matrix[1,1]
  s2 = matrix[1,2]
  result = s1-s2/max(matrix) #add abs back here to get rid of sign
  return(result)
}



plot1to1 = function (type, controldata, experimentdata, maxpolygon, lolim, uplim) {
  #plots 1:1 plot of data with correct polygon and legend
  #where maxpolygon gives the maximum limit of the polygon to be drawn,
    #lolim gives the lower limit to be used for xlim and ylim,
    #uplim gives the upper limit to be used for xlim and ylim, and
    #type is a character string used for the title of the graph
  plot(NA, NA, log = 'xy', pch = 19, xlim = c(lolim, uplim), ylim = c(lolim, uplim), xlab = '', ylab = '', cex.axis = 0.75)
  #coordinates to draw polygon fro 1.5 x's difference and 2x's difference
  x = c(0:maxpolygon, maxpolygon:0)
  y = c(2*(0:maxpolygon), 0.5*(maxpolygon:0))
  polygon(x, y, border = NA, col = "lightpink")
  y2 = c(1.5*(0:maxpolygon), (maxpolygon:0)/1.5)
  polygon(x, y2, border = NA, col = "lightgoldenrod1")
  abline(0, 1, lty = 2, lwd = 2, col = 'red')
  points(controldata, experimentdata, pch = 19, xlab = 'c', ylab = '', main = '')
  S_r2 = rsquare(controldata, experimentdata)
  rmse = rmse(experimentdata, controldata)
  mae = mae(experimentdata, controldata)
  legend('topleft', c(paste('r2 = ', round(S_r2,3), sep = ''), paste('rmse = ', round(rmse,3), sep = ''), paste('mae = ', round(mae,3), sep = '')), bty = 'n', cex = 0.75)
  mtext('Experiment', side = 2, line = 2, cex = 0.75)
  mtext('Control', side = 1, line = 2, cex = 0.75)
  mtext(type, side = 3, line = .25, cex = 0.75)
}



spp_dat = function(genus, species, dat) { 
  #get ranks of top 5 species, to use inside function "top_control_spp"
  if(nrow(dat)>0){
    sum_abun = sum(dat[,7])
    rank = which(dat$species==species & dat$genus == genus)
    if(length(rank)>0){
      rel_abun = round(dat[rank,7]/sum_abun,4)
      spp_data = list(rank, rel_abun)}
    if(length(rank)==0){
      rank = 0
      rel_abun = 0 
      spp_data = list(rank, rel_abun)}
  }
  return(spp_data)}


top_control_spp = function(control_data, experiment_data) {
  #finds top 5 species in control community, looks at how these species reorganize in experimental communities
  ranks_data = c("co_rank"=1,"exp_rank"=1, "co_relabun" = 1, "exp_relabun" = 1)

  control_data = control_data[order(-control_data[,7]),]
  experiment_data = experiment_data[order(-experiment_data[,7]),]
  #identify top 3 species
  spp = control_data[c(1:5),c(4,5)]  
  for (s in 1:nrow(spp)) {
    #returns list with [1]rank and [2]relative abundance
    c_dat = spp_dat(spp[s,1], spp[s,2], control_data)    
    e_dat = spp_dat(spp[s,1], spp[s,2], experiment_data)
    ranks_data = rbind(ranks_data,c(as.numeric(c_dat[1]), as.numeric(e_dat[1]), as.numeric(c_dat[2]), as.numeric(e_dat[2])))
  }
  ranks_data = ranks_data[-1,]
return(ranks_data)
}


count_RAD_shapes = function (cID, eID, Cshape, Eshape){
# return the number of unique shapes across the unique IDs, still has a few redundancies
  IDs = unique(c(unique(cID), unique(eID)))
  allids = c(cID, eID)
  allshapes = c(Cshape, Eshape)
  shapes = c()
  for (i in 1:length(IDs)){
    index = match(IDs[i], allids)
    shape = allshapes[index]
    shapes = append(shapes, shape)
  }
  nums = table(shapes)
  return (nums)
}


#----- Null Modeling Functions, credit to B. Weinstein (SBU)

nullN<-function(siteXspp){
# The total N (summed across both sites) remains the same, but total N observed at each site
# is allowed to differ due to the species-level randomizations.
  
  #Create an output matrix
  out<-matrix(nrow=nrow(siteXspp),ncol=ncol(siteXspp))
  
  #Draw new abundance distribution
  for (x in 1:ncol(siteXspp)){
    totalN<-sum(siteXspp[,x])
    N1<-sample(0:totalN,1)    #changed 0:totalN, instead of 1:totalN
    N2<-totalN - N1
    out[,x]<-c(N1,N2)
  }
  
  #Compute difference in abundances
  Nboth<-rowSums(out)
  Tstar<-Nboth[[1]] - Nboth[[2]]
  return(Tstar)
}

NullCommunityN<-function(siteXspp){
# input paired communities in site x species matrix, run randomization test
# to determine if difference in N is > than expected by random
# randomizes the abundance of each species in the paired communities, while still assuming
# that the total number of individuals within each species was observed. Note that this may also 
# change observed S for each site and/or form of abundance distribution (not analyzed here).

  #create an output matrix
  Nboth_obs<-rowSums(siteXspp)
  Tstar_obs<-Nboth_obs[[1]] - Nboth_obs[[2]] #OBSERVED: control N - manipulated N
  
  #replicate null distribution, decide the number of randomizations n=X
  nullDistribution<-replicate(n=100,expr=nullN(siteXspp))
  
  #Find quantile of the null distribution for the observed test statistic
  quant<-ecdf(nullDistribution) (Tstar_obs)      
  
  #output the quantile
  if(quant > .95 | quant < .05) {decision<-"Sign."}
  else { decision <- "Random"}
  #if(quant < .95 & quant > .05) {decision<-"Random"}
  
  return(as.list(c(decision, as.numeric(quant))))
}


