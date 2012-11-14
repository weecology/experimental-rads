# Looks at rodent movement
wd = "C://Documents and Settings//sarah//My Documents//Active Research Projects//Rodent Movement"
setwd(wd)
source = ("movement_fxns.R")

PB = read.csv("PB_2000-2009.csv")
PP = read.csv("PP_2000-2009.csv")
small_rodents = rbind(PB, PP)

#add treatment
small_rodents = trmt_to_rats(small_rodents)

# make sure text is text and not atomic
small_rodents$tag = as.character(small_rodents$tag) 
small_rodents$sex = as.character(small_rodents$sex)
small_rodents$species = as.character(small_rodents$species)

# give untagged indivs unique tag numbers (7 digits)
small_rodents = id_unknowns(small_rodents)

# get list of unique tags
tags = unique(small_rodents$tag)   

# output list of flagged data
flags = find_bad_data(small_rodents, tags)

# get list of unique "bad tags"
badtags=unique(flags$tag)

# delete bad tags from dataset for analysis
for (i in 1:length(badtags)) {
  small_rodents = subset(small_rodents, tag != badtags[i])
  }
  
# don't use periods with negative period numbers and periods with only one day of trapping
small_rodents = subset(small_rodents, period != 267 & period != 277 & period != 278 &
                                      period != 283 & period != 284 & period != 300 &
                                      period != 311 & period != 313 & period != 314 &
                                      period != 318 & period != 321 & period != 323 &
                                      period != 337 & period != 339 & period != 344 &
                                      period != 351)
  

### Create a set of capture histories by treatment and by plot
tags = unique(sort(small_rodents$tag))
periods = unique(sort(small_rodents$period))

mark_trmt = create_trmt_hist(small_rodents, tags, periods)
mark_plot = create_plot_hist(small_rodents, tags, periods)  

# get list of indivs that moved plots or treatment, species is included
moving_rats = find_rats_that_move(small_rodents, tags)

# get rodents that move plots
moving_tags = unique(moving_rats$tag)
outcount = 0
MARK_movers = data.frame("ch"=1, "censored"=1, "tag"=1, "spp"=1, "sex"=1, "mass"=1)
for (i in 1:length(moving_tags)) {
  mover = subset(mark_plot, tags == moving_tags[i])
  if (nrow(mover) > 0) {
  outcount = outcount + 1
  MARK_movers[outcount,] <- mover
  }}   

#count non-movers
moving_tags = unique(moving_plot$tag)
nonmover = mark_plot
for (i in 1:length(moving_tags)) {
  nonmover = subset(nonmover, tags != moving_tags[i])
}   

# find num captures/rat
rat_catches=num_captures(small_rodents, tags)

# plot captures to see if MARK is a good way to look at these
pdf(file="recaptures.pdf",10,10)
par(mfrow=c(1,1))

hist(rat_catches$captures, labels=TRUE, xlab = "# captures / rodent", ylab = "frequency",
      main = "Rodent recaptures", col = "orange3", breaks=31)

hist(rat_catches[which(rat_catches$spp=="PB"),4], labels=TRUE, xlab = "# captures / rodent", ylab = "frequency",
      main = "PB recaptures", col="plum4", breaks=31)
      
hist(rat_catches[which(rat_catches$spp=="PP"),4], labels=TRUE, xlab = "# captures / rodent", ylab = "frequency",
      main = "PP recaptures", col="cadetblue", breaks=18, xlim=c(0,20), ylim=c(0,2000))

dev.off()

# look at juveniles?

# how often are movers moving trmt?
trmt_moves = examine_trmt_moves(small_rodents, moving_tags)
  trmt_moves = trmt_moves[which(trmt_moves$num_moves > 0),]
  trmt_moves$num_moves = as.numeric(trmt_moves$num_moves)
  trmt_moves$c2r = as.numeric(trmt_moves$c2r)
  trmt_moves$r2c = as.numeric(trmt_moves$r2c)
  trmt_moves$c2e = as.numeric(trmt_moves$c2e)
  trmt_moves$r2e = as.numeric(trmt_moves$r2e)
plot_moves = examine_plot_moves(small_rodents, moving_tags)


hist(trmt_moves[,7], labels = TRUE)
hist(trmt_moves[,8], labels = TRUE)
hist(trmt_moves[,9], labels = TRUE)
hist(trmt_moves[,10], labels = TRUE)

M = trmt_moves[which(trmt_moves[,3]=="M"),]
F = trmt_moves[which(trmt_moves[,3]=="F"),]

PB = trmt_moves[which(trmt_moves[,2]=="PB"),]
PP = trmt_moves[which(trmt_moves[,2]=="PP"),]

PBmale = subset(M, sp == "PB")
PBfemale = subset(F, sp == "PB")
PPmale = subset(M, sp == "PP")
PPfemale = subset(F, sp == "PP")

hist(PB[,7], labels=T, main = "C - R")
hist(PB[,8], labels=T, main = "R - C")
hist(PP[,7], labels=T, main = "C - R")
hist(PP[,8], labels=T, main = "R - C")