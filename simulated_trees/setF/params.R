#--------------------------------------------------
# parameters for sim.coal.tree
#--------------------------------------------------

# Number of trees to generate
nTrees <- 50

# Number of samples from each patient
nSamples <- c(20, 20)

# Infection times, in years
t_inf <- c(0, 1.5/12)

# Sample times, in years
t_sam <- c(1+1/12, 2+1/12)

# Transmission window for recipient, in years
#   time after original infection of R where R could be infected again
#   Inf means no limit, including after sampling
tr_window_1 <- 0.01 # 0 causes failure

# Transmission from donor to recipient (in forward time)
#   if not 0, transmission is not necessarily a single event
qDR <- 0

# Starting effective population size (alpha) for each patient
a <- c(20, 20)
# Linear population growth rates (beta) for each patient, per day
b <- c(3, 3)

#--------------------------------------------------#

# Transmission tree
names <- c("D", "R") # donor and recipient
donor <- c(0, 1)     # 0 means donor is the index case; 1 means the recipient is infected by the donor
tt <- data.frame(names, donor, t_inf, t_sam, stringsAsFactors=F) # make transmission tree into data frame

nInd <- 2 # two individuals if no tree is supplied

tr_window <- rep(tr_window_1, nInd)
rm(tr_window_1)

# Transmission from donor to recipient (in forward time)
#   (reverse-time migration from R to D)
rhoR <- rep(qDR, dim(tt)[1])
rm(qDR)
# Transmission from recipient to donor (in forward time)
#   (reverse-time migration from D to R)
qRD <- 0
rhoD <- rep(qRD, dim(tt)[1])
rm(qRD)
# rhoD values for patients that did not infect anyone else have no effect
# the rhoR value for the index case (the one whose infector is not sampled) has no effect

# set random number generator seed (unless otherwise supplied)
set.seed(444)

# number of digits in tree names (so that output files sort more easily)
tree_name_digits <- nchar(as.character(nTrees))

# whether or not to save the tree files
save_trees <- TRUE

# whether or not to make pdf plots of the trees (defaults to FALSE if not overridden by CL or param file)
plot_trees <- TRUE
