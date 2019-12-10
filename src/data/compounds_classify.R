# This script takes in the transformed output of compounds_tx.R
# and classifies ICM cells in each embryo using two alternative approaches:
# * K-means clustering, as previously done in Saiz et al., (2016) and 
#   Morgani et al., (2018).
# * Hierarchical clustering, which we favor in the current study.

# Both methods are used for the sake of comparison. They produce similar results
# however, K-means tends to overestimate the size of the DP and DN populations,
# due to standard behavior of the algorithm - it tries to make clusters of 
# equal size. This is problematic, as the size of the DP compartment decreases
# with embryo development, and the DN is mostly present in very late stage 
# embryos, where it captures epiblast cells that have lost NANOG expression.

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Check if data is already loaded 
data.exsts <- exists('compos')
if(data.exsts == F) { 
  compos <- read.csv('./data/interim/compounds-tx.csv')
  compos.ref <- rbind(read.csv('./references/compounds_exp_ref.csv'))
}
rm(data.exsts)

# If raw data doesn't exist, run script to generate it
data.exsts <- exists('compos')
if(data.exsts == F) { 
  source('./src/data/compounds_tx.R')
  compos.ref <- rbind(read.csv('./references/compounds_exp_ref.csv'))
}
rm(data.exsts)

# Load functions that will be used in the script
source('./src/functions/icm_spread.R')
source('./src/functions/mkay.R')

# Slim down data by selecting only values for Channel == 'CH2'
# (all data is repeated x5, once per channel)
compos <- subset(compos, Channel == 'CH2')

# Read in unicos table previously generated
unicos <- read.csv('./references/compounds_unicos.csv')

# Create a vector with unique combinations of antibodies and stainings
# and a list with the corresponding experiments stained that way
if(exists('stains') == F) { 
  stains <- unique(paste(unicos$CH2, unicos$CH3, unicos$CH5))
  exps <- list()
  for(s in 1:length(stains)) { 
    exps[[s]] <- unique(unicos$Experiment[which(paste(unicos$CH2, 
                                                      unicos$CH3, 
                                                      unicos$CH5) == stains[s])])
  }
}

################################################################################
# K-means clustering
################################################################################

# Embryos stained for each NANOG antibody (rb and rat) cluster differently
# To process them separately but in parallel, introduce a new variable and split
# the data frame into two, one for each subgroup
compos$nanog.ab <- 'nada'
# The first element of the exps list should contain experiments
# that were stained for NANOG.rat (corresponding to stains[1])
compos$nanog.ab[which(compos$Experiment %in% exps[[1]])] <- 'ratty'
compos$nanog.ab[which(compos$Experiment %in% exps[[2]])] <- 'rabbity'

# Split dataset by nanog.ab
compos <- split(compos, as.factor(compos$nanog.ab))

# Create vectors of identities to assign based on the clusters K-means finds
# Note these may need editing if adding/removing data, as the order of the 
# corresponding clusters would change
id.vectors <- list(c('EPI', 'DN', 'PRE', 'DP'), 
                   c('PRE', 'EPI', 'DN'))
# Create a vector of k values for each subset of data (3 for both)
# For the rabbit antibody we'll need to add a DP cluster using mkay()
ks <- c(3, 3)

# Run mkay() on the compos list, which runs k-means and 
# automatically assigns identities using id.vectors based on the distance
# to cluster centers.
for(i in 1:length(compos)) { 
  cc <- compos[[i]] %>% filter(TE_ICM != 'TE') %>% 
    select(Cellcount, CH5.ebLogCor, CH3.ebLogCor.x)
  icm.spread(cc)
  if(unique(compos[[i]]$nanog.ab == 'ratty')) {
    compos[[i]] <- mkay(dataset = compos[[i]], miniset = cc, 
                        k = ks[i], ids = id.vectors[[i]], 
                        x.var = 'CH5.ebLogCor', y.var = 'CH3.ebLogCor.x')
  }
  if(unique(compos[[i]]$nanog.ab == 'rabbity')) { 
    compos[[i]] <- mkay(dataset = compos[[i]], miniset = cc, 
                        k = ks[i], ids = id.vectors[[i]], DP = T, 
                        x.var = 'CH5.ebLogCor', y.var = 'CH3.ebLogCor.x')
    }
}
compos <- do.call(rbind, compos)

# Uncomment below to visualize the result
# qplot(CH5.ebLogCor,  CH3.ebLogCor.x,
#       data = subset(compos, TE_ICM == 'ICM'),
#       color = Identity.km) + scale_color_manual(values = idcols) +
#   facet_grid(Treatment ~ Stage) + looks + theme(aspect.ratio = 1)

# Note the large size of the DN cluster, and how it encompasses 
# cells closer to the EPI as well as cells from the PrE cluster

################################################################################
# Hierarchical clustering
################################################################################

# Standardize NANOG and GATA6 values for each litter against litter maxima
# for the purpose of cell classification using Hierarchical clustering only.
# This is a rather noisy dataset and this step creates tighter clusters
compos <- split(compos, as.factor(compos$Litter))
for(l in 1:length(compos)) { 
  compos[[l]]$CH3.ebLogCor.xl <- 
    compos[[l]]$CH3.ebLogCor.x / max(compos[[l]]$CH3.ebLogCor.x)
  compos[[l]]$CH5.ebLogCor.xl <- 
    compos[[l]]$CH5.ebLogCor / max(compos[[l]]$CH5.ebLogCor)
}
compos <- do.call(rbind, compos)

# Separate TE and ICM cells, as in other datasets
# We will only perform Hierarchical clustering with the icm subset, 
# a much smaller dataset, which will be faster, and is the only data  
# we need to classify anyway
te <- subset(compos, TE_ICM == 'TE')
icm <- subset(compos, TE_ICM == 'ICM')

# Assign TE cells to a made up cluster 0
te$id.cluster <- 0
te$Identity.hc <- 'TE'

# Split icm cells into the two staining groups, as above
# The classification is much cleaner, as there are large differences
# between both groups, despite the transformation of NANOG values
ee <- list(subset(icm, Experiment %in% exps[[1]]), 
           subset(icm, Experiment %in% exps[[2]]))

# Vector of k values to use below
ks <- c(5, 6)

# Perform Hierarchical clustering on each subset of data
for(e in 1:length(ee)) {
  my.clusters <- hclust(dist(data.frame(ee[[e]]$CH5.ebLogCor.xl,
                                        ee[[e]]$CH3.ebLogCor.xl)),
                        method = 'average')
  # Uncomment below to visualize dendogram
  # plot(my.clusters)
  k <- ks[e]
  ee[[e]]$id.cluster <- cutree(my.clusters, k)
  my.table <- table(ee[[e]]$id.cluster, ee[[e]]$Identity.km)
  # Print table comparing H-clustering vs K-means clustering
  print(my.table)
  # Uncomment below to see output
  # my.plot <- qplot(CH5.ebLogCor.s,  CH3.ebLogCor.xs,
  #                  data = ee[[e]], color = id.cluster) +
  #   looks + scale_color_gradient2(low = 'black', mid = 'green',
  #                                 high = 'yellow', midpoint = (k+1)/2) +
  #   facet_grid(Genotype1 ~ Stage) + theme(aspect.ratio = 1)
  # print(my.plot)
}

# Generate vectors of identities for each cluster
idxclust <- list(data.frame(id.cluster = 1:ks[1], 
                            Identity.hc = c('PRE', 'EPI', 'PRE', 
                                            'EPI.lo', 'EPI.lo')), 
                 data.frame(id.cluster = 1:ks[2], 
                            Identity.hc = c('EPI', 'DP', 'PRE', 
                                            'EPI', 'EPI.lo', 'DN')))

# Assign identity to each cluster as per idxclust
for(e in 1:length(ee)) { 
  ee[[e]] <- merge(ee[[e]], idxclust[[e]])
}

# Combine list again into icm and combine TE and ICM cells
icm <- do.call(rbind, ee)
compos <- rbind(te, icm)
rm(te, icm)

compos$Identity.hc <- factor(compos$Identity.hc, 
                              levels = c('TE', 'PRE', 'DP', 'EPI', 
                                         'EPI.lo', 'DN'))

# Uncomment below to see the outcome
qplot(CH5.ebLogCor,  CH3.ebLogCor.x,
      data = subset(compos, TE_ICM == 'ICM'), color = Identity.hc) +
  looks + scale_color_manual(values = idcols) +
  facet_grid(Treatment ~ Stage) + theme(aspect.ratio = 1)

################################################################################
# Write out data to the ./data/processed folder
write.csv(compos, file = './data/processed/compounds-processed.csv', 
          row.names = F)

##