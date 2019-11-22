# This script takes in the transformed output of ablat_tx.R
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

# Check if data is already loaded and read it in if not 
data.exsts <- exists('ablat')
if(data.exsts == F) { 
  ablat <- read.csv('./data/interim/ablat-fixedall-tx.csv')
}
rm(data.exsts)

# If there is no interim data table, run the script to generate it
data.exsts <- exists('ablat')
if(data.exsts == F) { 
  source('./src/data/ablat_tx.R')
}
rm(data.exsts)

# Load functions that will be used in the script
source('./src/functions/icm_spread.R')
source('./src/functions/mkay.R')

# Data is represented x5, for the 5 channels. Slim down to one
ablat <- subset(ablat, Channel == 'CH2')

# Read in unicos table previously generated
unicos <- read.csv('./references/ablat_unicos.csv')

# Extract unique staining patterns and corresponding experiments, 
# as done in ablat_tx.R
if(exists('stains') == F) { 
  stains <- unique(paste(unicos$CH2, unicos$CH3, unicos$CH5))
  exps <- list()
  for(s in 1:length(stains)) { 
    exps[[s]] <- unique(unicos$Experiment[which(paste(unicos$CH2, 
                                                      unicos$CH3, 
                                                      unicos$CH5) == stains[s])])
  }
}

# Split embryos between 
# * those stained with GATA6(gt) and either NANOG antibody, 
# which is rather homogeneous and contain all controls and ablated embryos
# with some littermates
# * those stained with any other combination, which are mixed littermates
# and get classified in new-littermates_classify.R

# Embryos stained with NANOG(rb)
ng.rb <- unicos$Experiment[which(unicos$CH3 == 'NANOG.rb')]
# Embryos stained with NANOG(rat) - only Intact and Random Controls
ng.rat <- unicos$Experiment[which(unicos$CH2 == 'NANOG.rat' & 
                                    unicos$Treatment != 'Littermate')]
## and embryos stained with GATA6(gt)
g6.gt <- unicos$Experiment[which(unicos$CH5 == 'GATA6.gt')]

ablat.a <- subset(ablat, Experiment %in% c(ng.rb, ng.rat) & 
                    Experiment %in% g6.gt)
ablat.b <- subset(ablat, !Experiment %in% unique(ablat.a$Experiment))

# Standardize NANOG and GATA6 values for each litter against litter maxima
# for the purpose of cell classification with both K-means and H-clustering
# This is a rather noisy dataset and this step creates tighter clusters
ablat.a <- split(ablat.a, as.factor(ablat.a$Litter))
for(l in 1:length(ablat.a)) { 
  ablat.a[[l]]$CH3.ebLogCor.xl <- 
    ablat.a[[l]]$CH3.ebLogCor.x / max(ablat.a[[l]]$CH3.ebLogCor.x)
  ablat.a[[l]]$CH5.ebLogCor.xl <- 
    ablat.a[[l]]$CH5.ebLogCor / max(ablat.a[[l]]$CH5.ebLogCor)
}
ablat.a <- do.call(rbind, ablat.a)

################################################################################
# K-means clustering
################################################################################

# Embryos stained for each NANOG antibody (rb and rat) cluster a bit differently
# To process them separately but at once, introduce a new variable and split
# the data frame into two, one for each subgroup
ablat.a$nanog.ab <- 'nada'
ablat.a$nanog.ab[which(ablat.a$Experiment %in% ng.rb)] <- 'rabbity'
ablat.a$nanog.ab[which(ablat.a$Experiment %in% ng.rat)] <- 'ratty'

# Create vectors of identities to assign based on the clusters K-means finds
# Note these may need editing if adding/removing data, as the order of the 
# corresponding clusters would change
id.vectors <- list(c('EPI', 'PRE', 'DP'), 
                   c('PRE', 'EPI'))
# Create a vector of k values for each subset of data (3 for NANOG(rb) and 
# 2 for NANOG(rat), see the graphs generated below)
ks <- c(3, 2)

# Split dataset by nanog.ab and run mkay(), which runs k-means and 
# automatically assigns identities using id.vectors based on the distance
# to cluster centers.
ablat.a <- split(ablat.a, as.factor(ablat.a$nanog.ab))
for(i in 1:length(ablat.a)) { 
  aa <- ablat.a[[i]] %>% filter(TE_ICM != 'TE', 
                            Treatment %in% c('Control', 'Littermate')) %>% 
    select(Cellcount, CH5.ebLogCor.xl, CH3.ebLogCor.xl)
  icm.spread(aa)
  ablat.a[[i]] <- mkay(dataset = ablat.a[[i]], miniset = aa, 
                       k = ks[i], ids = id.vectors[[i]], 
                       x.var = 'CH5.ebLogCor.xl', y.var = 'CH3.ebLogCor.xl')
  }
ablat.a <- do.call(rbind, ablat.a)

# Uncomment below to visualize the result
# qplot(CH5.ebLogCor.xl,  CH3.ebLogCor.xl,
#       data = subset(ablat.a, TE_ICM == 'ICM'),
#       color = Identity.km) + scale_color_manual(values = idcols) +
#   facet_grid(Treatment ~ Stage.t0) + looks + theme(aspect.ratio = 1)

# Note that in this subset of embyros, because they have been grown in culture
# and are thus a bit more delayed, the epiblast maintains fairly robust
# NANOG expression and thus there really isn't a distinct EPI.lo cluster

################################################################################
# Hierarchical clustering
################################################################################

# Separate TE and ICM cells, as in new-lms_tx.R
# We will only perform Hierarchical clustering with the icm subset, 
# a much smaller dataset, which will be faster, and is the only data  
# we need to classify anyway
te <- subset(ablat.a, TE_ICM == 'TE')
icm <- subset(ablat.a, TE_ICM == 'ICM')

# Assign TE cells to a made up cluster 0
te$id.cluster <- 0
te$Identity.hc <- 'TE'

# Vector of k values to use below
k <- c(4, 3)

# Perform hierarchical clustering on all ICM cells using average linkage
# and cut tree at 4 clusters (k[1])
my.clusters <- hclust(dist(data.frame(icm$CH5.ebLogCor.xl, 
                                      icm$CH3.ebLogCor.xl)), 
                      method = 'average')
# Uncomment to see clustering tree
# plot(my.clusters)

# Cut tree at k = 4
icm$id.cluster <- cutree(my.clusters, k[1])

# Repeat only for those embryos stained for GATA6(gt) and NANOG(rat)
# but finding 3 clusters only. 

# Including these embryos in the seeding for the other group generates
# more robust clusters, but the clusters for this second iteration 
# are cleaner when classified separately
my.clusters <- hclust(dist(data.frame(
  icm$CH5.ebLogCor.xl[which(icm$nanog.ab == 'ratty')], 
  icm$CH3.ebLogCor.xl[which(icm$nanog.ab == 'ratty')])), 
  method = 'average')

# Uncomment to see clustering tree
# plot(my.clusters)

# Cut tree at k = 3
icm$id.cluster[which(icm$nanog.ab == 'ratty')] <- cutree(my.clusters, k[2])

## Uncomment below to see outcome
# qplot(CH5.ebLogCor,  CH3.ebLogCor.x,
#       data = icm, color = id.cluster) +
#   looks + scale_color_gradient2(low = 'black', mid = 'green',
#                                 high = 'yellow', midpoint = (k[1]+1)/2) +
#   facet_grid(Treatment ~ Stage.t0) + theme(aspect.ratio = 1)

# Make data frames with correspondence between clusters and identity
# and merge with main data frame to assign identities based on cluster
icm <- split(icm, as.factor(icm$nanog.ab))

idxclust <- list(data.frame(id.cluster = c(1:k[1]), 
                            Identity.hc = c('PRE', 'DP', 'EPI', 'EPI.lo')),
                 data.frame(id.cluster = c(1:k[2]), 
                            Identity.hc = c('PRE', 'EPI', 'DN')))
icm[[1]] <- merge(icm[[1]], idxclust[[1]])
icm[[2]] <- merge(icm[[2]], idxclust[[2]])
icm <- do.call(rbind, icm)

# Re generate ablat.a from TE and ICM cells
ablat.a <- rbind(te, icm)
rm(te, icm)
# Drop NANOG antibody variable
ablat.a$nanog.ab <- NULL

# Assign levels to Identity.hcf
# Note that I scored DN cells as EPI.lo for embryo stained with NANOG(rb), 
# but as DN for those stained with NANOG(rat). This is purely a judgment call. 
# In the first group, cluster 4 seems quite tight and close to the EPI, so
# I am considering them low-NANOG EPI cells, whereas in the second group, 
# cluster 3 is small and more spread out than in the previous case, 
# and some of the cells fall near PrE cells. Therefore I'm being conservative 
# and assuming they may be a mix of EPI low and DN (thus DN overall)
ablat.a$Identity.hc <- factor(ablat.a$Identity.hc, 
                              levels = c('TE', 'PRE', 'DP', 'EPI', 
                                         'EPI.lo', 'DN'))

# Uncomment below to see the outcome
## Visualize the result
# qplot(CH5.ebLogCor,  CH3.ebLogCor.x,
#       data = subset(ablat.a, TE_ICM == 'ICM'), color = Identity.hc) +
#   looks + scale_color_manual(values = idcols) +
#   facet_grid(Treatment ~ Stage.t0) + theme(aspect.ratio = 1)

################################################################################
# Incorporate data for littermates that were processed with new-littermates
# when running new-lms_tx.R and were stained for NANOG and GATA6
################################################################################

new.lms.red <- read.csv('./data/processed/new-lms-processed.csv')
# Select only embryos in new-littermates.A which are in ablat.b
new.lms.red <- subset(new.lms.red, Embryo_ID %in% unique(ablat.b$Embryo_ID))
# Rename column litter.median back to group.median
new.lms.red <- rename(new.lms.red, group.median = litter.median)
# Slim down the data and merge with ablat.b to incorporate identity info
# as well as transformed CH3 and CH5 data
new.lms.red <- new.lms.red %>% filter(Channel == 'CH2') %>% 
  select(Embryo_ID, Cell_ID, Identity.km, 
         Identity.hc, id.cluster, 
         CH1.ebLogCor.s, CH2.ebLogCor.s, 
         CH3.ebLogCor.s, CH5.ebLogCor.s, 
         CH3.ebLogCor.x,  CH3.ebLogCor.xs, 
         CH3.ebLogCor.xl, CH5.ebLogCor.xl)
ablat.b[which(colnames(ablat.b) %in% c('CH3.ebLogCor.x', 
                                    'CH1.ebLogCor.s', 
                                    'CH2.ebLogCor.s', 
                                    'CH3.ebLogCor.s', 
                                    'CH5.ebLogCor.s', 
                                    'CH3.ebLogCor.xs'))] <- NULL
ablat.b <- merge(ablat.b, new.lms.red)

# Combine ablat.a and ablat.b back into ablat
ablat <- rbind(ablat.a, ablat.b)

# Record number of embryos at end
n.end <- unique(ablat$Embryo_ID)

print(paste('N embryos loaded at beginning', length(n.start), sep = ": "))
print(paste('N embryos after processing', length(n.end), sep = ": "))

################################################################################
# Write out data to the ./data/processed folder
# All embryos
write.csv(ablat, file = './data/processed/ablat-processed.csv', row.names = F)
# processed littermates only
write.csv(subset(ablat, Treatment == 'Littermate'), 
          file = './data/processed/lms-ablat-processed.csv', row.names = F)
