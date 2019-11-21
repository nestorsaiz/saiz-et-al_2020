# This script takes in the transformed output of new-lms_tx.R
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
data.exsts <- exists('new.lms')
if(data.exsts == F) { 
  new.lms <- read.csv('./data/interim/new-lms-tx.csv')
  new.lms.ref <- rbind(read.csv('./references/new-littermates_exp_ref.csv'), 
                       read.csv('./references/fgf4-lms_exp_ref.csv'))
}
rm(data.exsts)

# Load functions that will be used in the script
source('./src/functions/icm_spread.R')
source('./src/functions/mkay.R')

# Select only embryos with 32 or more cells (blastocysts)
# and only data for Channel == 'CH2' (all rows are present x5, one per channel)
new.lms <- subset(new.lms, Cellcount > 31 & Channel == 'CH2')

# Read in unicos table previously
unicos <- read.csv('./references/new-lms_unicos.csv')

# Extract unique staining patterns and corresponding experiments, 
# as done in new-lms_tx.R
if(exists('stains') == F) { 
  stains <- unique(paste(unicos$CH2, unicos$CH3, unicos$CH5))
  exps <- list()
  for(s in 1:length(stains)) { 
    exps[[s]] <- unique(unicos$Experiment[which(paste(unicos$CH2, 
                                                      unicos$CH3, 
                                                      unicos$CH5) == stains[s])])
  }
}

# Split dataset dataset into:
# * Embryos stained with both NANOG and GATA6 (either antibody)
#   which will be classified using K-means and Hierarchical clustering
nanogata <- subset(new.lms, 
                   Experiment %in% c(exps[[1]], exps[[2]], exps[[3]], 
                                     exps[[4]], exps[[7]]))
# * Embryos stained with NANOG and SOX17 (no GATA6) and
# * embryos stained with NANOG and GATA4 (no GATA6)
# Embryos in these two last subsets also express Pdgfra:H2B-GFP. 
# Since they represent a rather small set of observations, 
# they will be classified using thresholding, instead of K-means clustering
# and also using Hierarchical clustering
s17 <- subset(new.lms, Experiment %in% exps[[6]])
g4 <- subset(new.lms, Experiment %in% exps[[5]])

################################################################################
# K-means clustering
################################################################################

# Define the miniset of the data that will be used to perform K-means
# (cell counts, NANOG and GATA6 intensity values)
ng <- nanogata %>% filter(TE_ICM == 'ICM') %>% 
  select(Cellcount, CH5.ebLogCor.xs, CH3.ebLogCor.xs)

# Visualize distribution of NANOG and GATA6 values in ICM cells
icm.spread(ng)

# Define a vector of identities - this should run on the provided dataset 
# but may need to be changed, if new data is added. 
# Assess the order of centers printed out, re-order the id.vector and re-run
# lines 61-86
id.vector <- c('EPI', 'DP', 'PRE', 'DN')

# Run mkay() to perform k-means and classify cells according to id.vector
nanogata <- mkay(dataset = nanogata, miniset = ng, 
                  k = 4, ids = id.vector, 
                  x.var = 'CH5.ebLogCor.xs', y.var = 'CH3.ebLogCor.xs')

# Uncomment below to visualize the result
# qplot(CH5.ebLogCor.xs,  CH3.ebLogCor.xs,
#       data = subset(nanogata, TE_ICM == 'ICM'),
#       color = Identity.km) + scale_color_manual(values = idcols) +
#   facet_grid(Genotype1 ~ Stage) + looks + theme(aspect.ratio = 1)

################################################################################
# Embryos not stained for both NANOG and GATA6
################################################################################

# Assign lineage identity to s17 using thresholding
# based on linear-scale values of GFP (CH2) and NANOG (rb, CH3) 
# note the inversion of the log-corrected variable here and below, for g4
s17$Identity.km <-
  ifelse(s17$TE_ICM == 'TE', 'TE',
         ifelse(exp(s17$CH3.ebLogCor) > 250 &
                  exp(s17$CH2.ebLogCor) < 130,
                'EPI',
                ifelse(exp(s17$CH3.ebLogCor) > 250 &
                         exp(s17$CH2.ebLogCor) > 130,
                       'DP',
                       ifelse(exp(s17$CH3.ebLogCor) < 250 &
                                exp(s17$CH2.ebLogCor) > 130,
                              'PRE', 'DN'))))

# Likewise, assign identity to g4 using thresholding
# based on linear-scale values of GFP (CH2) and NANOG (rat, CH5) 
g4$Identity.km <-
  ifelse(g4$TE_ICM == 'TE', 'TE',
         ifelse(exp(g4$CH5.ebLogCor) > 250 &
                  exp(g4$CH2.ebLogCor) < 130,
                'EPI',
                ifelse(exp(g4$CH5.ebLogCor) > 250 &
                         exp(g4$CH2.ebLogCor) > 130,
                       'DP',
                       ifelse(exp(g4$CH5.ebLogCor) < 250 &
                                exp(g4$CH2.ebLogCor) > 130,
                              'PRE', 'DN'))))

# Note the variable is still named Identity.km, even though it is 
# calculated using thresholding. This is purely semantic, I am just
# trying to create too many columns

# Regenerate new.lms from the three subsets of data
new.lms <- rbind(nanogata, s17, g4)
rm(nanogata, s17, g4)

################################################################################
# Hierarchical clustering
################################################################################

# Standardize NANOG and GATA6 values for each litter against litter maxima
# for the purpose of cell classification using Hierarchical clustering only.
# This is a rather noisy dataset and this step creates tighter clusters
new.lms <- split(new.lms, as.factor(interaction(new.lms$Background, 
                                                new.lms$Litter, 
                                                new.lms$Gene1)))
for(l in 1:length(new.lms)) { 
  new.lms[[l]]$CH3.ebLogCor.xl <- 
    new.lms[[l]]$CH3.ebLogCor.x / max(new.lms[[l]]$CH3.ebLogCor.x)
  new.lms[[l]]$CH5.ebLogCor.xl <- 
    new.lms[[l]]$CH5.ebLogCor.x / max(new.lms[[l]]$CH5.ebLogCor.x)
}
new.lms <- do.call(rbind, new.lms)

# Separate TE and ICM cells. TE cells are manually assigned to cluster 0,
# with Identity.hc == 'TE'
te <- subset(new.lms, TE_ICM == 'TE')
icm <- subset(new.lms, TE_ICM == 'ICM')

te$id.cluster <- 0
te$Identity.hc <- 'TE'

# Split dataset per staining groups, as I found the clustering result is cleaner
# that way (see notebooks). Somehow the population structure is highly affected
# by this, in ways that perhaps K-means is unsensitive to
# Even though groups 1 and 4 (stains[1] and stains[4]) are both 
# stained for NANOG.rb and GATA6.gt, somehow the clustering is cleaner
# if done separately. I also separate embryos from the Fgf4 dataset
ee <- list(subset(icm, Experiment %in% exps[[1]] & Gene1 != 'Fgf4'),
           subset(icm, Experiment %in% c(exps[[2]], 
                                         exps[[7]]) & Gene1 != 'Fgf4'), 
           subset(icm, Experiment %in% exps[[3]] & Gene1 != 'Fgf4'), 
           subset(icm, Experiment %in% exps[[4]] & Gene1 != 'Fgf4'), 
           subset(icm, Experiment %in% exps[[5]] & Gene1 != 'Fgf4'), 
           subset(icm, Experiment %in% exps[[6]] & Gene1 != 'Fgf4'), 
           subset(icm, Gene1 == 'Fgf4'))

# Set the k values to cut the trees (these have been determined empirically, 
# see notebooks for details)
ks <- c(5, 6, 5, 5, 4, 6)

# Loop through the groups, running hclust() on the corresponding channels
# (different for ee[[1:4]], e[[5]] and ee[[6]]), cutting the tree at the
# corresponding k value and printing a scatterplot for each to visualize
# the outcome.
for(e in 1:length(ee)) {
  if(e %in% 1:4) {
    my.clusters <- hclust(dist(data.frame(ee[[e]]$CH5.ebLogCor.xl,
                                          ee[[e]]$CH3.ebLogCor.xl)),
                          method = 'average')
    k <- ks[e]
    ee[[e]]$id.cluster <- cutree(my.clusters, k)
    my.table <- table(ee[[e]]$id.cluster, ee[[e]]$Identity.km)
    # Print table comparing H-clustering vs K-means clustering
    print(stains[e])
    print(my.table)
    # plot(my.clusters)
    # Uncomment below to see output
    # Note that I am plotting the .xs channel (on a 0-1 scale, for all embryos),
    # not the channel standardized for each litter (.xl)
    # my.plot <- qplot(CH5.ebLogCor.xs,  CH3.ebLogCor.xs,
    #                  data = ee[[e]], color = id.cluster) +
    #   looks + scale_color_gradient2(low = 'black', mid = 'green',
    #                                 high = 'yellow', midpoint = (k+1)/2) +
    #   facet_grid(Genotype1 ~ Stage) + theme(aspect.ratio = 1)
    # print(my.plot)
  }
  if(e == 5) {
    my.clusters <- hclust(dist(data.frame(ee[[e]]$CH2.ebLogCor.s,
                                          ee[[e]]$CH5.ebLogCor.s)),
                          method = 'average')
    k <- ks[e]
    ee[[e]]$id.cluster <- cutree(my.clusters, k)
    my.table <- table(ee[[e]]$id.cluster, ee[[e]]$Identity.km)
    ## Print table comparing H-clustering vs K-means clustering
    print(stains[e])
    print(my.table)
    # plot(my.clusters)
    # Uncomment below to see output
    # my.plot <- qplot(CH2.ebLogCor.s,  CH5.ebLogCor.s,
    #                  data = ee[[e]], color = id.cluster) +
    #   looks + scale_color_gradient2(low = 'black', mid = 'green',
    #                                 high = 'yellow', midpoint = (k+1)/2) +
    #   facet_grid(Genotype1 ~ Stage) + theme(aspect.ratio = 1)
    # print(my.plot)
  }
  if(e == 6) {
    my.clusters <- hclust(dist(data.frame(ee[[e]]$CH2.ebLogCor.s,
                                          ee[[e]]$CH3.ebLogCor.s)),
                          method = 'average')
    k <- ks[e]
    ee[[e]]$id.cluster <- cutree(my.clusters, k)
    my.table <- table(ee[[e]]$id.cluster, ee[[e]]$Identity.km)
    ## Print table comparing H-clustering vs K-means clustering
    print(stains[e])
    print(my.table)
    # plot(my.clusters)
    # Uncomment below to see output
    # my.plot <- qplot(CH2.ebLogCor.s,  CH3.ebLogCor.s,
    #                  data = ee[[e]], color = id.cluster) +
    #   looks + scale_color_gradient2(low = 'black', mid = 'green',
    #                                 high = 'yellow', midpoint = (k+1)/2) +
    #   facet_grid(Genotype1 ~ Stage) + theme(aspect.ratio = 1)
    # print(my.plot)
  }
  if(e == 7) { ee[[e]]$Identity.hc <- 'n/a'}
}

# Assign identities based on resulting clusters
idxclust <- list(data.frame(id.cluster = 1:ks[1], 
                            Identity.hc = c('DP', 'PRE', 'EPI', 
                                            'EPI.lo', 'DN')), 
                 data.frame(id.cluster = 1:ks[2],
                            Identity.hc = c('EPI', 'PRE', 'PRE',
                                            'EPI.lo', 'DP', 'DN')),
                 data.frame(id.cluster = 1:ks[3], 
                            Identity.hc = c('EPI', 'PRE', 'DP', 
                                            'DN', 'DN')), 
                 data.frame(id.cluster = 1:ks[4], 
                            Identity.hc = c('PRE', 'EPI.lo', 'EPI', 
                                            'DP', 'DN')), 
                 data.frame(id.cluster = 1:ks[5], 
                            Identity.hc = c('EPI.lo', 'EPI', 'PRE', 'EPI')), 
                 data.frame(id.cluster = 1:ks[6], 
                            Identity.hc = c('EPI', 'PRE', 'DN', 
                                            'DP', 'DP', 'EPI.lo')), 
                 data.frame(id.cluster = 0, 
                            Identity.hc = 'n/a'))

for(e in 1:length(ee)) { 
  ee[[e]] <- merge(ee[[e]], idxclust[[e]])
}

# Combine all four groups again
icm <- do.call(rbind, ee)
rm(idxclust, ee, my.clusters)

# Combine TE and ICM back into the main table
new.lms <- rbind(te, icm)

new.lms$Identity.hc <- factor(new.lms$Identity.hc, 
                              levels = c('TE', 'PRE', 'DP', 'EPI', 
                                         'EPI.lo', 'DN', 'n/a'))

