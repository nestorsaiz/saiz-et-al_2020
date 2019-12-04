# This script takes in the transformed output of emb-xim_tx.R
# it classifies GFP+ and GFP- cells automatically,
# and it classifies ICM cells in each embryo using two alternative approaches, 
# as in the other *_classify.R scripts:
# * K-means clustering, as previously done in Saiz et al., (2016) and 
#   Morgani et al., (2018).
# * Hierarchical clustering, which we favor in the current study.

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Check if data is already loaded and read it in if not
data.exsts <- exists('g6.chimeras')
if(data.exsts == F) { 
  g6.chimeras <- read.csv('./data/interim/emb-xim-tx.csv')
  g6.ref <- read.csv('./references/emb-xim_exp_ref.csv')
}
rm(data.exsts)

# If there is no interim data, run script to generate it
data.exsts <- exists('g6.chimeras')
if(data.exsts == F) { 
  source('./src/data/emb-xim_tx.R')
  g6.ref <- read.csv('./references/emb-xim_exp_ref.csv')
}
rm(data.exsts)

# Load functions that will be used in the script
source('./src/functions/icm_spread.R')
source('./src/functions/mkay.R')

# Load some extra packages
library('Ckmeans.1d.dp')

# Slim down dataset, as we don't need replicates for each channel
g6.chimeras <- subset(g6.chimeras, Channel == 'CH2' & 
                        H_genotype != 'unknown')

################################################################################
# Classify GFP+ and GFP- cells automatically
################################################################################

# Cells will be classified as host and donor, which are somewhat arbitrary
# categories I do not use in the paper.
# * Host cells generally refer to unlabelled (GFP-) cells, either wt or Gata6-/-
#   however, in non-chimeric control embryos, I classify all cells as 'host'
#   regardless of whether they have GFP or not, to facilitate data management.
# * Donor cells refer to GFP+ wild type cells

# Split dataset into Control embryos (100% host) and chimeras (mix)
cc <- subset(g6.chimeras, Treatment == 'Control')
xs <- subset(g6.chimeras, Treatment == 'Chimera')

# Introduce control to check that all embryos are accounted for
check <- length(g6.chimeras$Embryo_ID) - 
  (length(cc$Embryo_ID) + length(xs$Embryo_ID))
if(check != 0) {
  print('SOMETHING WENT WRONG')
}

# Classify all control cells as 'host'
cc$cell_type <- 'host'

## Split xs by litter and convert into a list
xs <- split(xs, as.factor(xs$Litter))

# Perform hierarchical clustering to classify cells in chimeras
# using the levels of GFP and the levels of Hoechst as a second variable.
# I do it for each litter individually to account for technical variation
# between litters, while dampening the noise found in individual embryos.
# I tried doing a global classification too (all embryos at once) and
# per-embryo classification, but both performed worse than doing per-litter. 
for(l in 1:length(xs)) { 
  
  my.clusters <- hclust(dist(data.frame(xs[[l]]$CH1.ebLogCor,
                                        xs[[l]]$CH2.logCor)),
                        method = 'average')
  xs[[l]]$cluster <- cutree(my.clusters, 2)
  
  bigger <- mean(subset(xs[[l]], cluster == 2)$CH2.logCor) - 
    mean(subset(xs[[l]], cluster == 1)$CH2.logCor) > 0
  if(bigger == T) {
    xs[[l]]$cell_type[xs[[l]]$cluster == 2] <- 'donor'
    xs[[l]]$cell_type[xs[[l]]$cluster == 1] <- 'host'
  }
  else {
    xs[[l]]$cell_type[xs[[l]]$cluster == 2] <- 'host'
    xs[[l]]$cell_type[xs[[l]]$cluster == 1] <- 'donor'
  }
  # To plot the dendrogram for each litter, uncomment the following
  # plot(my.clusters)
}

# Below is an alternative, method, using 1D K-means clustering
# using the Ckmeans.1d.dp package. It does a good job, 
# but Hierarchical clustering performs better

# # For each litter, extracth the values of CH2.logCor,
# # apply Ckmeans.1d.dp with a fixed k = 2 clusters (GFP+ and GFP-),
# # and classify cells based on the clustering
# for(l in 1:length(xs)) {
#         # Extract CH2 levels, which show a bimodal distribution
#         ss <- xs[[l]]$CH2.logCor
#         # Classify CH2 values in ss 
#         ssc <- Ckmeans.1d.dp(ss, 2)
#         # Assign cell types based on ssc
#         xs[[l]]$cell_type <- c('host', 'donor')[ssc$cluster]
# }

# Rbind all litters again into xs
xs <- do.call(rbind, xs)
# Drop the cluster variable
xs$cluster <- NULL

# # I also tried k-means clustering using Hoechst and GFP to classify cells
# # but does not perform as well as Hierarchical clustering
# bb <- xs %>% select(Cellcount, 
#                     CH1.ebLogCor, CH2.logCor)
# kk <- kmeans(bb[, 2:3], 2)
# kk$centers
# ssq <- matrix(0, length(xs$TE_ICM), 2)
# for(i in 1:2) {
#         ssq[,i] <- (xs$CH1.ebLogCor - kk$centers[i,1])^2 +
#                 (xs$CH2.logCor - kk$centers[i,2])^2
# }
# min.ssq <- apply(ssq, 1, which.min)
# xs$cell_type.km <- rep(NA, nrow(xs))
# xs$cell_type.km <- c('donor', 'host')[min.ssq]

# Combine both datasets (cc and xs) into g6.chimeras
g6.chimeras <- rbind(cc, xs)
rm(cc, xs)
# Make cell_type a factor and fix its levels
g6.chimeras$cell_type <- factor(g6.chimeras$cell_type, 
                                levels = c('host', 'donor'))

################################################################################
# Assign lineage identity using K-means clustering
################################################################################

# Define the miniset of the data that will be used to perform K-means
# (cell counts, NANOG and GATA6 intensity values).
# I exclude values from Gata6-/- cells, as they have a lot of NANOG-low
# cells and pull the EPI cluster away from its center in wild types
gg <- g6.chimeras %>% filter(TE_ICM == 'ICM', 
                             !H_genotype %in% c('ko', 'mosaic')) %>% 
  select(Cellcount, CH5.ebLogCor, CH3.ebLogCor)

# Visualize distribution of NANOG and GATA6 values in ICM cells
icm.spread(gg)

# Define a vector of identities - this should run on the provided dataset 
# but may need to be changed, if new data is added. 
id.vector <- c('EPI', 'PRE', 'DP', 'DN')

# Run mkay() to perform k-means and classify cells according to id.vector
g6.chimeras <- mkay(dataset = g6.chimeras, miniset = gg, 
                    k = 2, ids = id.vector,
                    DP = T, x.var = 'CH5.ebLogCor', y.var = 'CH3.ebLogCor')

# Uncomment below to visualize the result
# qplot(CH5.ebLogCor,  CH3.ebLogCor,
#       data = subset(g6.chimeras, TE_ICM == 'ICM'),
#       color = Identity.km) + scale_color_manual(values = idcols) +
#   facet_grid(Treatment ~ H_genotype) + looks + theme(aspect.ratio = 1)

################################################################################
# Assign lineage identity using Hierarchical clustering
################################################################################

# Separate TE and ICM cells. TE cells are manually assigned to cluster 0,
# with Identity.hc == 'TE'
te <- subset(g6.chimeras, TE_ICM == 'TE')
icm <- subset(g6.chimeras, TE_ICM == 'ICM')

te$id.cluster <- 0
te$Identity.hc <- 'TE'

# Perform hierarchical clustering with ICM cells using average linkage (UPGMA)
# Classify separately KO cells and all other genotypes, 
# as doing all together results in the non-ko genotypes pulling too many
# boundary ko cells towards the PrE cluster (Gata6-/- don't have PrE)
icm$grouping <- ifelse(icm$H_genotype != 'ko', 'NK', 'K')
icm <- split(icm, as.factor(icm$grouping))

# Loop through the list of two elements doing the same operation on each.
# Break into 4 and 5 clusters. No DP appear in KO cells though
ks <- c(4, 5)
for(i in 1:2) { 
  my.clusters <- hclust(dist(data.frame(icm[[i]]$CH5.ebLogCor, 
                                        icm[[i]]$CH3.ebLogCor)), 
                        method = 'average')
  # Uncomment below to plot dendrogram
  # plot(my.clusters)
  icm[[i]]$id.cluster <- cutree(my.clusters, ks[i])
  my.table <- table(icm[[i]]$id.cluster, icm[[i]]$Identity.km)
  print(my.table)
  my.plot <- qplot(CH5.ebLogCor,  CH3.ebLogCor,
                   data = icm[[i]], color = id.cluster) +
    looks + scale_color_gradient2(low = 'black', mid = 'green',
                                  high = 'yellow',
                                  midpoint = (ks[i]+1)/2) +
    facet_wrap(H_genotype ~ cell_type) + theme(aspect.ratio = 1)
  # Uncomment below to print the outcome as scatter plots
  # print(my.plot)
}

# Assign identities to each cluster
idxclust <- list(data.frame(id.cluster = 1:ks[1], 
                            Identity.hc = c('EPI', 'EPI.lo', 'PRE', 'EPI')), 
                 data.frame(id.cluster = 1:ks[2], 
                            Identity.hc = c('DP', 'EPI', 'PRE', 'PRE', 
                                            'EPI.lo')))

# Merge identities with icm dataset
for(e in 1:length(icm)) { 
  icm[[e]] <- merge(icm[[e]], idxclust[[e]])
}
icm <- do.call(rbind, icm)
icm$grouping <- NULL

# Combine TE and ICM cells again
g6.chimeras <- rbind(te, icm)
rm(te, icm)

# Uncomment below to visualize the result
# qplot(CH5.ebLogCor,  CH3.ebLogCor,
#       data = subset(g6.chimeras, TE_ICM == 'ICM'), color = Identity.hc) +
#   looks + scale_color_manual(values = idcols) +
#   facet_wrap(H_genotype ~ cell_type) + theme(aspect.ratio = 1)

# A few Gata6-/- cells (18) are classified as PrE cells. 
# Some of these cells were checked visually and they're actually wild type
# GFP+ cells with very faint GFP.
# As can be seen on the plot below, these cells show higher GATA6 than the rest
# of GFP- mutant EPI cells, but their GFP levels are comparable
# to the GFP- mutant cells and were mistakenly classified as GFP-.

# Uncomment below to generate plot
# qplot(CH5.ebLogCor, CH2.logCor, 
#       data = subset(g6.chimeras, H_genotype == 'ko' & 
#                       cell_type == 'host' & TE_ICM == 'ICM'), 
#       geom = 'jitter', color = Identity.hc) + looks + 
#   scale_color_manual(values = idcols) + theme(aspect.ratio = 2)

# Therefore, Gata6 ko (GFP-, host) cells classified as PrE, with GFP (CH2.logCor)
# values > 5.25 will be re-classified as GFP+ (donor), those with values <5.25
# will be left as GFP- but assigned to EPI - those cells are also on the edge
# of the EPI cluster when looking at log[GATA6] vs log[NANOG]
aa <- subset(g6.chimeras, H_genotype == 'ko' & 
               cell_type == 'host' & Identity.hc == 'PRE')
bb <- subset(g6.chimeras, H_genotype != 'ko' | 
               cell_type != 'host' | Identity.hc != 'PRE')
aa$cell_type <- ifelse(aa$CH2.logCor > 5.25, 'donor', 'host')
aa$Identity.hc <- ifelse(aa$CH2.logCor < 5.25, 'EPI', 'PRE')
g6.chimeras <- rbind(aa, bb)
rm(aa, bb)

# Make Identity.hc a factor and order it
g6.chimeras$Identity.hc <- factor(g6.chimeras$Identity.hc, 
                                  levels = c('TE', 'PRE', 'DP', 
                                             'EPI', 'EPI.lo', 'DN'))

################################################################################
# Write out data to disk
write.csv(g6.chimeras, file = './data/processed/emb-xim-processed.csv', 
          row.names = F)

