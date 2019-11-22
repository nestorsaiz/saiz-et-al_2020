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
## Assign lineage identity #####################################################
################################################################################


