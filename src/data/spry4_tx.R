# This scripts re-classifies ICM cells in the Spry4 dataset using
# Hierarchical clustering, instead of K-means clustering, as we have
# used for data in the present study

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# If data is not present, read in from interim data folder
if(exists('spry.lms') == F) { 
  spry.lms <- read.csv('./data/interim/spry4-lms-tidy.csv')
}

# If still not present, run all-lms_read.R to generate it
if(exists('spry.lms') == F) { 
  source('./src/data/all-lms_read.R')
  }

# Standardize NANOG and GATA6 levels per litter for classification purposes, 
# as done for other datasets
spry.lms <- split(spry.lms, as.factor(spry.lms$Experiment))
for(l in 1:length(spry.lms)) { 
  spry.lms[[l]]$CH3.ebLogCor.l <- spry.lms[[l]]$CH3.ebLogCor / 
    max(spry.lms[[l]]$CH3.ebLogCor)
  spry.lms[[l]]$CH5.ebLogCor.l <- spry.lms[[l]]$CH5.ebLogCor / 
    max(spry.lms[[l]]$CH5.ebLogCor)
}
spry.lms <- do.call(rbind, spry.lms)

## Impose identity on TE cells
te <- subset(spry.lms, TE_ICM == 'TE')
te$id.cluster <- 0
te$Identity.hc <- 'TE'

# and on cells classified as 'in' or 'out', which are too early (morulas)
rest <- subset(spry.lms, Identity.km == 'morula')
rest$id.cluster <- 6
rest$Identity.hc <- 'morula'

# Classify ICM cells in blastocysts using Hierarchical clustering
icm <- subset(spry.lms, TE_ICM == 'ICM' & Identity.km != 'morula')
my.clusters <- hclust(dist(data.frame(icm$CH5.ebLogCor.l, 
                                      icm$CH3.ebLogCor.l)), 
                      method = 'average')
# Uncomment to print dendrogram
# plot(my.clusters)

# Split in five clusters
k = 5
icm$id.cluster <- cutree(my.clusters, k)
table(icm$id.cluster, icm$Identity.km)

# Uncomment to visualize the result
# my.plot <- qplot(CH5.ebLogCor.l,  CH3.ebLogCor.l,
#                  data = icm, color = id.cluster) +
#         looks + scale_color_gradient2(low = 'black', mid = 'green',
#                                       high = 'yellow', midpoint = (k+1)/2) +
#         facet_grid(Genotype1 ~ Stage) + theme(aspect.ratio = 1)
# print(my.plot)

# Assing identity to each cluster and combine with icm
idxclust <- data.frame(id.cluster = 1:k, 
                       Identity.hc = c('DP', 'PRE', 'EPI', 'EPI.lo', 'DN'))
icm <- merge(icm, idxclust)

# Renerate spry.lms dataset and remove other objects
spry.lms <- rbind(icm, rbind(te, rest))
rm(icm, te, rest, my.clusters, idxclust)

spry.lms$Identity.hc <- factor(spry.lms$Identity.hc, 
                               levels = c('TE', 'PRE', 'DP', 'EPI', 'EPI.lo', 
                                          'DN', 'morula'))

# Write transformed data to file
write.csv(spry.lms, file = './data/processed/spry4-lms-processed.csv', 
          row.names = F)

