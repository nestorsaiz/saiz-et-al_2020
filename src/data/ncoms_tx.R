# This scripts re-classifies ICM cells in the Nature Communications dataset 
# using Hierarchical clustering, instead of K-means clustering, as we have
# used for data in the present study

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# If data is not present, read in from interim data folder
if(exists('ncoms.lms') == F) { 
  ncoms.lms <- read.csv('./data/interim/ncoms-lms-tidy.csv')
}

# If still not present, run all-lms_read.R to generate it
if(exists('ncoms.lms') == F) {
  source('./src/data/all-lms_read.R')
}

# Order factor levels for plotting
ncoms.lms$Stage <- factor(ncoms.lms$Stage, 
                          levels = c('32_64', '64_90', '90_120', 
                                     '120_150', '>150'))

# Standardize NANOG and GATA6 levels per litter for classification purposes, 
# as done for other datasets
ncoms.lms <- split(ncoms.lms, as.factor(ncoms.lms$Experiment))
for(l in 1:length(ncoms.lms)) { 
  ncoms.lms[[l]]$CH4.ebLogCor.l <- ncoms.lms[[l]]$CH4.ebLogCor / 
    max(ncoms.lms[[l]]$CH4.ebLogCor)
  ncoms.lms[[l]]$CH5.ebLogCor.l <- ncoms.lms[[l]]$CH5.ebLogCor / 
    max(ncoms.lms[[l]]$CH5.ebLogCor)
}
ncoms.lms <- do.call(rbind, ncoms.lms)

# Impose identity on TE cells
te <- subset(ncoms.lms, TE_ICM == 'TE')
te$id.cluster <- 0
te$Identity.hc <- 'TE'

# Classify ICM cells in blastocysts using Hierarchical clustering
icm <- subset(ncoms.lms, TE_ICM == 'ICM')
my.clusters <- hclust(dist(data.frame(icm$CH4.ebLogCor.l, 
                                      icm$CH5.ebLogCor.l)), 
                      method = 'average')
# Uncomment to print dendrogram
# plot(my.clusters)

k = 15
icm$id.cluster <- cutree(my.clusters, k)
table(icm$id.cluster, icm$Identity.km)

# Uncomment to visualize the result
my.plot <- qplot(CH4.ebLogCor,  CH5.ebLogCor,
                 data = icm, color = id.cluster) +
        looks + scale_color_gradient2(low = 'black', mid = 'green',
                                      high = 'yellow', midpoint = (k+1)/2) +
        facet_grid(Genotype1 ~ Stage) + theme(aspect.ratio = 1)
print(my.plot)

# Assing identity to each cluster and combine with icm
idxclust <- data.frame(id.cluster = 1:k, 
                       Identity.hc = c('PRE', 'DP', 'EPI', 
                                       'PRE', 'EPI', 'EPI.lo', 
                                       'DP', 'DN', 'DN', 'PRE', 
                                       'EPI', 'PRE', 'DN', 
                                       'EPI', 'DN'))
icm <- merge(icm, idxclust)

# Regenerate ncoms.lms dataset and remove te and icm
ncoms.lms <- rbind(te, icm)
rm(te, icm, my.clusters, idxclust)

ncoms.lms$Identity.hc <- factor(ncoms.lms$Identity.hc, 
                               levels = c('TE', 'PRE', 'DP', 'EPI', 
                                          'EPI.lo', 'DN'))

