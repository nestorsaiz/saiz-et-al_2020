# This script takes in the transformed output of esc-xim_tx.R
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
data.exsts <- exists('esc.chimeras')
if(data.exsts == F) { 
  esc.chimeras <- read.csv('./data/interim/esc-xim-tx.csv')
  esc.ref <- read.csv('./references/esc-xim_exp_ref.csv')
}
rm(data.exsts)

# If there is no raw data, run script to generate it
data.exsts <- exists('esc.chimeras')
if(data.exsts == F) { 
  source('./src/data/esc-xim_tx.R')
  esc.ref <- read.csv('./references/esc-xim_exp_ref.csv')
}
rm(data.exsts)

# Load functions that will be used in the script
source('./src/functions/icm_spread.R')
source('./src/functions/mkay.R')

# Load some extra packages
library('Ckmeans.1d.dp')

# Read in unicos table previously generated
unicos <- read.csv('./references/esc-xim_unicos.csv')

# Data is represented x5, for the 5 channels. Slim down to one
# and drop embryos of unknown genotype
esc.chimeras <- subset(esc.chimeras, Channel == 'CH2' & 
                         Genotype2 != 'unknown')

# Extract unique staining patterns and corresponding experiments, 
# as done esc-xim_tx.R
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

# This dataset is very heterogeneous in terms of the combinations of markers 
# used for staining populations (see stains). To accurately classify ICM 
# cell types, we will split the data into four groups, based on markers used:

# 1. Embryos stained for GATA6.gt and NANOG (either .rb or .rat)
#    using our more standard scheme of GATA6.gt in CH5 (red), 
#    which have NANOG.rb in CH3 (far red) or NANOG.rat in CH2 (green)
#    These patterns correspond to stains[c(1, 2, 4, 5, 6, 10)]

# 2. Embryos in which GATA6.gt is on CH2 (green) and NANOG.rb in CH3 (far red),
#    a more unusual combination used only for a handful of experiments
#    with H2B-tdTomato ESCs - in stains[c(7, 8)] or c(exps[[7]], exps[[8]])

# 3. Embryos stained for GATA6.rb and NANOG.rat, which includes many
#    chimeras/controls with H2B-tdTomato cells (stains[c(3, 9)])

# 4. Embryos stained for GATA6.rb alone (have both H2B-tdTomato ESCs and 
#    carry the Pdgfra:H2B-GFP allele, (heterozygous for Pdgfra))
#    They correspond to stains[c(11, 12)] and Pdgfra hets embryos in stains[3]
#    NOTE: There is a mistake in the IF reference sheet for Litter CV
#          (contained in exps[[3]]), it states all embryos have NANOG.rat
#          in CH2, but that only applies to wt embryos, het embryos for Pdgfra
#          express GFP in CH2. This is fixed with the subsetting below.

list.chimeras <- list(subset(esc.chimeras, Experiment %in% 
                               c(exps[[1]], exps[[2]], exps[[4]], 
                                 exps[[5]], exps[[6]], exps[[10]])), 
                      subset(esc.chimeras, Experiment %in% 
                               c(exps[[7]], exps[[8]])), 
                      subset(esc.chimeras, Experiment %in% 
                               c(exps[[3]], exps[[9]]) & 
                               Genotype2 == 'wt'), 
                      subset(esc.chimeras, Experiment %in% 
                               c(exps[[11]], exps[[12]]) &
                               (Genotype2 != 'wt' | Gene2 == 'Pdgfra') | 
                               Experiment %in% exps[[3]] & 
                               Genotype2 != 'wt'))

# Create a list containing the 'minisets' of data to use for K-means clustering
# for each of the first 3 elements of list.chimeras
minilist <- list()
minilist[[1]] <- list.chimeras[[1]] %>% filter(TE_ICM != 'TE') %>% 
  select(Cellcount, CH5.ebLogCor.x, CH3.ebLogCor.x)
minilist[[2]] <- list.chimeras[[2]] %>% filter(TE_ICM != 'TE') %>% 
  select(Cellcount, CH2.ebLogCor, CH3.ebLogCor)
minilist[[3]] <- list.chimeras[[3]] %>% filter(TE_ICM != 'TE', 
                                               Identity != 'ESC') %>%
  select(Cellcount, CH5.ebLogCor.x, CH3.ebLogCor.x)

# Define the number of clusters to look for in each of the minisets above
ks <- c(2, 3, 2)

# Define the vectors of identities that will be used to assign identities
id.vectors <- list(c('PRE', 'EPI', 'DP'), 
                   c('DP', 'EPI', 'PRE'), 
                   c('PRE', 'EPI', 'DP'))

# Loop through list.chimeras[[1]] to [[3]] and perform K-means
# on the corresponding miniset, finding the number of clusters defined in ks
for(i in 1:3) { 
  icm.spread(minilist[[i]])
  if(i == 2) {
    list.chimeras[[i]] <- mkay(dataset = list.chimeras[[i]],
                               miniset = minilist[[i]], k = ks[i],
                               ids = id.vectors[[i]], x.var = 'CH2.ebLogCor',
                               y.var = 'CH3.ebLogCor')
  }
  else {
    list.chimeras[[i]] <- mkay(dataset = list.chimeras[[i]],
                               miniset = minilist[[i]], k = ks[i], DP = T,
                               ids = id.vectors[[i]], x.var = 'CH5.ebLogCor.x',
                               y.var = 'CH3.ebLogCor.x')
  }
  # Uncomment below to visualize the result for each subset
  # if(i == 2) {
  #   my.plot <- qplot(CH2.ebLogCor,  CH3.ebLogCor,
  #         data = subset(list.chimeras[[i]], TE_ICM == 'ICM'),
  #         color = Identity.km) + scale_color_manual(values = idcols) +
  #     facet_grid(Stage.t0 ~ Treatment) + looks + theme(aspect.ratio = 1)
  # }
  # else {
  #   my.plot <- qplot(CH5.ebLogCor.x,  CH3.ebLogCor.x,
  #         data = subset(list.chimeras[[i]], TE_ICM == 'ICM'),
  #         color = Identity.km) + scale_color_manual(values = idcols) +
  #     facet_grid(Stage.t0 ~ Treatment) + looks + theme(aspect.ratio = 1)
  # }
  # print(my.plot)
}

# For the subset #4 of list.chimeras (see above), we use the Ckmeans.1d package
# to calculate the centers of the quasi-bimodal distribution of GATA6.rb levels 
# to classify PrE and epiblast fates based on that alone

# Create conditions to test for ICM (is.icm), ESC (is.esc) and TE (is.te)
is.icm <- list.chimeras[[4]]$Identity == 'ICM'
is.esc <- list.chimeras[[4]]$Identity == 'ESC'
is.te <- list.chimeras[[4]]$Identity == 'TE'

# Extract values for transformed GATA6.rb from embryos in list.chimeras[[4]] 
# and embryos from those same litters, which are contained in list.chimeras[[3]], 
# to get a more robust result
zz <- c(subset(list.chimeras[[4]], TE_ICM != 'TE')$CH5.ebLogCor.x, 
        subset(list.chimeras[[3]], TE_ICM != 'TE' & 
                 Litter %in% unique(list.chimeras[[4]]$Litter))$CH5.ebLogCor.x)

## Use package Ckmeans.1d.dp to calculate the two centers of the distribution
zzc <- Ckmeans.1d.dp(zz, 2)

# Make a matrix to hold sum of squares (rows = icm rows, 2 columns)
ssq <- matrix(0, length(list.chimeras[[4]]$Identity[is.icm]), 2)
# and populate with min sum of squares for each cell
# (min distance to each center for each cell's [GATA6] value)
for(i in 1:2) {
  ssq[,i] <- (list.chimeras[[4]]$CH5.ebLogCor.x[is.icm] - zzc$centers[i])^2
}
# calculate what center each cell is closest to 
# (which sum of squares is smallest)
min.ssq <- apply(ssq, 1, which.min)

# Create new variable to hold the lineage identity
# assigned using k-means clustering (Identity.km)
list.chimeras[[4]]$Identity.km <- rep(NA, nrow(list.chimeras[[4]]))
# TE cells remain unchaged
list.chimeras[[4]]$Identity.km[is.te] <- 'TE'
# ESC remain unchanged
list.chimeras[[4]]$Identity.km[is.esc] <- 'ESC'
# Assign identity to ICM cells based on min.ssq values
list.chimeras[[4]]$Identity.km[is.icm] <- c('EPI', 'PRE')[min.ssq]

# Uncomment below to visualize outcome
# qplot(CH5.ebLogCor.x,  CH1.ebLogCor,
#               data = subset(list.chimeras[[4]], TE_ICM == 'ICM'),
#               color = Identity.km) + scale_color_manual(values = idcols) +
#           facet_grid(Stage.t0 ~ Treatment) + looks + theme(aspect.ratio = 1)

# Unlist data into a data frame
esc.chimeras <- do.call(rbind, list.chimeras)
rm(list.chimeras)

# Order levels in Identity.km factor
esc.chimeras$Identity.km <- factor(esc.chimeras$Identity.km, 
                                   levels = c('TE', 'PRE', 'DP', 
                                              'EPI', 'DN', 'ESC'))

################################################################################
# Hierarchical clustering
################################################################################






