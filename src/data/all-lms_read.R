# This script reads in Littermates data from multiple sources
# cleans up the data and combines them into lists for plotting.
# It generates some output tables along the way, but does not write out
# the clean datasets. These are all combined in a list, that is 
# meant to be used for plotting.
# It generates some summary statistics tables that are written out to 
# ./results

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Source necessary functions
source('./src/functions/do_counts.R')
source('./src/functions/stage.R')

################################################################################
# Read in littermates from Saiz et al (2016) Nature Communications 
################################################################################

# If data is not present, read it in from ./data/raw
if(exists('ncoms.lms') == F) { 
  ncoms.lms <- read.csv('./data/raw/ncoms-lms-raw.csv')
}

# If raw table has not been generated, source gatherer.R
if(exists('ncoms.lms') == F) { 
  source('./src/data/gatherer.R')
  }

# Add missing columns
ncoms.lms$Gene1 <- 'any'
ncoms.lms$Gene2 <- 'any'
ncoms.lms$Genotype1 <- 'wt'
ncoms.lms$Genotype2 <- 'wt'
ncoms.lms$Background <- 'CD1/mixed'
# Create Litter column from Experiment (they're the same in this case)
ncoms.lms$Litter <- ncoms.lms$Experiment
# And remove unnecessary columns
ncoms.lms[which(colnames(ncoms.lms) %in% 
                  c('Markers', 'Identity.lin', 'Regime', 
                    'Tt_length', 'Xpoint', 'CH1.logCor', 'CH2.logCor', 
                    'CH3.logCor', 'CH4.logCor', 'CH5.logCor', 
                    'CH1.Sum', 'CH2.Sum', 'CH3.Sum', 'CH4.Sum', 
                    'CH5.Sum'))] <- NULL

# Drop Cellcount variable and calculate total and ICM counts 
# as well as litter average for each embryo
ncoms.lms$Cellcount <- NULL
ncoms.lms <- do.counts(ncoms.lms, sep.treatment = F)

# Order factors in Identity.km
ncoms.lms$Identity.km <- factor(ncoms.lms$Identity.km, 
                                levels = c('TE', 'PRE', 'DP', 'EPI', 'DN'))

# Rename Identity variable (manually assigned) as Identity.man
ncoms.lms <- rename(ncoms.lms, Identity.man = Identity)

# Calculate cell counts for each ICM lineage
ncoms.counts <- ncoms.lms %>% 
  group_by(Experiment, Litter, 
           Exp_date, Img_date, 
           Embryo_ID, Background, 
           Treatment, TE_ICM, 
           Cellcount, Stage, 
           litter.median, icm.count, 
           Gene1, Genotype1, 
           Gene2, Genotype2, 
           Identity.km) %>% 
  summarize(count = n())
# Calculate the % of ICM that each ICM lineage represents
# (the result will be > 100% for TE cells, but it's irrelevant)
ncoms.counts$pc.icm <- ncoms.counts$count / ncoms.counts$icm.count * 100

# Force ncoms.counts into a data.frame
ncoms.counts <- data.frame(ncoms.counts)

# Write out tidy datasets to file
write.csv(ncoms.lms, file = './data/interim/ncoms-lms-tidy.csv', row.names = F)
write.csv(ncoms.counts, file = './data/interim/ncoms-counts-tidy.csv', 
          row.names = F)

################################################################################
# Read in littermates from Morgani et al (2018) Developmental Biology
################################################################################

# If data is not present, read it in from ./data/raw
if(exists('spry.lms') == F) { 
  spry.lms <- read.csv('./data/raw/spry4-lms-raw.csv')
}

# If raw table has not been generated, source gatherer.R
if(exists('spry.lms') == F) { 
  source('./src/data/gatherer.R')
}

# Rename Spry4:H2B-Venus/H2B-Venus (homozygous, 'homo') to 'ko'
# for consistency with other datasets (it is a knockin knockout)
spry.lms$Genotype1[which(spry.lms$Genotype1 == 'homo')] <- 'ko'

# Remove unnecessary variables
spry.lms[which(colnames(spry.lms) %in% 
                 c('farred.marker', 'red.marker', 
                   'farred.ab2', 'red.ab2', 'Id_man', 
                   'Tt_stage', 'Tt_length'))] <- NULL

# Rename variables to match rest of datasets
spry.lms <- rename(spry.lms, CH1.Avg = hoechst.Avg, 
                   CH2.Avg = green.Avg, CH3.Avg = farred.Avg,
                   CH4.Avg = bf.Avg, CH5.Avg = red.Avg, 
                   CH1.ebLogCor = hoechst.ebLogCor, 
                   CH2.ebLogCor = green.ebLogCor, 
                   CH3.ebLogCor = farred.ebLogCor,
                   CH4.ebLogCor = bf.ebLogCor, 
                   CH5.ebLogCor = red.ebLogCor, 
                   Marker = green.marker, 
                   ab.2 = green.ab2, 
                   litter.median = litter.mean)

# Add missing variables
spry.lms$Channel <- 'CH2'
spry.lms$Background <- 'CD1/mixed'

# Order factor levels
spry.lms$Identity.km <- factor(spry.lms$Identity.km, 
                               levels = c('TE', 'PRE', 'DP', 'EPI', 'DN', 
                                          'morula'))
spry.lms$Stage <- factor(spry.lms$Stage, 
                         levels = c('<8', '8_16', '16_32', 
                                    '32_64', '64_90', 
                                    '90_120', '120_150', '>150'))

# Write out tidy dataset to file
write.csv(spry.lms, file = './data/interim/spry4-lms-tidy.csv', row.names = F)

# Re-classify cells using Hierarchical Clustering, as we've done for
# other data in this study, using all-lms_tx.R

# Delete existing data table
rm(spry.lms)
# Read in output of all-lms_tx.R
spry.lms <- read.csv('./data/processed/spry4-lms-processed.csv')

# If the processed table does not exist, run all-lms_tx.R
if(exists('spry.lms') == F) { 
  source('./src/data/all-lms_tx.R')
}

# Reorder factor in case data was loaded from file
spry.lms$Identity.hc <- factor(spry.lms$Identity.hc, 
                               levels = c('TE', 'PRE', 'DP', 'EPI', 'EPI.lo', 
                                          'DN', 'morula'))

# Calculate cell counts for each ICM lineage in blastocysts only
spry.counts <- spry.lms %>% 
  filter(Cellcount > 31) %>% 
  group_by(Experiment, Litter, Embryo_ID, 
           TE_ICM, Exp_date, Img_date, 
           Treatment, Gene1, Gene2, 
           Genotype1, Genotype2, 
           Cellcount, Stage, Background, 
           litter.median, Identity.hc, icm.count) %>% 
  summarize(count = n())

# Account for zeroes
# Split ICM from TE cells
m.spy <- subset(spry.counts, TE_ICM == 'ICM')
t.spy <- subset(spry.counts, TE_ICM != 'ICM')
# Cast ICM cells to wide format
m.spy <- dcast(m.spy, Exp_date + Img_date + Experiment + Litter + Embryo_ID + 
                 Cellcount + TE_ICM + Stage + Treatment + Background + 
                 icm.count + litter.median + Genotype1 + Gene1 + 
                 Genotype2 + Gene2 ~ Identity.hc, value.var = 'count')
# Replace NAs in each ICM lineage denomination with zeros
m.spy[is.na(m.spy)] <- 0
# Melt ICM back to long format
m.spy <- melt(m.spy, id.vars = c('Exp_date', 'Img_date', 'Experiment', 'Litter', 
                                 'Embryo_ID', 'Cellcount', 'TE_ICM', 'Treatment', 
                                 'Stage', 'Background', 'Gene1', 'Genotype1',  
                                 'Gene2', 'Genotype2', 'icm.count', 
                                 'litter.median'), 
              variable.name = 'Identity.hc', value.name = 'count')
## Combine TE and ICM cells again
spry.counts <- rbind.fill(m.spy, t.spy)
rm(m.spy, t.spy)

## Calculate the % of the ICM that each ICM lineage represents
# (the result will be > 100% for TE cells, but it's irrelevant)
spry.counts$pc.icm <- spry.counts$count / spry.counts$icm.count * 100

# Write counts file out to ./data/processed
write.csv(spry.counts, file = './data/processed/spry4-lms-counts.csv', 
          row.names = F)

################################################################################
# Read in littermates generated in this study
################################################################################

# Read in processed new litttermates data
new.lms <- read.csv('./data/processed/new-lms-processed.csv')

# If data doesn't exist, generate it
if(exists('new.lms') == F) {
  source('./src/new-lms_runall.R')
}

# Remove some unnecessary variables
new.lms[which(colnames(new.lms) %in% 
                    c('MINS_correct', 'CH1.Sum', 'CH2.Sum', 'CH3.Sum', 
                      'CH4.Sum', 'CH5.Sum', 'id.cluster'))] <- NULL

# Order factors
new.lms$Identity.km <- factor(new.lms$Identity.km,
                                  levels = c('TE', 'PRE', 'DP', 'EPI', 'DN'))
new.lms$Identity.hc <- factor(new.lms$Identity.hc, 
                                  levels = c('TE', 'PRE', 'DP', 'EPI',
                                             'EPI.lo', 'DN', 'n/a'))
new.lms$Identity <- factor(new.lms$Identity, 
                           levels = c('TE', 'PRE', 'DP', 'EPI',
                                      'EPI.lo', 'DN'))
new.lms$Stage <- factor(new.lms$Stage,
                            levels = c('<8', '8_16', '16_32', '32_64', 
                                       '64_90', '90_120', '120_150', '>150'))

# Read in counts table
if(exists('new.lms.counts') == F) {
  new.lms.counts <- read.csv('./data/processed/new-lms-counts.csv')
}

# If not present, run script to generate it
if(exists('new.lms.counts') == F) {
  source('./src/data/new-lms_counter.R')
}

new.lms.counts$Identity <- factor(new.lms.counts$Identity, 
                                  levels = c('TE', 'PRE', 'DP', 'EPI', 
                                             'EPI.lo', 'DN'))

# Read in processed ablations littermates data, some of which
# are contained within new.lms
ablat.lms <- read.csv('./data/processed/ablat-lms-processed.csv')

# If table is not present, generate it from scratch
if(exists('ablat.lms') == F) {
  source('./src/ablat_runall.R')
}

# Rename group.median as litter.median, 
# since the only group here are littermates
ablat.lms <- rename(ablat.lms, litter.median = group.median)
# Remove unnecessary variables
ablat.lms[which(colnames(ablat.lms) %in% 
                      c('MINS_correct', 'Cellcount_t0', 'id.cluster', 
                        'Recovery', 'Cell_diff', 'target', 
                        'target_t0', 'target_killed', 
                        'CH1.Sum', 'CH2.Sum', 'CH3.Sum', 
                        'CH4.Sum', 'CH5.Sum', 
                        'Gene3', 'Genotype3'))] <- NULL

# Order factors and duplicate Identity.hc as Identity
ablat.lms$Identity.hc <- factor(ablat.lms$Identity.hc,
                                    levels = c('TE', 'PRE', 'DP', 'EPI', 
                                               'EPI.lo', 'DN'))
ablat.lms$Identity.km <- factor(ablat.lms$Identity.km,
                                    levels = c('TE', 'PRE', 'DP', 'EPI', 'DN'))
ablat.lms$Identity <- ablat.lms$Identity.hc

# Read in ICM lineage counts
if(exists('ablat.lms.counts') == F) {
  ablat.lms.counts <- read.csv('./data/processed/ablat-counts.csv')
}

# If it doesn't exist, run script to generate it
if(exists('ablat.lms.counts') == F) {
  source('./src/data/ablat_counter.R')
}

ablat.lms.counts <- subset(ablat.lms.counts, Treatment == 'Littermate')

# Drop unnecessary columns
ablat.lms.counts[which(colnames(ablat.lms.counts) %in% 
                         c('Recovery', 'Cell_diff', 'target', 
                           'target_t0', 'Gene3', 'Genotype3', 
                           'Cellcount_t0', 'litter.median.t0'))] <- NULL
# Rename variables to match other datasets
# Call Identity.hc Identity, group.median litter.median, etc, 
# to match other .counts datasets
ablat.lms.counts <- rename(ablat.lms.counts, 
                           litter.median = group.median, 
                           Stage = Stage.t0, 
                           Identity = Identity.hc)

# Order Identity levels
ablat.lms.counts$Identity <- factor(ablat.lms.counts$Identity,
                                       levels = c('TE', 'PRE', 'DP', 'EPI', 
                                                  'EPI.lo', 'DN'))

# Re-organize these data to separate:
# - Ablation littermates
# - New CD1 littermates, excluding ablation littermates and Fgf4 embryos
# - FVB littermates
# - Fgf4 littermates
fvb.raw <- subset(new.lms, Background == 'FVB')
fvb.counts <- subset(new.lms.counts, Background == 'FVB')

f4.lms.raw <- subset(new.lms, Gene1 == 'Fgf4')
f4.lms.raw$Genotype1[which(f4.lms.raw$Genotype1 == 'wt/fl')] <- 'wt'
f4.lms.raw$Genotype1[which(f4.lms.raw$Genotype1 == 'fl/ko')] <- 'het'

f4.lms.counts <- subset(new.lms.counts, Gene1 == 'Fgf4')
f4.lms.counts$Genotype1[which(f4.lms.counts$Genotype1 == 'wt/fl')] <- 'wt'
f4.lms.counts$Genotype1[which(f4.lms.counts$Genotype1 == 'fl/ko')] <- 'het'

new.lms <- subset(new.lms, Background != 'FVB' & 
                    !Embryo_ID %in% unique(ablat.lms$Embryo_ID) & 
                    Gene1 != 'Fgf4')
new.lms.counts <- subset(new.lms.counts, 
                         Background != 'FVB' & 
                           !Embryo_ID %in% unique(ablat.lms.counts$Embryo_ID) & 
                           Gene1 != 'Fgf4')

################################################################################
## Read in FGF4 lineage counts data from *Kang et al.*, (Development) 2014
## and those acquired in this study
################################################################################

# Load table with embryo data from a series of studies:
# * Kang et al., (2013) *Development*
# * Artus et al., (2010) *Development*
# * Artus et al., (2011) *Developmental Biology*
f4.counts <- read.csv('./data/mined-data/mining-lms-counts.csv')
# Extract data for Fgf4 embryos only, from Kang et al, (2013)
f4.counts <- subset(f4.counts, Gene == 'FGF4')
# Rename variables and add missing ones
f4.counts <- rename(f4.counts, Gene1 = Gene, 
                    Genotype1 = Genotype, 
                    count = Count)
f4.counts$Gene1 <- 'Fgf4'
f4.counts$Gene2 <- 'any'
f4.counts$Genotype2 <- 'wt'
f4.counts$Litter <- 'unknown'
f4.counts$Background <- "CD1/mixed"
# Create Experiment variable from Litter, to match other datasets
f4.counts$Experiment <- f4.counts$Litter

# Calculate number of ICM cells as the difference of Cellcount - TE
# Extract number of TE cells
f4.icm <- f4.counts %>% filter(TE_ICM == 'TE') %>% 
  group_by(Embryo_ID, count, Cellcount) %>% 
  summarize()
# Substract from total and reduce dataset to icm count to merge with main table
f4.icm$icm.count <- f4.icm$Cellcount - f4.icm$count
f4.icm <- f4.icm %>% group_by(Embryo_ID, icm.count) %>% 
  summarize()
f4.counts <- merge(f4.counts, f4.icm)
# Calculate the % of ICM that each lineage represents
f4.counts$pc.icm <- f4.counts$count / f4.counts$icm.count * 100

## Combine with Fgf4 littermates acquired in this study
f4.counts <- rbind.fill(f4.counts, f4.lms.counts)

# Drop markers variable, which isn't informative
f4.counts$Markers <- NULL

################################################################################
# Read in embryos for the Gata6; Gata4 and the Gata6; Nanog allelic series
# generated in this study
################################################################################

compos <- read.csv('./data/processed/compounds-processed.csv')

# If data doesn't exist, generate it
if(exists('compos') == F) { 
  source('./src/compounds_runall.R')
}

# Remove unnecessary variables
compos[which(colnames(compos) %in% 
               c('MINS_correct', 'CH1.Sum', 'CH2.Sum', 'CH3.Sum', 
                 'CH4.Sum', 'CH5.Sum'))] <- NULL


compos$Identity.km <- factor(compos$Identity.km, 
                             levels = c('TE', 'PRE', 'DP', 'EPI', 
                                        'EPI.lo', 'DN'))
compos$Identity.hc <- factor(compos$Identity.hc, 
                             levels = c('TE', 'PRE', 'DP', 'EPI', 
                                        'EPI.lo', 'DN'))
compos$Stage <- factor(compos$Stage,
                       levels = c('<8', '8_16', '16_32', '32_64', 
                                  '64_90', '90_120', '120_150', '>150'))

# If it is not loaded, read in lineage counts
if(exists('compos.counts') == F) { 
  compos.counts <- read.csv('./data/processed/compounds-counts.csv')
}

# If processed data does not exist, run script to generate it
if(exists('compos.counts') == F) {
  source('./src/data/compounds_counter.R')
}

# Order factors
compos.counts$Identity.hc <- factor(compos.counts$Identity.hc, 
                                    levels = levels(compos$Identity.hc))

################################################################################
# Add embryos from the Gata6 allelic series from Schrode et al (2014) Dev Cell
# and for a Gata4 allelic series generated by me (unpublished)
################################################################################



################################################################################
## Combine all counts datasets 
################################################################################

# Rename Identity variables to Identity in Nat Comms and Spry4 datasets
ncoms.counts <- rename(ncoms.counts, Identity = Identity.km)
spry.counts <- rename(spry.counts, Identity = Identity.hc)

# Make a list with desired datasets
allcounts.list <- list(ncoms.counts, spry.counts, 
                       ablat.lms.counts, new.lms.counts, 
                       fvb.counts, f4.counts #,
                       #compos.counts, fgfr.counts, 
                       #gata4.counts, gata6.counts
                       )

# Re-stage embryos
for(e in 1:length(allcounts.list)) { 
  allcounts.list[[e]] <- stage(allcounts.list[[e]])
}

# Exclude possible mosaic embryos and morulas 
# and slim down the data frame by dropping unnecessary variables
droppers <- c('Markers', 'Treatment', 'Exp_date', 
              'Img_date', 'Litter',  'pc.icm')

for(e in 1:length(allcounts.list)) {
  allcounts.list[[e]] <- filter(allcounts.list[[e]], 
                                Genotype1 != 'mosaic', 
                                Genotype2 != 'mosaic', 
                                Cellcount > 29)
  allcounts.list[[e]][which(colnames(allcounts.list[[e]]) %in% 
                              droppers)] <- NULL
}
rm(droppers)

# Add up EPI and EPI.lo cells in each dataset
# and re-calculate the % of ICM each compartment represents (including all.EPI)
for(e in 1:length(allcounts.list)) {
  # Split TE and ICM cells
  te <- subset(allcounts.list[[e]], TE_ICM == 'TE')
  icm <- subset(allcounts.list[[e]], TE_ICM == 'ICM')
  
  # Cast data frame to wide format
  te <- dcast(te, Experiment + Embryo_ID + Cellcount + 
                TE_ICM + Stage + 
                icm.count + litter.median +
                Genotype1 + Genotype2 + 
                Gene1 + Gene2 + 
                Background ~ Identity, 
              value.var = 'count')
  te[which(colnames(te) == 'NA')] <- NULL
  
  icm <- dcast(icm, Experiment + Embryo_ID + Cellcount + 
                 TE_ICM + Stage + 
                 icm.count + litter.median + 
                 Genotype1 + Genotype2 + 
                 Gene1 + Gene2 + 
                 Background ~ Identity, 
               value.var = 'count')
  icm[which(colnames(icm) == 'NA')] <- NULL
  
  # Turn NAs into zeroes
  icm[is.na(icm)] <- 0
  icm$all.EPI <- 0
  # Add the missing DN variable to Gata6 dataset, if present
  # and fill with zeroes
  if(e == 10) {
    icm$DN <- 0
  }
  
  # In Nat Comms dataset, where there is no EPI.lo category
  # split embryos into late and early and 
  # add EPI and DN cells as all.EPI in late embryos only
  if(e == 1) { 
    # Split into early and late blastocysts
    early <- subset(icm, Cellcount < 100)
    late <- subset(icm, Cellcount >= 100)
    # Add EPI and DN cell numbers in late embryos only
    early$all.EPI <- early$EPI
    late$all.EPI <- late$EPI + late$DN
    # Combine early and late embryos again
    icm <- rbind(early, late)
    }
  
  # For the rest of datasets, add EPI and EPI.lo cell numbers 
  else {
    icm$all.EPI <- icm$EPI + icm$EPI.lo
  }
  
  # Split data frame into embryos with both EPI and PrE
  # and embryos that are missing either lineage
  conepi <- subset(icm, all.EPI > 0 & PRE > 0)
  sinepi <- subset(icm, all.EPI == 0 | PRE == 0)
  # Calculate the PrE:EPI ratio in embryos with both EPI & PrE
  conepi$prepi.ratio <- conepi$PRE / conepi$all.EPI
  icm <- rbind.fill(conepi, sinepi)
  
  # Melt data frame back to long format
  te <- melt(te, id.vars = colnames(te)[
    which(!colnames(te) %in%
            c('TE'))],
    variable.name = 'Identity',
    value.name = 'count')

  icm <- melt(icm, id.vars = colnames(icm)[
    which(!colnames(icm) %in%
            c('PRE', 'DP', 'EPI', 'EPI.lo',
              'DN', 'all.EPI', 'morula'))],
    variable.name = 'Identity',
    value.name = 'count')
  # Rbind TE and ICM again
  allcounts.list[[e]] <- rbind.fill(te, icm)

  # Calculate the percentage of the ICM that each lineage represents
  allcounts.list[[e]]$pc.icm <-
    allcounts.list[[e]]$count / allcounts.list[[e]]$icm.count * 100

  # Order Genotypes & Identity
  allcounts.list[[e]]$Genotype1 <- factor(allcounts.list[[e]]$Genotype1,
                                          levels = c('wt', 'het', 'ko'))
  allcounts.list[[e]]$Genotype2 <- factor(allcounts.list[[e]]$Genotype2,
                                          levels = c('wt', 'het', 'mKate/+',
                                                     'ko', 'mKate/mKate'))
  allcounts.list[[e]]$Identity <- factor(allcounts.list[[e]]$Identity,
                                         levels = c('TE', 'PRE', 'DP', 'EPI',
                                                    'EPI.lo', 'all.EPI', 'DN'))
}

