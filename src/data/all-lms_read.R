# This script reads in Littermates data from multiple sources
# cleans up the data and combines them into lists for plotting

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Source necessary functions
source('./src/functions/do_counts.R')

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

# Calculate cell counts for each ICM lineage
spry.counts <- spry.lms %>% 
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
write.csv(spry.counts, file = './data/processed/spry-lms-counts.csv', 
          row.names = F)




# spry.counts$Identity <- spry.counts$Identity.hc
# spry.raw$Identity <- spry.raw$Identity.hc
