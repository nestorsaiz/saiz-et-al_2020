# This script counts cells of each type in the various ablation datasets
# generated by ablat_tx.R and ablat_classify.R.
# Note I am using the ICM cell identity determined by Hierarchical clustering

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Check if data is already loaded and read it in if not
data.exsts <- exists('ablat')
if(data.exsts == F) { 
  ablat <- read.csv('./data/processed/ablat-processed.csv')
  ablat.t0 <- read.csv('./data/processed/ablat-t0-processed.csv')
  ablat.ref <- rbind(read.csv('./references/ablat_exp_ref.csv'))
}
rm(data.exsts)

# If there is no interim data table, run the script to generate it
data.exsts <- exists('ablat')
if(data.exsts == F) { 
  source('./src/data/ablat_classify.R')
  ablat.t0 <- read.csv('./data/processed/ablat-t0-processed.csv')
  ablat.ref <- rbind(read.csv('./references/ablat_exp_ref.csv'))
}
rm(data.exsts)

################################################################################
# Calculate the number of cells per lineage for each embryo 
# for fixed end point data
################################################################################

ablat.lincounts <- ablat %>% filter(Channel == 'CH2') %>% 
  group_by(Experiment, Litter, 
           Exp_date, Img_date, 
           Embryo_ID, TE_ICM, 
           Cellcount, litter.median.t0, 
           Background, Treatment, 
           Gene1, Genotype1, 
           Gene2, Genotype2, 
           Gene3, Genotype3, 
           target, Identity.hc, icm.count, 
           group.median, Cell_diff, 
           Recovery, Stage.t0, 
           Cellcount_t0, target_t0) %>%
  summarize(count = n())

# In embryos that lack a certain cell type, that isn't registered as a 0
# by default summarizing the data in the step above

# Split ICM from TE cells
m.abl <- subset(ablat.lincounts, TE_ICM == 'ICM')
t.abl <- subset(ablat.lincounts, TE_ICM != 'ICM')

# Cast ICM cells to wide format so that each ICM denomination
# becomes a variable (PRE, DP, etc)
m.abl <- dcast(m.abl, Experiment + Exp_date + Img_date + Litter + Embryo_ID + 
                 Cellcount + litter.median.t0 + TE_ICM + Treatment + 
                 Background + target + icm.count + group.median + 
                 Genotype1 + Gene1 + Genotype2 + Gene2 + 
                 Cell_diff + Recovery + Stage.t0 + Cellcount_t0 ~ 
                 Identity.hc, value.var = 'count')
# Replace the resulting NAs in each ICM lineage denomination with zeros
m.abl[is.na(m.abl)] <- 0

# Melt ICM back to long format where all cell types are under Identity.hc
m.abl <- melt(m.abl, id.vars = c('Experiment', 'Litter', 'Embryo_ID', 
                                 'Exp_date', 'Img_date', 'Cellcount', 
                                 'TE_ICM', 'Treatment', 'target', 'Recovery', 
                                 'Background', 'Gene1', 'Genotype1',  
                                 'Gene2', 'Genotype2', 'icm.count', 
                                 'group.median', 'litter.median.t0', 
                                 'Cell_diff', 'Stage.t0', 'Cellcount_t0'), 
              variable.name = 'Identity.hc', value.name = 'count')

# Combine TE and ICM cells again
ablat.lincounts <- rbind.fill(m.abl, t.abl)
rm(m.abl, t.abl)

# Calculate the % of the ICM each lineage represents
ablat.lincounts$pc.icm <- ablat.lincounts$count / 
  ablat.lincounts$icm.count * 100

################################################################################
# Calculate number of cells per lineage for live data.
################################################################################

# Subset non-littermates for which I have cell counts of the whole embryo
# or the whole ICM AND which have identity assignments on live images
# (embryos carrying the mKate2 allele AND the Pdgfra allele)
ablat.t0.new <- subset(ablat.t0, Genotype2 %in% c('mKate/+', 'mKate/mKate') & 
                         Treatment != 'Littermate' & 
                         is.na(Identity) == F & 
                         is.na(Stage.t0) == F)

# Calculate the number of ICM cells per embryo for live data
icm.counts <- ablat.t0.new %>% filter(TE_ICM == 'ICM') %>% 
  group_by(Embryo_ID) %>% summarize(icm.count = n())
ablat.t0.new <- merge(ablat.t0.new, icm.counts)
rm(icm.counts)

# Calculate the number of cells per ICM lineage per embryo for live data
live.lincounts <- ablat.t0.new %>% filter(TE_ICM == 'ICM') %>% 
  group_by(Embryo_ID, Experiment, Litter, Stage.t0, timepoint, 
           Gene1, Genotype1, Gene2, Genotype2, 
           Gene3, Genotype3, icm.count, 
           Cellcount_t0, Identity) %>% 
  summarize(count = n())
live.lincounts <- rename(live.lincounts, Cellcount = Cellcount_t0)
live.lincounts$dataset <- 'live'

# And for fixed data
fixed.lincounts <- ablat %>% filter(Treatment == 'Littermate', 
                                    TE_ICM == 'ICM', 
                                    Channel == 'CH2') %>% 
  group_by(Embryo_ID, Experiment, 
           Litter, Stage.t0,  
           Gene1, Genotype1, 
           Gene2, Genotype2, 
           Gene3, Genotype3, 
           icm.count, Cellcount, 
           Identity.hc) %>% 
  summarize(count = n())
fixed.lincounts <- rename(fixed.lincounts, Identity = Identity.hc)
fixed.lincounts$timepoint <- 't0'
fixed.lincounts$dataset <- 'fixed'

# Account for zeroes in any lineage on both datasets
all.lincounts <- list(live.lincounts, fixed.lincounts)
for(i in 1:2) { 
  m.lincounts <- dcast(all.lincounts[[i]], Embryo_ID + Experiment + 
                         Litter + Stage.t0 + timepoint + 
                         Gene1 + Genotype1 + Gene2 + Genotype2 + 
                         Gene3 + Genotype3 + dataset + icm.count + 
                         Cellcount ~ Identity, value.var = 'count')
  m.lincounts[is.na(m.lincounts)] <- 0
  m.lincounts <- melt(m.lincounts, 
                      id.vars = c('Embryo_ID', 'Experiment', 'Litter', 
                                  'Stage.t0', 'timepoint', 'Gene1', 
                                  'Genotype1', 'Gene2', 'Genotype2', 
                                  'Gene3', 'Genotype3',
                                  'dataset', 'icm.count', 'Cellcount'), 
                      variable.name = 'Identity', value.name = 'count')
  all.lincounts[[i]] <- m.lincounts
}

live.lincounts <- all.lincounts[[1]]
fixed.lincounts <- all.lincounts[[2]]

live.lincounts$Stage.t0 <- factor(live.lincounts$Stage.t0, 
                                  levels = c('[30,50)', '[50,70)', '[70,90)', 
                                             '[90,110)', '[110,130)', '[130,150)', 
                                             '[150,170)', '[170,190)', '[190,210]'))
fixed.lincounts$Stage.t0 <- factor(fixed.lincounts$Stage.t0, 
                                   levels = levels(live.lincounts$Stage.t0))

# Extract the number of cells of each denomination at t=0 for each embryo
ablat.t0.counts <- live.lincounts %>% 
  group_by(Embryo_ID, Litter, timepoint, dataset, 
           icm.count, Cellcount, Identity, count) %>% 
  summarize() 
ablat.t0.counts <- dcast(ablat.t0.counts, Embryo_ID + Litter + 
                           timepoint + dataset + icm.count + 
                           Cellcount ~ Identity, value.var = 'count')
ablat.t0.counts <- rename(ablat.t0.counts, PRE.t0 = PRE, 
                          DP.t0 = DP, EPI.t0 = EPI)

################################################################################
# Write out data to the ./data/processed folder
write.csv(ablat.lincounts, './data/processed/ablat-counts.csv', row.names = F)
write.csv(ablat.t0.counts, './data/processed/ablat-t0-counts.csv', 
          row.names = F)

# Write out live.lincounts and fixed.lincounts  to ./data/interim folder
# for plotting later
write.csv(live.lincounts, './data/interim/live-lincounts.csv', row.names = F)
write.csv(fixed.lincounts, './data/interim/fixed-lincounts.csv', 
          row.names = F)


##