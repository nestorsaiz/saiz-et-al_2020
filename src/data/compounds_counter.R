# This script counts cells of each type in the compounds dataset
# generated by compounds_tx.R and compounds_classify.R.
# Note I am using the ICM cell identity determined by Hierarchical clustering

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Check if data is already loaded 
data.exsts <- exists('compos')
if(data.exsts == F) { 
  compos <- read.csv('./data/processed/compounds-processed.csv')
  compos.ref <- rbind(read.csv('./references/compounds_exp_ref.csv'))
}
rm(data.exsts)

# If raw data doesn't exist, run script to generate it
data.exsts <- exists('compos')
if(data.exsts == F) { 
  source('./src/data/compounds_classify.R')
  compos.ref <- rbind(read.csv('./references/compounds_exp_ref.csv'))
}
rm(data.exsts)

################################################################################
# Count the number of cells per lineage in each embryo
# using Identity as determined using Hierarchical Clustering (Identity.hc)
################################################################################

compos.counts <- compos %>% 
  group_by(Experiment, Litter, 
           Embryo_ID, TE_ICM, 
           Cellcount, Stage, 
           Exp_date, Img_date, 
           Background, 
           Gene1, Genotype1, 
           Gene2, Genotype2, 
           Treatment, Identity.hc, 
           icm.count, litter.median) %>%
  summarize(count = n())

# Account for zeroes in mutants (or elsewhere)
## Split ICM from TE cells
m.comp <- subset(compos.counts, TE_ICM == 'ICM')
t.comp <- subset(compos.counts, TE_ICM != 'ICM')
## Cast ICM cells to wide format
m.comp <- dcast(m.comp, Experiment + Exp_date + Img_date + Litter + Embryo_ID + 
                  Cellcount + TE_ICM + Stage + Treatment + Background + 
                  icm.count + litter.median + Genotype1 + Gene1 + 
                  Genotype2 + Gene2 ~ Identity.hc, value.var = 'count')
## Replace NAs in each ICM lineage denomination with zeros
m.comp[is.na(m.comp)] <- 0
## Melt ICM back to long format
m.comp <- melt(m.comp, id.vars = c('Experiment', 'Litter', 'Embryo_ID', 
                                   'Exp_date', 'Img_date', 
                                   'Cellcount', 'TE_ICM', 'Treatment', 
                                   'Stage', 'Background', 'Gene1', 'Genotype1',  
                                   'Gene2', 'Genotype2', 
                                   'icm.count', 'litter.median'), 
               variable.name = 'Identity.hc', value.name = 'count')
## Combine TE and ICM cells again
compos.counts <- rbind.fill(m.comp, t.comp)
rm(m.comp, t.comp)

# Calculate the percentage of the ICM each lineage represents
compos.counts$pc.icm <- compos.counts$count / compos.counts$icm.count * 100

################################################################################
# Write out counts to file
write.csv(compos.counts, file = './data/processed/compounds-counts.csv', 
          row.names = F)


##