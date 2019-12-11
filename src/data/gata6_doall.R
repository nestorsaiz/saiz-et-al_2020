# This script reads in the raw data from the study characterizing the
# Gata6-/- phenotype, published by Schrode *et al* (2014) in Developmental Cell

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load necessary functions 
source('./src/functions/eb_cor.R')

# Read in raw data table
gata6.data <- read.csv('./data/mined-data/ALLtable-G6-NG.csv')

################################################################################
# Clean up data and rename varibles
################################################################################

# Clean up to make variable names match current structure
gata6.data <- rename(gata6.data, Embryo_ID = embryo, 
                     Cell_ID = Cell.ID, Genotype1 = genotype, 
                     Identity.man = identity, Cellcount = count, 
                     E.Stage = stage, Stage = countstage)
gata6.data$number <- NULL
gata6.data$Gene1 <- 'Gata6'
gata6.data$Gene2 <- 'any'
gata6.data$Genotype2 <- 'wt'
gata6.data$Treatment <- 'Littermate'

# Rename mutant genotype to ko
gata6.data$Genotype1[which(gata6.data$Genotype1 == 'mut')] <- 'ko'

# Rename ICM identity as DP
gata6.data$Identity.man[which(gata6.data$Identity.man == 'ICM')] <- 'DP'

# Generate TE vs ICM variable from Identity
gata6.data$TE_ICM <- ifelse(gata6.data$Identity.man == 'TE', 'TE', 'ICM')

# Re-stage embryos using stage()
gata6.data <- stage(gata6.data)

# Calculate number of ICM cells
icmcounts <- gata6.data %>% filter(TE_ICM == 'ICM') %>% 
  group_by(Embryo_ID) %>% 
  summarize(icm.count = n())
gata6.data <- merge(gata6.data, icmcounts)
rm(icmcounts)

################################################################################
# Correct Z-associated fluorescence decay
################################################################################

# Define the possible channels to go through
channels <- c('CDX2', 'GATA6', 'NANOG')
# Create an empty matrix to hold the corrected values
# with the length of the number of cells in the dataset
ebLogCor <- matrix(0, nrow = length(gata6.data$Embryo_ID), 
                   # and as many columns as channels
                   ncol = length(channels),
                   dimnames = list(c(), channels))
# For each channel in 'channels', perform EB correction 
# and store in the corresponding column in 'ebLogCor'
for (c in channels) {
  ebLogCor[, c] <- ebcor(gata6.data, c)
}
# Convert ebLogCor to data frame and rename columns to MARKER.ebLogCor
ebLogCor <- data.frame(ebLogCor)
ebLogCor <- rename(ebLogCor, CDX2.ebLogCor = CDX2, 
                   GATA6.ebLogCor = GATA6, 
                   NANOG.ebLogCor = NANOG)

# Incorporate ebLogCor values into main table
gata6.data <- cbind(gata6.data, ebLogCor)
rm(ebLogCor)

################################################################################
# Count cell types
################################################################################

# Calculate the number of cells per lineage for each embryo
gata6.counts <- gata6.data %>%  
  group_by(Embryo_ID, TE_ICM, 
           Cellcount, Stage, 
           Gene1, Genotype1, 
           Gene2, Genotype2, 
           Treatment, Background, 
           Identity.man, icm.count) %>%
  summarize(count = n())
gata6.counts <- as.data.frame(gata6.counts)

# Account for zeroes in mutants (or else)
## Split ICM from TE cells
m.g6 <- subset(gata6.counts, TE_ICM == 'ICM')
t.g6 <- subset(gata6.counts, TE_ICM != 'ICM')
## Cast ICM cells to wide format
m.g6 <- dcast(m.g6, Embryo_ID + Cellcount + TE_ICM + Stage + Treatment + 
                Background + icm.count + Genotype1 + Gene1 + 
                Genotype2 + Gene2 ~ Identity.man, value.var = 'count')
## Replace NAs in each ICM lineage denomination with zeros
m.g6[is.na(m.g6)] <- 0
## Melt ICM back to long format
m.g6 <- melt(m.g6, id.vars = c('Embryo_ID', 'Cellcount', 'TE_ICM', 'Treatment', 
                               'Stage', 'Gene1', 'Genotype1', 'Background', 
                               'Gene2', 'Genotype2', 'icm.count'), 
             variable.name = 'Identity.man', value.name = 'count')
## Combine TE and ICM cells again
gata6.counts <- rbind(m.g6, t.g6)
rm(m.g6, t.g6)

# Calculate the percentage of the ICM each lineage represents
gata6.counts$pc.icm <- gata6.counts$count / gata6.counts$icm.count * 100

################################################################################
# Write data out to file
################################################################################

write.csv(gata6.data, file = './data/interim/gata6-tidy.csv', 
          row.names = F)
write.csv(gata6.counts, file = './data/interim/gata6-counts-tidy.csv', 
          row.names = F)


##