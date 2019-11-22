# This file reads in curated files for the dataset 
# corresponding to embryo-embryo chimeras
# from the corresponding folder (data/corfiles),
# combines them into one table and cleans up the data

# Check that setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load functions that will be used in the script
source('./src/functions/do_counts.R')

# Load some extra packages
library('Ckmeans.1d.dp')

# Define the location of files to be loaded (relative path to working directory)
my.dir <- './data/corfiles/emb-xim_corf/'

# Read in each of the files in the folder and combine in one data frame
my.files <- dir(my.dir)
g6.chimeras <- do.call(rbind.fill, 
                       lapply(my.files, 
                              function(x) read.csv(paste(my.dir, x, sep = ''))))

# Record the number of embryos loaded
n.start <- unique(g6.chimeras$Embryo.Id)

# Load the experimental reference file
g6.ref <- read.csv('./references/emb-xim_exp_ref.csv')

# Calculate the initial proportion of donor (H2B-GFP+) cells 
g6.ref$t0.donor <- g6.ref$D_cells/(g6.ref$H_cells + 
                                     g6.ref$D_cells) * 100

# Rename columns to match exp.ref
# and remove empty and unnecessary columns
g6.chimeras <- rename(g6.chimeras, Embryo_ID = Embryo.Id,
                      Cell_ID = Cell.ID, TE_ICM = Identity)
g6.chimeras$X.1 <- NULL
g6.chimeras$X.2 <- NULL
g6.chimeras$Cell_diff <- NULL

# Report missing embryos (mismatch between experimental reference file
# and combination of corrected MINS output files (g6.chimeras))
ee <- g6.ref$Embryo_ID[which(!g6.ref$Embryo_ID %in% n.start)]
print(paste('unaccounted for embryos:', length(ee), sep = ' '))

# Combine main table and experimental reference
g6.chimeras <- merge(g6.chimeras, g6.ref)

# Rename Tt_stage to Stage.t0 as in other datasets
g6.chimeras <- rename(g6.chimeras, Stage.t0 = Tt_stage)

# Count the number of cells per embryo, icm and litter/group average
g6.chimeras <- do.counts(g6.chimeras, sep.treatment = T)

# Remove Litters Y and Z, as they were not stained for GFP as all the others
# and will just cause trouble when classifying cells
g6.chimeras <- subset(g6.chimeras, !Litter %in% c('Y', 'Z'))

# Order the levels of factors as desired
g6.chimeras$Treatment <- factor(g6.chimeras$Treatment, levels = c('Control', 
                                                                  'Chimera'))
g6.chimeras$Embryo_size <- factor(g6.chimeras$Embryo_size, 
                                  levels = c('Single', 'Double'))
g6.chimeras$TE_ICM <- factor(g6.chimeras$TE_ICM, 
                             levels = c('TE', 'ICM'))
g6.chimeras$H_genotype <- factor(g6.chimeras$H_genotype, 
                                 levels = c('wt', 'het', 'ko', 'mosaic', 
                                            'CAG:H2B-GFP', 'unknown', 'tbd'))

################################################################################
# Write file out to disk
write.csv(g6.chimeras, file = './data/raw/emb-xim-raw.csv', row.names = F)
