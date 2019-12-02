# This file reads in curated files for the ESC chimeras dataset (esc-xim)
# from the corresponding folder (data/corfiles), combines them together 
# and performs basic counts of cell numbers

# Check that setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load functions that will be used in the script
source('./src/functions/do_counts.R')
source('./src/functions/stage.R')

# Define the location of files to be loaded (relative path to working directory)
my.dir <- './data/corfiles/esc-xim_corf/'

# Read in each of the files in the esc-xim folder and combine in one
# data frame
my.files <- dir(my.dir)
esc.chimeras <- do.call(rbind.fill, 
                   lapply(my.files, 
                          function(x) read.csv(paste(my.dir, x, sep = ''))))

# Remove spurious variables
esc.chimeras[which(colnames(esc.chimeras) %in% c('Cellcount', 
                                                 'X.1', 'X.2'))] <- NULL

# Record the number of embryos loaded
n.start <- unique(esc.chimeras$Embryo_ID)

# Convert manually assigned Identity to TE/ICM/ESC in all embryos
esc.chimeras$Identity <- ifelse(esc.chimeras$Identity == 'TE', 'TE', 
                                ifelse(esc.chimeras$Identity == 'ESC', 'ESC',
                                       'ICM'))

# Load the experimental reference file
esc.ref <- read.csv('./references/esc-xim_exp_ref.csv')
esc.ref$Background <- 'CD1/mixed'

# Incorporate genotype info for ESCs
esc.genotype <- read.csv('./references/esc-genotype.csv')
esc.ref <- merge(esc.ref, esc.genotype)
rm(esc.genotype)

# Reformat embryo genotype names to be consistent with other datasets
viejos <- c("WT", "F4FL_KO", "F4_KO", "F4FL_FL", "F4FL_KO?")
nuevos <- c('wt', 'fl/ko', 'ko', 'fl/fl', 'unclear')
for(i in 1:length(viejos)) {
  esc.ref$Genotype1[grepl(viejos[i], esc.ref$Genotype1)] <- 
    nuevos[i]
} 
rm(viejos, nuevos)

# Apply ESC culture and ESC line information for each experiment
# to all Control, Chimeras and Littermates, for plotting purposes
# (even though Controls and Littermates do not actually have ESCs)
esc.info <- esc.ref %>% filter(Treatment == 'Chimera') %>%
  group_by(Litter, ESC_line, ES_culture) %>%
  summarize()
esc.ref$ESC_line <- NULL
esc.ref$ES_culture <- NULL
# separate from esc.ref litters where there are only controls (would be dropped)
uno <- subset(esc.ref, Litter %in% esc.info$Litter)
dos <- subset(esc.ref, !Litter %in% esc.info$Litter)
# fill in ESC genotype and culture with 'no.esc'
dos$ESC_line <- 'no.esc'
dos$ES_culture <- 'no.esc'
# Combine all subsets again
esc.ref <- merge(uno, esc.info)
esc.ref <- rbind(esc.ref, dos)
rm(esc.info, uno, dos)

# Drop duplicate rows generated for Litter DI & CY
esc.ref <- filter(esc.ref, 
                      !(Embryo_ID %in% c('010919Xim_F4_DI3', 
                                         '010919Xim_F4_DI4') & 
                          ES_culture == '2i/LIF'), 
                      !(Embryo_ID %in% c('010919Xim_F4_DI5', 
                                         '010919Xim_F4_DI6', 
                                         '010919Xim_F4_DI7', 
                                         '010919Xim_F4_DI8') & 
                          ES_culture == 'S/LIF'), 
                      !(Embryo_ID %in% c('111418Xim_F4_CY2', 
                                         '111418Xim_F4_CY3') & 
                          ES_culture == '2i/LIF'), 
                      !(Embryo_ID %in% c('111418Xim_F4_CY5', 
                                         '111418Xim_F4_CY6') & 
                          ES_culture == 'S/LIF'))

# Combine main table and experimental reference
esc.chimeras <- merge(esc.chimeras, esc.ref)

# Create TE_ICM variable based on Identity (TE, ICM and ESCs)
esc.chimeras$TE_ICM <- ifelse(esc.chimeras$Identity == 'TE', 'TE', 'ICM')

# Run do.counts to calculate cell count per embryo, 
# number of ICM cells per embryo and average cell count per litter
esc.chimeras <- do.counts(esc.chimeras)

# Stage embryos
esc.chimeras <- stage(esc.chimeras)

# Rename Tt_stage to Stage.t0
esc.chimeras <- rename(esc.chimeras, Stage.t0 = Tt_stage)

# The experiments (litters) done between Feb 14 2018 and Feb 22 2018 did not
# grow beyond 85 cells, unlike most of the other experiments.
# Around those dates we discovered an issue with the mineral oil 
# under which we grow embryos (Sigma BioReagent light mineral oil (neat), 
# cat# M8410). We had this issue in the past, and we had determined the oil
# was killing embryos. Other labs had the same problem. 
# Therefore we exclude these embryos from later analyses
small.odd <- unique(subset(esc.chimeras, 
                           Exp_date >= 20180214 & 
                             Exp_date <= 20180222)$Litter)

################################################################################
# Write this raw, unprocessed data frame out to file
write.csv(esc.chimeras, file = './data/raw/esc-xim-raw.csv')
