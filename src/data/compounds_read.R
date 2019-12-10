# This file reads in curated files for the compounds dataset,
# which contains embryos for two allelic series:
# * Gata6; Gata4
# * Gata6; Nanog

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
my.dir <- './data/corfiles/compounds_corf/'

# Read in each of the files in the compounds folder and combine in one
# data frame
my.files <- dir(my.dir)
compos <- do.call(rbind.fill, 
                   lapply(my.files, 
                          function(x) read.csv(paste(my.dir, x, sep = ''))))
# Remove spurious variable that MS Excel generated
compos$X.1 <- NULL

# Log the number of embryos loaded in
n.start <- unique(compos$Embryo_ID)

# Replace NAs in Cell_cycle with 'unknown' - cell cycle was not scored by HG
compos$Cell_cycle[!compos$Cell_cycle %in% c('I', 'M')] <- 'unknown'

# Read experimental reference file and combine with main table
compos.ref <- read.csv('./references/compounds_exp_ref.csv')
compos <- merge(compos, compos.ref)

# Calculate total cell count, ICM cell count and litter median with do_counts
compos <- do.counts(compos, sep.treatment = F)
# and stage embryos based on their cell count
compos <- stage(compos)

# Write this raw, unprocessed data frame out to file
write.csv(compos, file = './data/raw/compounds-raw.csv', row.names = F)


##