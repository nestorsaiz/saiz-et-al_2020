# This is the driver script for the reference littermates dataset
# (new-littermates). It runs the necessary scripts to:
# * read in the manually curated files (new-lms_read.R)
# * correct and transform data (new-lms_tx.R)
# * classify cells into lineages (new-lms_classify.R)
# * generate associated graphs

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Run read script to load in curated files
source('./src/data/new-lms_read.R')

# Check if ablat.lms has been generated
if(exists('ablat.lms') == F) {
  ablat.lms <- read.csv('./data/interim/ablat-lms-tx.csv')
}

# If the interim *-tx.csv file does not exist, source the scripts to generate it
if(exists('ablat.lms') == F) { 
  source('./src/data/ablat_read.R')
  source('./src/data/ablat_tx.R')
}

# Run transform script to correct for fluorescence decay along the Z-axis
# and re-scale NANOG and GATA6 values
source('./src/data/new-lms_tx.R')

# Run classify script to assign identity to ICM cells
source('./src/data/new-lms_classify.R')