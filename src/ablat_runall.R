# This is the driver script for the cell ablation dataset
# (ablations). It runs the necessary scripts to:
# * read in the manually curated files (ablat_read.R)
# * correct and transform data (ablat_tx.R)
# * classify cells into lineages (ablat_classify.R)
# * generate associated graphs

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Run read script to load in curated files
source('./src/data/ablat_read.R')

# Run transform script to correct for fluorescence decay along the Z-axis
source('./src/data/ablat_tx.R')

# Run classify script to assign identity to ICM cells
source('./src/data/ablat_classify.R')

# Run script to count cells in each lineage and write them out to disk
source('./src/data/ablat_counter.R')
