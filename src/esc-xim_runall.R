# This is the driver script for the ESC chimeras dataset
# (esc-xim). It runs the necessary scripts to:
# * read in the manually curated files (esc-xim_read.R)
# * correct and transform data (esc-xim_tx.R)
# * classify cells into lineages (esc-xim_classify.R)
# * count cell populations (esc-xim_counter.R)
# * generate associated graphs

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Run read script to load in curated files
source('./src/data/esc-xim_read.R')

# Run transform script to correct for fluorescence decay along the Z-axis
# and re-scale NANOG and GATA6 values
source('./src/data/esc-xim_tx.R')

# Run classify script to assign identity to ICM cells
source('./src/data/esc-xim_classify.R')

# Run script to count cells in each lineage and write them out to disk
source('./src/data/esc-xim_counter.R')
