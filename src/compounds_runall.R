# This is the driver script for the compounds dataset
# (compounds). It runs the necessary scripts to:
# * read in the manually curated files (compounds_read.R)
# * correct and transform data (compounds_tx.R)
# * classify cells into lineages (compounds_classify.R)
# * count number of cells for each lineage (compounds_counter.R)
# * generate associated graphs

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Run read script to load in curated files
source('./src/data/compounds_read.R')

# Run transform script to correct for fluorescence decay along the Z-axis
# and re-scale NANOG values
source('./src/data/compounds_tx.R')

# Run classify script to assign identity to ICM cells
source('./src/data/compounds_classify.R')

# Run script to count cells in each lineage and write them out to disk
source('./src/data/compounds_counter.R')