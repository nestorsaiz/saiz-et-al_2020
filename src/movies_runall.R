# This is the driver script for the movies dataset 
# (time-lapse data). It runs the necessary scripts to:
# * read in the manually curated files (movies_read.R)
# * correct and transform data (movies_tx.R)
# * classify cells into lineages (movies_classify.R)
# * generate associated graphs

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Run read script to load in files containing tracking information
source('./src/data/movies_read.R')

# Run transform script to correct for fluorescence decay along the Z-axis
source('./src/data/movies_tx.R')

# Run classify script to assign identity to ICM cells
source('./src/data/movies_classify.R')
