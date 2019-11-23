# This is the driver script for the embryo-embryo chimeras dataset
# (emb-xim). It runs the necessary scripts to:
# * read in the manually curated files (emb-xim_read.R)
# * correct and transform data (emb-xim_tx.R)
# * classify cells into lineages (emb-xim_classify.R)
# * count cell populations (emb-xim_counter.R)
# * generate associated graphs

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Run read script to load in curated files
source('./src/data/emb-xim_read.R')

# Run transform script to correct for fluorescence decay along the Z-axis
source('./src/data/emb-xim_tx.R')

# Run classify script to find GFP+ and GFP- cells and 
# to assign identity to ICM cells
source('./src/data/emb-xim_classify.R')

# Run script to count cells in each lineage and write them out to disk
source('./src/data/emb-xim_counter.R')