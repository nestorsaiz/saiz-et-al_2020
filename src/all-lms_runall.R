# This is the driver script to find and combine the multiple
# Littermates datasets analyzed in the study:
# * those published in Saiz et al., (2016) Nat Comms
# * those published in Morgani et al., (2018) Dev Biol
# * those generated in this study

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Run script to read all littermates from the raw files
# tidy up data, do cell counts and combine them all in a list
source('./src/data/all-lms_read.R')

# Run script to calculate summary statistics of ICM composition
# for each subset of the data
source('./src/data/all-lms_stats.R')