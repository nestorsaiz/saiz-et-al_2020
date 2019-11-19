# This script uses data from embryos stained with two different NANOG
# or two different GATA6 antibodies to model the relationship between
# the results obtained with both antibodies for each marker 
# and use it to transform IF data of embryos stained with one 
# into the other's equivalent values.

# Our typical staining uses a rabbit anti-NANOG and a goat anti-GATA6
# to identify cell types in the ICM of blastocysts (see article for catalog #). 
# In this experiment, embryos were stained with either: 
# * rabbit anti-NANOG and rat anti-NANOG + goat anti-GATA6 as counterstain
# * goat anti-GATA6 and rabbit anti-GATA6 + rat anti-NANOG as counterstain
# Each primary antibody was always combined with the same secondary 
# (see abtest.if for combinations)

# The data is loaded and corrected for fluorescence decay along the Z-axis
# as we have done throughout the study. Therafter, a linear regression model
# is fitted to the values of each pair of antibodies.
# The intercept and slope can be used to transform the corresponding values
# of fluorescence detected with the rat anti-NANOG or the rabbit anti-GATA6.

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load functions that will be used in the script
source('./src/functions/do_counts.R')
source('./src/functions/stage.R')
source('./src/functions/eb_cor.R')
source('./src/functions/find-pairs.R')

# Define the location of files to be loaded (relative path to working directory)
my.dir <- './data/corfiles/ab-test_corf/'

# Read in each of the files in the ab-test folder and combine in one data frame
my.files <- dir(my.dir)
abtest <- do.call(rbind.fill, 
                   lapply(my.files, 
                          function(x) read.csv(paste(my.dir, x, sep = ''))))
# Remove spurious variable that MS Excel generated
abtest$X.1  <- NULL

################################################################################
# Read in and tidy the data 
################################################################################

# Read in experimental reference file
abtest.ref <- read.csv('./references/ab-test_exp_ref.csv')

# Remove spurious column and merge abtest.ref with main file
abtest$X.1 <- NULL
abtest <- merge(abtest, abtest.ref)

# Remove experiment from October 5, 2105, as those embyros were acquired 
# in our old Zeiss LSM 510 Meta, with different parameters and channel setup
abtest <- subset(abtest, Litter != '100515_CD1')

# Run do counts
abtest <- do.counts(abtest, sep.treatment = F)

# Stage embryos
abtest <- stage(abtest)

# Read in immunofluorescence reference file
abtest.if <- read.csv("./references/ab-test_if.csv")

## Incorporate IF data into main table
abtest <- merge(abtest, abtest.if)

## Order factors
abtest$TE_ICM <- factor(abtest$TE_ICM, levels = c('TE', 'ICM'))

################################################################################
# Correct Z-associated fluorescence decay
################################################################################

# Define the possible channels to go through
channels <- c('CH1.Avg', 'CH2.Avg', 'CH3.Avg', 'CH5.Avg')
# Create an empty matrix to hold the corrected values
# with the length of the number of cells in the dataset
ebLogCor <- matrix(0, nrow = length(abtest$Embryo_ID), 
                   # and as many columns as channels
                   ncol = length(channels),
                   dimnames = list(c(), channels))
# For each channel in 'channels', perform EB correction 
# and store in the corresponding column in 'ebLogCor'
for (c in channels) {
        ebLogCor[, c] <- ebcor(abtest, c)
}
# Convert ebLogCor to data frame and rename columns to CHX.ebLogCor
ebLogCor <- data.frame(ebLogCor)
ebLogCor <- rename(ebLogCor, CH1.ebLogCor = CH1.Avg, 
                   CH2.ebLogCor = CH2.Avg, 
                   CH3.ebLogCor = CH3.Avg,
                   CH5.ebLogCor = CH5.Avg)

# Incorporate ebLogCor values into main table
abtest <- cbind(abtest, ebLogCor)
rm(ebLogCor)

################################################################################
# Write out data to a file in ./data/interim
################################################################################

write.csv(abtest, file = './data/interim/abtest.csv', row.names = F)

################################################################################
# Fit a linear regression model to the relationship 
# between NANOG and GATA6 antibodies 
################################################################################

# Run the function find-pairs.R to establish the 
# pairs of antibodies and embryos to use in the models.
# Its otput is ab.pairs
find.pairs(abtest)

# Subset embryos stained with both NANOG antibodies and exclude 
# data from 050617, as it deviates notably from the other five experiments
ng.all <- subset(abtest, Experiment %in% ab.pairs[[1]][[2]][
        which(ab.pairs[[1]][[2]] != "050617NGvNG")])
# ng.new <- subset(abtest, Experiment %in% ab.pairs[[1]][[2]][1:2])
# ng.old <- subset(abtest, Experiment %in% ab.pairs[[1]][[2]][4:6])

# Subset embryos stained with both GATA6 antibodies
gatas <- subset(abtest, Experiment %in% ab.pairs[[2]][[2]])

# Fit a linear model to the NANOG stainings
ng.model <- lm(ng.all$CH2.ebLogCor ~ ng.all$CH3.ebLogCor)
# new.model <- lm(ng.new$CH2.ebLogCor ~ ng.new$CH3.ebLogCor)
# old.model <- lm(ng.old$CH2.ebLogCor ~ ng.old$CH3.ebLogCor)
# Fit a linear model to the GATA6 stainings
gata.model <- lm(gatas$CH3.ebLogCor ~ gatas$CH5.ebLogCor)
