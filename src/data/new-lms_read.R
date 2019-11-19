# This file reads in curated files for the new-lms (new littermates) dataset
# from the corresponding folder (data/corfiles)

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
my.dir <- './data/corfiles/new-lms_corf/'

# Read in each of the files in the new-littermates folder and combine in one
# data frame
my.files <- dir(my.dir)
new.lms <- do.call(rbind.fill, 
                   lapply(my.files, 
                          function(x) read.csv(paste(my.dir, x, sep = ''))))
# Remove spurious variable that MS Excel generated
new.lms$X.1 <- NULL

# And again for Fgf4 littermates
my.dir <- './data/corfiles/fgf-lms_corf/'

my.files <- dir(my.dir)
f4.lms <- do.call(rbind.fill, 
                  lapply(my.files, 
                         function(x) read.csv(paste(my.dir, x, sep = ''))))
# Remove spurious columns and rename (manual) Identity variable
f4.lms[which(colnames(f4.lms) %in% c('X.1', 'Cellcount'))] <- NULL
f4.lms <- rename(f4.lms, Identity.man = Identity)

# Fill in empty values in TE_ICM variable using (manual) Identity variable
f4.lms$TE_ICM[which(is.na(f4.lms$TE_ICM) == T)] <- 
  ifelse(f4.lms$Identity.man[
    which(is.na(f4.lms$TE_ICM) == T)] == 'TE', 'TE', 'ICM')

# Combine both datasets
new.lms <- rbind.fill(new.lms, f4.lms)
rm(f4.lms)

# Log embryos loaded
n.start <- unique(new.lms$Embryo_ID)

# Read in experimental reference files and incorporate into main data frame
new.lms.ref <- rbind(read.csv('./references/new-littermates_exp_ref.csv'), 
                     read.csv('./references/fgf4-lms_exp_ref.csv'))
new.lms <- merge(new.lms, new.lms.ref)

# Calculate total cell count, ICM cell count and litter median with do_counts
new.lms <- do.counts(new.lms, sep.treatment = F)
# and stage embryos based on their cell count
new.lms <- stage(new.lms)

# Write this raw, unprocessed data frame out to file
write.csv(new.lms, file = './data/raw/new-lms-raw.csv', row.names = F)
