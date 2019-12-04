# This script reads in data from the studies by
# * Saiz *et al* (2016) in *Nature Communications* and by
# * Morgani *et al* (2018) in *Developmental Biology*
# and writes out Littermates data as files to ./data/raw

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load extra packages
library('RCurl')

# Read .csv of transformed data for Saiz et al. (2016) from its GitHub repo
my.path <- getURL('https://raw.githubusercontent.com/nestorsaiz/saiz-et-al_2016/master/FGF_all_pooled_trans.csv')
ncoms.all <- read.csv(text = my.path)

# Extract Littermates only from the dataset
ncoms.lms <- subset(ncoms.all, Treatment == 'Littermate')

# Incorporate experiment and imaging date, from experimental reference file, 
# which are missing for some reason.

# Read file from GitHub
my.path <- getURL('https://raw.githubusercontent.com/nestorsaiz/saiz-et-al_2016/master/FGFonCD1_exp_ref.csv')
ncoms.ref <- read.csv(text = my.path)
# Extract Exp_date and Img_date for each experiment present in main table
ncoms.ref <- ncoms.ref %>% 
  filter(Experiment %in% unique(ncoms.all$Experiment), 
         Treatment == 'Littermate') %>% 
  group_by(Experiment, Exp_date, Img_date) %>% 
  summarize()

# Merge into main data frame
ncoms.lms <- merge(ncoms.lms, ncoms.ref)
rm(ncoms.ref)

# Write file out to disk
write.csv(ncoms.lms, file = './data/raw/ncoms-lms-raw.csv', row.names = F)


# Read .csv of transformed data for Morgani et al. (2018) from is GitHub repo
my.path <- getURL('https://raw.githubusercontent.com/nestorsaiz/morgani-et-al_2018/master/spry-all-processed.csv')
spry.all <- read.csv(text = my.path)

# Extract Littermates only stained for NANOG and GATA6 and drop
# embryos which are not wild type for other alleles than Spry4
spry.lms <- subset(spry.all, Treatment == 'Littermate' & 
                     red.marker == 'GATA6.gt' & 
                     farred.marker == 'NANOG.rb' & 
                     Genotype1 != 'unknown' & 
                     Genotype2 == 'wt')

# Write Spry4 littermates out to disk
write.csv(spry.lms, file = './data/raw/spry4-lms-raw.csv', row.names = F)


##