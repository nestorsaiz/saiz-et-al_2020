# This script reads in the files for time lapse movies 
# from a collection of folders, binds them together, 
# cleans the data, extracts information for each cell from a 'labels' string
# and writes out a tidy .csv with all movie data

# Check that setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load some extra packages
library('zoo')
library('tidyquant')

# Source functions that will be used in the script
source('./src/functions/eb_cor.R')
source('./src/functions/el-classifier.R')

# Define the location of files to be loaded (relative path to working directory
my.dir <- './data/moviefiles/'

# Record the existing folders (one per movie)
folders <- dir(my.dir)

# Create an empty list to store the movie data frames and initiate a counter
movies <- list()
cycles <- 0

# Modify read.table() with the specifics to read these files
my.reader <- function(dataset) { 
  read.table(dataset, header = T, sep = ',', skip = 3)
}

# Read in the relevant files for each movie and clean up the data
for(m in 1:length(folders)) {
  # Set the relative path for each folder in moviefiles
  subdir <- paste(my.dir, folders[m], sep = '')
  # Read in the file names in the folder and select those of interest
  files <- dir(subdir)
  goners <- c(grep('Overall', files), grep('Time', files))
  files <- files[-goners]
  
  # Record the movie name
  mov.name <- paste(strsplit(folders[m], split = "_")[[1]][1], 
                    strsplit(folders[m], split = "_")[[1]][2], 
                    strsplit(folders[m], split = "_")[[1]][3], 
                    sep = '_')
  
  # Make an empty list to store the individual data frames
  mov.files <- list()
  # Read the data frames corresponding to selected files 
  # and store in mov.files
  for(f in 1:length(files)) { 
    da.file <- paste(subdir, '/', files[f], sep = '')
    mov.files[[f]] <- my.reader(da.file)
    mov.files[[f]][which(colnames(mov.files[[f]]) %in%
                           c('Unit', 'X', 'Category'))] <- NULL
  }
  
  ## merge the files corresponding to Min, Max, Median and Mean 
  ## fluorescence values for each channel and rbind all channels
  if(length(mov.files) == 13) {
    ch1 <- join_all(mov.files[c(1, 4, 7, 10)])
    ch2 <- join_all(mov.files[c(2, 5, 8, 11)])
    ch3 <- join_all(mov.files[c(3, 6, 9, 12)])
    channels <- rbind(ch1, ch2, ch3)
  }
  if(length(mov.files) == 9) {
    ch1 <- join_all(mov.files[c(1, 3, 5, 7)])
    ch2 <- join_all(mov.files[c(2, 4, 6, 8)])
    channels <- rbind(ch1, ch2)
  }
  
  # Merge the fluorescence values with the position file
  my.mov <- merge(channels, mov.files[length(mov.files)])
  my.mov$Image <- NULL
  my.mov$Collection <- NULL
  my.mov$Embryo_ID <- mov.name
  my.mov <- rename(my.mov,
                   X = Position.X,
                   Y = Position.Y,
                   Z = Position.Z,
                   Mean = Intensity.Mean,
                   Median = Intensity.Median,
                   Min = Intensity.Min,
                   Max = Intensity.Max, 
                   timeframe = Time)
  
  movies[[m]] <- my.mov
  cycles <- cycles + 1
}

# Print out the cycle counter
print(paste('Number of files loaded:', cycles, sep = ' '))

# Cells which only have one timeframe because they die right away
# have no TrackID assigned. This can cause problems later when counting cells
# Correct this by giving them a TrackID that consists of the number
# 100000 + the corresponding spot ID number
for(m in 1:length(movies)) {
  bien <- subset(movies[[m]], is.na(TrackID) == F)
  mal <- subset(movies[[m]], is.na(TrackID) == T)
  mal$TrackID <- mal$ID + 100000
  movies[[m]] <- rbind(bien, mal)
}

# Extract information stored in the labels variable (NS.labs), such as cell
# type, generation, division, etc., and make new variables for each label
for(m in 1:length(movies)) { 
  movies[[m]]$identity_t0 <- ifelse(grepl('PRE', movies[[m]]$NS.labs), 'PRE', 
                                    ifelse(grepl('EPI', movies[[m]]$NS.labs), 'EPI', 
                                           ifelse(grepl('DP', movies[[m]]$NS.labs), 'DP', 
                                                  'unknown')))
  movies[[m]]$cell_treatment <- ifelse(grepl('targeted', movies[[m]]$NS.labs), 
                                       'targeted', 'intact')
  movies[[m]]$TrackID <- paste(movies[[m]]$TrackID, 
                               movies[[m]]$identity_t0,
                               movies[[m]]$cell_treatment, sep = '.')
  movies[[m]]$genx <- ifelse(grepl('Granddaughter', movies[[m]]$NS.labs), 'gd', 
                             ifelse(grepl('d1', movies[[m]]$NS.labs) | 
                                      grepl('d2', movies[[m]]$NS.labs), 'd', 'm'))
  movies[[m]]$gd.branch <- ifelse(grepl('Granddaughter1', movies[[m]]$NS.labs), 0.1, 
                                  ifelse(grepl('Granddaughter2', movies[[m]]$NS.labs), 0.2, 0))
  movies[[m]]$d.branch <- ifelse(grepl('d1', movies[[m]]$NS.labs), 1, 
                                 ifelse(grepl('d2', movies[[m]]$NS.labs), 2, 0))
  movies[[m]]$branch <- movies[[m]]$d.branch + movies[[m]]$gd.branch
  movies[[m]]$division <- ifelse(grepl('Mito', movies[[m]]$NS.labs), TRUE, FALSE)
  movies[[m]]$death <- ifelse(grepl('Apoptosis', movies[[m]]$NS.labs), TRUE, FALSE)
  movies[[m]]$in.frame <- ifelse(grepl('off', movies[[m]]$NS.labs), FALSE, TRUE)
  movies[[m]]$Cell_ID <- paste(movies[[m]]$TrackID, movies[[m]]$genx, 
                               movies[[m]]$branch, sep = '.')
  movies[[m]]$d.branch <- NULL
  movies[[m]]$gd.branch <- NULL
}

# Rename embryos, which currently have the file name, 
# to match the format in the main ablation experimental reference file
for(m in 1:length(movies)) { 
  # Break Embryo_ID content into components
  parts <- strsplit(movies[[m]]$Embryo_ID, split = '_')
  # Reformat date from the existing DDMonYY to MMDDYY
  my.date <- as.Date(parts[[1]][1], '%d%b%y')
  my.date <- format(my.date, '%m%d%y')
  # Reformat as MMDDYY + Abl + _ + Litter + Embryo number
  ee <- paste(paste(my.date, 'Abl', sep = ''), parts[[1]][3], sep = '_')
  movies[[m]]$Embryo_ID <- ee
}

## Rbind all elements of movies into a data frame
movies <- do.call(rbind, movies)

# Make TrackID and Cell_ID unique by incorporating the corresponding embryo name
movies <- rename(movies, spot_ID = ID)
movies$TrackID <- paste(movies$Embryo_ID, movies$TrackID, sep = '.')
movies$Cell_ID <- paste(movies$Embryo_ID, movies$Cell_ID, sep = '.')

# Add TE vs ICM variable
movies$TE_ICM <- 'ICM'

# Log and print out the actual number of movies loaded
# (there are more than one file per movie/embryo)
n.movies <- unique(movies$Embryo_ID)
print(paste('Number of movies loaded:', length(n.movies), sep = ' '))

################################################################################
# Add metadata
################################################################################

# Read in experimental reference
ablat.ref <- read.csv('./references/ablat_exp_ref.csv')
ablat.t0.stage <- read.csv('./references//ablat_t0-stage.csv')

# Extract data from ablat.ref only for embryos contained in ablat.movs dataset
# and extract just some of the variables in it
ablat.ref <- subset(ablat.ref, 
                        Embryo_ID %in% unique(movies$Embryo_ID))
ablat.ref <- ablat.ref %>% 
  select(Experiment, Litter, Embryo_ID, 
         Treatment, Cellcount_t0, 
         Gene1, Genotype1, Gene2, Genotype2, 
         target, Cell_diff, Recovery)
ablat.ref <- merge(ablat.ref, ablat.t0.stage)
rm(ablat.t0.stage)

# Combine into one set and remove exp.ref
movies <- merge(movies, ablat.ref)
rm(ablat.ref)

# Rename targets so controls are assigned 'none' as target cell
movies$target[which(movies$Treatment == 'Control')] <- 'none'

# Calculate number of embryos for each category and write to file
n.embryos <- movies %>% 
  group_by(Embryo_ID, target, Stage.t0, Treatment, Cell_diff) %>% 
  summarize() %>% 
  group_by(target, Stage.t0, Treatment, Cell_diff) %>% 
  summarize(N = n())
write.csv(n.embryos, file = './results/movies_N-embryos.csv', row.names = F)

################################################################################
# Write out tidy table to use in other applications
write.csv(movies, file = './data/raw/movies-all.csv', row.names = F)

