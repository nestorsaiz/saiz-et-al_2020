# This script reads in the files from time lapse movies


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

