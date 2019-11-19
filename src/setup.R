# This script loads a handful of required packages for all scripts
library('packrat')
library('reshape2')
library("plyr")
library("dplyr")
library("ggplot2")
library('RColorBrewer')

# the 'plotting-aes.R' file with some objects I use for plotting
source('./src/plotting-aes.R')

# and sets a seed for reproducibility
set.seed(21)
