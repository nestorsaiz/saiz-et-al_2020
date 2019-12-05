## This script creates a series of named character vectors containing 
## my colors of choice to give to values of some variables

## It also creates an object, looks, with some preset theme configuration
## for ggplot2, to make plotting code simpler

## Load required packages
library('ggplot2')
library('colorspace')

## Create vector to define cell lineage identity colors
idcols <- c('EPI' = '#e5331a', 'all.EPI' = '#e5331a', 'EPI.lo' = '#ff6666',  
            'PRE' = '#1a33e5', 'DP' = '#993399', 'DN' = '#acacac', 
            'TE' = '#a1d99b', 'ICM' = '#993399', 'morula' = '#ce1256', 
            'ESC' = '#ccff33', 'out' = '#339933', 'in' = '#993399')

## Create vector to define cell type color ("donor", GFP- vs "host", GFP+)
cellcols <- c('host' = '#cccccc', 'donor' = '#33cc33')

## Create a generic vector for genotype colors
gencols <- c('wt' = '#D6F9DD', 'het' = '#379A54', 
             'homo' = '#007300', 'unknown' = 'black')

## Create a vector for colors to use for different ESC lines if desired
escols <- c("CAG:H2B-GFP" = '#33cc33', 
            "H2B-tdTomato" = '#cc3333', 
            "F103_CAG:H2B-GFP-B2" = '#33cc33', 
            "F105_GFP" = '#007300', 'no.esc' = 'black')

## Make object containing aesthetics for the plots (font size, etc)
looks <- theme_bw() + theme(panel.grid = element_blank(), 
                            strip.background = element_blank(), 
                            panel.border = element_rect(color = 'black', 
                                                        size = 1), 
                            axis.ticks = element_line(color = 'black', 
                                                      size = 0.5), 
                            axis.text = element_text(size = 10, 
                                                     color = 'black'), 
                            axis.title = element_text(size = 12, 
                                                      color = 'black'), 
                            legend.text = element_text(size = 10, 
                                                       color = 'black'), 
                            legend.title = element_text(size = 12, 
                                                        color = 'black'), 
                            strip.text = element_text(size = 14, 
                                                      color = 'black'))