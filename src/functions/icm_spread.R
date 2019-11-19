# This function makes a standard scatterplot using ggplot2
# to visualize the distribution of values for two given channels in ICM cells.

# It's a merely a shortcut to simplify scripts, as this is a good way
# to visualize the data before setting up K values for K-means clustering

# It requires ggplot2
library('ggplot2')

icm.spread <- function(dataset) {
  ol <- ggplot(dataset, aes(x = dataset[, 2], y = dataset[, 3]))
  ol <- ol + geom_jitter(aes(color = Cellcount), size = 2, alpha = 0.6) + 
    geom_density_2d()
  ol <- ol + theme_bw() + 
    scale_color_distiller(direction = 1, palette = 'Blues') + 
    theme(aspect.ratio = 1, 
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 20), 
          legend.text = element_text(size = 15), 
          legend.title = element_text(size = 20))
  print(ol)
}