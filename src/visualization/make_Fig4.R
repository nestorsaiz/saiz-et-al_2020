# This script generates the plots in Figure 4 of the paper

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load necessary data
if(exists('movies') == F) { 
  movies <- read.csv('./data/processed/movies-all-processed.csv')
  }

# Load some extra packages
library('zoo')
library('tidyquant')

################################################################################
# Generate endings table as in movies_classify.R
################################################################################

# Split data by Cell_ID (segment of a track between mitosis)
movies <- split(movies, as.factor(movies$Cell_ID))
endings <- list()
# Select the last spot of each track
for(c in 1:length(movies)) { 
  my.sum <- movies[[c]] %>% filter(timeframe == max(timeframe)) %>% 
    group_by(Cell_ID, identity_t0, timeframe, division, 
             death, in.frame, end.point) %>% summarize()
  endings[[c]] <- my.sum
}
movies <- do.call(rbind, movies)
endings <- do.call(rbind, endings)

# Label the end points based on label information gathered earlier
# to indicate wether that cell died, divided, went off frame, was lost, 
# or simply whether that is the last time frame on the movie
endings$ending <- ifelse(endings$death == T, 'death', 
                         ifelse(endings$division == T, 'division', 
                                ifelse(endings$end.point == T, 'end point', 
                                       ifelse(endings$in.frame == F, 'off frame', 
                                              'lost'))))
endings <- endings %>% group_by(Cell_ID, identity_t0, ending) %>% 
  select(Cell_ID, identity_t0, ending)

################################################################################
# Figure 4f
# Survival of intact vs targeted cells, as total number of cells over time 
################################################################################

# Extract Cell_ID for cells that go off frame and cells thad die eventually
off.frame <- unique(endings$Cell_ID[which(endings$ending == 'off frame')])
ded.cells <- unique(endings$Cell_ID[which(endings$ending == 'death')])

# Count the number of cells of each combination of 
# Treatment (Control or Ablated) x cell_treatment (intact or targeted cells)
# for each time point up to 15h (most movies are 15h or less).
# Exclude cells that go off frame
timecounts <- movies %>% filter(!Cell_ID %in% off.frame & 
                                  hours < 15) %>%
  group_by(Treatment, cell_treatment, hours) %>%
  summarize(count = n())

# Calculate the percentage of the total that number represents at each point
timecounts$group <- paste(timecounts$Treatment, 
                          timecounts$cell_treatment, sep = '.')
timecounts <- split(timecounts, as.factor(timecounts$group))
for(t in 1:length(timecounts)){
  timecounts[[t]]$pc.count <- timecounts[[t]]$count / 
    max(timecounts[[t]]$count)
}
timecounts <- do.call(rbind, timecounts)

# Calculate the value of time (x) that corresponds to near 50% survival
my.x <- median(timecounts$hours[which(timecounts$pc.count > 0.45 & 
                                        timecounts$pc.count < 0.55)])
my.y <- median(timecounts$pc.count[which(timecounts$pc.count > 0.45 & 
                                    timecounts$pc.count < 0.55)])
my.y <- round(my.y, digits = 2)

# Generate plot
panel.f <- ggplot(data = timecounts,
                aes(x = hours, y = pc.count))
panel.f <- panel.f + geom_line(aes(color = interaction(Treatment, 
                                                       cell_treatment)), 
                           size = 0.5)
panel.f <- panel.f + geom_segment(aes(x = 0, xend = my.x,
                                  y = my.y, yend = my.y),
                              size = 0.35, linetype = 5)
panel.f <- panel.f + geom_segment(aes(x = my.x, xend = my.x,
                                 y = 0, yend = my.y),
                             size = 0.5, linetype = 5)
panel.f <- panel.f + scale_color_manual(values = c('Control.intact' = 'green', 
                                               'Ablated.intact' = 'darkgreen', 
                                               'Ablated.targeted' = 'gray60'))
panel.f <- panel.f + scale_x_continuous(breaks = seq(0, 15, by = 1))
panel.f <- panel.f + looks + theme(aspect.ratio = 0.5)
panel.f <- panel.f + labs(y = 'Survival (% cells)', 
                      color = 'Embryo & treatment', 
                      title = 'Figure 4f', x = 'Time (hours)')
# Uncomment print() below to generate plot
# print(panel.f)

################################################################################
# Figure 4g
# Survival of example subset of cells surrounding a targeted cell
# Plot Y coordinates over time, offset by initial X position
################################################################################

# Extract a TrackID for a subset of cells for the embryo shown in the image,
# which are located around two given targeted cells
my.embryo <- "021119Abl_EK1"
my.cells <- unique(subset(movies, Embryo_ID == my.embryo &  
                            timeframe == 1 &
                            X < 28 & X > 15 & 
                            Y > 15 & Y < 35)$TrackID)

p.data <- movies %>% filter(Channel == 1, 
                            Embryo_ID == my.embryo, 
                            branch %in% c(0, 1))

# Calculate the relative X and Y for each cell over time
p.data <- split(p.data, as.factor(p.data$TrackID))
for(c in 1:length(p.data)) { 
  # Make X and Y relative to X0 and Y0 for each track
  x0 <- p.data[[c]]$X[which(p.data[[c]]$timeframe == 1)]
  y0 <- p.data[[c]]$Y[which(p.data[[c]]$timeframe == 1)]
  p.data[[c]]$x0 <- x0
  p.data[[c]]$y0 <- y0
  p.data[[c]]$Xrel <- p.data[[c]]$X - x0
  p.data[[c]]$Yrel <- p.data[[c]]$Y - y0
}

# Calculate the moving average of X and Y for each cell 
# over 1h bins (4x time frames)
k <- 4
for(c in 1:length(p.data)){ 
  my.check <- length(p.data[[c]]$Embryo_ID) < k+1
  if(my.check == T) {
    p.data[[c]]$mavgX <- p.data[[c]]$X
    p.data[[c]]$mavgY <- p.data[[c]]$Y
  }
  else { 
    x.zoo <- zoo(x = p.data[[c]]$X, 
                 order.by = p.data[[c]]$timeframe)
    y.zoo <- zoo(x = p.data[[c]]$Y, 
                 order.by = p.data[[c]]$timeframe)
    mavgX <- rollmean(x.zoo, k = k, fill = 'extend', align = 'left')
    mavgY <- rollmean(y.zoo, k = k, fill = 'extend', align = 'left')
    timeframe <- p.data[[c]]$timeframe[order(p.data[[c]]$timeframe)]
    mavgX <- data.frame(mavgX, timeframe)
    mavgY <- data.frame(mavgY, timeframe)
    p.data[[c]] <- merge(p.data[[c]], mavgX)
    p.data[[c]] <- merge(p.data[[c]], mavgY)
  }
}
p.data <- do.call(rbind, p.data)

# Generate plot
panel.g <- ggplot(data = p.data %>% 
                   filter(timeframe == 1, 
                          TrackID %in% my.cells), 
                 aes(x = x0 + timeframe, y = mavgY * 0.277))
panel.g <- panel.g + geom_point(aes(color = cell_treatment), 
                              size = 15) + 
  geom_line(data = p.data %>% 
              filter(TrackID %in% my.cells), 
            aes(color = cell_treatment, 
                x = x0 + timeframe, y = mavgY * 0.277, 
                group = TrackID), size = 1)
panel.g <- panel.g + looks + 
  scale_color_manual(values = c('intact' = 'darkgreen',
                                'targeted' = 'gray60')) + 
  theme(aspect.ratio = 1/2) + 
  labs(title = 'Figure 4g', color = 'Cell type', 
       y = 'Y coordinate', x = 'Time (hours)') +
  scale_x_continuous(breaks = seq(15, 85, by = 5)) + 
  scale_y_continuous(breaks = seq(2, 14, by = 2))
# Uncomment print() below to visualize plot
# print(panel.g)

################################################################################
# Generate PDF
################################################################################

pdf(file = './figures/fig4all_NS.pdf', width = 11, paper = 'a4r')
print(panel.f)
print(panel.g)
dev.off()


##