# This script generates the plots in Figure 2 
# corresponding to experimental data (a, b and d)
# and in associated Supplementary Figure 3

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load necessary data
# Read in or generate reference littermates datasets
if (exists('allcounts.list') == F) { 
  source('./src/data/all-lms_read.R')
}

# Make vectors of datasets
datasets <- c('NComms', 'Spry4', 'new-lms')
# and colors to use for each of them
deez.cols <- c('#a6cee3', '#1f78b4', '#b2df8a')

################################################################################
# Figure 2a
# Box plot showing PrE:EPI ratio over time (at sequential stages of development)
# for data generated in this study
################################################################################

# Define data to be plotted
my.data <- rbind(allcounts.list[[3]], allcounts.list[[4]])
# Slim down data to select wild type embryos > 31 cells with both PrE and EPI
my.data <- my.data %>% 
  filter(!Stage %in% c('16_32'), 
         Genotype1 == 'wt', 
         Genotype2 == 'wt', 
         is.na(prepi.ratio) == F) %>% 
  group_by(Embryo_ID, Experiment, 
           Stage, prepi.ratio) %>% 
  summarize()
my.data$dataset <- datasets[3]

# Generate plot
panel.a <- ggplot(data = my.data, 
                  aes(x = Stage, y = log2(prepi.ratio)))
panel.a <- panel.a + 
  geom_rect(ymin = log2(1.2), ymax = log2(1.5), 
                               xmin = -Inf, xmax = Inf,
                               fill = 'gray75') + 
  geom_boxplot(aes(fill = dataset), 
               color = 'black', 
               outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F) + 
  geom_jitter(color = 'black', 
              shape = 20, width = 0.2)
panel.a <- panel.a + 
  looks + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = deez.cols[3]) + 
  labs(title = unique(my.data$dataset))
# Uncomment print() below to visualize plot
# print(panel.a)

################################################################################
# Supplementary Figure 3 a and b
# Box plot showing PrE:EPI ratio over time (at sequential stages of development)
# for data published in Saiz et al (2016) Nat Comms and 
# Morgani et al (2018) Dev Bio
################################################################################

# Define data to be plotted
my.data <- list(allcounts.list[[1]], 
                allcounts.list[[2]])
# Slim down data to select wild type embryos > 31 cells with both PrE and EPI
for(i in 1:2) { 
  my.data[[i]] <- my.data[[i]] %>% 
    filter(!Stage %in% c('16_32'), 
           Genotype1 == 'wt', 
           Genotype2 == 'wt', 
           is.na(prepi.ratio) == F) %>% 
    group_by(Embryo_ID, Experiment, 
             Stage, prepi.ratio) %>% 
    summarize()
  my.data[[i]]$dataset <- datasets[i]
}

# Generate plots
# Panel a
panel.S3a <- ggplot(data = my.data[[1]], 
                    aes(x = Stage, y = log2(prepi.ratio)))
panel.S3a <- panel.S3a + 
  geom_rect(ymin = log2(1.2), ymax = log2(1.5), 
            xmin = -Inf, xmax = Inf,
            fill = 'gray75') + 
  geom_boxplot(aes(fill = dataset), 
               color = 'black', 
               outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F) + 
  geom_jitter(color = 'black', 
              shape = 20, width = 0.2)
panel.S3a <- panel.S3a + 
  looks + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = deez.cols[1]) + 
  labs(title = unique(my.data[[1]]$dataset))
# Uncomment print() below to visualize plot
# print(panel.S3a)

# Panel b
panel.S3b <- ggplot(data = my.data[[2]], 
                    aes(x = Stage, y = log2(prepi.ratio)))
panel.S3b <- panel.S3b + 
  geom_rect(ymin = log2(1.2), ymax = log2(1.5), 
            xmin = -Inf, xmax = Inf,
            fill = 'gray75') + 
  geom_boxplot(aes(fill = dataset), 
               color = 'black', 
               outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F) + 
  geom_jitter(color = 'black', 
              shape = 20, width = 0.2)
panel.S3b <- panel.S3b + 
  looks + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = deez.cols[2]) + 
  labs(title = unique(my.data[[2]]$dataset))
# Uncomment print() below to visualize plot
# print(panel.S3b)

################################################################################
# Figure 2b 
# Line plot showing coefficient of variation over time for all three datasets
# This study, Saiz et al (2016) and Morgani et al (2018)
################################################################################

# Define data to be plotted 
my.data <- list(allcounts.list[[1]], allcounts.list[[2]], 
                rbind(allcounts.list[[3]], 
                      allcounts.list[[4]]))
# Slim down data to select wild type embryos > 31 cells with both PrE and EPI
for(i in 1:3) { 
  my.data[[i]] <- my.data[[i]] %>% 
    filter(!Stage %in% c('16_32'), 
           Genotype1 == 'wt', 
           Genotype2 == 'wt', 
           is.na(prepi.ratio) == F) %>% 
    group_by(Embryo_ID, Experiment, 
             Stage, prepi.ratio) %>% 
    summarize()
  my.data[[i]]$dataset <- datasets[i]
}

# Calculate mean(ratio), SD(ratio) and coefficient of variation
# for each stage and dataset
my.var <- list()
for(i in 1:3) { 
  my.var[[i]] <- my.data[[i]] %>% group_by(Stage, dataset) %>%
    summarize(ratio.sd = sd(prepi.ratio), 
              ratio.avg = mean(prepi.ratio))
  my.var[[i]]$ratio.cv <- my.var[[i]]$ratio.sd / my.var[[i]]$ratio.avg
}

# Generate plot
cv.ratio <- ggplot(data = rbind(my.var[[1]], my.var[[2]], my.var[[3]]), 
                   aes(x = Stage, y = ratio.cv)) +
  # geom_point(aes(color = dataset)) + 
  geom_path(aes(color = dataset, group = dataset)) + 
  looks + theme(aspect.ratio = 1) + 
  scale_color_manual(values = deez.cols)
# Uncomment print() below to visualize plot
# print(cv.ratio)


