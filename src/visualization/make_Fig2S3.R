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

# Source necssary functions
source('./src/functions/stage.R')

# Make vectors of datasets
datasets <- c('NComms', 'Spry4', 'new-lms')
# and colors to use for each of them
deez.cols <- c('#a6cee3', '#1f78b4', '#b2df8a')

################################################################################
# Figure 2c
# Growth curves for each ICM lineage in three datasets
# Saiz et al (2016) Nat Comms
# Morgani et al (2018) Dev Bio
# This study
################################################################################

# Define data to be plotted 
my.data <- list(allcounts.list[[1]], allcounts.list[[2]], 
                rbind(allcounts.list[[3]], 
                      allcounts.list[[4]]))


## allcounts.list[[1]] are Nat Comms data (Saiz et al)
## allcounts.list[[2]] are Spry4 data (Morgani et al)
## allcounts.list[[3]] are littermates from ablations experiments
## allcounts.list[[4]] are newly generated littermates, including embryos
## for a Pdgfra allelic series 

## Re-stage embryos in each dataset to generate 10 cell bins
## and select only ICM cells from wild type embryos of up to 160 cells

for(e in 1:length(my.data)) { 
  my.data[[e]] <- stage(my.data[[e]], bin = 15)
  my.data[[e]] <- subset(my.data[[e]], TE_ICM == 'ICM' & 
                               Genotype1 == 'wt' & 
                               Genotype2 == 'wt' & 
                               Cellcount < 165 & 
                           Stage != '16_32')
  my.data[[e]] <- dcast(my.data[[e]], Experiment +
                              Embryo_ID + Cellcount +
                              TE_ICM + Stage + icm.count +
                              Genotype1 + Genotype2 +
                              Gene1 + Gene2 + Background +
                              prepi.ratio +
                              Stage2 ~ Identity,
                            value.var = 'pc.icm')
}

## Plot growth curves for each ICM cell type over time
## For each dataset separately, as labelled

my.grocs <- ggplot() +
  geom_line(data = subset(my.data[[1]]), 
            aes(x = Stage2, y = PRE), stat = 'summary', fun.y = 'median', 
            color = "#1a33e5", linetype = 'solid', group = 1) + 
  geom_line(data = subset(my.data[[1]]), 
            aes(x = Stage2, y = all.EPI), stat = 'summary', fun.y = 'median', 
            color = "#e5331a", linetype = 'solid', group = 1) + 
  geom_line(data = subset(my.data[[1]]), 
            aes(x = Stage2, y = DP), stat = 'summary', fun.y = 'median', 
            color = "#803380", linetype = 'solid', group = 1)
my.grocs <- my.grocs + labs(x = 'Total cell count', y = '% of ICM')
my.grocs <- my.grocs + looks + theme(aspect.ratio = 0.75, 
                                     axis.text.x = element_text(angle = 45, 
                                                                hjust = 1))
my.grocs <- my.grocs + scale_color_manual(values = c(idcols, 'all.EPI' = 'red'))
my.grocs <- my.grocs + labs(title = 'Saiz et al. (2016) Nat Comms')
my.grocs <- my.grocs + ylim(0, 100)
print(my.grocs)

my.grocs <- ggplot() +
  geom_line(data = subset(my.data[[2]]), 
            aes(x = Stage2, y = PRE), stat = 'summary', fun.y = 'median', 
            color = "#1a33e5", linetype = 'dashed', group = 1) + 
  geom_line(data = subset(my.data[[2]]), 
            aes(x = Stage2, y = all.EPI), stat = 'summary', fun.y = 'median', 
            color = "#e5331a", linetype = 'dashed', group = 1) + 
  geom_line(data = subset(my.data[[2]]), 
            aes(x = Stage2, y = DP), stat = 'summary', fun.y = 'median', 
            color = "#803380", linetype = 'dashed', group = 1)
my.grocs <- my.grocs + labs(x = 'Total cell count', y = '% of ICM')
my.grocs <- my.grocs + looks + theme(aspect.ratio = 0.75, 
                                     axis.text.x = element_text(angle = 45, 
                                                                hjust = 1))
my.grocs <- my.grocs + scale_color_manual(values = c(idcols, 'all.EPI' = 'red'))
my.grocs <- my.grocs + labs(title = 'Morgani et al (2018) Dev Bio')
my.grocs <- my.grocs + ylim(0, 100)
print(my.grocs)

my.grocs <- ggplot() +
  geom_line(data = subset(my.data[[3]]), 
            aes(x = Stage2, y = PRE), stat = 'summary', fun.y = 'median', 
            color = "#1a33e5", linetype = 'dotted', group = 1) + 
  geom_line(data = subset(my.data[[3]]), 
            aes(x = Stage2, y = all.EPI), stat = 'summary', fun.y = 'median', 
            color = "#e5331a", linetype = 'dotted', group = 1) + 
  geom_line(data = subset(my.data[[3]]), 
            aes(x = Stage2, y = DP), stat = 'summary', fun.y = 'median', 
            color = "#803380", linetype = 'dotted', group = 1)
my.grocs <- my.grocs + labs(x = 'Total cell count', y = '% of ICM')
my.grocs <- my.grocs + looks + theme(aspect.ratio = 0.75, 
                                     axis.text.x = element_text(angle = 45, 
                                                                hjust = 1))
my.grocs <- my.grocs + scale_color_manual(values = c(idcols, 'all.EPI' = 'red'))
my.grocs <- my.grocs + labs(title = 'This study')
my.grocs <- my.grocs + ylim(0, 100)
print(my.grocs)

# Generate plot
fig.2d <- ggplot(data = my.data[[1]] %>% 
                   filter(Identity %in% c('DP', 'PRE', 'all.EPI')), 
                 aes(x = Cellcount, y = pc.icm))
fig.2d <- fig.2d + geom_smooth(aes(color = Identity), 
                               size = 1, span = 0.8, 
                               linetype = 'solid', se = F) + 
  geom_smooth(data = my.data[[2]] %>% 
                filter(Identity %in% c('DP', 'PRE', 'all.EPI')), 
              aes(x = Cellcount, y = pc.icm, color = Identity), 
              size = 1, span = 0.8, linetype = 'dashed', se = F) + 
  geom_smooth(data = my.data[[3]] %>% 
                filter(Identity %in% c('DP', 'PRE', 'all.EPI')), 
              aes(x = Cellcount, y = pc.icm, color = Identity), 
              size = 1, span = 0.8, linetype = 'dotted', se = F)
fig.2d <- fig.2d + 
  looks + 
  theme(aspect.ratio = 0.75) + 
  scale_color_manual(values = idcols)
print(fig.2d)

################################################################################
# Supplementary Figure 3a
# Box plot showing PrE:EPI ratio over time (at sequential stages of development)
# for data generated in this study
################################################################################

# Define data to be plotted
my.data <- rbind(allcounts.list[[3]], 
                 allcounts.list[[4]])
# Slim down data to select wild type embryos > 31 cells with both PrE and EPI
my.data <- my.data %>% 
  filter(!Stage %in% c('16_32'), 
         Cellcount < 165, 
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
# Supplementary Figure 3 b and c
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
           Cellcount < 165, 
           Genotype1 == 'wt', 
           Genotype2 == 'wt', 
           is.na(prepi.ratio) == F) %>% 
    group_by(Embryo_ID, Experiment, 
             Stage, prepi.ratio) %>% 
    summarize()
  my.data[[i]]$dataset <- datasets[i]
}

# Generate plots
# Panel b
panel.S3b <- ggplot(data = my.data[[1]], 
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
  scale_fill_manual(values = deez.cols[1]) + 
  labs(title = unique(my.data[[1]]$dataset))
# Uncomment print() below to visualize plot
# print(panel.S3b)

# Panel c
panel.S3c <- ggplot(data = my.data[[2]], 
                    aes(x = Stage, y = log2(prepi.ratio)))
panel.S3c <- panel.S3c + 
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
panel.S3c <- panel.S3c + 
  looks + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = deez.cols[2]) + 
  labs(title = unique(my.data[[2]]$dataset))
# Uncomment print() below to visualize plot
# print(panel.S3c)

################################################################################
# Figure S3d 
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
           Cellcount < 165, 
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
    # filter(log2(prepi.ratio) < 4 & log2(prepi.ratio) > -4) %>% 
    summarize(ratio.sd = sd(prepi.ratio), 
              ratio.avg = mean(prepi.ratio))
  my.var[[i]]$ratio.cv <- my.var[[i]]$ratio.sd / my.var[[i]]$ratio.avg * 100
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

# Generate PDF with plots
pdf(file = './figures/FigS3_NS.pdf')
print(panel.a)
print(panel.S3b)
print(panel.S3c)
print(cv.ratio)
dev.off()

##