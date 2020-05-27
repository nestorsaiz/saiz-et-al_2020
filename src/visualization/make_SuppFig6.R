# This script generates the plots in Supplementary Figure 6 

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load necessary data
if(exists('ablat') == F) {
  ablat <- read.csv('./data/processed/ablat-processed.csv')
}

if(exists('movies') == F) {
  movies <- read.csv('./data/processed/movies-all-processed.csv')
}

# Order factors for plotting purposes
ablat$Identity.hc <- factor(ablat$Identity.hc, 
                            levels = c('TE', 'PRE', 'DP', 'EPI',
                                       'EPI.lo', 'DN'))
ablat.lincounts$Identity.hc <- factor(ablat.lincounts$Identity.hc, 
                                      levels = c('TE', 'PRE', 'DP', 'EPI',
                                                 'EPI.lo', 'DN'))
ablat$Stage.t0 <- factor(ablat$Stage.t0, levels = c('[30,50)', '[50,70)', 
                                                    '[70,90)', '[90,110)', 
                                                    '[110,130)'))
ablat$Treatment <- factor(ablat$Treatment, 
                          levels = c('Littermate', 'Control', 
                                     'Ablated', 'Random'))
ablat$target <- factor(ablat$target, levels = c('none', 'PRE', 'EPI', 
                                                'PRE+DP', 'EPI+DP'))

movies$identity.t <- factor(movies$identity.t, 
                            levels = c('PRE', 'DP', 'EPI'))
movies$Treatment <- factor(movies$Treatment, 
                           levels = c('Control', 'Ablated', 'Random'))
movies$Stage.t0 <- factor(movies$Stage.t0, 
                          levels = c('[30,50)', '[50,70)', '[70,90)', 
                                     '[90,110)', '[110,130)'))
movies$Cell_diff <- factor(movies$Cell_diff, 
                           levels = c('zero', 'kill_0.25', 'kill_0.5', 
                                      'kill_0.75', 'kill_all'))
movies$identity_t0.th <- factor(movies$identity_t0.th, 
                                levels = levels(movies$identity.t))
movies$target <- factor(movies$target, levels = c('none', 'PRE', 'EPI'))

################################################################################
# Figure S6a
# Box plots showing the absolute number of cells killed at each stage
# and experimental group
################################################################################

# Define the data to be plotted
# by merging the number of cells hit per embryo from the exp reference
# and the metadata for each embryo from live.lincounts.
# Select only embryos for the desired experimental groups ant targets
my.data <- merge(
  select(ablat.ref, Embryo_ID, target, 
         target_t0, target_killed, Cell_diff), 
  live.lincounts
) %>% 
  filter(Stage.t0 %in% c('[50,70)', '[70,90)', '[90,110)'), 
         target %in% c('EPI', 'PRE'))

# Make Cell_diff a factor to plot in the right order
my.data$Cell_diff <- factor(my.data$Cell_diff, 
                            levels = c('zero', 'kill_0.5', 
                                       'kill_0.75', 'kill_all'))

# Generate plot
panel.a <- ggplot(data = my.data, 
                 aes(x = Cell_diff, y = target_killed))
panel.a <- panel.a + 
  geom_boxplot(aes(fill = target), color = 'black',
               outlier.shape = 1, outlier.size = 2) + 
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F) + 
  geom_jitter(color = 'black', shape = 20, width = 0.2)
panel.a <- panel.a + 
  facet_grid( ~ target + Stage.t0) +
  looks + 
  scale_fill_manual(values = idcols) + 
  theme(aspect.ratio = 1, 
        strip.text.x = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = 'No of cells targeted', 
       x = 'Experimental group', 
       title = 'Figure S6a')
# Uncomment print() below to visualize plot  
print(panel.a)

################################################################################
# Figure S6b
# Box plots showing the total cell counts at the end of the experiment
# for each experimental group
################################################################################

# Extract total cell number for each experimental group and experiment
cellcounts <- ablat %>% filter(Recovery %in% c('24h', '0min'), 
                               !Stage.t0 %in% c('[30,50)', '[110,130)'), 
                               Cell_diff != 'kill_0.25', 
                               Exp_date > 20170601, 
                               target %in% c('none', 'PRE', 'EPI')) %>% 
  group_by(Embryo_ID, Cellcount, Exp_date, 
           Treatment, target, target_killed,
           Litter, group.median, 
           Experiment, Gene1, Genotype1, 
           Stage.t0, icm.count, Cell_diff) %>% 
  summarize()
cellcounts <- data.frame(cellcounts)

# Arrange Litters in ascending order of Littermate average cellcount
cellcounts$Litter <- factor(cellcounts$Litter)
# Order factors for plotting
cellcounts$Cell_diff <- factor(cellcounts$Cell_diff, 
                               levels = c('zero', 'kill_0.25', 'kill_0.5',
                                          'kill_0.75', 'kill_all'))
cellcounts$target <- factor(cellcounts$target, 
                            levels = c('none', 'PRE', 'EPI'))

# Generate plot
panel.b <- ggplot(data = cellcounts %>% 
                    filter(target %in% c('none', 'PRE') | 
                             target == 'EPI' & 
                             Cell_diff != 'kill_0.75'), 
                  aes(x = Stage.t0, 
                      y = Cellcount))
panel.b <- panel.b + geom_boxplot(aes(fill = Treatment), color = I('black'), 
                                  outlier.shape = 1, outlier.size = 2, 
                                  show.legend = F) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)+ 
  geom_jitter(aes(color = target), 
              shape = 20, size = 1.5, width = 0.2)
panel.b <- panel.b + 
  looks + 
  facet_grid( ~ target + Treatment + Cell_diff) + 
  scale_fill_manual(values = c('Littermate' = 'white', 
                               'Control' = 'gray90', 
                               'Random' = 'gray75', 
                               'Ablated' = 'gray60')) + 
  scale_color_manual(values = c('none' = 'black', idcols)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        aspect.ratio = 9/3) + 
  labs(y = 'Total cell count', x = NULL, 
       title = 'Figure S6b') + 
  scale_y_continuous(breaks = seq(30, 150, by = 15), 
                     limits = c(30, 140))
# Uncomment print() below to visualize plot
# print(panel.b)

################################################################################
# Figure S6c
# Box plots showing the number of ICM cells at the end of the experiment
# for each experimental group
################################################################################

# Generate plot
panel.c <- ggplot(data = cellcounts %>% 
                    filter(target %in% c('none', 'PRE') | 
                             target == 'EPI' & 
                             Cell_diff != 'kill_0.75'), 
                  aes(x = Stage.t0, 
                      y = icm.count))
panel.c <- panel.c + geom_boxplot(aes(fill = Treatment), color = I('black'), 
                                  outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F) + 
  geom_jitter(aes(color = target), 
              shape = 20, size = 1.5, width = 0.2)
panel.c <- panel.c + 
  looks + 
  facet_grid( ~ target + Treatment + Cell_diff) + 
  scale_fill_manual(values = c('Littermate' = 'white', 
                               'Control' = 'gray90', 
                               'Random' = 'gray75', 
                               'Ablated' = 'gray60')) + 
  scale_color_manual(values = c('none' = 'black', idcols)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        aspect.ratio = 9/3) + 
  scale_y_continuous(breaks = seq(0, 50, by = 10), 
                     limits = c(0, 50)) + 
  labs(y = 'ICM cells', x = NULL, 
       title = 'Figure S6c')
# Uncomment print() below to visualize plot
# print(panel.c)

################################################################################
# Figure S6d
# Stacked bar plots showing average relative ICM composition at the end
# of the experiment for each experimental group
################################################################################

# Define data to plot
my.data <- ablat %>% filter(TE_ICM == 'ICM', 
                            Treatment != 'Littermate',
                            Recovery %in% c('24h', '0min'),  
                            !Stage.t0 %in% c('[30,50)', '[110,130)'),
                            Cell_diff != 'kill_0.25', 
                            Exp_date > 20170601, 
                            target %in% c('none', 'PRE') | 
                              target == 'EPI' & 
                              Cell_diff != 'kill_0.75')

## Average, grouped by stage at time of ablation, average per treatment group
panel.d <- ggplot(data = my.data, 
                  aes(x = Stage.t0, 
                      fill = Identity.hc))
panel.d <- panel.d + geom_bar(position = 'fill') + 
  geom_hline(yintercept = 0.40, linetype = 2)
panel.d <- panel.d + 
  looks + 
  facet_grid( ~ target + Treatment + Cell_diff) + 
  scale_fill_manual(values = idcols) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        aspect.ratio = 9/3) + 
  labs(y = '% of ICM', fill = 'Identity', 
       x = 'Stage at time of cell ablation (total cell count)', 
       title = 'Figure S6d')
# Uncomment print() below to visualize plot
# print(panel.d)

################################################################################
# Figure S6e
# Stacked bar plots showing average relative ICM composition
# at the end of the movie in mKate/mKate or mKate/+ embryos 
# targeted at the 70-90 cell stage only
# In mKate/mKate often I was able to follow all/most ICM cells
################################################################################

# Define possible stages to plot
stages <- c('[50,70)', '[70,90)', '[90,110)')

# Index for the stage to be plotted for the figure ([70-90) cells at t0)
s <- 2

# Extract end point of mKate/+ or mKate/mKate embryos 
# targeted at the 70-90 cell stage
end.data <- movies %>% filter(Channel == 1, 
                              end.point == T, 
                              Genotype2 != 'wt', 
                              Stage.t0 == stages[s])

# Generate plot
panel.e <- ggplot(data = end.data, 
                  aes(x = Cell_diff, fill = identity.t))
panel.e <- panel.e + geom_bar(position = 'fill') + 
  geom_hline(yintercept = 0.40, linetype = 2)
panel.e <- panel.e + 
  scale_fill_manual(values = idcols) + 
  looks + 
  facet_grid(cols = vars(target), scales = 'free') + 
  theme(aspect.ratio = 9/3) + 
  labs(y = '% of ICM at end of movie', fill = 'Identity', 
       title = 'Figure S6e')
# Uncomment print() below to visualize plot
# print(panel.e)

################################################################################
# Generate PDF
################################################################################

pdf(file = './figures/figS6all_NS.pdf', width = 11, paper = 'a4r')
print(panel.a)
print(panel.b)
print(panel.c)
print(panel.d)
print(panel.e)
dev.off()


##