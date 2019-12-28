# This script generates the plots in Figure 6 of the paper

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

if(exists('ablat.lincounts') == F) { 
  ablat.lincounts <- read.csv('./data/processed/ablat-counts.csv')
}

if(exists('movies') == F) {
  movies <- read.csv('./data/processed/movies-all-processed.csv')
}

# Load extra packages
library('RColorBrewer')

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
ablat.lincounts$Stage.t0 <- factor(ablat.lincounts$Stage.t0, 
                                   levels = c('[30,50)', '[50,70)', 
                                              '[70,90)', '[90,110)', 
                                              '[110,130)'))
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
# Figure 6a
# Relative ICM composition of reference fixed embryos 
# at sequential stages of development, at which cell ablation was performed
################################################################################

# Define data to plot
my.data <- ablat %>% filter(Treatment == 'Littermate', 
                            TE_ICM == 'ICM')

# Generate table with N numbers for data in this panel
fig.6a.N <- my.data %>% group_by(Embryo_ID, Stage.t0) %>% summarize() %>%
  group_by(Stage.t0) %>% summarize(N = n())
# And write out to disk
write.csv(fig.6a.N, file = './results/fig6a_N-numbers.csv', row.names = F)

# Generate plot
panel.a <- ggplot(data = my.data, 
                  aes(x = Stage.t0, fill = Identity.hc))
panel.a <- panel.a + geom_bar(position = 'fill')
panel.a <- panel.a + looks + scale_fill_manual(values = idcols)
panel.a <- panel.a + theme(aspect.ratio = 9/5, 
                             axis.text.x = element_text(angle = 45, 
                                                        hjust = 1))
panel.a <- panel.a + labs(x = 'Stage', y = '% of ICM',
                          title = 'Figure 6a')
# Uncomment print() below to visualize plot
# print(panel.a)

################################################################################
# Figure 6c
# PrE:EPI log(ratio) over developmental time, for each % of ablated PrE cells.
# Best fitting lines (linear model) for each group overlaid in the same color
# Gray box spans the log(1.2)-log(1.5) region (corresponding to composition 
# between 60%PrE:40% EPI and 55%PrE:45%EPI)
# Random ablation data is overlaid as another layer, in dark green (no legend)
################################################################################

# Calculate PrE:EPI ratio for each experimental group

## Extract PrE and EPI counts from the lineage counts table
ablat.ratio <- subset(ablat.lincounts, Identity.hc %in% c('EPI', 'PRE') & 
                        Recovery %in% c('24h', '0min') &  
                        # Stage.t0 != '[110,130)' & 
                        Exp_date > 20170601)
## Restructure to wide format
ablat.ratio <- dcast(ablat.ratio, Experiment + Litter + Embryo_ID + Cellcount + 
                       Treatment + target + icm.count + group.median + 
                       Genotype1 + Gene1 + Genotype2 + Gene2 + 
                       Cell_diff + Recovery + Stage.t0 + Cellcount_t0 + 
                       litter.median.t0 ~ Identity.hc, 
                     value.var = 'pc.icm')
## Divide the number of PrE by the number of EPI cells to obtain the ratio
## 60% PrE / 40% EPI = 1.5; 55% PrE / 45% EPI ~ 1.2
ablat.ratio$ratio <- ablat.ratio$PRE / ablat.ratio$EPI

ablat.ratio$Treatment <- factor(ablat.ratio$Treatment, 
                                levels = c('Littermate', 'Control', 
                                           'Ablated', 'Random'))
ablat.ratio$Cell_diff <- factor(ablat.ratio$Cell_diff, 
                                levels = c('zero', 'kill_0.25', 'kill_0.5', 
                                           'kill_0.75', 'kill_all'))

# Define data to plot
ratiodata <- ablat.ratio %>% 
  filter(Treatment != 'Littermate', 
         litter.median.t0 %in% seq(50, 110, by = 1), 
         Cell_diff %in% c('zero', 'kill_0.5', 
                          'kill_0.75', 'kill_all'), 
         target %in% c('none', 'PRE', 'EPI'))

# Generate plot
panel.c <- ggplot(data = subset(ratiodata, Treatment != 'Random' & 
                                    target %in% c('none', 'PRE')), 
                    aes(x = litter.median.t0, 
                        y = log(ratio)))
panel.c <- panel.c + geom_rect(ymin = log(1.2), ymax = log(1.5), 
                                   xmin = -Inf, xmax = Inf,
                                   fill = 'gray85')
panel.c <- panel.c + geom_jitter(aes(color = Cell_diff), 
                                     shape = 20, size = 2, width = 0.2)
panel.c <- panel.c + geom_smooth(aes(color = Cell_diff), 
                                     size = 1, method = 'lm', se = F)
panel.c <- panel.c + 
  geom_jitter(data = subset(ratiodata,
                            Treatment == 'Random' & 
                              target != 'EPI'), 
              aes(x = litter.median.t0, 
                  y = log(ratio)),
              color = 'gray55', shape = 20, 
              size = 2, width = 0.2)
panel.c <- panel.c + 
  geom_smooth(data = subset(ratiodata,
                            Treatment == 'Random' &
                              target != 'EPI'), 
              aes(x = litter.median.t0, 
                  y = log(ratio)),
              color = 'gray55', size = 1, 
              method = 'lm', se = F)
panel.c <- panel.c + 
  looks + 
  scale_color_manual(values = brewer.pal(9, 'Blues')[c(3, 5, 7, 9)]) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'Figure 6c - PrE ablation', y = 'log(PrE:EPI ratio)', 
       x = 'Median cell count of litter at start') + 
  ylim(-3, 3)
# Uncomment print() below to visualize plot
# print(panel.c)

################################################################################
# Figure 6d
# Same as above, for EPI ablations
################################################################################

# Geneate plot
panel.d <- ggplot(data = subset(ratiodata, Treatment != 'Random' &
                                    Cell_diff != 'kill_0.75' & 
                                    target %in% c('none', 'EPI')), 
                    aes(x = litter.median.t0, 
                        y = log(ratio)))
panel.d <- panel.d + geom_rect(ymin = log(1.2), ymax = log(1.5), 
                                   xmin = -Inf, xmax = Inf,
                                   fill = 'gray85')
panel.d <- panel.d + geom_jitter(aes(color = Cell_diff), 
                                     shape = 20, size = 2, width = 0.2)
panel.d <- panel.d + geom_smooth(aes(color = Cell_diff), 
                                     size = 1, method = 'lm', se = F)
panel.d <- panel.d + 
  geom_jitter(data = subset(ratiodata,
                            Treatment == 'Random' &
                              target != 'PRE'), 
              aes(x = litter.median.t0, 
                  y = log(ratio)),
              color = 'gray55', shape = 20, 
              size = 2, width = 0.2)
panel.d <- panel.d + 
  geom_smooth(data = subset(ratiodata, 
                            Treatment == 'Random' & 
                              target != 'PRE'), 
              aes(x = litter.median.t0, 
                  y = log(ratio)),
              color = 'gray55', size = 1, 
              method = 'lm', se = F)
panel.d <- panel.d + 
  looks + 
  scale_color_manual(values = brewer.pal(9, 'Reds')[c(3, 6, 9)]) +
  theme(aspect.ratio = 1) + 
  labs(title = 'Figure 6d - EPI ablation', y = 'log(PrE:EPI ratio)', 
       x = 'Median cell count of litter at start') + 
  ylim(-3, 3)
# Uncomment print() below to visualize plot
# print(panel.d)

################################################################################
# Figure 6g
# Pdgfra:H2B-GFP levels of intact DP cells over time 
# for different experimental groups targeted at the 70-90 cell stage
################################################################################

# Define possible stages to plot
stages <- c('[50,70)', '[70,90)', '[90,110)')

# Index for the stage to be plotted for the figure ([70-90) cells at t0)
s <- 2

# Define data to plot
my.data <- movies %>% filter(Channel == 1, 
                             cell_treatment == 'intact', 
                             target %in% c('none', 'PRE', 'EPI'), 
                             # Cell_diff %in% c('zero', 'kill_0.5', 
                             #                  'kill_all'), 
                             Stage.t0 == stages[s], 
                             # Genotype2 == 'mKate/mKate'' &', 
                             hours < 16)

# Generate plot
panel.g <- ggplot(data = my.data %>% 
                    filter(identity_t0.th == 'DP'), 
                aes(x = hours, y = mavg))
panel.g <- panel.g + geom_line(aes(color = identity.t,
                               group = interaction(TrackID, Cell_ID)),
                           size = 0.75, alpha = 0.4)
panel.g <- panel.g + geom_smooth(aes(color = identity_t0.th), size = 2)
# Plot intact PrE and EPI cells as faint background lines, for reference
panel.g <- panel.g +
  geom_line(data = my.data %>% 
              filter(identity_t0.th != 'DP'),
            aes(x = hours, y = mavg, color = identity.t,
                group = interaction(TrackID, Cell_ID)),
            size = 0.75, alpha = 0.1)
# panel.g <- panel.g + 
# geom_smooth(data = my.data %>%
#               filter(identity_t0.th != 'DP'),
#             aes(x = hours, y = mavg, color = identity_t0.th),
#             alpha = 0.5)
panel.g <- panel.g + facet_grid(target ~ Cell_diff)
panel.g <- panel.g + 
  scale_color_manual(values = idcols) + 
  looks + 
  theme(aspect.ratio = 0.75, 
        strip.text.x = element_text(size = 6)) + 
  labs(title = paste('Figure 6g -', 
                     'Ablation at', 
                     stages[s], 
                     'cells', sep = ' '))
# Uncomment print() below to visualize plot
# print(panel.g)

################################################################################
# Figure 6h
# Stacked bar plot showing the end fate of intact cells
# that were DP at the beginning of the experiment
################################################################################

# Define data to plot
my.data <- movies %>% filter(Channel == 1, 
                    end.point == T, 
                    Stage.t0 == stages[s], 
                    identity_t0.th == 'DP')

# Generate plot
panel.h <- ggplot(data = my.data, 
              aes(x = Embryo_ID, fill = identity.t))
panel.h <- panel.h + geom_bar(position = 'fill')
panel.h <- panel.h + 
  scale_fill_manual(values = idcols)+ 
  looks + 
  theme(aspect.ratio = 2, 
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  facet_grid(Stage.t0 ~ target + Cell_diff, 
             scales = 'free') + 
  labs(title = 'Figure 6h', y = 'DP cells end fate (%)', 
       x = 'Individual embryos (as in Fig. S7a-c)')
# Uncomment print() below to visualize plot
# print(panel.h)

################################################################################
# Figure 6i
# Box plot of PrE:EPI log(ratio) at the end of the movie
# in mKate/mKate or mKate/+ embryos targeted at the 70-90 cell stage 
# In mKate/mKate often I was able to follow all/most ICM cells
################################################################################

# Extract end point of mKate/+ or mKate/mKate embryos 
# targeted at the 70-90 cell stage
end.data <- movies %>% filter(Channel == 1, 
                              end.point == T, 
                              Genotype2 != 'wt', 
                              Stage.t0 == stages[s])

# Define data to plot
## Calculate counts for each ICM cell type per embryo
end.counts <- end.data %>% 
  group_by(Embryo_ID, Litter, 
           Experiment, Genotype2, 
           Treatment, target, 
           Cell_diff, Stage.t0, 
           identity.t) %>% 
  summarize(count = n())
## Cast to wide format to calculate PrE/EPI ratio
end.counts <- dcast(end.counts, Experiment + Litter + Embryo_ID + 
                      Treatment + target + 
                      Genotype2 + Cell_diff + 
                      Stage.t0 ~ identity.t, 
                    value.var = 'count')
end.counts$prepi.ratio <- end.counts$PRE / end.counts$EPI

# Order factors for plotting
end.counts$Treatment <- factor(end.counts$Treatment, 
                               levels = c('Littermate', 'Control', 
                                          'Ablated', 'Random'))
end.counts$Cell_diff <- factor(end.counts$Cell_diff, 
                               levels = c('zero', 'kill_0.25', 'kill_0.5', 
                                          'kill_0.75', 'kill_all'))

# Generate plot
panel.i <- ggplot(data = end.counts, 
                 aes(x = Cell_diff, y = prepi.ratio))
panel.i <- panel.i + geom_boxplot(aes(fill = Treatment), color = I('black'), 
                                outlier.shape = 1, outlier.size = 2, 
                                show.legend = F) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)  +
  geom_jitter(aes(color = target), 
              shape = 20, size = 1.5, 
              width = 0.2) + 
  geom_hline(yintercept = 1.5, linetype = 2)
panel.i <- panel.i + 
  looks + 
  scale_color_manual(values = c('none' = 'black', 
                                'PRE' = '#1a33e5', 
                                'EPI' = "#e5331a")) + 
  scale_fill_manual(values = c('Littermate' = 'white', 
                               'Control' = 'gray90', 
                               'Random' = 'gray75', 
                               'Ablated' = 'gray60')) + 
  facet_grid(cols = vars(target), scales = 'free') + 
  theme(aspect.ratio = 9/5) + 
  labs(y = 'log(PrE:EPI ratio)')
# Uncomment print() below to visualize plot
# print(panel.i)

################################################################################
# Generate PDF
################################################################################

pdf(file = './figures/fig6all_NS.pdf', width = 11, paper = 'a4r')
print(panel.a)
print(panel.c)
print(panel.d)
print(panel.g)
print(panel.h)
print(panel.i)
dev.off()


##