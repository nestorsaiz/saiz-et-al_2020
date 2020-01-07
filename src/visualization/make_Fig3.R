# This script generates the plots contained in Figure 3 of the paper

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Read in ESC chimeras data
if (exists('esc.chimeras') ==  F) { 
  esc.chimeras <- read.csv('./data/processed/esc-xim-processed.csv')
  esc.ref <- read.csv('./references/esc-xim_exp_ref.csv')
  esc.xim.lincounts <- read.csv('./data/processed/esc-xim-counts.csv')
  esc.xim.lincounts2 <- read.csv('./data/processed/esc-xim-widecounts.csv')
}

# If there is no raw data loaded, generate it from scratch
if (exists('esc.chimeras') ==  F) { 
  source('./src/esc-xim_runall.R')
}

# Order Identity levels for plotting
esc.chimeras$Identity.hc <- factor(esc.chimeras$Identity.hc, 
                                   levels = c('TE', 'PRE', 'DP', 
                                              'EPI', 'EPI.lo', 'DN', 'ESC'))
esc.xim.lincounts$Identity.hc <- factor(esc.xim.lincounts$Identity.hc, 
                                        levels = c('TE', 'PRE', 'DP', 
                                                   'EPI', 'EPI.lo', 'DN', 'ESC'))

################################################################################
# Generate table with embryo counts for reference
################################################################################

# List litters to be excluded from analysis because embryos are too small
# (see annotation in esc-xim_read.R)
small.odd <- unique(subset(esc.chimeras, 
                           Exp_date >= 20180214 & 
                             Exp_date <= 20180222)$Litter)

# Extract final #ESCs per embryo from esc.xim.lincounts2
esc.end <- select(esc.xim.lincounts2, Embryo_ID, esc.end)

# Generate table with embryos included in Fig 3 and their cellcounts
fig.3_nums <- merge(esc.chimeras, esc.end) %>% 
  filter(TE_ICM == 'ICM', Genotype1 == 'wt', 
         Stage.t0 == '8cell', ES_culture == 'S/LIF', 
         !Litter %in% small.odd)  %>%
  group_by(Embryo_ID, Litter, Treatment, ESC_line, 
           esc.end, Cellcount, icm.count) %>%
  summarize()

# Generate table with N numbers for embryos in Fig 3
fig.3N <- fig.3_nums %>% 
  group_by(Treatment, esc.end) %>%
  summarize(N = n())

# Write it out to file
write.csv(fig.3N, file = './results/fig3_N-numbers.csv', row.names = F)

################################################################################
# Figure 3c
# Total number of ICM cells in series of chimeras (host-derived + ESCs)
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts %>%
  filter(!Litter %in% small.odd, 
         Genotype1 == 'wt', 
         ES_culture == 'S/LIF', 
         Stage.t0 == '8cell', 
         ESC_line %in% c('CAG:H2B-GFP', 
                         'H2B-tdTomato', 
                         'no.esc'), 
         Exp_date > 20170101) %>% 
  group_by(Embryo_ID, Treatment, 
           esc.end, icm.count, 
           ESC_line) %>%
  summarize()

# Generate plot
fig.3c <- ggplot(data = my.data, 
                 aes(x = esc.end, y = icm.count))
fig.3c <- fig.3c + geom_boxplot(color = 'black',
                                outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F) + 
  geom_jitter(color = 'black', shape = 20, width = 0.2, 
                               show.legend = F)
fig.3c <- fig.3c + looks + 
  theme(aspect.ratio = 9/6) + 
  ylim(0, 100) + 
  scale_color_manual(values = escols) + 
  labs(x = '#ESCs at 48h (as xEPI)', 
       y = 'Total ICM cells at 48h', 
       title = 'Figure 3c') 
# Uncomment print() below to visualize plot
print(fig.3c)

################################################################################
# Figure 3d
# Number of host-derived ICM cells only 
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts %>%
  filter(!Litter %in% small.odd, 
         Genotype1 == 'wt', 
         ES_culture == 'S/LIF',
         Stage.t0 == '8cell', 
         ESC_line %in% c('CAG:H2B-GFP', 
                         'H2B-tdTomato', 
                         'no.esc')) %>% 
  group_by(Embryo_ID, Treatment, 
           esc.end, host.icm, 
           ES_culture, ESC_line) %>%
  summarize()
  
# Generate plot
fig.3d <- ggplot(data = my.data, 
                 aes(x = esc.end, y = host.icm))
fig.3d <- fig.3d + geom_boxplot(color = 'black',
                                outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F) + 
  geom_jitter(color = 'black', shape = 20, width = 0.2, 
              show.legend = F)
fig.3d <- fig.3d + looks + 
  theme(aspect.ratio = 9/6) + 
  ylim(0, 100) + 
  #facet_wrap( ~ ES_culture, nrow = 1) + 
  scale_color_manual(values = escols) + 
  labs(x = '#ESCs at 48h (as xEPI)', 
       y = '#host-derived ICM cells at 48h', 
       title = 'Figure 3d') 
# Uncomment print() below to visualize plot
# print(fig.3d)

################################################################################
# Figure 3e
# ICM composition per embryo. Each bar is the ICM of one embryo, 
# arrangedby increasing absolute number of ESCs
################################################################################

# Extract the absolute number of ESCs per embryo from counts table
esc.counts <- select(esc.xim.lincounts2, Embryo_ID, ESC, esc.end)

# Define data to plot
my.data <- merge(esc.chimeras, esc.counts) %>% 
  filter(TE_ICM == 'ICM', 
         Genotype1 == 'wt', 
         Cellcount > 80, 
         Stage.t0 == '8cell', 
         ES_culture == 'S/LIF', 
         !Litter %in% small.odd)

# Generate plot
fig.3e <- ggplot(data = my.data, 
                  aes(x = reorder(interaction(Embryo_ID, esc.end), ESC), 
                      fill = Identity.hc))
fig.3e <- fig.3e + geom_bar(position = 'fill') + 
  geom_hline(yintercept = 0.4, linetype = 'dashed')
fig.3e <- fig.3e + looks + 
  scale_fill_manual(values = idcols) + 
  #facet_wrap( ~ Treatment + esc.end, scales = 'free', nrow = 1) + 
  labs(y = '% of ICM', title = 'Figure 3e', 
       x = 'Embryo (by final number of ESCs, after 48h)') + 
  theme(aspect.ratio = 0.33, 
        axis.text.x = element_text(angle = 30, size = 6, hjust = 1)) + 
  annotate('text', y = 0.05, 
           x = length(unique(my.data$Embryo_ID)) * 0.975, 
           label = paste('N =', print(length(unique(my.data$Embryo_ID))), 
                         sep = ' '))
# Uncomment print() below to visualize plot
# print(fig.3e)

################################################################################
# Figure 3f
# Total number of PrE cells in each group of the chimeric series
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts2 %>%
  filter(!Litter %in% small.odd, 
         Genotype1 == 'wt', 
         Cellcount > 80, 
         ES_culture == 'S/LIF',
         Stage.t0 == '8cell')

# Generate plot
fig.3f <- ggplot(data = my.data, 
                 aes(x = esc.end, y = PRE))
fig.3f <- fig.3f + geom_boxplot(color = 'black',
                                outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F) + 
  geom_jitter(color = idcols['PRE'], 
              shape = 20, width = 0.2, 
              show.legend = F)
fig.3f <- fig.3f + looks + 
  theme(aspect.ratio = 9/6) + 
  # facet_wrap( ~ ES_culture, nrow = 1) + 
  scale_color_manual(values = escols) + 
  labs(x = '#ESCs at 48h (as xEPI)', 
       y = '#PrE cells at 48h', 
       title = 'Figure 3f')
# Uncomment print() below to visualize plot
# print(fig.3f)

################################################################################
# Generate PDF
################################################################################

pdf(file = './figures/fig3all_NS.pdf', width = 11, paper = 'a4r')
print(fig.3c)
print(fig.3d)
print(fig.3e)
print(fig.3f)
dev.off()

##