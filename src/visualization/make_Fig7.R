# This script generates the plots contained in Figure 7 of the paper

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

# Order genotyhpes for plotting
esc.xim.lincounts2$Genotype1 <- factor(esc.xim.lincounts2$Genotype1, 
                                       levels = c('wt', 'ko'))

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


# Generate table with embryos included in Fig 6 and their cellcounts
fig.6_nums <- merge(esc.chimeras, esc.end) %>% 
  filter(TE_ICM == 'ICM', 
         Genotype1 == 'ko', 
         Stage.t0 == '8cell', 
         !Litter %in% small.odd, 
         ES_culture == 'S/LIF', 
         ESC_genotype %in% c('no.esc', 'wt', 'ko'), 
         ESC_line != "F103_CAG:H2B-GFP-B2")  %>%
  group_by(Embryo_ID, Litter, Treatment, ESC_line, ESC_genotype, 
           ES_culture, esc.end, Cellcount, icm.count) %>%
  summarize()

# Generate table with N numbers for embryos in Fig 6
fig.6N <- fig.6_nums %>% 
  group_by(Treatment, ESC_line, ES_culture, 
           ESC_genotype, esc.end) %>%
  summarize(N = n())

# Write it out to file
write.csv(fig.6N, file = './results/fig7_N-numbers.csv', row.names = F)

################################################################################
# Figure 7c
# Box plots showing the number of PrE cells for each group in the chimeric
# series compared to the number of PrE cells in wild type controls
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts2 %>%
  filter(!Litter %in% small.odd, 
         (Genotype1 == 'ko' & 
            Treatment %in% c('Control', 
                             'Chimera') | 
            Genotype1 == 'wt' & 
            Treatment == 'Control' & 
            esc.end == '0'), 
         ES_culture == 'S/LIF', 
         ESC_line %in% c('CAG:H2B-GFP', 'H2B-tdTomato'), 
         Stage.t0 == '8cell', 
         Exp_date > 20170101)

# Generate plot
fig.6c <- ggplot(data = my.data, 
                 aes(x = interaction(Genotype1, esc.end), 
                     y = PRE))
fig.6c <- fig.6c + geom_boxplot(color = 'black',
                                outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
fig.6c <- fig.6c + geom_jitter(color = idcols['PRE'], 
                               shape = 20, width = 0.2, 
                               show.legend = F)
fig.6c <- fig.6c + looks + facet_wrap( ~ ES_culture, nrow = 1)
fig.6c <- fig.6c + theme(aspect.ratio = 9/7, 
                         axis.text.x = element_text(angle = 30, hjust = 1))
fig.6c <- fig.6c + scale_color_manual(values = escols)
fig.6c <- fig.6c + labs(x = 'Final #ESCs (as xEPI)', 
                        y = 'Final number of PrE cells', 
                        title = 'Figure 7c')
# Uncomment print() below to visualize plot
# print(fig.6c)

################################################################################
# Figure 7d
# ICM composition per embryo. Each bar is the ICM of one embryo, 
# arrangedby increasing absolute number of ESCs
################################################################################

# Extract absolute number of ESCs per embryo
esc.counts <- select(esc.xim.lincounts2, Embryo_ID, ESC, esc.end)

# Define data to plot
my.data <- merge(esc.chimeras, esc.counts) %>% 
  filter(TE_ICM == 'ICM', 
         Genotype1 == 'ko', 
         Cellcount > 80, 
         Treatment %in% c('Control', 'Chimera'), 
         ESC_genotype %in% c('wt', 'no.esc'), 
         ESC_line %in% c('CAG:H2B-GFP', 
                         'H2B-tdTomato'), 
         Stage.t0 == '8cell', 
         ES_culture %in% c('S/LIF', 'no.esc'), 
         !Litter %in% small.odd, 
         Exp_date > 20170101)

# Generate plot
fig.6d <- ggplot(data = my.data, 
                 aes(x = reorder(Embryo_ID, ESC), fill = Identity.hc))
fig.6d <- fig.6d + geom_bar(position = 'fill')
fig.6d <- fig.6d + geom_hline(yintercept = 0.4, linetype = 'dashed')
fig.6d <- fig.6d + looks + scale_fill_manual(values = idcols)
# fig.6d <- fig.6d + facet_wrap( ~ ESC_line, scales = 'free', nrow = 1)
fig.6d <- fig.6d + labs(title = 'Figure 7d', y = '% of ICM', 
                        x = 'Embryo (by final number of ESCs at 48h)')
fig.6d <- fig.6d + theme(aspect.ratio = 0.5, 
                         axis.text.x = element_text(angle = 30, hjust = 1))
fig.6d <- fig.6d + annotate('text', y = 0.05, 
                            x = length(unique(my.data$Embryo_ID)) * 0.975, 
                            label = paste('N =', 
                                          print(length(
                                            unique(my.data$Embryo_ID))), 
                                          sep = ' '))
# Uncomment print() below to visualize plot
# print(fig.6d)

################################################################################
# Generate PDF
################################################################################

pdf(file = './figures/fig7all_NS.pdf', width = 11, paper = 'a4r')
print(fig.6c)
print(fig.6d)
dev.off()



