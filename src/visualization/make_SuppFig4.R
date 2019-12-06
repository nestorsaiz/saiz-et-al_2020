# This script generates the plots contained in Supplementary Figure 4

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

################################################################################
# Figure S4a
# Boxplot showing the number of ESCs after 48h in each chimera group
# as defined by number of cells at t=0 (D_cells) and ESC genotype in S/LIF
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts2 %>%
  filter(!Litter %in% small.odd, 
         Genotype1 == 'wt', 
         ES_culture == 'S/LIF',
         Stage.t0 == '8cell')

# Generate plot
fig.s4a <- ggplot(data = my.data, 
                  aes(x = D_cells, y = ESC))
fig.s4a <- fig.s4a + geom_boxplot(color = 'black',
                                  outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
fig.s4a <- fig.s4a + geom_jitter(color = 'black', 
                                 shape = 20, width = 0.2, 
                                 show.legend = F)
fig.s4a <- fig.s4a + looks # + facet_wrap( ~ ESC_line, scales = 'free_x')
# fig.s4a <- fig.s4a + looks + facet_wrap( ~ ES_culture, scales = 'free_x')
fig.s4a <- fig.s4a + theme(aspect.ratio = 9/5)
fig.s4a <- fig.s4a + scale_color_manual(values = escols)
fig.s4a <- fig.s4a + labs(x = 'Initial number of ESCs', 
                          y = 'Final number of ESCs', 
                          title = 'Figure S4a')
# Uncomment print() below to visualize plot
# print(fig.s4a)

################################################################################
# Figure S4c
# Stacked bar plots showing relative ICM composition for each group of chimeras
# Same data as in Figure 3e but binned and averaged by group
################################################################################

# Define data to plot
my.data <- merge(esc.chimeras, esc.counts) %>% 
  filter(TE_ICM == 'ICM', 
         Genotype1 == 'wt', 
         Cellcount > 80, 
         Stage.t0 == '8cell', 
         ES_culture == 'S/LIF', 
         !Litter %in% small.odd)

# Generate plot
fig.s4c <- ggplot(data = my.data, 
                 aes(x = esc.end, fill = Identity.hc))
fig.s4c <- fig.s4c + geom_bar(position = 'fill')
fig.s4c <- fig.s4c + geom_hline(yintercept = 0.4, linetype = 'dashed')
fig.s4c <- fig.s4c + looks + scale_fill_manual(values = idcols)
fig.s4c <- fig.s4c + labs(y = '% of ICM', x = '#ESCs at 48h (as xEPI)', 
                        title = 'Figure S4c')
fig.s4c <- fig.s4c + theme(aspect.ratio = 9/6)
# Uncomment print() below to visualize plot
# print(fig.s4c)

################################################################################
# Figure S4d
# Stacked bar plots showing average ICM composition for each group of chimeras
# Same as S4c, but showing average absolute numbers
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts %>% 
  filter(TE_ICM == 'ICM', 
         Genotype1 == 'wt', 
         ES_culture == 'S/LIF', 
         Stage.t0 == '8cell', 
         !Litter %in% small.odd) %>% 
  group_by(Genotype1, #ESC_line, ES_culture,
           esc.end, Identity.hc) %>% 
  summarize(meancount = mean(count), 
            N = n())

# Generate plot
fig.s4d <- ggplot(data = my.data, 
                  aes(x = esc.end, y = meancount, 
                      fill = Identity.hc))
fig.s4d <- fig.s4d + geom_bar(stat = 'identity')
fig.s4d <- fig.s4d + geom_text(aes(label = N), check_overlap = T,
                               vjust = 1, y = 5, size = 8)
fig.s4d <- fig.s4d + looks + scale_fill_manual(values = idcols)
# fig.s4d <- fig.s4d + facet_wrap( ~ ESC_line)
fig.s4d <- fig.s4d + labs(y = 'Number of cells', 
                          x = '#ESCs at 48h (as xEPI)', 
                          title = 'Figure S4d')
fig.s4d <- fig.s4d + theme(aspect.ratio = 9/6)
# Uncomment print() below to visualize plot
# print(fig.s4d)

################################################################################
# Figure S4e
# Stacked bar plots showing average relative ICM composition for host-derived
# ICM cells only in each group of chimeras
################################################################################

# Define data to plot
my.data <- merge(esc.chimeras, esc.end) %>% 
  filter(TE_ICM == 'ICM', 
         Identity.hc != 'ESC', 
         Genotype1 == 'wt', 
         Stage.t0 == '8cell', 
         ES_culture == 'S/LIF', 
         !Litter %in% small.odd)

# Generate plot
fig.s4e <- ggplot(data = my.data, 
                  aes(x = esc.end, fill = Identity.hc))
fig.s4e <- fig.s4e + geom_bar(position = 'fill')
fig.s4e <- fig.s4e + looks + scale_fill_manual(values = idcols)
fig.s4e <- fig.s4e + labs(y = '% of ICM', x = '#ESCs at 48h (as xEPI)', 
                          title = 'Figure S4e')
fig.s4e <- fig.s4e + theme(aspect.ratio = 9/6)
# Uncomment print() below to visualize plot
# print(fig.s4e)

################################################################################
# Figure S4h
# Stacked bar plots with ICM compostion for chimeras, control embryos
# and reference littermates for chimeras where ESCs were injected 
# at the blastocyst stage
################################################################################

# Generate table with embryos included in Fig S4f-h and their cellcounts
fig.s4g_nums <- esc.chimeras %>% 
  filter(TE_ICM == 'ICM', 
         Genotype1 == 'wt', 
         Litter != 'CH', 
         Stage.t0 == 'blastocyst', 
         ES_culture == 'S/LIF') %>%
  group_by(Embryo_ID, Litter, Treatment, ESC_line, 
           group.median, D_cells, Cellcount, icm.count) %>%
  summarize()

# Write it out to disk
write.csv(fig.s4g_nums, file = './results/figS4h_N-numbers.csv', row.names = F)

# Define data to plot
esc.chimeras$Treatment <- factor(esc.chimeras$Treatment, 
                                 levels = c('Littermate', 'Control', 'Chimera'))
my.data <- esc.chimeras %>% 
  filter(TE_ICM == 'ICM', 
         Genotype1 == 'wt', 
         Stage.t0 == 'blastocyst', 
         ES_culture == 'S/LIF')

# Generate plot
fig.s4g <- ggplot(data = my.data, 
                  aes(x = Treatment, fill = Identity.hc))
fig.s4g <- fig.s4g + geom_bar(position = 'fill')
fig.s4g <- fig.s4g + looks + scale_fill_manual(values = idcols)
fig.s4g <- fig.s4g + facet_wrap( ~ ESC_line) 
fig.s4g <- fig.s4g + labs(y = '% of ICM', title = 'Figure S4h')
fig.s4g <- fig.s4g + theme(aspect.ratio = 9/4, 
                           axis.text.x = element_text(angle = 30, hjust = 1))
# Uncomment print() below to visualize plot
# print(fig.s4g)

################################################################################
# Generate PDF
################################################################################

pdf(file = './figures/figS4all_NS.pdf', width = 11, paper = 'a4r')
print(fig.s4a)
print(fig.s4c)
print(fig.s4d)
print(fig.s4e)
print(fig.s4h)
dev.off()

##