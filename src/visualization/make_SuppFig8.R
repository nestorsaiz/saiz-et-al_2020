# This script generates the plots contained in Supplementary Figure 8

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
esc.chimeras$Identity.km <- factor(esc.chimeras$Identity.km, 
                                   levels = c('TE', 'PRE', 'DP', 
                                              'EPI', 'EPI.lo', 'DN', 'ESC'))
esc.xim.lincounts$Identity.hc <- factor(esc.xim.lincounts$Identity.hc, 
                                        levels = c('TE', 'PRE', 'DP', 
                                                   'EPI', 'EPI.lo', 'DN', 'ESC'))

# Order factors for plotting
esc.xim.lincounts2$Genotype1 <- factor(esc.xim.lincounts2$Genotype1, 
                                       levels = c('wt', 'fl/fl', 'fl/ko', 
                                                  'ko', 'tbd', 'mosaic', 
                                                  'unknown'))
esc.xim.lincounts2$ESC_genotype <- factor(esc.xim.lincounts2$ESC_genotype, 
                                          c('no.esc', 'wt', 'het', 'ko'))
esc.xim.lincounts2$ESC_line <- factor(esc.xim.lincounts2$ESC_line, 
                                      levels = c('no.esc', 'CAG:H2B-GFP', 
                                                 'H2B-tdTomato', 
                                                 'F103_CAG:H2B-GFP-B2', 
                                                 'F105_GFP'))
esc.chimeras$Treatment <- factor(esc.chimeras$Treatment, 
                                       levels = c('Littermate', 'Control', 
                                                  'Chimera'))
esc.xim.lincounts2$Treatment <- factor(esc.xim.lincounts2$Treatment, 
                                       levels = c('Littermate', 'Control', 
                                                  'Chimera'))

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
write.csv(fig.6N, file = './results/fig6_N-numbers.csv', row.names = F)

################################################################################
# Figure S8a
# Box plots showing the number of PrE cells in wild type controls, 
# Fgf4-/- controls, chimeras with wt ESCs and Fgf4-/- hosts and
# chimeras with Fgf4-/- ESCs and Fgf4-/- hosts
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts2 %>%
  filter(!Litter %in% small.odd, 
         (Genotype1 == 'ko' & 
            Treatment %in% c('Control', 
                             'Chimera') & 
            ESC_genotype %in% c('wt', 'ko', 
                                'no.esc') & 
            esc.end %in% c('0', '0.5x', '1x', '2x') | 
            Genotype1 == 'wt' & 
            Treatment == 'Control' & 
            esc.end == '0'), 
         ES_culture == 'S/LIF', 
         ESC_line %in% c('CAG:H2B-GFP', 'H2B-tdTomato', 
                         'F105_GFP', 'no.esc'), 
         Stage.t0 == '8cell', 
         Exp_date > 20170101)

# Generate plot
fig.s8a <- ggplot(data = my.data, 
                  aes(x = interaction(Treatment, ESC_genotype, Genotype1), 
                      y = PRE))
fig.s8a <- fig.s8a + geom_boxplot(color = 'black',
                                  outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
fig.s8a <- fig.s8a + geom_jitter(aes(color = ESC_genotype), 
                                 shape = 20, width = 0.2)
fig.s8a <- fig.s8a + looks + facet_wrap( ~ ES_culture, nrow = 1)
fig.s8a <- fig.s8a + theme(aspect.ratio = 9/4, 
                           axis.text.x = element_text(angle = 30, hjust = 1))
fig.s8a <- fig.s8a + scale_color_manual(values = c('wt' = '#ccff33', 
                                                   'ko' = '#007300', 
                                                   'no.esc' = 'black'))
fig.s8a <- fig.s8a + labs(x = 'Host genotype, ESC genotype', 
                          y = 'Final number of PrE cells', 
                          title = 'Figure S8a')
# Uncomment print() below to visualize plot
# print(fig.s8a)

################################################################################
# Figure S8b
# Boxplot showing the number of ESCs after 48h in each chimera group
# as defined by number of cells at t=0 (D_cells) and ESC genotype in S/LIF
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts2 %>%
  filter(!Litter %in% small.odd, 
         Genotype1 == 'ko', 
         ES_culture == 'S/LIF', 
         ESC_line %in%  c("CAG:H2B-GFP", 
                          "H2B-tdTomato"),
         Stage.t0 == '8cell', 
         D_cells != '5',
         Exp_date > 20170101)

# Generate plot
fig.s8b <- ggplot(data = my.data, 
                  aes(x = D_cells, y = ESC))
fig.s8b <- fig.s8b + geom_boxplot(color = 'black',
                                  outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
fig.s8b <- fig.s8b + geom_jitter(color = 'black', 
                                 shape = 20, width = 0.2, 
                                 show.legend = F)
fig.s8b <- fig.s8b + looks + facet_wrap( ~ ES_culture, 
                                         nrow = 1, scales = 'free_x')
fig.s8b <- fig.s8b + theme(aspect.ratio = 9/5)
fig.s8b <- fig.s8b + scale_color_manual(values = escols)
fig.s8b <- fig.s8b + labs(x = 'Initial number of ESCs', 
                          y = 'Final number of ESCs', 
                          title = 'Figure S8b')
# Uncomment print() below to visualize plot
# print(fig.s8b)

################################################################################
# Figure S8c
# Box plots showing the number of host-derived ICM cells in each group
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts %>%
  filter(!Litter %in% small.odd, 
         TE_ICM == 'ICM', 
         Genotype1 == 'ko', 
         ES_culture == 'S/LIF', 
         Stage.t0 == '8cell', 
         ESC_line %in% c('CAG:H2B-GFP', 
                         'H2B-tdTomato'), 
         Exp_date > 20170101) %>% 
  group_by(Embryo_ID, Treatment, 
           esc.end, host.icm, 
           ESC_line) %>%
  summarize()

# Generate plot
fig.s8c <- ggplot(data = my.data, 
                  aes(x = esc.end, y = host.icm))
fig.s8c <- fig.s8c + geom_boxplot(color = 'black',
                                  outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
fig.s8c <- fig.s8c + geom_jitter(color = 'black', 
                                 shape = 20, width = 0.2, 
                                 show.legend = F)
fig.s8c <- fig.s8c + looks + theme(aspect.ratio = 9/6) + ylim(0, 100)
fig.s8c <- fig.s8c + scale_color_manual(values = escols)
fig.s8c <- fig.s8c + labs(x = 'Final#ESCs (as xEPI)', 
                          y = 'Host-derived ICM cells', 
                          title = 'Figure S8c') 
# Uncomment print() below to generate plot
# print(fig.s8c)

################################################################################
# Figure S8d
# Stacked bar plot with average relative ICM composition for each group
# Same data as in Figure 6d but binned by #ESCs
################################################################################

# Define data to plot
my.data <- merge(esc.chimeras, esc.end) %>% 
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
fig.s8d <- ggplot(data = my.data, 
                  aes(x = esc.end, fill = Identity.hc))
fig.s8d <- fig.s8d + geom_bar(position = 'fill') + 
  geom_hline(yintercept = 0.4, linetype = 'dashed')
fig.s8d <- fig.s8d + looks + scale_fill_manual(values = idcols)
fig.s8d <- fig.s8d + facet_wrap( ~ ES_culture, nrow = 1)
fig.s8d <- fig.s8d + labs(y = '% of ICM', x = 'Final #ESCs (as xEPI)', 
                          title = 'Figure S8d')
fig.s8d <- fig.s8d + theme(aspect.ratio = 9/6)
# Uncomment print() below to visualize plot
# print(fig.s8d)

################################################################################
# Figure S8e
# Stacked bar plots showing average ICM composition for each group of chimeras
# Same as S8d, but showing average absolute numbers
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts %>% 
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
         Exp_date > 20170101) %>% 
  group_by(ES_culture, #ESC_line, 
           esc.end, Identity.hc) %>% 
  summarize(meancount = mean(count), 
            N = n())

# Generate plot
fig.s8e <- ggplot(data = my.data, 
                  aes(x = esc.end, y = meancount, 
                      fill = Identity.hc))
fig.s8e <- fig.s8e + geom_bar(stat = 'identity')
fig.s8e <- fig.s8e + geom_text(aes(label = N), check_overlap = T,
                               vjust = 1, y = 7, size = 8)
fig.s8e <- fig.s8e + looks + scale_fill_manual(values = idcols)
fig.s8e <- fig.s8e + facet_wrap( ~ ES_culture, nrow = 1)
fig.s8e <- fig.s8e + labs(y = 'Number of cells', 
                          x = 'Final #ESCs (as xEPI)', 
                          title = 'Figure S8e')
fig.s8e <- fig.s8e + theme(aspect.ratio = 9/6) + ylim(0, 100)
# Uncomment print() below to visualize plot
# print(fig.s8e)

################################################################################
# Figure S8f
# Box plot showing the size of the ICM in wt vs Fgf4-/- control blastocysts
################################################################################

# Define data to plot
my.data <- esc.xim.lincounts2 %>%
  filter(!Litter %in% small.odd, 
         Genotype1 %in% c('wt', 'ko'), 
         Treatment == 'Control', 
         Stage.t0 == '8cell')

# Generate plot
fig.s8f <- ggplot(data = my.data, 
                  aes(x = Genotype1, y = icm.count))
fig.s8f <- fig.s8f + geom_boxplot(color = 'black',
                                  outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
fig.s8f <- fig.s8f + geom_jitter(color = 'black', 
                                 shape = 20, width = 0.2)
fig.s8f <- fig.s8f + looks + theme(aspect.ratio = 9/3) + ylim(0, 70)
fig.s8f <- fig.s8f + labs(x = 'Genotype', 
                          y = 'Number of ICM cells', 
                          title = 'Figure S8f')
# Uncomment print() below to visualize plot
# print(fig.s8f)

################################################################################
# Figure S8h
# Stacked bar plots with ICM compostion for chimeras, control embryos
# and reference littermates for chimeras where ESCs were injected 
# at the blastocyst stage
################################################################################

# Extract the absolute number of ESCs per embryo from counts table
esc.counts <- select(esc.xim.lincounts2, Embryo_ID, ESC, esc.end)

# Define data to plot
my.data <- merge(esc.chimeras, esc.counts) %>% 
  filter(TE_ICM == 'ICM', 
         Genotype1 == 'ko', 
         ESC_genotype %in% c('wt', 'no.esc'), 
         ESC_line %in% c('CAG:H2B-GFP', 
                         'H2B-tdTomato'), 
         Stage.t0 == 'blastocyst', 
         ES_culture %in% c('2i/LIF', 'no.esc'), 
         !Litter %in% small.odd, 
         Exp_date > 20170101)

# Generate plot
fig.s8h <- ggplot(data = my.data, 
                  aes(x = reorder(Embryo_ID, ESC), fill = Identity.hc))
fig.s8h <- fig.s8h + geom_bar(position = 'fill')
fig.s8h <- fig.s8h + geom_hline(yintercept = 0.4, linetype = 'dashed')
fig.s8h <- fig.s8h + looks + scale_fill_manual(values = idcols)
fig.s8h <- fig.s8h + facet_wrap( ~ ES_culture + Treatment, 
                                 scales = 'free', nrow = 1)
fig.s8h <- fig.s8h + labs(y = '% of ICM', title = 'Figure S8h', 
                          x = 'Embryo (by final number of ESCs, after 24h')
fig.s8h <- fig.s8h + theme(aspect.ratio = 1, 
                           axis.text.x = element_text(angle = 30, hjust = 1))
# Uncomment print() below to visualize plot
# print(fig.s8h)

################################################################################
# Generate PDF
################################################################################



