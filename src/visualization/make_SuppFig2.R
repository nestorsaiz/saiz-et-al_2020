# This script generates the plots contained in Supplementary Figure 2

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Read in embryo-embryo chimeras data
if (exists('g6.chimeras') ==  F) { 
  g6.chimeras <- read.csv('./data/processed/emb-xim-processed.csv')
  g6.ref <- read.csv('./references/emb-xim_exp_ref.csv')
  g6.lincounts <- read.csv('./data/processed/emb-xim-counts.csv')
  g6.lintype <- read.csv('./data/processed/emb-xim-typecounts.csv')
}

# Run emb-xim_counter.R to generate counts tables
source('./src/data/emb-xim_counter.R')

# If tables aren't present in disk, generate them from scratch
if (exists('g6.chimeras') == F) { 
  source('./src/emb-xim_runall.R')
}

# Create vector to define cell type color (GFP+ vs GFP- (wt or Gata6-/-))
cellcols <- c('host' = '#cccccc', 'donor' = '#33cc33', 
              'host.wt' = '#cccccc', 'host.ko' = '#ff3399', 
              'host.CAG:H2B-GFP' = '#33cc33', 
              'donor.wt' = '#33cc33', 'donor.ko' = '#33cc33')

################################################################################
# Make function to generate plots in panels b, c, h and i, which are similar,
# by just changing the data and the filters used in each
################################################################################

g6.bp <- function(dataset, my.filters, y.var, y.label) { 
  plot.0 <- ggplot(data = dataset %>% filter(H_genotype %in% 
                                               my.filters, 
                                             Embryo_size == 'Single'), 
                   aes(x = as.factor(t0.donor), 
                       y = !! as.name(y.var)))
  plot.0 <- plot.0 + 
    geom_boxplot(color = 'black',
                 outlier.shape = 1, outlier.size = 2) +
    stat_summary(fun.y = mean, colour = "black", geom = "point", 
                 shape = 4, size = 3, show.legend = F)
  plot.0 <- plot.0 + 
    geom_jitter(data = dataset %>% filter(H_genotype %in% 
                                            my.filters, 
                                          Embryo_size == 'Single'), 
                aes(color = interaction(cell_type, H_genotype)), 
                shape = 20, width = 0.2)
  plot.0 <- plot.0 + 
    geom_jitter(data = dataset %>% filter(H_genotype %in% 
                                            my.filters, 
                                          Embryo_size == 'Double'), 
                color = '#007300', shape = 20, width = 0.2)
  plot.0 <- plot.0 + scale_color_manual(values = cellcols) + looks
  plot.0 <- plot.0 + facet_wrap( ~ cell_type)
  plot.0 <- plot.0 + theme(aspect.ratio = 1)
  plot.0 <- plot.0 + labs(y = y.label, x = 'Initial % H2B-GFP cells', 
                          title = paste('Genotypes: ', 
                                        my.filters, sep = ''))
  return(plot.0)
}

################################################################################
# Supplementary Figure 2b
# Box plots showing contribution of wt GFP- and wt GFP+ cells to the chimera
################################################################################

# Use function to generate plot
figS2b <- g6.bp(dataset = g6.typecount, my.filters = c('wt', 'CAG:H2B-GFP'), 
      y.var = 'pc.type', y.label = '% of total')
# Uncomment below to visualize plot
# print(figS2b)

################################################################################
# Supplementary Figure 2c
# Box plots showing contribution of wt GFP- and wt GFP+ cells to the ICM only
################################################################################

# Use function to generate plot
figS2c <- g6.bp(dataset = g6.icmtype, my.filters = c('wt', 'CAG:H2B-GFP'), 
                y.var = 'pc.icm', y.label = '% of ICM')
# Uncomment below to visualize plot
# print(figS2c)

################################################################################
# Same as Figure 2c, but showing contribution to TE - not shown in paper
################################################################################

# Use function to generate plot
# Can be used to plot number of any other lineage by simply 
# changing the subsetting condition (Identity.hc)
figS2x <- g6.bp(dataset = subset(g6.lintype, 
                                    Identity.hc == 'TE'), 
                   my.filters = c('wt', 'CAG:H2B-GFP'), 
                y.var = 'pc.lin', y.label = '% of TE')
# Uncomment below to visualize plot
# print(figS2x)

################################################################################
# Supplementary Figure 2e
# Box plots showing contribution of wt GFP- and wt GFP+ cells 
# to PrE and EPI of the resulting chimera
################################################################################

# Generate plot. Similar as above, but x variable and facetting change
figS2e <- ggplot(data = g6.lintype %>% 
                   filter(H_genotype %in% c('wt'), 
                          Identity.hc %in% c('PRE', 'EPI')), 
                 aes(x = donor.end, 
                     y = pc.lin)) 
figS2e <- figS2e + geom_boxplot(color = 'black',
                                outlier.shape = 1, 
                                outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
figS2e <- figS2e + geom_jitter(data = g6.lintype %>% 
                                 filter(H_genotype %in% c('wt'), 
                                        Identity.hc %in% c('PRE', 'EPI'), 
                                        Embryo_size == 'Single'), 
                               aes(color = interaction(cell_type, H_genotype)), 
                               shape = 20, width = 0.2)
figS2e <- figS2e + geom_jitter(data = g6.lintype %>% 
                                 filter(H_genotype %in% c('wt'), 
                                        Identity.hc %in% c('PRE', 'EPI'), 
                                        Embryo_size == 'Double'), 
                               color = '#007300', 
                               shape = 20, width = 0.2)
figS2e <- figS2e + facet_grid(Identity.hc ~ interaction(cell_type, H_genotype))
figS2e <- figS2e + scale_color_manual(values = cellcols)
figS2e <- figS2e + looks + theme(aspect.ratio = 1, 
                                 axis.text.x = element_text(angle = 45, 
                                                            hjust = 1))
figS2e <- figS2e + labs(x = 'Final % of H2B-GFP ICM cells', 
                        y = '% of lineage', 
                        color = 'Cell type, host genotype', 
                        title = 'Figure S2e')
# Uncomment below to visualize plot
# print(figS2e)

################################################################################
# Supplementary Figure 2f
# Stacked bar plots showing the global lineage composition of the chimeric ICM
# in chimeras containing wt GFP- and wt GFP+ cells
################################################################################

# Generate plot
figS2f <- ggplot(subset(g6.chimeras, TE_ICM == 'ICM' & 
                          H_genotype %in% c('wt')), 
                 aes(x = donor.end, fill = Identity.hc))
figS2f <- figS2f + geom_bar(position = 'fill')
figS2f <- figS2f + scale_fill_manual(values = idcols) 
figS2f <- figS2f + facet_grid(H_genotype ~ .)
figS2f <- figS2f + looks + theme(aspect.ratio = 1, 
                                 axis.text.x = element_text(angle = 45, 
                                                            hjust = 1))
figS2f <- figS2f + labs(x = 'Final % of H2B-GFP ICM cells', 
                        y = '% ICM', 
                        fill = 'Identity', 
                        title = 'Figure S2f')
# Uncomment below to visualize plot
# print(figS2f)

################################################################################
# Supplementary Figure 2h
# Box plots showing contribution of Gata6-/- GFP- and wt GFP+ cells to chimera
################################################################################

# Use function to generate plot
figS2h <- g6.bp(dataset = g6.typecount, my.filters = c('ko'), 
                y.var = 'pc.type', y.label = '% of total')
# Uncomment below to visualize plot
# print(figS2h)

################################################################################
# Supplementary Figure 2i
# Box plots showing contribution of Gata6-/- GFP- and wt GFP+ cells to the ICM 
################################################################################

# Use function to generate plot
figS2i <- g6.bp(dataset = g6.icmtype, my.filters = c('ko'), 
                y.var = 'pc.icm', y.label = '% of ICM')
# Uncomment below to visualize plot
# print(figS2i)

################################################################################
# Same as Figure 2i, but showing contribution to TE - not shown in paper
################################################################################

# Use function to generate plot
# Can be used to plot number of any other lineage by simply 
# changing the subsetting condition (Identity.hc)
figS2xx <- g6.bp(dataset = subset(g6.lintype, 
                                    Identity.hc == 'TE'), 
                   my.filters = c('ko'), 
                   y.var = 'pc.lin', y.label = '% of TE')
# Uncomment below to visualize plot
# print(figS2xx)

################################################################################
# Supplementary Figure 2k
# Box plots showing contribution of Gata6-/- GFP- and wt GFP+ cells 
# to PrE and EPI of the resulting chimera
################################################################################

# Generate plot. Similar as above, but x variable and facetting change
figS2k <- ggplot(data = g6.lintype %>% 
                   filter(H_genotype %in% c('ko'), 
                          Identity.hc %in% c('PRE', 'EPI')), 
                 aes(x = donor.end, 
                     y = pc.lin)) 
figS2k <- figS2k + geom_boxplot(color = 'black',
                                outlier.shape = 1, 
                                outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
figS2k <- figS2k + geom_jitter(data = g6.lintype %>% 
                                 filter(H_genotype %in% c('ko'), 
                                        Identity.hc %in% c('PRE', 'EPI'), 
                                        Embryo_size == 'Single'), 
                               aes(color = interaction(cell_type, H_genotype)), 
                               shape = 20, width = 0.2)
figS2k <- figS2k + geom_jitter(data = g6.lintype %>% 
                                 filter(H_genotype %in% c('ko'), 
                                        Identity.hc %in% c('PRE', 'EPI'), 
                                        Embryo_size == 'Double'), 
                               color = '#007300', 
                               shape = 20, width = 0.2)
figS2k <- figS2k + facet_grid(Identity.hc ~ interaction(cell_type, H_genotype))
figS2k <- figS2k + scale_color_manual(values = cellcols)
figS2k <- figS2k + looks + theme(aspect.ratio = 1, 
                                 axis.text.x = element_text(angle = 45, 
                                                            hjust = 1))
figS2k <- figS2k + labs(x = 'Final % of H2B-GFP ICM cells', 
                        y = '% of lineage', 
                        color = 'Cell type, host genotype', 
                        title = 'Figure S2k')
# Uncomment below to visualize plot
# print(figS2k)

################################################################################
# Supplementary Figure 2l
# Stacked bar plots showing the global lineage composition of the chimeric ICM
# in chimeras containing Gata6-/- GFP- and wt GFP+ cells
################################################################################

# Generate plot
figS2l <- ggplot(subset(g6.chimeras, TE_ICM == 'ICM' & 
                          H_genotype %in% c('ko')), 
                 aes(x = donor.end, fill = Identity.hc))
figS2l <- figS2l + geom_bar(position = 'fill')
figS2l <- figS2l + scale_fill_manual(values = idcols) 
figS2l <- figS2l + facet_grid(H_genotype ~ .)
figS2l <- figS2l + looks + theme(aspect.ratio = 1, 
                                 axis.text.x = element_text(angle = 45, 
                                                            hjust = 1))
figS2l <- figS2l + labs(x = 'Final % of H2B-GFP ICM cells', 
                        y = '% ICM', 
                        fill = 'Identity', 
                        title = 'Figure S2l')
# Uncomment below to visualize plot
# print(figS2l)

################################################################################
# Calculate the % of each subpopulation (GFP-, GFP+) that adopted EPI fate
# in each embryo
################################################################################

# Extract same info as for g6.lincounts, but without total icm counts
g6.lincounts2 <- g6.chimeras %>% 
  filter(Channel == 'CH2', TE_ICM == 'ICM') %>% 
  group_by(Experiment, Embryo_ID, Treatment,
           Litter, TE_ICM, Cellcount, 
           H_genotype, D_genotype, 
           Embryo_size, t0.donor, donor.end, 
           cell_type, Identity.hc) %>%
  summarise(lin.count = n())
# Account for zeroes
m.g6 <- dcast(g6.lincounts2, Experiment + Litter + Embryo_ID + Treatment + 
                Cellcount + TE_ICM + H_genotype + D_genotype + 
                Embryo_size + t0.donor + cell_type + 
                donor.end ~ Identity.hc, 
              value.var = 'lin.count')
# Replace NAs in each ICM lineage denomination with zeros
m.g6[is.na(m.g6)] <- 0
m.g6$all.EPI <- m.g6$EPI + m.g6$EPI.lo
# Melt ICM back to long format
g6.lincounts2 <- melt(m.g6, id.vars = c('Experiment', 'Litter', 'Embryo_ID', 
                                        'Treatment', 'Cellcount', 'TE_ICM', 
                                        'H_genotype', 'D_genotype', 'Embryo_size', 
                                        't0.donor', 'donor.end', 'cell_type'), 
                      variable.name = 'Identity.hc', value.name = 'lin.count')
# remove wide format table
rm(m.g6)

# Extract ICM counts for each cell type and incorporate to above table
gg <- g6.icmtype %>% 
  filter(Embryo_ID %in% unique(g6.lincounts2$Embryo_ID)) %>% 
  group_by(Embryo_ID, cell_type, icm.type) %>% 
  summarize()
g6.lincounts2 <- merge(g6.lincounts2, gg)
rm(gg)
# Calculate the % of each cell type in each subpopulation
g6.lincounts2$pc.icmtype <- g6.lincounts2$lin.count / 
  g6.lincounts2$icm.type * 100

################################################################################
# Supplementary Figure 2m
# Box plots showing the % of each compartment (GFP+ or GFP-) that 
# adopted EPI fate (the probability of those cells becoming EPI)
# in chimeras where both compartments are wild type
################################################################################

figS2m <- ggplot(data = g6.lincounts2 %>% 
                   filter(Identity.hc == 'all.EPI',
                          H_genotype %in% c('wt')), 
                 aes(x = donor.end, y = pc.icmtype))
figS2m <- figS2m + geom_boxplot(color = 'black', outlier.shape = 1, 
                                outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
figS2m <- figS2m + geom_jitter(aes(color = Identity.hc), 
                               shape = 20, width = 0.2)
figS2m <- figS2m + scale_color_manual(values = idcols)
figS2m <- figS2m + facet_grid(H_genotype ~ cell_type)
figS2m <- figS2m + looks + theme(aspect.ratio = 1, 
                                 axis.text.x = element_text(angle = 45, 
                                                            hjust = 1))
figS2m <- figS2m + labs(x = 'Final % of H2B-GFP ICM cells', 
                        y = '% of cells with EPI identity', 
                        color = 'Identity', 
                        title = 'Figure S2m')
# Uncomment to visualize plot
# print(figS2m)

################################################################################
# Supplementary Figure 2m
# Box plots showing the % of each compartment (GFP+ or GFP-) that 
# adopted EPI fate (the probability of those cells becoming EPI)
# in chimeras in which GFP- cells are Gata6-/-
################################################################################

figS2n <- ggplot(data = g6.lincounts2 %>% 
                   filter(Identity.hc == 'all.EPI',
                          H_genotype %in% c('ko')), 
                 aes(x = donor.end, y = pc.icmtype))
figS2n <- figS2n + geom_boxplot(color = 'black', outlier.shape = 1, 
                                outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
figS2n <- figS2n + geom_jitter(aes(color = Identity.hc), 
                               shape = 20, width = 0.2)
figS2n <- figS2n + scale_color_manual(values = idcols)
figS2n <- figS2n + facet_grid(H_genotype ~ cell_type)
figS2n <- figS2n + looks + theme(aspect.ratio = 1, 
                                 axis.text.x = element_text(angle = 45, 
                                                            hjust = 1))
figS2n <- figS2n + labs(x = 'Final % of H2B-GFP ICM cells', 
                        y = '% of cells with EPI identity', 
                        color = 'Identity', 
                        title = 'Figure S2n')
# Uncomment to visualize plot
# print(figS2n)

# Generate PDF
pdf(file = './figures/figS2all_NS.pdf', width = 11, paper = 'a4r')
print(figS2b)
print(figS2c)
print(figS2e)
print(figS2f)
print(figS2h)
print(figS2i)
print(figS2k)
print(figS2l)
print(figS2m)
print(figS2n)
dev.off()


