# This script generates the plots in Figure 1 of the paper

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

# If tables aren't present in disk, generate them from scratch
if (exists('g6.chimeras') == F) { 
  source('./src/emb-xim_runall.R')
}

################################################################################
# Figure 1c
# Stacked bar plots showing ICM composition in chimeras made with
# wt GFP+ and wt GFP- blastomeres
################################################################################

# Define levels of Identity.hc for plotting
g6.chimeras$Identity.hc <- factor(g6.chimeras$Identity.hc, 
                                  levels = c('TE', 'PRE', 'DP', 
                                             'EPI', 'EPI.lo', 'DN'))

# Generate plot (uncomment print() to visualize plot)
fig1c <- ggplot(subset(g6.chimeras, TE_ICM == 'ICM' & 
                          H_genotype %in% c('wt')), 
                 aes(x = donor.end, fill = Identity.hc))
fig1c <- fig1c + geom_bar(position = 'fill')
fig1c <- fig1c + scale_fill_manual(values = idcols) 
fig1c <- fig1c + facet_grid(H_genotype ~ cell_type)
fig1c <- fig1c + looks + theme(aspect.ratio = 1, 
                               axis.text.x = element_text(angle = 45, 
                                                          hjust = 1))
fig1c <- fig1c + labs(x = 'Final % of H2B-GFP ICM cells', 
                      y = '% population', 
                      fill = 'Identity')
# print(fig1c)

# Generate table with number of embryos in each group in panel C
fig1c.N <- g6.chimeras %>% 
  filter(TE_ICM == 'ICM', 
         H_genotype %in% c('wt')) %>%
  group_by(cell_type, Embryo_ID, H_genotype, donor.end) %>%
  summarize() %>% 
  group_by(cell_type, donor.end) %>%
  summarize(N = n())

################################################################################
# Figure 1f
# Stacked bar plots showing ICM composition in chimeras made with
# wt GFP+ and Gata6-/- GFP- blastomeres
################################################################################

# Generate plot (uncomment print() to visualize plot)
fig1f <- ggplot(subset(g6.chimeras, TE_ICM == 'ICM' & 
                         H_genotype %in% c('ko')), 
                aes(x = donor.end, fill = Identity.hc))
fig1f <- fig1f + geom_bar(position = 'fill')
fig1f <- fig1f + scale_fill_manual(values = idcols) 
fig1f <- fig1f + facet_grid(H_genotype ~ cell_type)
fig1f <- fig1f + looks + theme(aspect.ratio = 1, 
                               axis.text.x = element_text(angle = 45, 
                                                          hjust = 1))
fig1f <- fig1f + labs(x = 'Final % of H2B-GFP ICM cells', 
                      y = '% population', 
                      fill = 'Identity')
# print(fig1f)

# Generate table with number of embryos in each group in panel F
fig1f.N <- g6.chimeras %>% 
  filter(TE_ICM == 'ICM', 
         H_genotype %in% c('ko')) %>%
  group_by(cell_type, Embryo_ID, H_genotype, donor.end) %>%
  summarize() %>% 
  group_by(cell_type, donor.end) %>%
  summarize(N = n())

# Generate PDF with figures
pdf(file = './figures/fig1cf_NS_lin-distrib.pdf', width = 11, paper = 'a4r')
print(fig1c)
print(fig1f)
dev.off()

