# This script generates the plots in Supplementary Figure 1

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

# Load extra packages
library('plotly')

################################################################################
# Supplementary Figure 1a, b
# Stacked bar plots with the % of ICM that each lineage represents over time
################################################################################

# Define data to be plotted
my.data <- list(rbind.fill(subset(new.lms, TE_ICM == 'ICM' & 
                                    Genotype1 == 'wt' & 
                                    Genotype2 == 'wt' & 
                                    Cellcount > 30 & 
                                    Cellcount < 161), 
                           subset(ablat.lms, TE_ICM == 'ICM' & 
                                    Genotype1 == 'wt' & 
                                    Cellcount > 30 & 
                                    Cellcount < 161) #, 
                           # subset(compos, TE_ICM == 'ICM' & 
                           #          Genotype1 == 'wt' & 
                           #          Genotype2 == 'wt' & 
                           #          Cellcount > 30 & 
                           #          Cellcount < 161)
                           ), 
                subset(spry.lms, TE_ICM == 'ICM' & 
                         Genotype1 %in% c('wt', 'het') & 
                         Genotype2 == 'wt' & 
                         Cellcount > 30 & 
                         Cellcount < 161),
                subset(ncoms.lms, TE_ICM == 'ICM' & 
                         Cellcount > 30 & 
                         Cellcount < 161))

# Define Y variable for each plot
my.y <- c('Identity.hc', 'Identity.hc', 'Identity.km')

# List titles for each plot
my.titles <- c('This study', 'Morgani et al. (2018)', 'Saiz et al. (2016)')

# Generate plot
pdf(file = './figures/figS1ab_NS_ICM-comp.pdf', width = 11, paper = 'a4r')
for(d in c(1:3)) { 
  pc.icm <- ggplot(data = my.data[[d]], 
                   aes(x = reorder(interaction(Embryo_ID, Cellcount), 
                                   Cellcount), 
                       fill = !! as.name(my.y[d])))
  pc.icm <- pc.icm + geom_bar(position = 'fill')
  pc.icm <- pc.icm + looks + scale_fill_manual(values = idcols)
  pc.icm <- pc.icm + annotate('text', y = 0.10, 
                              x = length(unique(my.data[[d]]$Embryo_ID)) * 0.975, 
                              label = print(length(unique(my.data[[d]]$Embryo_ID))))
  pc.icm <- pc.icm + theme(axis.text.x = element_text(angle = 45, 
                                                      hjust = 1, 
                                                      size = 6), 
                           aspect.ratio = 0.3)
  pc.icm <- pc.icm + labs(title = my.titles[d], 
                          x = 'Embryo (by total cell number)', 
                          y = '% of ICM', fill = 'Identity')
  print(pc.icm)
}
dev.off()

################################################################################
# Supplementary Figure 1d
# Box plots comparing SOX17 and GATA4 levels between ICM lineages
################################################################################

# Load IF chart for 'new littermates' and 'ablation'
unicos.nl <- read.csv('./references/new-lms_unicos.csv')
unicos.ab <- read.csv('./references/ablat_unicos.csv')
unicos <- rbind.fill(unicos.nl, unicos.ab)

# Make a set of colors for the PrE markers
precols <- c('GATA6.gt' = '#eff3ff', 'GATA6.rb' = '#bdd7e7', 
             'PDGFRa' = '#6baed6', 'SOX17' = '#3182bd', 'GATA4' = '#08519c')

# Extract the experiments stained for each marker combination
s17 <- unicos$Experiment[which(unicos$CH5 == 'SOX17.gt' & 
                                 unicos$CH3 == 'GATA6.rb' & 
                                 unicos$CH2 == 'NANOG.rat')]
g4 <- unicos$Experiment[which(unicos$CH3 == 'GATA4.rb' & 
                                unicos$CH5 == 'GATA6.gt' & 
                                unicos$CH2 == 'NANOG.rat')]

# Generate plots for each marker (uncomment print() line to view plot)
## SOX17
s17.box <- ggplot(data = rbind.fill(subset(new.lms, TE_ICM == 'ICM' & 
                                             Experiment %in% s17), 
                                    subset(ablat.lms, TE_ICM == 'ICM' & 
                                             Experiment %in% s17)), 
                  aes(x = Identity.hc, y = CH5.ebLogCor.s))
s17.box <- s17.box + geom_boxplot(fill = precols[4], color = 'black', 
                                  outlier.shape = 1, outlier.stroke = 0.5) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
# s17.box <- s17.box + geom_jitter(shape = 20, width = 0.2, 
#                                  show.legend = F)
s17.box <- s17.box + looks + theme(aspect.ratio = 1)
s17.box <- s17.box + labs(x = 'Identity', y = 'log[SOX17]', 
                          title = 'Fig S1d')
# print(s17.box)

## GATA4
g4.box <- ggplot(data = rbind.fill(subset(new.lms, TE_ICM == 'ICM' & 
                                            Experiment %in% g4), 
                                   subset(ablat.lms, TE_ICM == 'ICM' & 
                                            Experiment %in% g4)), 
                 aes(x = Identity.hc, y = CH3.ebLogCor.s))
g4.box <- g4.box + geom_boxplot(fill = precols[5], color = 'black', 
                                outlier.shape = 1, outlier.stroke = 0.5) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)
# g4.box <- g4.box + geom_jitter(shape = 20, width = 0.2, 
#                                  show.legend = F)
g4.box <- g4.box + looks + theme(aspect.ratio = 1)
g4.box <- g4.box + labs(x = 'Identity', y = 'log[GATA4]', 
                        title = 'Fig S1d')
# print(g4.box)

# Generate PDF
pdf(file = './figures/figS1d_NS_S17G4boxplots.pdf', width = 11, paper = 'a4r')
print(s17.box)
print(g4.box)
dev.off()

################################################################################
# Supplementary Figure 1e
# 3D scatter plots of NANOG, GATA6 and either SOX17 or GATA4 levels 
# for embryos in Figure S1d above. 
# Cells are color coded for lineage as assigned using NANOG vs GATA6 only
################################################################################

# Define aesthetics for the plotly plot
my.axistitle <- list(family = 'Helvetica', size = 18, color = 'black')
my.ticks <- list(family = 'Helvetica', size = 14, color = 'black')

# Generate plots
## NANOG vs GATA6 vs SOX17
s17.plot <- plot_ly(rbind.fill(subset(new.lms, TE_ICM == 'ICM' & 
                                        Experiment %in% s17), 
                               subset(ablat.lms, TE_ICM == 'ICM' & 
                                        Experiment %in% s17)), 
                    y = ~CH3.ebLogCor.xs, z = ~CH5.ebLogCor.s, 
                    x = ~CH5.ebLogCor.xs, color = ~Identity.hc, 
                    colors = idcols, 
                    marker = list(size = 4)) %>% add_markers() %>% 
  layout(scene = list(xaxis = list(title = '[GATA6] (AU)', 
                                   gridcolor = 'rgb(0, 0, 0)', 
                                   titlefont = my.axistitle, 
                                   tickfont = my.ticks), 
                      yaxis = list(title = '[NANOG] (AU)',
                                   gridcolor = 'rgb(0, 0, 0)', 
                                   titlefont = my.axistitle, 
                                   tickfont = my.ticks), 
                      zaxis = list(title = '[SOX17] (AU)', 
                                   gridcolor = 'rgb(0, 0, 0)', 
                                   titlefont = my.axistitle, 
                                   tickfont = my.ticks)), 
         paper_bgcolor = 'rgb(255, 255, 255)',
         plot_bgcolor = 'rgb(255, 255, 255)')
print(s17.plot)

## NANOG vs GATA6 vs GATA4
g4.plot <- plot_ly(rbind.fill(subset(new.lms, TE_ICM == 'ICM' & 
                                       Experiment %in% g4), 
                              subset(ablat.lms, TE_ICM == 'ICM' & 
                                       Experiment %in% g4)), 
                   y = ~CH3.ebLogCor.xs, z = ~CH3.ebLogCor.s, 
                   x = ~CH5.ebLogCor.xs, color = ~Identity.hc, 
                   colors = idcols, 
                   marker = list(size = 4)) %>% add_markers() %>% 
  layout(scene = list(xaxis = list(title = '[GATA6] (AU)', 
                                   gridcolor = 'rgb(0, 0, 0)', 
                                   titlefont = my.axistitle, 
                                   tickfont = my.ticks), 
                      yaxis = list(title = '[NANOG] (AU)',
                                   gridcolor = 'rgb(0, 0, 0)', 
                                   titlefont = my.axistitle, 
                                   tickfont = my.ticks), 
                      zaxis = list(title = '[GATA4] (AU)', 
                                   gridcolor = 'rgb(0, 0, 0)', 
                                   titlefont = my.axistitle, 
                                   tickfont = my.ticks)), 
         paper_bgcolor = 'rgb(255, 255, 255)',
         plot_bgcolor = 'rgb(255, 255, 255)')
print(g4.plot)

# RStudio cannot export these dynamic plots as PDFs, so I just screenshot them





