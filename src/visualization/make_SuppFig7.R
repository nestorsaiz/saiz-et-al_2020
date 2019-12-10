# This script generates the plots in Supplementary Figure 7

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load necessary data
if(exists('movies') == F) {
  movies <- read.csv('./data/processed/movies-all-processed.csv')
}

# Generate tables with cell and events counts 
if(exists('count.events') == F) { 
  source('./src/data/movies_counter.R')
}

# Order factors for plotting purposes
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
count.events$identity.t <- factor(count.events$identity.t, 
                                  levels = levels(movies$identity.t))

################################################################################
# Figure S7a-c
# Three sets of plots for each target group ('none', 'PRE' and 'EPI'),
# for movies taken after perturbation was done at the 70-90 cell stage
# (same embryos as in Figure 5, panels g, h, i):
# * Growth curves for each ICM lineage in individual embryos 
# * Dot plots showing PrE deaths, divisions and new PrE from DP differentiation
# * Dot plots showing EPI deaths, divisions and new EPI from DP differentiation
################################################################################

# Create indexes for each target group and for each final ICM lineage (PRE, EPI)
targets <- c('none', 'PRE', 'EPI')
ids <- c('PRE', 'EPI')

# Generate all 9 plots printed onto a PDF
pdf(file = './figures/figS7a_NS.pdf', width = 11, paper = 'a4r')
for(t in targets) { 
  ## Define data to go in growth plot
  gro.data <- gro.movs %>% filter(target ==  t, 
                                  !Cell_diff %in% c('kill_0.25'), 
                                  Stage.t0 == '[70,90)', 
                                  hours < 16)
  # Generate lineage growth curves
  grop <- ggplot(data = gro.data, 
                 aes(x = hours, y = count.t, 
                     color = identity.t))
  grop <- grop + geom_line(size = 0.75)
  grop <- grop + 
    scale_color_manual(values = idcols) + 
    looks + 
    labs(title = paste('Target', t, sep = ' ')) + 
    theme(aspect.ratio = 1) + 
    facet_wrap( ~ Cell_diff + Embryo_ID, nrow = 1)
  print(grop) 
  
  # Define events data to go in dot plots (death, division, switches)
  evt.data <- events.t %>%
    filter(target == t,
           Cell_diff != 'kill_0.25',
           Stage.t0 == '[70,90)',
           hours < 16)
  
  # Generate two equivalent dot plots, one for PrE, one for EPI
  for(i in ids) {
    evtplot <- ggplot(data = evt.data,
                      aes(x = hours))
    evtplot <- evtplot +
      geom_dotplot(data = evt.data %>%
                     filter(death == T,
                            identity.t == i),
                   color = 'black', fill = 'black',
                   binwidth = 0.5, dotsize = 1,
                   stackgroups = T, method = 'histodot') +
      geom_dotplot(data = evt.data %>%
                     filter(division == T,
                            identity.t == i),
                   color = 'black', fill = 'cyan2',
                   binwidth = 0.5, dotsize = 1,
                   stackgroups = T, method = 'histodot') +
      geom_dotplot(data = evt.data %>%
                     filter(switch == T,
                            identity.t == i),
                   color = 'black', fill = 'orange2',
                   binwidth = 0.5, dotsize = 1,
                   stackgroups = T, method = 'histodot')
    evtplot <- evtplot + looks +
      labs(title = paste('Target', t,
                         'Cell type', i,
                         sep = ' ')) +
      theme(aspect.ratio = 1) +
      facet_wrap( ~ Cell_diff + Embryo_ID, nrow = 1)
    print(evtplot)
  }
}
dev.off()

################################################################################
# Figure S7d
# Pdgfra:H2B-GFP dynamics in each cell of one embryo from b,
# in which 75% of PrE was killed
################################################################################

# Set index to embryo "010719Abl_EH5" (-75% PrE)
my.pre <- unique(subset(movies, Treatment == 'Ablated' & 
                          target == 'PRE')$Embryo_ID)
i <- 3

# Generate plot
panel.d <- ggplot(data = subset(movies, Channel == 1 & 
                                  Embryo_ID == my.pre[i]), 
                  aes(x = hours, y = mavg))
panel.d <- panel.d + geom_line(aes(color = identity.t, 
                                   group = interaction(TrackID, Cell_ID)), 
                               size = 0.75) + 
  geom_point(data = subset(movies, Channel == 1 & 
                             Embryo_ID == my.pre[i] & 
                             switch == T), 
             aes(x = hours, y = mavg), 
             shape = 18, size = 2) +
  geom_point(data = subset(movies, Channel == 1 & 
                             Embryo_ID == my.pre[i] & 
                             death == T), 
             aes(x = hours, y = mavg), 
             shape = 4, size = 2)
panel.d <- panel.d + facet_wrap( ~ interaction(TrackID, identity_t0.th) + 
                                   cell_treatment)
panel.d <- panel.d + looks + 
  scale_color_manual(values = idcols) + 
  labs(title = paste(my.pre[i]), y = 'log(PdgfraH2B-GFP)') + 
  theme(aspect.ratio = 0.75, 
        # strip.text.x = element_text(size = 8))
        strip.text.x = element_blank())
# Uncomment print() below to visualize plot
# print(panel.d)

################################################################################
# Figure S7e
# Pdgfra:H2B-GFP dynamics in each cell of one embryo from c,
# in which 75% of EPI was killed
################################################################################

# Set index to embryo "102218Abl_DV3" (-75% EPI)
my.epi <- unique(subset(movies, Treatment == 'Ablated' & 
                          target == 'EPI')$Embryo_ID)
i <- 5

# Generate plot
panel.e <- ggplot(data = subset(movies, Channel == 1 & 
                                  Embryo_ID == my.epi[i]), 
                  aes(x = hours, y = mavg))
panel.e <- panel.e + geom_line(aes(color = identity.t, 
                                   group = interaction(TrackID, Cell_ID)), 
                               size = 0.75) + 
  geom_point(data = subset(movies, Channel == 1 & 
                             Embryo_ID == my.epi[i] & 
                             switch == T), 
             aes(x = hours, y = mavg), 
             shape = 18, size = 2) +
  geom_point(data = subset(movies, Channel == 1 & 
                             Embryo_ID == my.epi[i] & 
                             death == T), 
             aes(x = hours, y = mavg), 
             shape = 4, size = 2)
panel.e <- panel.e + facet_wrap( ~ interaction(TrackID, identity_t0.th) + 
                                   cell_treatment)
panel.e <- panel.e + looks + 
  scale_color_manual(values = idcols) + 
  labs(title = paste(my.epi[i]), y = 'log(PdgfraH2B-GFP)') + 
  theme(aspect.ratio = 0.75, 
        # strip.text.x = element_text(size = 8))
        strip.text.x = element_blank())
# Uncomment print() below to visualize plot
# print(panel.e)

################################################################################
# Generate PDF
################################################################################

pdf(file = './figures/figS7d-e_NS.pdf', width = 11, paper = 'a4r')
print(panel.d)
print(panel.e)
dev.off()

################################################################################
# Figure 7f
# Box plots showing the % of cell death for each lineage
# in each experimental group
################################################################################

# Define data to plot
my.events <- count.events %>% 
  filter(Stage.t0 %in% c('[70,90)'), 
         Cell_diff != 'kill_0.25', 
         cell_treatment == 'intact')
my.target <- c('PRE', 'EPI')

# Generate plot for each target (PRE and EPI)
pdf(file = './figures/figS7f-g_NS.pdf', width = 11, paper = 'a4r')
for(i in 1:length(my.target)) {
  panel.f <- ggplot(data = my.events %>% 
                      filter(target %in% c('none', my.target[i])), 
                    aes(x = Cell_diff, 
                        y = pc.ded))
  panel.f <- panel.f + geom_boxplot(aes(fill = identity.t), color = 'black',
                                    outlier.shape = 1, outlier.size = 2) +
    stat_summary(fun.y = mean, colour = "black", geom = "point", 
                 shape = 4, size = 3, show.legend = F) + 
    geom_jitter(color = 'black', shape = 20, width = 0.2)
  panel.f <- panel.f + 
    facet_grid(Stage.t0 ~ identity.t) + 
    labs(title = paste('Figure S7f, ', 'Target:', my.target[i]), 
         x = 'Experimental group', y = '% cell death') + 
    scale_fill_manual(values = idcols) + 
    looks + theme(aspect.ratio = 1, 
                  axis.text.x = element_text(angle = 30, hjust = 1), 
                  strip.text.x = element_text(size = 10))
  print(panel.f)
}

################################################################################
# Figure 7g
# Box plots showing the % of cell division foreach lineage 
# in each experimental group
################################################################################

# Generate plot for each target (PRE and EPI)
for(i in 1:length(my.target)) {
  panel.g <- ggplot(data = my.events %>% 
                     filter(target %in% c('none', my.target[i])),
                   aes(x = Cell_diff, 
                       y = pc.mito))
  panel.g <- panel.g + geom_boxplot(aes(fill = identity.t), color = 'black', 
                                  outlier.shape = 1, outlier.size = 2) + 
    stat_summary(fun.y = mean, colour = "black", geom = "point", 
                 shape = 4, size = 3, show.legend = F) + 
    geom_jitter(color = 'black', shape = 20, width = 0.2)
  panel.g <- panel.g + 
    facet_grid(Stage.t0 ~ identity.t) + 
    labs(title = paste('Figure S7f, ', 'Target:', my.target[i]), 
         x = 'Experimental group', y = '% cell division') + 
    scale_fill_manual(values = idcols) + 
    looks + theme(aspect.ratio = 1, 
                  axis.text.x = element_text(angle = 30, hjust = 1), 
                  strip.text.x = element_text(size = 10))
  print(panel.g)
}
dev.off()

##