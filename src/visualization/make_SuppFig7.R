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
# Figure S7d-e
#
################################################################################
