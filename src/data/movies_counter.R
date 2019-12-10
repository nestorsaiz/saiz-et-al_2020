# This script performs calculations on movies dataset and generates
# dataframes that will be used for plotting

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load necessary data
if(exists('movies') == F) {
  movies <- read.csv('./data/processed/movies-all-processed.csv')
  endings <- read.csv('./references/movies-endings.csv')
}

if(exists('endings') == F) {
  source('./src/data/movies_classify.R')
}

################################################################################
# Calculate lineage size for each embryo
################################################################################

# At t=0 (time of experiment)
count.t0 <- movies %>% filter(Channel == 1, timeframe == 1) %>% 
  group_by(Embryo_ID, identity.t, cell_treatment, .drop = F) %>% 
  summarize(t0.count = n())

# Over time
gro.movs <- movies %>% filter(Channel == 1) %>% 
  group_by(Embryo_ID, identity.t, target, 
           Treatment, Cell_diff, Stage.t0, 
           # division, death, 
           timeframe, hours) %>% 
  summarize(count.t = n())

# Order Cell difference factors
gro.movs$Cell_diff <- factor(gro.movs$Cell_diff, 
                             levels = c('zero', 'kill_0.25', 
                                        'kill_0.5', 'kill_0.75', 
                                        'kill_all'))

################################################################################
# Calculate the frequency of each cell behavior, namely:
# * division 
# * death
# * fate choice (switch)
################################################################################

# Obtain timing for each event
events.t <- movies %>% 
  filter(Channel == 1) %>% 
  group_by(Embryo_ID, Stage.t0, Cell_diff, 
           Treatment, cell_treatment, target, 
           division, death, switch, 
           timeframe, hours, identity.t) %>% 
  summarize()

# Count the maximum number of cells per lineage in each embryo
count.max <- movies %>% filter(Channel == 1) %>% 
  group_by(Embryo_ID, target, Stage.t0, 
           Treatment, cell_treatment, Cell_diff, 
           timeframe, hours, identity.t) %>% 
  summarize(count.t = n()) %>% 
  group_by(Embryo_ID, target, Stage.t0, 
           Treatment, Cell_diff, cell_treatment, 
           identity.t) %>% 
  summarize(max.count = max(count.t))

# Calculate the number of cells present at t = 0 that divide, for each embryo
count.mito <- movies %>% filter(Channel == 1, 
                                timeframe == 1) %>% 
  group_by(Embryo_ID, identity.t, Cell_ID, .drop = F) %>% 
  filter(Cell_ID %in% 
           endings$Cell_ID[which(endings$ending == 'division')]) %>% 
  summarize() %>% 
  group_by(Embryo_ID, identity.t, .drop = F) %>% 
  summarize(n.mito = n())

# Calculate the total number of cell deaths for each compartment, per embryo
count.ded <- movies %>% filter(Channel == 1, death == T) %>% 
  group_by(Embryo_ID, identity.t, cell_treatment, .drop = F) %>% 
  summarize(n.deads = n())

# Combine all into one table of events 
# (dplyr::join_all() causes a fatal error, don't use here)
## Divisions
count.mito <- merge(count.t0, count.mito, all.x = T)
count.mito$n.mito[is.na(count.mito$n.mito)] <- 0
## Deaths
count.ded <- merge(count.max, count.ded, all.x = T)
count.ded$n.deads[is.na(count.ded$n.deads)] <- 0
## Combine both
count.events <- merge(count.mito, count.ded)

# Calculate the fraction of each initial compartment that divide and
# the fraction of max size of each compartment that die
count.events$pc.mito <- count.events$n.mito / count.events$t0.count * 100
count.events$pc.ded <- count.events$n.deads / count.events$max.count * 100

# Extract relevant metadata only for each embryo
mov.meta <- movies %>% 
  group_by(Embryo_ID, Treatment, 
           target, Cell_diff, Stage.t0) %>% 
  summarize()
# and combine with count.events
count.events <- merge(count.events, mov.meta)

# Order factors in count.events
count.events$Cell_diff <- factor(count.events$Cell_diff, 
                                 levels = c('zero', 'kill_0.25', 'kill_0.5', 
                                            'kill_0.75', 'kill_all'))
count.events$Treatment <- factor(count.events$Treatment, 
                                 levels = c('Control', 'Ablated'))



