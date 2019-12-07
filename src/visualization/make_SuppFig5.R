# This script generates the plots in Supplementary Figure 5

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Load necessary data
if(exists('ablat.t0') == F) {
  ablat.t0 <- read.csv('./data/processed/ablat-t0-processed.csv')
}

if(exists('new.lms') == F) { 
  new.lms <- read.csv('./data/processed/new-lms-processed.csv')
  unicos.nl <- read.csv('./references/new-lms_unicos.csv')
}

if(exists('live.lincounts') == F) {
  live.lincounts <- read.csv('./data/interim/live-lincounts.csv')
  fixed.lincounts <- read.csv('./data/interim/fixed-lincounts.csv')
}

if(exists('movies') == F) {
  movies <- read.csv('./data/processed/movies-all-processed.csv')
}

# Order levels of factors for plotting
ablat.t0$Identity <- factor(ablat.t0$Identity, levels = c('PRE', 'DP', 'EPI'))
new.lms$Identity.hc <- factor(new.lms$Identity.hc, 
                                  levels = c('TE', 'PRE', 'DP', 'EPI', 'DN'))
new.lms$Stage <- factor(new.lms$Stage,
                            levels = c('<8', '8_16', '16_32', '32_64', 
                                       '64_90', '90_120', '120_150', '>150'))
live.lincounts$Identity <- factor(live.lincounts$Identity, 
                                  levels = levels(ablat.t0$Identity))
fixed.lincounts$Identity <- factor(fixed.lincounts$Identity, 
                                   levels = levels(live.lincounts$Identity))
live.lincounts$Stage.t0 <- factor(live.lincounts$Stage.t0, 
                                  levels = c('[30,50)', '[50,70)', '[70,90)', 
                                             '[90,110)', '[110,130)', 
                                             '[130,150)', '[150,170)', 
                                             '[170,190)', '[190,210]'))
fixed.lincounts$Stage.t0 <- factor(fixed.lincounts$Stage.t0, 
                                   levels = levels(live.lincounts$Stage.t0))
movies$identity_t0.th <- factor(movies$identity_t0.th, 
                                levels = c('PRE', 'DP', 'EPI'))

################################################################################
# Figure S5b
# Box plots comparing H2B-GFP levels in each ICM lineage 
# of live Pdgfra:H2B-GFP/+ embryos prior to cell ablation (t=0).
# Lineages in this dataset were assigned manually at the time of experiment
################################################################################

# Define data to plot
pd.t0 <- ablat.t0 %>% 
  group_by(Embryo_ID, Identity, 
           Stage.t0, Genotype2) %>% 
  summarize(pd.mean = mean(CH1.ebLogCor))  %>% 
  filter(Identity %in% c('EPI', 'PRE', 'DP'), 
         Stage.t0 %in% c('[30,50)', '[50,70)', 
                         '[70,90)', '[90,110)'), 
         Genotype2 %in% c('mKate/+', 
                          'mKate/mKate'))

# Generate plot
panel.b <- ggplot(data = pd.t0, 
                aes(x = Stage.t0, y = pd.mean))
panel.b <- panel.b + geom_boxplot(color = 'black', 
                              outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)  
panel.b <- panel.b + geom_jitter(aes(color = Identity), 
                             shape = 20, size = 1.5, width = 0.2)
panel.b <- panel.b + scale_color_manual(values = idcols) + looks 
panel.b <- panel.b + facet_wrap( ~ Identity, nrow = 1)
panel.b <- panel.b + theme(aspect.ratio = 2) 
panel.b <- panel.b + labs(y = 'log[Pdgfra:H2B-GFP] (mean per embryo)', 
                      x = 'Stage (number of cells)', 
                      title = 'Figure S5b')
# Uncomment print() below to visualize plot
# print(panel.b)

################################################################################
# Figure S5d
# Box plots comparing H2B-GFP levels in each ICM lineage
# of fixed Pdgfra:H2B-GFP/+ embryos after staining for NANOG and GATA6.
# Lineages in these dataset were assigned automatically 
# using Hierarchical clustering based on NANOG and GATA6 levels.
################################################################################

# Extract only fixed littermates that are Pdgfra:H2B-GFP/+
pd.lms <- subset(new.lms, Experiment %in% 
                   unicos.nl$Experiment[which(unicos.nl$CH2 == 'GFP')] & 
                   Gene1 == 'Pdgfra' & 
                   Genotype1 == 'het')

# Re-stage using the same 20-cell bins as for t0 embryos
cutz <- seq(30, 210, by = 20)
pd.lms$Stage.t0 <- cut(pd.lms$Cellcount, 
                       breaks = cutz, right = F, include.lowest = T)
rm(new.lms, unicos.nl)

# Calculate the average level of GFP in each ICM cell type per embryo
# and define data to plot
pd.lms <- pd.lms %>% 
  group_by(Embryo_ID, Identity.hc, Stage.t0) %>% 
  summarize(pd.mean = mean(CH2.ebLogCor)) %>% 
  filter(Identity.hc %in% c('EPI', 'PRE', 'DP'), 
         Stage.t0 %in% c('[30,50)', '[50,70)', 
                         '[70,90)', '[90,110)'))

# Generate plot
panel.d <- ggplot(data = pd.lms, 
                 aes(x = Stage.t0, y = pd.mean))
panel.d <- panel.d + geom_boxplot(color = 'black', 
                                outlier.shape = 1, outlier.size = 2) +
  stat_summary(fun.y = mean, colour = "black", geom = "point", 
               shape = 4, size = 3, show.legend = F)  
panel.d <- panel.d + geom_jitter(aes(color = Identity.hc), 
                               shape = 20, size = 1.5, width = 0.2)
panel.d <- panel.d + scale_color_manual(values = idcols) + looks 
panel.d <- panel.d + facet_wrap( ~ Identity.hc, nrow = 1)
panel.d <- panel.d + theme(aspect.ratio = 2) 
panel.d <- panel.d + labs(y = 'log[Pdgfra:H2B-GFP] (mean per embryo)', 
                        x = 'Stage (number of cells)', 
                        title = 'Figure S5d')
# Uncomment print() below to visualize plot
# print(panel.d)

################################################################################
# Figure S5e
# Compare lineage growth over time in live and fixed embryos
################################################################################

# Define data to plot
my.data <- rbind(live.lincounts, fixed.lincounts) %>% 
  filter(Stage.t0 %in% c('[30,50)', '[50,70)', 
                         '[70,90)', '[90,110)', 
                         '[110,130)'), 
         Identity != 'DN', 
         Cellcount %in% seq(50, 110, by = 1))

# Generate plot
panel.e <- ggplot(data = my.data, 
                     aes(x = Cellcount, 
                         y = count))
panel.e <- panel.e + geom_jitter(aes(color = dataset), 
                                       shape = 20, size = 1.5, 
                                       width = 0.2)
panel.e <- panel.e + geom_smooth(aes(color = dataset), 
                                       size = 1.2, span = 1)
panel.e <- panel.e + scale_color_manual(values = c('live' = 'green', 
                                                         'fixed' = 'black')) 
panel.e <- panel.e + looks + theme(aspect.ratio = 1)
panel.e <- panel.e + facet_wrap( ~ Identity, nrow = 1) 
panel.e <- panel.e + scale_y_continuous(breaks = seq(0, 30, by = 5)) + 
  scale_x_continuous(breaks = seq(50, 130, by = 10))
panel.e <- panel.e + labs(x = 'Cell count', y = 'Number of cells', 
                          title = 'Figure S5e')
# Uncomment print() below to visualize plot
# print(panel.e)

################################################################################
# Figure S5f
# Plot dynamics of Pdgfra expression and cell identity 
# in a subset of control embryos (2x embryos for each stage)
################################################################################

# List selected embryos
my.controls <- c("012119Abl_EI6", "013118Abl_CM4", 
                 "012318Abl_CJ4", "121018Abl_EB5", 
                 "020618Abl_CQ5", "020618Abl_CQ7")

# Define data to plot
my.data <- subset(movies, Channel == 1 & 
         Embryo_ID %in% my.controls &
         hours < 16)
# Separate data into non-DP and DP cells
my.data$panel <- 'top'
my.data$panel[which(my.data$identity_t0.th == 'DP')] <- 'bottom'
my.data$panel <- factor(my.data$panel, levels = c('top', 'bottom'))

# Generate top plot
panel.f <- ggplot(data = my.data, 
                aes(x = hours, y = mavg))
panel.f <- panel.f + geom_line(aes(color = identity.t,
                               group = interaction(TrackID, Cell_ID)),
                           size = 0.75, alpha = 0.5)
panel.f <- panel.f + geom_smooth(data = my.data %>% 
                                   filter(identity_t0.th != 'DP'), 
                                 aes(color = identity_t0.th), size = 1.5)
panel.f <- panel.f + facet_grid(panel ~ Stage.t0)
panel.f <- panel.f + 
  scale_color_manual(values = idcols) + 
  looks + 
  theme(aspect.ratio = 0.75, 
        strip.text.x = element_text(size = 10)) + 
  labs(x = 'Time (hours)', y = 'log[Pdgfra:H2B-GFP]', 
       title = 'Figure S5f')
# Uncomment print() below to visualize plot
# print(panel.f)

################################################################################
# Figure S5g
# Dynamics of Pdgfra expression for one control embryo ("020618Abl_CQ5")
# at the 90-100 cell stage at t=0, for each cell individually,
# highlighting apoptosis and cell fate switches
################################################################################

# Define data to plot
i <- 5
my.data <- subset(movies, Channel == 1 & 
                    Embryo_ID == my.controls[i])

# Generate plot
panel.g <- ggplot(data = my.data, 
                  aes(x = hours, y = mavg))
panel.g <- panel.g + geom_line(aes(color = identity.t, 
                                   group = interaction(TrackID, Cell_ID)), 
                               size = 0.75) + 
  geom_point(data = subset(movies, Channel == 1 & 
                             Embryo_ID == my.controls[i] & 
                             switch == T), 
             aes(x = hours, y = mavg), 
             shape = 18, size = 2) +
  geom_point(data = subset(movies, Channel == 1 & 
                             Embryo_ID == my.controls[i] & 
                             death == T), 
             aes(x = hours, y = mavg), 
             shape = 4, size = 2)
panel.g <- panel.g + facet_wrap( ~ interaction(TrackID, identity_t0.th))
panel.g <- panel.g + looks + 
  scale_color_manual(values = idcols) + 
  labs(title = paste('Figure S5g. Embryo:', 
                     my.controls[i], sep = ' '), 
       y = 'log(PdgfraH2B-GFP)', x = 'Time (hours)') + 
  theme(aspect.ratio = 0.75, 
        # strip.text.x = element_text(size = 8))
        strip.text.x = element_blank())
# Uncomment print() below to visualize plot
# print(panel.g)

################################################################################
# Generate PDF
################################################################################

pdf(file = './figures/figS5all_NS.pdf', width = 11, paper = 'a4r')
print(panel.b)
print(panel.d)
print(panel.e)
print(panel.f)
print(panel.g)
dev.off()


##