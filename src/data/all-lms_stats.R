# This script calculates summary statistics for each stage
# for the lineage counts of each set of littermates:
# * mean ICM size for each stage
# * N embryos for each set and stage
# * mean % of ICM that each lineage represents per stage
# * median % of ICM that each lineage represents per stage
# * SD for the % of ICM that each lineage represents per stage

# Check if setup.R has been ran
setup.ran <- exists('looks')
if (setup.ran == F) { 
  source('./src/setup.R')
}
rm(setup.ran)

# Check if data has been loaded, and run all-lms_read.R otherwise
if(exists('allcounts.list') == F) { 
  source('./src/data/all-lms_read.R')
}

# Calculate stats for the Nat Comms (Saiz et al., (2016)) dataset
ncoms.sum <- allcounts.list[[1]] %>% filter(TE_ICM == 'ICM', 
                                            !Stage %in% c('16_32', '>150'), 
                                            !Identity %in% c('EPI', 'DN')) %>% 
  group_by(Stage, Gene1, Genotype1, Identity) %>% 
  summarize(mean.icm.size = mean(icm.count), 
            N.embryos = n(), 
            mean.pc.icm = mean(pc.icm), 
            median.pc.icm = median(pc.icm), 
            sd.pc.icm = sd(pc.icm))

# Calculate stats for wild type embryos in the 
# Spry4 (Morgani et al., (2018)) dataset
spry.sum.wt <- allcounts.list[[2]] %>% filter(TE_ICM == 'ICM', 
                                              !Stage %in% c('16_32', '>150'), 
                                              !Identity %in% c('morula', 'EPI', 
                                                               'DN'), 
                                              Genotype1 == 'wt') %>% 
  group_by(Stage, Gene1, Genotype1, Identity) %>% 
  summarize(mean.icm.size = mean(icm.count), 
            N.embryos = n(), 
            mean.pc.icm = mean(pc.icm), 
            median.pc.icm = median(pc.icm), 
            sd.pc.icm = sd(pc.icm))

# Calculate stats for Spry4:H2B-Venus/+ embryos in the 
# Spry4 (Morgani et al., (2018)) dataset
spry.sum.het <- allcounts.list[[2]] %>% filter(TE_ICM == 'ICM', 
                                               !Stage %in% c('16_32', '>150'), 
                                               !Identity %in% c('morula', 'EPI', 
                                                                'DN'), 
                                               Genotype1 == 'het') %>% 
  group_by(Stage, Gene1, Genotype1, Identity) %>% 
  summarize(mean.icm.size = mean(icm.count), 
            N.embryos = n(), 
            mean.pc.icm = mean(pc.icm), 
            median.pc.icm = median(pc.icm), 
            sd.pc.icm = sd(pc.icm))

# Calculate stats for the ablation littermates dataset
ablat.lms.sum <- allcounts.list[[3]] %>% filter(TE_ICM == 'ICM', 
                                             !Stage %in% c('16_32', '>150'), 
                                             !Identity %in% c('EPI', 'DN'), 
                                             Genotype1 == 'wt') %>% 
  group_by(Stage, Genotype1, Identity) %>% 
  summarize(mean.icm.size = mean(icm.count), 
            N.embryos = n(), 
            mean.pc.icm = mean(pc.icm), 
            median.pc.icm = median(pc.icm), 
            sd.pc.icm = sd(pc.icm))

# Calculate stats for the new littermates dataset
new.lms.sum <- allcounts.list[[4]] %>% filter(TE_ICM == 'ICM', 
                                             !Stage %in% c('16_32', '>150'), 
                                             !Identity %in% c('EPI', 'DN'), 
                                             Genotype1 == 'wt') %>% 
  group_by(Stage, Genotype1, Identity) %>% 
  summarize(mean.icm.size = mean(icm.count), 
            N.embryos = n(), 
            mean.pc.icm = mean(pc.icm), 
            median.pc.icm = median(pc.icm), 
            sd.pc.icm = sd(pc.icm))

# Calculate stats for the Fgf4+/- littermates dataset
fgf4.sum <- allcounts.list[[6]] %>% filter(TE_ICM == 'ICM', 
                                           !Stage %in% c('16_32', '>150'), 
                                           !Identity %in% c('EPI', 'DN'), 
                                           Genotype1 == 'het') %>% 
  group_by(Stage, Gene1, Genotype1, Identity) %>% 
  summarize(mean.icm.size = mean(icm.count), 
            N.embryos = n(), 
            mean.pc.icm = mean(pc.icm), 
            median.pc.icm = median(pc.icm), 
            sd.pc.icm = sd(pc.icm))

# r1.sum <- allcounts.list[[8]] %>%
#   filter(TE_ICM == 'ICM', 
#          !Stage %in% c('16_32', '>150'), 
#          !Identity %in% c('EPI', 'DN'), 
#          Gene1 == 'Fgfr1', 
#          Genotype1 == 'ko', 
#          Genotype2 == 'wt') %>% 
#   group_by(Stage, Gene1, Genotype1, 
#            Gene2, Genotype2, Identity) %>% 
#   summarize(mean.icm.size = mean(icm.count), 
#             N.embryos = n(), 
#             mean.pc.icm = mean(pc.icm), 
#             median.pc.icm = median(pc.icm), 
#             sd.pc.icm = sd(pc.icm))
# 
# r2.sum <- allcounts.list[[8]] %>%
#   filter(TE_ICM == 'ICM', 
#          !Stage %in% c('16_32', '>150'), 
#          !Identity %in% c('EPI', 'DN'), 
#          Gene2 == 'Fgfr2', 
#          Genotype2 == 'ko', 
#          Genotype1 == 'wt') %>% 
#   group_by(Stage, Gene1, Genotype1, 
#            Gene2, Genotype2, Identity) %>% 
#   summarize(mean.icm.size = mean(icm.count), 
#             N.embryos = n(), 
#             mean.pc.icm = mean(pc.icm), 
#             median.pc.icm = median(pc.icm), 
#             sd.pc.icm = sd(pc.icm))
# 
# ng.sum <- allcounts.list[[4]] %>% 
#   filter(TE_ICM == 'ICM', 
#          !Identity %in% c('EPI', 'DN'), 
#          Gene1 == 'Gata6', 
#          Gene2 == 'Nanog') %>% 
#   group_by(Stage, Gene1, Genotype1, 
#            Gene2, Genotype2, Identity) %>% 
#   summarize(mean.icm.size = mean(icm.count), 
#             N.embryos = n(), 
#             mean.pc.icm = mean(pc.icm), 
#             median.pc.icm = median(pc.icm), 
#             sd.pc.icm = sd(pc.icm))

write.csv(ncoms.sum, file = './results/stats-NComms.csv', row.names = F)
write.csv(spry.sum.wt, file = './results/stats-Spry4-wt.csv', row.names = F)
write.csv(spry.sum.het, file = './results/stats-Spry4-het.csv', row.names = F)
write.csv(ablat.lms.sum, file = './results/stats-ablat-lms.csv', row.names = F)
write.csv(new.lms.sum, file = './results/stats-new-lms.csv', row.names = F)
write.csv(fgf4.sum, file = './results/stats-Fgf4-het.csv', row.names = F)

