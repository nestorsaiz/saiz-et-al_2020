## This function applies the dplyr package to perform 
## three basic calculations on the raw dataset:
## 1. Total cell count per embryo
## 2. Number of ICM cells per embryo
## 3. Average cell count per litter or experimental group

library("plyr")
library("dplyr")

do.counts <- function(dataset) {
        ## Count the number of cells per embryo and add to main table
        counts <- dataset %>% group_by(Embryo_ID) %>% 
                summarize(Cellcount = n())
        dataset <- merge(dataset, counts)
        ## Calculate the number of ICM cells per embryo
        tecounts <- dataset %>% filter(TE_ICM == 'TE') %>% 
                group_by(Embryo_ID) %>% 
                summarize(te.count = n())
        dataset <- merge(dataset, tecounts)
        dataset$icm.count <- dataset$Cellcount - dataset$te.count
        dataset$te.count <- NULL
        ## Calculate the median cellcount per litter
        med.litter <- dataset %>% group_by(Experiment, Treatment) %>% 
                summarize(litter.median = median(Cellcount))
        ## Combine with main table and remove avg table
        dataset <- merge(dataset, med.litter)
        return(dataset)
}
