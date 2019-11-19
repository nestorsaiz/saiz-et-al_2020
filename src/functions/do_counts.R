## This function applies the dplyr package to perform 
## some basic calculations on the raw dataset and return, as new variables:
## 1. Total cell count per embryo (Cellcount)
## 2. Number of ICM cells per embryo (icm.count)
## 3. Average cell count per litter or experimental group 
##    (litter.median or group.median)

## It requires dplyr to run
library("plyr")
library("dplyr")

## This function takes a dataset as argument, 
## as well as the names given to the variables for:
## * embryo identifier (embryo.var, by default 'Embryo_ID'), 
## * experimental unit (expunit.var, by default 'Experiment') 
## * treatment (tt.var, by default 'Treatment')
## It also takes a logical option to separate the experimental unit 
## by experimental treatment (if TRUE, the default) or not (if FALSE) - 
## this option allows the user to calculate the median cell count 
## for the entire experimental unit or for each of the treatment groups
## given by tt.var

do.counts <- function(dataset, embryo.var = 'Embryo_ID', 
                      TI.var = 'TE_ICM', expunit.var = 'Experiment', 
                      tt.var = 'Treatment', sep.treatment = T) {
        ## Count the number of cells per embryo and add to main table
        counts <- dataset %>% group_by(!! as.name(embryo.var)) %>% 
                summarize(Cellcount = n())
        dataset <- merge(dataset, counts)
        ## Extract unique values of TI.var in dataset
        TI.vals <- as.data.frame(unique(dataset[which(colnames(dataset) == 
                                                              TI.var)]))
        ## Find which one may correspond to trophectoderm
        te <- TI.vals[, 1][which(TI.vals[, 1] %in% c('TE', 'Te', 'te', 'T', 
                                                     'Trophectoderm', 'TB', 
                                                     'tb', 'Tb', 
                                                     'Trophoblast'))]
        ## Calculate the number of ICM cells per embryo
        tecounts <- dataset %>% filter(!! as.name(TI.var) == te) %>% 
                group_by(!! as.name(embryo.var)) %>% 
                summarize(te.count = n())
        dataset <- merge(dataset, tecounts)
        dataset$icm.count <- dataset$Cellcount - dataset$te.count
        dataset$te.count <- NULL
        ## Calculate the median cellcount per either experimental group 
        ## (each litter, or whatever has been defined as such) 
        ## or per treatment group within the experimental group
        if (sep.treatment == T) { 
                med.litter <- dataset %>% 
                        group_by(!! as.name(expunit.var), 
                                 !! as.name(tt.var)) %>% 
                        summarize(group.median = median(Cellcount))
        }
        else { 
                med.litter <- dataset %>% 
                        group_by(!! as.name(expunit.var)) %>% 
                        summarize(litter.median = median(Cellcount))
        }
        ## Combine with main table and return
        dataset <- merge(dataset, med.litter)
        return(dataset)
}
