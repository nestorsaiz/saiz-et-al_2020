## This function applies the dplyr package to perform 
## three basic calculations on the raw dataset:
## 1. Total cell count per embryo
## 2. Number of ICM cells per embryo
## 3. Average cell count per litter or experimental group

## This function requires dplyr to run
library("plyr")
library("dplyr")

## This function requires a dataset, as well as the names of the variables for
## * embryo identifier (embryo), 
## * total cell count per embryo (cell.count), 
## * experimental unit (exp.unit) 
## It also takes an option to separate the experimental unit
## by experimental treatment or not - this option allows the user to calculate 
## the median cell count for the entire experimental unit or 
## for each of the treatment groups (median for controls vs each treatment)

do.counts <- function(dataset, embryo.var = 'Embryo_ID', 
                      cellcount.var = 'Cellcount', TI.var = 'TE_ICM', 
                      expunit.var = 'Experiment', treatment.var = 'Treatment', 
                      sep.treatment = T) {
        ## Count the number of cells per embryo and add to main table
        counts <- dataset %>% group_by(!! as.name(embryo.var)) %>% 
                summarize(!! as.name(cellcount.var) = n())
        dataset <- merge(dataset, counts)
        ## Extract unique values of TI.var in dataset
        TI.vals <- as.data.frame(unique(dataset[which(colnames(dataset) == 
                                                              TI.var)]))
        ## Find which one may correspond to trophectoderm
        te <- TI.vals[, 1][which(TI.vals[, 1] %in% c('TE', 'Te', 'te', 'T', 
                                                     'Trophectoderm'))]
        ## Calculate the number of ICM cells per embryo
        tecounts <- dataset %>% filter(!! as.name(TI.var) == te) %>% 
                group_by(!! as.name(embryo.var)) %>% 
                summarize(te.count = n())
        dataset <- merge(dataset, tecounts)
        dataset$icm.count <- dataset[which(colnames(dataset) == 
                                                   cellcount.var)] - 
                dataset$te.count
        dataset$te.count <- NULL
        ## Calculate the median cellcount per experimental grp (litter or else) 
        ## or per treatment group
        if (sep.treatment == T) { 
                med.litter <- dataset %>% 
                        group_by(!! as.name(expunit.var), 
                                 !! as.name(treatment.var)) %>% 
                        summarize(litter.median = 
                                          median(!! as.name(cellcount.var)))
        }
        else { 
                med.litter <- dataset %>% 
                        group_by(!! as.name(expunit.var)) %>% 
                        summarize(litter.median = 
                                          median(!! as.name(cellcount.var)))
        }
        ## Combine with main table and remove avg table
        dataset <- merge(dataset, med.litter)
        return(dataset)
}
