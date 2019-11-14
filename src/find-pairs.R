## This small function takes in a dataset to be explored
## and a vector of marker gene names (capitalization doesn't matter) 
## to search pairs of.
## It returns a list of 2-element lists, each containing:
## a) the pairs of antibodies for each marker found in the dataset
## b) the names of the experiments where said antibody pairs were used

find.pairs <- function(dataset, markers = c('nanog', 'gata6')) {
        mms <- list()
        for(m in 1:length(markers)) {
                mms[[m]] <- unique(grep(paste(markers[m], '+', sep = ''), 
                                        ignore.case = T, dataset$Marker, 
                                        perl = T, value = T))
        }
        ab.pairs <- list()
        for(m in 1:length(mms)) {
                if(length(mms[[m]]) > 1) {
                        cond1 <- unique(dataset$Experiment[
                                which(dataset$Marker == mms[[m]][1])])
                        cond2 <- unique(dataset$Experiment[
                                which(dataset$Marker == mms[[m]][2])])
                        expts <- cond1[which(cond1 %in% cond2)]
                        ab.pairs[[m]] <- list(mms[[m]], expts)
                }
        }
        ab.pairs <<- ab.pairs
}