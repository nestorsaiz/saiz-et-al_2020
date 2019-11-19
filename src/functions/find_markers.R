find.markers <- function(dataset) { 
        unicos <- matrix(0, nrow = length(unique(dataset$Experiment)), 
                         ncol = 4, dimnames = list(NULL, c('Experiment', 
                                                           'CH2', 'CH3', 
                                                           'CH5')))
        unicos[, 1] <- unique(dataset$Experiment)
        for(j in 1:length(unicos[, 1])) {
                e <- unicos[j, 1]
                unicos[j, 2] <- unique(subset(dataset,
                                              Experiment == e &
                                                      Channel == 'CH2')$Marker)
                unicos[j, 3] <- unique(subset(dataset,
                                              Experiment == e &
                                                      Channel == 'CH3')$Marker)
                unicos[j, 4] <- unique(subset(dataset,
                                              Experiment == e &
                                                      Channel == 'CH5')$Marker)
        }
        unicos <<- data.frame(unicos)
}