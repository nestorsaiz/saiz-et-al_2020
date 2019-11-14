view.tx <- function(dataset) { 
        ## Visualize the outcome of the transformation
        ## Create three groups of experiments based on staining
        aa <<- unicos$Experiment[which(unicos$CH2 == 'NANOG.rat' & 
                                               unicos$CH5 == 'GATA6.gt')]
        bb <<- unicos$Experiment[which(unicos$CH2 == 'NANOG.rat' & 
                                               unicos$CH3 == 'GATA6.rb')]
        cc <<- unicos$Experiment[which(unicos$CH3 == 'NANOG.rb' & 
                                               unicos$CH5 == 'GATA6.gt')]
        
        ## Plot GATA6 vs NANOG values (log, corrected, transformed)
        ## and color code for each subset of experiments (aa, bb, cc)
        source('plotting-aes.R')
        
        ## Before transformation 
        before.x <- ggplot(data = subset(dataset, Experiment %in% cc & 
                                                 TE_ICM == 'ICM' & Cellcount > 15),
                           aes(x = CH5.ebLogCor, 
                               y = CH3.ebLogCor))
        before.x <- before.x + geom_jitter(color = 'gray70', size = 1)
        before.x <- before.x + geom_jitter(data = subset(dataset, 
                                                         Experiment %in% bb & 
                                                                 TE_ICM == 'ICM'),
                                           aes(x = CH3.ebLogCor, 
                                               y = CH2.ebLogCor), color = 'blue', 
                                           size = 1)
        before.x <- before.x + geom_jitter(data = subset(dataset, 
                                                         Experiment %in% aa & 
                                                                 TE_ICM == 'ICM'),
                                           aes(x = CH5.ebLogCor, 
                                               y = CH2.ebLogCor), color = 'green', 
                                           size = 1)
        before.x <- before.x + facet_grid(. ~ Stage) + looks +
                theme(aspect.ratio = 1)
        before.x <- before.x + coord_fixed(1) + ylim(3, 9)
        print(before.x)
        
        ## and after transformation
        after.x <- ggplot(data = subset(dataset, Experiment %in% cc & 
                                                TE_ICM == 'ICM' & Cellcount > 15),
                          aes(x = CH5.ebLogCor.x, 
                              y = CH3.ebLogCor.x))
        after.x <- after.x + geom_jitter(color = 'gray70', size = 1)
        after.x <- after.x + geom_jitter(data = subset(dataset, 
                                                       Experiment %in% bb & 
                                                               TE_ICM == 'ICM'),
                                         aes(x = CH5.ebLogCor.x, 
                                             y = CH3.ebLogCor.x), color = 'blue', 
                                         size = 1)
        after.x <- after.x + geom_jitter(data = subset(dataset, 
                                                       Experiment %in% aa & 
                                                               TE_ICM == 'ICM'),
                                         aes(x = CH5.ebLogCor.x, 
                                             y = CH3.ebLogCor.x), color = 'green', 
                                         size = 1)
        after.x <- after.x + facet_grid(. ~ Stage) + looks + 
                theme(aspect.ratio = 1)
        print(after.x)
        }
