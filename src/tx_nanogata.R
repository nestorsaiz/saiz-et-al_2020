tx.nanogata <- function(dataset) { 
        ## Source script modelling the relationship betweeen 
        ## NANOG.rat vs NANOG.rb and GATA6.rb vs GATA6.gt 
        source('abtest.R')
        
        ## Split dataset in embryos stained for NANOG.rat and embryos which aren't
        rat <- subset(dataset, Experiment %in% ng.rat)
        norat <- subset(dataset, !Experiment %in% ng.rat)
        
        ## Transform CH2.ebLogCor in rat subset using intercept and slope from ng.model
        ## obtained in abtest.R and make equivalent column in norat subset
        ## where the values are the same as CH3.ebLogCor (NANOG.rb).
        ## Call the transformed value CH3.ebLogCor.x, as it will be equivalent to 
        ## CH3.ebLogCor, where NANOG.rb values are found.
        rat$CH3.ebLogCor.x <- (rat$CH2.ebLogCor - ng.model$coefficients[1])/ 
                ng.model$coefficients[2]
        norat$CH3.ebLogCor.x <- norat$CH3.ebLogCor
        ## Combine both subsets into new.lms
        dataset <- rbind(rat, norat)
        rm(rat, norat)
        
        ## Split dataset in embryos stained for GATA6.rb and embryos which aren't
        rb <- subset(dataset, Experiment %in% g6.rb)
        norb <- subset(dataset, !Experiment %in% g6.rb)
        
        ## Transform CH3.ebLogCor in rb subset using intercept and slope from gata.model
        ## obtained in abtest.R and make equivalent column in norb subset
        ## where the values are the same as CH5.ebLogCor
        rb$CH5.ebLogCor.x <- (rb$CH3.ebLogCor - gata.model$coefficients[1])/ 
                gata.model$coefficients[2]
        norb$CH5.ebLogCor.x <- norb$CH5.ebLogCor
        ## Combine both subsets into new.lms
        dataset <- rbind.fill(rb, norb)
        rm(rb, norb)
        
        ## Return the transformed dataset
        return(dataset)
        }