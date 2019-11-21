# A classifier to assign lineage identity to cells over time
# based on Pdgfra expression levels and manual thresholds

# Rule #1: cells that become PrE or EPI do not change identity 
# for the remainder of the movie

# Rule #2: to become PrE or EPI, they need to express Pdgfra above or below
# the manually set threshold (see movie-analysis script) and remain
# there for at least 2h after (based on area between curve and threshold)

el.classifier <- function(dataset, u = 7) { 
        for(i in 1:length(dataset)) {
                ## Create an index (ndx) for the ascending timeframe
                ndx <- order(dataset[[i]]$timeframe)
                ## and apply the index to order all values of each track
                dataset[[i]] <- dataset[[i]][ndx,]
                ## Check whether the track has more than one value (length > 1)
                shorty <- length(dataset[[i]]$timeframe) < 2
                if(shorty != T) { 
                        for(t in 2:length(dataset[[i]]$timeframe)) {
                                ## Evaluate identity of previous time point
                                previo <- dataset[[i]]$identity.t[t-1]
                                ## Evaluate putative identity of current time point
                                este <- dataset[[i]]$put.id[t]
                                ## Once cells become PrE or EPI, they won't change
                                if(previo == 'PRE') {
                                        dataset[[i]]$identity.t[t] <- 'PRE'
                                } else if(previo == 'EPI') {
                                        dataset[[i]]$identity.t[t] <- 'EPI'
                                } else if(previo == 'DP') {
                                        ## Check if area under curve 
                                        ## for PrE and EPI thresholds is > 0 or < 0
                                        ## for the next 2h or until end of track
                                        fin <- length(dataset[[i]]$timeframe)
                                        a.check <- t + u <= fin 
                                        if(a.check == T) {
                                                pre.auc <- sum(dataset[[i]]$diff.p[t:t+u]) > 0
                                                epi.auc <- sum(dataset[[i]]$diff.e[t:t+u]) < 0
                                        } else {
                                                pre.auc <- sum(dataset[[i]]$diff.p[t:fin]) > 0
                                                epi.auc <- sum(dataset[[i]]$diff.e[t:fin]) < 0
                                        }
                                        
                                        ## Putative DP cells are DP regardless
                                        if(este == 'D') {
                                                dataset[[i]]$identity.t[t] <- 'DP'
                                        } 
                                        ## To become PrE, need to meet these conditions
                                        ## putative PrE + PrE.AUC for next 2h > 0
                                        else if(este == 'P') {
                                                if(pre.auc == T) {
                                                        dataset[[i]]$identity.t[t] <- 'PRE'
                                                        dataset[[i]]$switch[t] <- T
                                                } else {
                                                        dataset[[i]]$identity.t[t] <- 'DP'
                                                }
                                        } 
                                        ## To become EPI, need to meet these conditions
                                        ## putative EPI + EPI.AUC for next 2h < 0
                                        else if(este == 'E') {
                                                if(epi.auc == T) {
                                                        dataset[[i]]$identity.t[t] <- 'EPI'
                                                        dataset[[i]]$switch[t] <- T
                                                } else {
                                                        dataset[[i]]$identity.t[t] <- 'DP'
                                                }
                                        }
                                }
                        }
                }
                
        }
        return(dataset)
}
