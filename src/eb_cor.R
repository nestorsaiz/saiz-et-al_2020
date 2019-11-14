## This function performs Empirical Bayes correction on any given 
## vector, as performed for Saiz *et al* (2016) *Nature Comms*
## The function requires two arguments: x and channel
## x is the dataset where the vector to correct is located (matrix, data frame)
## channel is the vector itself (fluorescence channel). It can be given as an
## integer (1 - 5) or as a string (the full name of the vector, in quotations)
## The arguments group and logadd are optional and have defaults
## group is a grouping variable to separate embryos by
## logadd is the value to add to all fluorescence measurements to avoid NAs due 
## to zero values (by default logadd = 0.0001)

ebcor <- function(x, channel, group = NULL, logadd = 0.0001) {
    # get unique embryo IDs
    embryos <- unique(x$Embryo_ID)
    # if using integers (1 - 5) for channel, convert to standard names
    if (is.numeric(channel) == TRUE) {
            channels <- c('CH1.Avg', 'CH2.Avg', 'CH3.Avg', 'CH4.Avg', 'CH5.Avg')
            channel <- channels[channel]
    }
    # else, use whatever channel name given
    # data
    x$CH.Avg <- get(channel, x)
    # fitted regression coefficients and their standard errors
    coefs <- matrix(0, length(embryos), 2)
    for(i in 1:length(embryos)) {
        xi <- x[x$Embryo_ID==embryos[i],]
        coefs[i, 1:2] <- summary(lm(log(CH.Avg + logadd) ~ Z + (TE_ICM == "TE"), 
                                    data=xi))$coefficients[2, 1:2]
    }
    # if grouping variable is NULL create a dummy vector of 1s
    if (missing(group)) group <- rep(1, nrow(x))
    # group indicator of each embryo
    egrp <- tapply(group, x$Embryo_ID, function(y) {unique(y)[1]})
    # Emperical Bayes correction across the embryos in a group
    ebcoefs <- rep(0, length(embryos))
    for (i in unique(egrp)) {
        ebcoefs[egrp==i] <- mean(coefs[egrp==i, 1]) + 
                (1 - coefs[egrp==i, 2]^2/(coefs[egrp==i, 2]^2 + 
                                                  var(coefs[egrp==i, 1])))*(coefs[egrp==i, 1] - 
                                                                                    mean(coefs[egrp==i, 1]))
    }
    # EB corrected log signal
    CH.ebLogCor <- rep(NA, nrow(x))
    for(i in 1:length(embryos)) {
        ii <- x$Embryo_ID==embryos[i]
        CH.ebLogCor[ii] <- log(logadd + x$CH.Avg[ii]) - ebcoefs[i]*x$Z[ii]
    }
    # return EB corrected log signal
    CH.ebLogCor
}
