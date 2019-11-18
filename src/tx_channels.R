# This function transforms fluorescence values in a given dataset
# using a regression model previously calculated.
# The function takes a number of arguments (whith no defaults):
# * dataset: data with the channel to be transformed
# * what.subset: a vector of experiments which defines 
#   the subset of data to be transformed
# * what.model: the regression model to be used for the linear transformation
# * input.ch: the channel to be transformed, either as a number (1-5) or
#   in the form of a string "CHx", with x being the number of the channel (1-5)
# * end.ch: the channel to which we want to make input.ch equivalent,
#   in the same format as input.ch
#   (for instance, if we want to transform CH2 into CH3-equivalent values, 
#   input.ch is CH2, and end.ch is CH3)

# Check if model has already been loaded and if not, 
# run ab-test.R, which generates the models
model.exsts <- exists('ng.model')
if(model.exsts == F) { 
        source('./src/ab-test.R')
}
rm(model.exsts)


tx.channel <- function(dataset, what.subset, what.model, 
                        input.ch, end.ch) { 
        # Format input and end channel as "CHx"
        if (is.numeric(input.ch) == TRUE) {
                input.ch <- paste('CH', input.ch, sep = '')
        }
        if (is.numeric(end.ch) == TRUE) {
                end.ch <- paste('CH', end.ch, sep = '')
        }
        # Append the desired termination to the channels to match data tables
        input.ch <- paste(input.ch, '.ebLogCor', sep = '')
        intact.ch <- paste(end.ch, '.ebLogCor', sep = '')
        output.ch <- paste(end.ch, '.ebLogCor.x', sep = '')
        
        # Split dataset in data to be transformed and data not to be trasnformed
        tx <- subset(dataset, Experiment %in% what.subset)
        notx <- subset(dataset, !Experiment %in% what.subset)
        
        # Transform channel in tx using intercept and slope from given model
        # (what.model)
        tx.col <- (tx[input.ch] - what.model$coefficients[1]) / 
                what.model$coefficients[2]
        tx.col <- data.frame(tx.col)
        colnames(tx.col) <- output.ch
        
        # For values in notx, duplicate channel appending .x
        notx.col <- notx[intact.ch]
        notx.col <- data.frame(notx.col)
        colnames(notx.col) <- output.ch
        
        # Combine each transformed variable with corresponding data subset
        tx <- cbind(tx, tx.col)
        notx <- cbind(notx, notx.col)
        
        # Combine both subsets and return the transformed dataset
        dataset <- rbind(tx, notx)
        return(dataset)
}