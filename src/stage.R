# This function takes in a dataset and stages categorizes entries 
# (cells and embryos) according to their total cell count. 
# It defaults to the bins we typically use, namely 8-16, 16-32, 32-64, 64-90, 
# 90-120, 120-150, >150 cells, as in our previous papers
# but it can take in any integer indicating bin size, in number of cells.

# Optionally, it can also take in the name of the variable containing the
# total number of cells - defaults to 'Cellcount'.

stage <- function(dataset, count.var = 'Cellcount',
                  bin = NULL){
        dataset$Cellcount <- get(count.var, dataset)
        # 1. If no intervals are given, default to our standard
        if(missing(bin)) {
                cuts <- c(1, 8, 16, 32, 64, 90, 120, 150, Inf)
                ## 2. Label the stages based on default intervals
                stages <- c('<8', '8_16', '16_32', '32_64', '64_90', 
                            '90_120', '120_150', '>150')
                ## 3. Apply the cut function for the defined intervals 
                ## leaving intervals open on the right (excludes upper limit)
                dataset$Stage <- cut(dataset$Cellcount, breaks = cuts, 
                                     labels = stages, right = F)
                # Convert 'Stage' into a factor with the levels ordered
                # in increasing number of cells
                dataset$Stage <- factor(dataset$Stage, 
                                        levels = c('<8', '8_16', '16_32', 
                                                   '32_64', '64_90', 
                                                   '90_120', '120_150', 
                                                   '>150'))
        }
        ## If a bin size is given, use to break into bins and stage accordingly
        ## label as intervals and leave open on the right
        else {
                cuts <- c(seq(1, 250, by = bin), Inf)
                dataset$Stage2 <- cut(dataset$Cellcount, breaks = cuts, 
                                      right = F)
        }
        return(dataset)
}