# This function implements K-menans to classify ICM cells from a given dataset
# into PrE, EPI, double positive (DP) and double negative (DN) cells

mkay <- function(dataset, miniset, k, ids, 
                 DP = F, DN = F, x.var, y.var, 
                 TI.var = 'TE_ICM') {
  # Set a seed for reproducibility
  set.seed(21)
  
  # Extract the index for the columns of dataset containing the x and y vars
  x <- which(colnames(miniset) == x.var)
  y <- which(colnames(miniset) == y.var)
  
  # Define conditions to test for TE, ICM or ESC identity
  is.te <- dataset[TI.var] == 'TE'
  is.esc <- dataset[TI.var] == 'ESC'
  is.icm <- dataset[TI.var] == 'ICM'
  
  # Create new variable in dataset to hold the lineage identity
  # assigned using k-means clustering (Identity.km)
  dataset$Identity.km <- rep(NA, nrow(dataset))
  
  # Perform K-means using k clusters on the columns x and y of the miniset
  km <- kmeans(miniset[, x:y], k)
  # Extract the centers of the clusters found
  centers <- km$centers
  # If either DP or DN are required (TRUE), define DP and/or DN centers
  # using the maxima or minima of centers, respectively
  if(DP == T) {
    dp <- c(centers[, 1][which.max(centers[, 1])],
            centers[, 2][which.max(centers[, 2])])
    centers <- rbind(centers, dp)
  }
  if(DN == T) {
    dn <- c(centers[, 1][which.min(centers[, 1])],
            centers[, 2][which.min(centers[, 2])])
    centers <- rbind(centers, dn)
  }
  print(centers)
  
  # Make matrix to hold sum of squares (rows = icm rows, k columns)
  ssq <- matrix(0, length(dataset[TI.var][is.icm]), k)
  # and populate with min sum of squares for each cell
  # (min distance to each center for each cell's [GATA6] vs [NANOG] value)
  for(i in 1:k) {
    ssq[,i] <- (dataset[x.var][is.icm] - centers[i,1])^2 + 
      (dataset[y.var][is.icm] - centers[i,2])^2
  }
  # calculate what center each cell is closest to 
  # (i.e. which sum of squares is smallest)
  min.ssq <<- apply(ssq, 1, which.min)
  
  # Assign Identity
  # TE cells and ESCs (if present) remain unchaged
  dataset$Identity.km[is.te] <- 'TE'
  dataset$Identity.km[is.esc] <- 'ESC'
  # Determine ICM population based on min.ssq values
  dataset$Identity.km[is.icm] <- ids[min.ssq]
  dataset$Identity.km <- factor(dataset$Identity.km, 
                                levels = c('TE', 'PRE', 'DP', 'EPI', 'DN'))
  return(dataset)
}

