## Use Random Forests to learn a partition over feature space and
## then use bandit algorithms to choose the best partition

#' Get the partition of the feature space defined by a random forest
#'
#' @param data n x d sampled points
#' @param targets n sampled function values
#' @param bounds d x 2 array of upper and lower bounds for each dimension
#' @param n_tree The number of trees to fit the random forest with
#' @param maxnodes The maximum number of leaves a tree can have, limits grid size
rf_partition <- function(data, targets, bounds, n_tree, max_nodes) {
    # fit a random forest with ntree trees
    rf <- randomForest::randomForest(data, targets, ntree=n_tree,
                                     maxnodes=max_nodes)
    # get the split points for the trees
    splits <- plyr::ldply(1:n_tree,
                          function(x) randomForest::getTree(rf, x)[,c(3,4)])
    # add in the bounds on the space to splits 
    bounds_melt <- reshape2::melt(bounds)[c(1,3)]
    names(splits) <- c("split.var", "split.point")  
    names(bounds_melt) <-names(splits)
    splits <- rbind(splits, bounds_melt)
    #get the left and right end points for each box in each dimension
    one_dim_bounds <- lapply(with(splits, split(split.point, split.var)),
                             function(x) zoo::rollapply(sort(x), 2, c))[-1]
    # get all combinations of the boxes
    n_boxes_dim <- lapply(one_dim_bounds, function(x) 1:dim(x)[1])
    grid_points <- expand.grid(n_boxes_dim)
    # combine these combination of box edges to partition the space
    new_bounds <- array(Reduce(rbind,
                               mapply(function(x,y) x[y,],
                                      one_dim_bounds,
                                      as.list(grid_points),
                                      SIMPLIFY=F)),
                        dim=c(dim(grid_points),2))
    return(new_bounds)
}
