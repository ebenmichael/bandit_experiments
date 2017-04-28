
source("optimization.R")
## Use Random Forests to learn a partition over feature space and
## then use bandit algorithms to choose the best partition

#' Get the partition from a decision tree
#'
#' @param tree Tree data.frame from a random forest
#' @param i The index of the current node
#' @param bounds The upper and lower bounds for each dimension
#'
#' @return A d x 2 x n_leaf_nodes tensor of the upper and lower bounds for each
#'         dimension and each element of the partition
partition_recursive <- function(tree, i, bounds) {
    if(tree[i, 3] == 0) {
        return(bounds)
    } else {
        var <- tree[i, 3]
        value <- tree[i, 4]
        
                                        # move to the left 
        left_bound <- bounds
        left_bound[var, 2] <- value
        left_child <- tree[i, 1]
        from_left <- partition_recursive(tree, left_child, left_bound)
        
                                        # move to the right
        right_bound <- bounds
        right_bound[var, 1] <- value
        right_child <- tree[i, 2]
        from_right <- partition_recursive(tree, right_child, right_bound)

        return(abind::abind(from_left, from_right, along=3))
    }
}

#' Get the partition of the feature space defined by a random forest
#'
#' @param data n x d sampled points
#' @param targets n sampled function values
#' @param bounds d x 2 array of upper and lower bounds for each dimension
#' @param n_tree The number of trees to fit the random forest with
#' @param max_nodes The maximum number of leaves a tree can have, limits grid size
#'
#' @return A d x 2 x n_leaf_nodes tensor of the upper and lower bounds for each
#'         dimension and each element of the partition
rf_partition <- function(data, targets, bounds, n_tree, max_nodes=NULL) {

    # fit a random forest with ntree trees
    rf <- randomForest::randomForest(data, targets, ntree=n_tree,
                                     maxnodes=max_nodes)
    # get the split points for the trees
    splits <- plyr::ldply(1:n_tree,
                          function(x) randomForest::getTree(rf, x))
    # recursively get the partition
    return(partition_recursive(splits, 1, bounds))
}

f_arm <- function(fun, noise_model, bounds) {
    arm <- list(fun=fun, noise_model=noise_model, bounds=bounds)
    class(arm) <- c("f_arm", "arm")
    return(arm)
}

sample.f_arm <- function(arm, n_samples, with_values=FALSE) {
    # sample uniformly in the bounds
    values <- box_runif(n_samples, arm$bounds)
    f_vals <- apply(values, 1,
                    function(x) sample_function(arm$fun,
                                                arm$noise_model,
                                                x, 1))
    
    # sample from the function
    if(! with_values) {
        return(f_vals)
    } else {
        return(list(values, f_vals))
    }

}


f_partition_bandit <- function(fun, noise_model, partition) {

    arms <- apply(partition, 3, function(x) f_arm(fun, noise_model, x))
    class(arms) <- c("f_partition_bandit", "bandit")
    return(arms)
}


pull_arm.f_partition_bandit <- function(bandit, i, n_samples) {
    return(sample(bandit[[i]], n_samples))
}


#' Use a decision tree to partition the space and give the centers of the
#' partitions to sequential halving
tree_sequential_halving <- function(objective, noise_model, bounds,
                                    get_values, budget, max_nodes,
                                    proportion) {
    # query the function using proportion of the budget
    values <- get_values(floor(budget * proportion), bounds)
    targets <- apply(values, 1, function(x) sample_function(objective,
                                                            noise_model,
                                                            x,
                                                            1))
    # fit a decision tree and partition the space
    tree <- randomForest::randomForest(values,
                                       targets,
                                       mtry=dim(values)[2],
                                       n_tree=1,
                                       #replace=FALSE,
                                       #sampsize=dim(values)[1],
                                       maxnodes=max_nodes)
    tree <- randomForest::getTree(tree)
    partition <- partition_recursive(tree, 1, bounds)
    # make a bandit from the mid points of the partitions
    mid_points <- t(apply(partition, 3, rowMeans))
    band <- create_bandit(objective, noise_model, mid_points)
    # run sequential halving with the remaining budget

    best_mid <- sequential_halving(band, floor(budget * (1-proportion)))[[1]]
    return(mid_points[best_mid, ])
}

#' Use a Random forest to sequentially partition the search space and find the
#' arg max of a function
#'
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param get_values A function to propose a set of candidate values
#' @param budget The total number of samples allowed
#' @param eta The proportion to keep at each round
#' @param max_nodes The maximum number of leaf nodes per tree
#' @param n_tree The number of trees in the random forest
sequential_tree <- function(objective, noise_model, bounds, get_values, budget,
                            eta, max_nodes, n_tree) {
    n_arms <- max_nodes
    dimension <- dim(bounds)[1]
    # create a partition with only one element, the whole space
    partition <- array(bounds, dim=c(dim(bounds), 1))
    # keep track of partitions to look at later
    partitions <- list()
    for( r in 1:(floor(logb(n_arms, eta)))) {
        #print(r)
        #print(dim(partition))
        # number of times to pull each disjoint partition
        n_pulls <- floor(budget / (dim(partition)[3] * ceiling(logb(n_arms, eta))))
        #print(n_pulls)
        # for each element in the partiton, sample in that space

        part_band <- f_partition_bandit(objective, noise_model, partition)
        # train a decision tree on the data and get the partitions
        data <- lapply(part_band,
                       function(x) sample(x,
                                          n_pulls,
                                          with_values=TRUE))
        rfs <- lapply(data,
                      function(x)
                          randomForest::randomForest(x[[1]],
                                                     x[[2]],
                                                     maxnodes=max_nodes,
                                                     ntree=n_tree,
                                                     mtry=dimension))
        trees <- lapply(rfs, randomForest::getTree)
        
        tree_partitions <- lapply(1:length(trees),
                                  function(i)
                                      partition_recursive(trees[[i]],
                                                          1,
                                                          partition[,,i]))
        # take the n_arms/4 partitions with the best average
        mid_points <- lapply(tree_partitions,
                             function(x) t(apply(x, 3, rowMeans)))
        preds <- unlist(mapply(predict, rfs, mid_points))
        n_partitions <- length(preds)
        #print(n_partitions)
        if(n_partitions <= eta ^ 2) {
            best_partition <- which.max(preds)
            max_pred <- preds[best_partition]
            #print(best_partition)
            #print(partition)
            partition <- array(do.call(abind::abind,
                                       tree_partitions)[,, best_partition],
                               dim=c(dim(partition)[1:2], 1))
        } else {
            top <- sort(preds,
                        partial=n_partitions -
                            floor(n_partitions/(eta ^ 2)) +
                            1)[n_partitions - floor(n_partitions / (eta ^ 2)) + 1]
            #print(top)
            keep_bool <- preds >= top
                                        #print(keep_bool)
            partition <- do.call(abind::abind, tree_partitions)[,, keep_bool]
            if(length(dim(partition)) == 2) {
                partition <- array(partition, dim=c(dim(partition), 1))
            }
        }
        partitions[[r]] <- partition
        #print(partition)
        # set max_nodes to 2, always dividing the current partitions into eta
        max_nodes = eta
        #print("----")
    }
    return(list(rowMeans(partition[,,1]),
                max_pred,
                partition[,,1],
                partitions))
}


#' Use hyperband with the sequential tree algorithm
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param get_values A function to propose a set of candidate values
#' @param limit The maximum amount of pulls for any parameter setting
#' @param eta The proportion of configurations discarded at each round
#' @param n_tree The number of trees in the random forest
#'
#' @return The best parameter setting found and the number of times each
#'         setting was tried
hypertree <- function(objective, noise_model, bounds, get_values, limit, eta,
                      n_tree) {
    # set up
    s_max <-floor(logb(limit, eta))
    B <- (s_max +1) * limit
    #print((s_max + 1)^2 * limit)
    best_val <- -Inf
    best_arg_max <- NA
    for(s in s_max:0) {
        # get the number of arms and the number of samples for sequential tree
        n_values <- ceiling( B / limit / (s+1) * eta^s)
        resources <- floor(limit * eta^(-s))
        budget <- B
        #print(c(n_values, resources, budget))        
        # run sequential tree for this setting of n and r
        output <- sequential_tree(objective, noise_model, bounds, get_values,
                                  budget, eta, max_nodes=n_values, n_tree)
        arg_max <- output[[1]]
        value <- output[[2]]
        if(value >= best_val) {
            best_val <- value
            best_arg_max <- arg_max
        }
        #print(c(best_arg_max, best_val))
    }
    return(best_arg_max)
    
}


#' Use a Random forest to sequentially partition the search space and find the
#' arg max of a function
#'
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param get_values A function to propose a set of candidate values
#' @param budget The total number of samples allowed
#' @param eta The proportion to keep at each round
#' @param rounds The number of times to make the partition finer
#' @param n_tree The number of trees in the random forest
partition_tree <- function(objective, noise_model, bounds, get_values, budget,
                            eta, rounds, n_tree) {
    max_nodes <- eta^2
    dimension <- dim(bounds)[1]
    # create a partition with only one element, the whole space
    partition <- array(bounds, dim=c(dim(bounds), 1))
    for( r in 1:rounds) {
        #print(r)
        #print(dim(partition))
        # number of times to pull each disjoint partition
        n_pulls <- floor(budget / (dim(partition)[3] * rounds))
        #print(n_pulls)
        # for each element in the partiton, sample in that space
        part_band <- f_partition_bandit(objective, noise_model, partition)
        # train a decision tree on the data and get the partitions
        data <- lapply(part_band,
                       function(x) sample(x,
                                          n_pulls,
                                          with_values=TRUE))
        rfs <- lapply(data,
                      function(x)
                          randomForest::randomForest(x[[1]],
                                                     x[[2]],
                                                     maxnodes=max_nodes,
                                                     ntree=n_tree,
                                                     mtry=dimension))
        trees <- lapply(rfs, randomForest::getTree)
        
        tree_partitions <- lapply(1:length(trees),
                                  function(i)
                                      partition_recursive(trees[[i]],
                                                          1,
                                                          partition[,,i]))
        # take the n_arms/eta partitions with the best average
        mid_points <- lapply(tree_partitions,
                             function(x) t(apply(x, 3, rowMeans)))
        preds <- unlist(mapply(predict, rfs, mid_points))
        n_partitions <- length(preds)
        #print(n_partitions)

        top <- sort(preds,
                    partial=n_partitions -
                        floor(n_partitions/(eta)) +
                        1)[n_partitions - floor(n_partitions / (eta)) + 1]

        keep_bool <- preds >= top
        #print(keep_bool)
        partition <- do.call(abind::abind, tree_partitions)
        #print(partition)        
        best_partition_idx <- which.max(preds)
        best_partition <- partition[,, best_partition_idx]
        partition <- partition[,, keep_bool]
        #print("-")
        #print(partition)
        #print(sum(apply(apply(partition,1, diff), 1, prod)))
        # if the volume gets too small, return the middle of the box
        volume <- sum(apply(apply(partition,1, diff), 1, prod))
        if(volume < sqrt(.Machine$double.eps)) {
            return(rowMeans(best_partition))
        }
        # set max_nodes to 2, always dividing the current partitions into eta
        max_nodes <- eta
        #print("----")
    }
    return(rowMeans(best_partition))
}
