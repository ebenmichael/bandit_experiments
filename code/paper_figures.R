###Code to create the figures for the paper and poster
source("experiments.R")
source("random_forest_opt.R")
source("optimization.R")
library(ggplot2)
library(gridExtra)
### functions and bounds for testing

## Branin function for testing optimization
branin <- function(x) {
    a <- 1
    b <- 5.1 / (4 * pi^2)
    c <- 5 / pi
    r <- 6
    s <- 10
    t <- 1 / (8 * pi)
    return(1 * (x[2] - b * x[1]^2 + c * x[1] - r)^2 + s * (1-t) * cos(x[1]) + s)}

neg_branin <- function(x) -branin(x)
bran_bd <- rbind(c(-5,10), c(0,15))
bran_arg_max <- rbind(c(-pi, 12.275), c(pi, 2.275), c(9.42478, 2.475))



## Hartmann 3d function
hart3 <- function(xx)
{  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(3.0, 10, 30,
         0.1, 10, 35,
         3.0, 10, 30,
         0.1, 10, 35)
  A <- matrix(A, 4, 3, byrow=TRUE)
  P <- 10^(-4) * c(3689, 1170, 2673,
                   4699, 4387, 7470,
                   1091, 8732, 5547,
                   381, 5743, 8828)
  P <- matrix(P, 4, 3, byrow=TRUE)
	
  xxmat <- matrix(rep(xx,times=4), 4, 3, byrow=TRUE)
  inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
  outer <- sum(alpha * exp(-inner))
	
  y <- -outer
  return(y)
}
neg_hart3 <- function(x) -hart3(x)
hart3_bd <- rbind(c(0,1), c(0,1), c(0,1))
hart3_arg_max <- c(0.114614, 0.555649, 0.852547)

### More fine grained experiments for sequential tree

#' Count the proportion of times that the final partition in sequential
#' tree contains an arg max
#' 
#' @param resources Vector of resources levels to try
#' @param n_per_resources Number of times to optimize at each resource level
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_arg_max A matrix of true arg maxes to compare to
#' @param eta Amount to cut the number of boxes by
#' @param max_nodes The maximum number of leaf nodes per tree
#' @param n_tree Number of trees in forest
arg_max_in_partition_exp <- function(resources, n_per_resource, objective,
                                     noise_model, bounds, true_arg_max, eta,
                                     max_nodes, n_tree) {

    print(c(eta, n_tree, max_nodes))
    partial_algo <- pryr::partial(sequential_tree,
                                  objective=objective,
                                  noise_model=noise_model,
                                  bounds=bounds,
                                  eta=eta,
                                  max_nodes=max_nodes,
                                  n_tree=n_tree)
    
    results <- vapply(resources,
                      function(x)
                          rowMeans(
                              replicate(n_per_resource,
                                        contains_arg_max(
                                            true_arg_max,
                                            sequential_tree(objective,
                                                            noise_model,
                                                            bounds,
                                                            box_runif,
                                                            x,
                                                            eta,
                                                            max_nodes,
                                                            n_tree)[[3]]))),
                      double(2))
    #results <- data.frame(t(results))
    #names(results) <-  c("percent.true", "max.side.length")
    #results$n_resources <- resources
    #results$n_tree <- n_tree
    #results$eta <- eta
    #results$max_nodes <- max_nodes
    return(t(results))
}

contains_arg_max <- function(true_arg_max, box) {
    # get the maximum side length of the box
    max_len <- max(box[,2] - box[,1])
    return(c( 1* any(sapply(1:dim(true_arg_max)[1] ,
                      function(i)
                          any(sapply(1:dim(true_arg_max)[2],
                                     function(j)
                                         true_arg_max[i,j] >= box[j,1] &&
                                         true_arg_max[i,j] <= box[j,2])))),
             max_len))
}

#' Count the proportion of times that the final partition in sequential
#' tree contains an arg max
#' 
#' @param resources Vector of resources levels to try
#' @param n_per_resources Number of times to optimize at each resource level
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_arg_max A matrix of true arg maxes to compare to
#' @param etas Vector of amounts to cut the number of boxes by
#' @param max_nodes Vector of the maximum number of leaf nodes per tree
#' @param n_trees Vector of number of trees in forest
arg_max_in_partition_mult_exp <- function(resources, n_per_resource, objective,
                                     noise_model, bounds, true_arg_max, etas,
                                     max_nodes, n_trees) {

    grid <- expand.grid(etas, max_nodes, n_trees)
    dat <- as.data.frame(t(mapply(
        function(eta, max_node, n_tree) arg_max_in_partition_exp(resources,
                                               n_per_resource,
                                               objective,
                                               noise_model,
                                               bounds,
                                               true_arg_max,
                                               eta,
                                               max_node,
                                               n_tree),
        grid$Var1,
        grid$Var2,
        grid$Var3)))
    names(dat) <- c("percent.true", "max.side.length")
    dat$eta <- grid$Var1
    dat$max_nodes <- grid$Var2
    dat$n_tree <- grid$Var3
    return(dat)
}

#' Count the proportion of times that seq tree doesn't throw away the arg max
#' 
#' @param resources Vector of resources levels to try
#' @param n_per_resources Number of times to optimize at each resource level
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_arg_max A matrix of true arg maxes to compare to
#' @param eta Amount to cut the number of boxes by
#' @param n_tree Number of trees in forest
correct_per_round_exp <- function(resources, n_per_resource, objective,
                                     noise_model, bounds, true_arg_max, eta,
                                     n_tree) {

    print(c(eta, n_tree))
    partial_algo <- pryr::partial(seq_tree_eta_n_tree,
                                  objective=objective,
                                  noise_model=noise_model,
                                  bounds=bounds,
                                  eta=eta,
                                  n_tree=n_tree)
    
    results <- sapply(resources,
                      function(x)
                          rowMeans(
                              replicate(n_per_resource,
                                        sapply(partial_algo(x)[[4]],
                                               function(part)
                                                   max(
                                                       apply(
                                                           part,
                                                           3,
                                                           function(box)
                                                               contains_arg_max(
                                                                   true_arg_max,
                                                                   box)[1]))))))
    if(is.list(results)) {
        results <- plyr::ldply(1:length(resources),
                               function(i)
                                   data.frame(pct.contain=results[[i]],
                                              round=1:length(results[[i]]),
                                              n_resources=resources[i],
                                              eta=eta,
                                              n_tree=n_tree))
    }
    else {
        results <- data.frame(pct.contain = results,
                              round = 1:dim(results)[1],
                              n_resources = resources[1],
                              eta = eta,
                              n_tree = n_tree)
    }
    return(results)
}

#' Count the proportion of times that seq tree doesn't throw away the arg max,
#' for various etas and number of trees
#' 
#' @param resources Vector of resources levels to try
#' @param n_per_resources Number of times to optimize at each resource level
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_arg_max A matrix of true arg maxes to compare to
#' @param etas Vector of amount to cut the number of boxes by
#' @param n_trees Vector of number of trees in forest
correct_per_round_mult_exp <- function(resources, n_per_resource, objective,
                                       noise_model, bounds, true_arg_max, etas,
                                       n_trees) {
    grid <- expand.grid(etas, n_trees)
    return(plyr::ldply(
        1:dim(grid)[1],
        function(i)
            correct_per_round_exp(resources,
                                  n_per_resource,
                                  objective,
                                  noise_model,
                                  bounds,
                                  true_arg_max,
                                  grid[i,1],
                                  grid[i,2])))

}


### Empirical results and plots

seq_tree_arg_max <- function(budget, n_exps) {
    max_max_nodes <- floor(exp(emdbook::lambertW_base(2 * budget * log(2)) - .1))
    res <- arg_max_in_partition_mult_exp(c(budget),
                                         50,
                                         neg_branin,
                                         "gaussian",
                                         bran_bd,
                                         bran_arg_max,
                                         c(2,3),
                                         round(seq(10,max_max_nodes,
                                                   length.out=n_exps)),
                                         c(1))

    # make a plot
    melt_res <- reshape2::melt(res[,c("percent.true", "max.side.length", "eta",
                                      "max_nodes")], id.vars=c("eta","max_nodes"))
    plt <- ggplot(melt_res, aes(x=max_nodes, y=value, color=as.factor(eta))) +
        geom_point() + geom_line(size=1.5) +
        facet_grid(variable ~ . , scales="free_y",
                   labeller = as_labeller(
                       c("percent.true" = "Proportion Containing Arg Max",
                         "max.side.length" = "Average Maximum Side Length of Box"))) +
        theme_minimal() +
        xlab("Number of Elements in Initial Partition") +
        ylab("") +
        scale_color_discrete("Eta")
                       
    return(res)
}
