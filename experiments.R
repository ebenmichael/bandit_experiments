source("bandits.R")
source("optimization.R")
source("random_forest_opt.R")
library(ggplot2)
library(pryr)
library(dplyr)

## Run an experiment with a given bandit
run_experiment <- function(bandit, best_arm, n_trials, n_min) {
   algo_type <- rep("", 4 * n_trials)
    n_samples <- integer(4 * n_trials)
    correct <- integer(4 * n_trials)
    for(i in 1:n_trials) {
        print("se")
        algo_type[i] <- "successive_elimination"
        output <- successive_elimination(bandit, 0.1, n_min)
        n_samples[i] <- sum(output$arm_pull_mat[,2])
        correct[i] <- 1 * (output$best_idx == best_arm)

        print("lil ucb")
        algo_type[i + n_trials] <- "lil_ucb"
        output <- lil_ucb(bandit, 0.1, 0.01, 9, 1, n_min)
        n_samples[i + n_trials] <- sum(output$arm_pull_mat[,2])
        correct[i + n_trials] <- 1 * (output$best_idx == best_arm)

        print("lil lucb")
        algo_type[i + 2 * n_trials] <- "lil_lucb"
        output <- lil_lucb(bandit, 0.1, 0.01, n_min)
        n_samples[i + 2 * n_trials] <- sum(output$arm_pull_mat[,2])
        correct[i + 2 * n_trials] <- 1 * (output$best_idx == best_arm)

        print("batch racing")
        algo_type[i + 3 * n_trials] <- "batch_racing"
        output <- batch_racing(bandit, 0.1, 1, n_min, n_min)
        n_samples[i + 3 * n_trials] <- sum(output$arm_pull_mat[,2])
        correct[i + 3 * n_trials] <- 1 * (output$best_idxs == best_arm)        
    }
    return(data.frame("algo_type"=algo_type, "n_samples"=n_samples, "correct"=correct))
}

## Run the "sparse" experiment from Jamieson 2013
sparse_experiment <- function(n_arms, n_trials, n_min) {
    sparse_bandit <- normal_bandit(c(1/4, rep(0, n_arms - 1)), rep(1/4,n_arms))
    return(run_experiment(sparse_bandit, 1, n_trials, n_min))
}


## Run the "alpha" experiment from Jamieson 2013
alpha_experiment <- function(alpha, n_arms, n_trials, n_min) {
    alpha_bandit <- normal_bandit(1 - (1:n_arms)^alpha, rep(1/4,n_arms))
    return(run_experiment(alpha_bandit, 1, n_trials, n_min))
}


####Optimization experiments

## Branin function for testing optimization
branin <- function(x) {
    a <- 1
    b <- 5.1 / (4 * pi^2)
    c <- 5 / pi
    r <- 6
    s <- 10
    t <- 1 / (8 * pi)
    return(1 * (x[2] - b * x[1]^2 + c * x[1] - r)^2 + s * (1-t) * cos(x[1]) + s)}


#' Try to optimize a function many times at different levels of resource use
#'
#' @param resources Vector of resources levels to try
#' @param n_per_resources Number of times to optimize at each resource level
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param true_max The true maximum to compare the error to
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param algo The algorithm to use which only takes in (
#'             objective, noise_model, resources)
#' @param algo_name Name of the algorithm
run_opt_exp <- function(resources, n_per_resource, objective,
                        noise_model, bounds, true_max, algo,
                        algo_name) {
    print(algo_name)
    partial_algo <- pryr::partial(algo,
                                  objective=objective,
                                  noise_model=noise_model,
                                  bounds=bounds)
    
    results <- vapply(resources,
                      function(x) replicate(n_per_resource,
                                            true_max -
                                            objective(partial_algo(x))),
                      double(n_per_resource))
    results <- data.frame(t(results))
    results$n_resources = resources
    results$algorithm = algo_name
    return(results)
}

#' Run run_opt_exp for multiple algorithms
#'
#' @param resources Vector of resources levels to try
#' @param n_per_resources Number of times to optimize at each resource level
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_max The true maximum to compare the error to
#' @param algos Names of the algorithms to use which only take in (
#'             objective, noise_model, resources)
run_mult_opt_exp <- function(resources, n_per_resource, objective,
                             noise_model, bounds, true_max, algos) {
    return(plyr::ldply(algos,
          function(a) run_opt_exp(resources,
                                  n_per_resource,
                                  objective,
                                  noise_model,
                                  bounds,
                                  true_max,
                                  match.fun(a),
                                  a)))
}

summary_exp <- function(results) {
    melted_results <- reshape2::melt(results,
                                     id.vars=c("n_resources", "algorithm"))
    melted_results$log.val <- log10(melted_results$value)
    return(summarize_results(melted_results))
}

summarize_results <- function(melted_results) {
    return(plyr::ddply(melted_results, c("n_resources","algorithm"),
                       summarize,
                       N = length(!is.na(value)),
                       mean = mean(value[!is.na(value)]),                       
                       median = median(value[!is.na(value)]),
                       sd = sd(value[!is.na(value)]),
                       se = sd / sqrt(N),
                       log.mean = mean(log.val),
                       log.sd = sd(log.val),
                       log.se = log.sd / N))
}

plot_summary_exp <- function(results) {
    summary_res <- summary_exp(results)
    plot_summarized_results(summary_res)
}

plot_summarized_results <- function(summary_res) {
    return(ggplot(summary_res, aes(x=log10(n_resources), y=median, color=algorithm)) +
        geom_line() +
        geom_errorbar(aes(ymin=mean-1.96 * se,
                          ymax=mean+1.96*se),
                      width=.1) +
        geom_point() +
        scale_y_log10() + 
        xlab("Total Number of Samples (Log Scale)") +
        ylab("Function Error (Log Scale)"))

}

## A collection of partially applied functions to use with run_opt_exp

## Sequential halving random search optimization where limit / n_samples = n_values
seq_halving_n_rand  <- function(objective, noise_model, bounds,
                                limit, n_values) {
    b_max <- bandit_opt(objective,
               noise_model,
               function(y) box_runif(y, bounds),
               n_values,
               limit,
               sequential_halving)
    return(b_max[[1]])
}

## Sequential halving random search optimization where limit / n_samples = 100
seq_halving_100_rand  <- function(objective, noise_model, bounds, limit) {
    n_values <- floor(limit / 100)
    b_max <- seq_halving_n_rand(objective, noise_model, bounds, limit, n_values)
    return(b_max)
}

## Sequential halving random search optimization where we have full exploration
seq_halving_max_rand  <- function(objective, noise_model, bounds, limit) {
    n_values <-  floor(exp(emdbook::lambertW_base(limit * log(2)) - .1))
    b_max <- seq_halving_n_rand(objective, noise_model, bounds, limit, n_values)
    return(b_max)
}

## Sequential halving random search optimization with less exploration
seq_halving_less_rand  <- function(objective, noise_model, bounds, limit) {
    n_values <- floor(exp(emdbook::lambertW_base(limit/10 * log(2))))
    b_max <- seq_halving_n_rand(objective, noise_model, bounds, limit, n_values)
    return(b_max)
}

## Bayesian optimization with fixed number of experiments
bayes_opt_fixed <- function(objective, noise_model, bounds, limit) {
    n_values <- 50
    n_samples <- limit / n_values
    return(bayes_opt(objective, noise_model, n_samples, n_values, bounds))
}

## Bayesian optimization with a growing number of experiments and samples
bayes_opt_growing <- function(objective, noise_model, bounds, limit) {
    n_samples <- floor(20 * sqrt(limit / 20))
    n_values <- floor(n_samples / 20)
    return(bayes_opt(objective, noise_model, n_samples, n_values, bounds))
}

## Sequential Tree with a certain number of trees
seq_tree_eta_n_tree <- function(objective, noise_model, bounds, limit, eta,
                                n_tree) {
    max_nodes <- floor(exp(emdbook::lambertW_base(eta * limit * log(eta)) - .1))
    return(sequential_tree(objective,
                           noise_model,
                           bounds,
                           box_runif,
                           limit,
                           eta,
                           max_nodes,
                           n_tree))
}

seq_tree_1_tree <- function(objective, noise_model, bounds, limit) {
    return(seq_tree_eta_n_tree(objective, noise_model, bounds, limit, 2,
                               1)[[1]])
}

seq_tree_10_tree <- function(objective, noise_model, bounds, limit) {
    return(seq_tree_eta_n_tree(objective, noise_model, bounds, limit, 2,
                               10)[[1]])
}

seq_tree_3 <- function(objective, noise_model, bounds, limit) {
    return(seq_tree_eta_n_tree(objective, noise_model, bounds, limit, 3,
                               1)[[1]])
}

seq_tree_4 <- function(objective, noise_model, bounds, limit) {
    return(seq_tree_eta_n_tree(objective, noise_model, bounds, limit, 4,
                               1)[[1]])
}

hyperband_eta <- function(objective, noise_model, bounds, limit, eta) {
    # solve for the number of resources to give the required total budget
    r = floor(1 /eta * exp(2 * emdbook::lambertW_base(sqrt(eta) * sqrt(limit) * log(eta) / 2) +
                           0.1))
    return(hyperband(objective, noise_model, function(x) box_runif(x, bounds),
                     r, eta)[[1]])
}

hyperband_3 <- function(objective, noise_model, bounds, limit) {
    return(hyperband_eta(objective, noise_model, bounds, limit, 3))
}

hyperband_4 <- function(objective, noise_model, bounds, limit) {
    return(hyperband_eta(objective, noise_model, bounds, limit, 4))
}

hypertree_eta <- function(objective, noise_model, bounds, limit, eta) {
    # solve for the number of resources to give the required total budget
    r = floor(1 /eta * exp(2 * emdbook::lambertW_base(sqrt(eta) * sqrt(limit) * log(eta) / 2) +
                           .1))
    return(hypertree(objective, noise_model, bounds, box_runif, r, eta, 1))
}

hypertree_3 <- function(objective, noise_model, bounds, limit) {
    return(hypertree_eta(objective, noise_model, bounds, limit, 3))
}

hypertree_4 <- function(objective, noise_model, bounds, limit) {
    return(hypertree_eta(objective, noise_model, bounds, limit, 4))
}


## Partition tree with fixed number of rounds
partition_tree_fixed_eta <- function(objective, noise_model, bounds, limit, eta) {
    rounds <- 20
    return(partition_tree(objective,
                          "gaussian",
                          bounds,
                          box_runif,
                          limit,
                          eta,
                          rounds,
                          1))
}

partition_tree_fixed_2 <- function(objective, noise_model, bounds, limit) {
    return(partition_tree_fixed_eta(objective,
                          "gaussian",
                          bounds,
                          limit,
                          2))
}

partition_tree_fixed_3 <- function(objective, noise_model, bounds, limit) {
    return(partition_tree_fixed_eta(objective,
                          "gaussian",
                          bounds,
                          limit,
                          3))
}

## Partition tree with a growing number of rounds
partition_tree_growing_eta <- function(objective, noise_model, bounds, limit, eta) {
    rounds <- floor(sqrt(limit / 20))
    return(partition_tree(objective,
                          "gaussian",
                          bounds,
                          box_runif,
                          limit,
                          eta,
                          rounds,
                          1))
}

partition_tree_growing_2 <- function(objective, noise_model, bounds, limit) {
    return(partition_tree_growing_eta(objective,
                          "gaussian",
                          bounds,
                          limit,
                          2))
}

partition_tree_growing_3 <- function(objective, noise_model, bounds, limit) {
    return(partition_tree_growing_eta(objective,
                          "gaussian",
                          bounds,
                          limit,
                          3))
}


tree_seq_halving_prop <- function(objective, noise_model, bounds, limit,
                                  proportion) {
    max_nodes <-  floor(exp(emdbook::lambertW_base(limit *
                                                   log(2) *
                                                   (1 - proportion)) - .1))
    print(max_nodes)
    return(tree_sequential_halving(objective,
                                   noise_model,
                                   bounds,
                                   box_runif,
                                   limit,
                                   max_nodes,
                                   proportion))
                                   
}

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
#' @param n_tree Number of trees in forest
arg_max_in_partition_exp <- function(resources, n_per_resource, objective,
                                     noise_model, bounds, true_arg_max, eta,
                                     n_tree) {

    print(c(eta, n_tree))
    partial_algo <- pryr::partial(seq_tree_eta_n_tree,
                                  objective=objective,
                                  noise_model=noise_model,
                                  bounds=bounds,
                                  eta=eta,
                                  n_tree=n_tree)
    
    results <- vapply(resources,
                      function(x) rowMeans(replicate(n_per_resource,
                                            contains_arg_max(true_arg_max,
                                                             partial_algo(x)[[3]]))),
                      double(2))
    results <- data.frame(t(results))
    names(results) <-  c("percent.true", "max.L.inf")
    results$n_resources <- resources
    results$n_tree <- n_tree
    results$eta <- eta
    return(results)
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
#' @param n_trees Vector of number of trees in forest
arg_max_in_partition_mult_exp <- function(resources, n_per_resource, objective,
                                     noise_model, bounds, true_arg_max, etas,
                                     n_trees) {

    grid <- expand.grid(etas, n_trees)
    return(t(mapply(
        function(eta, n_tree) arg_max_in_partition_exp(resources,
                                               n_per_resource,
                                               objective,
                                               noise_model,
                                               bounds,
                                               true_arg_max,
                                               eta,
                                               n_tree),
        grid$Var1,
        grid$Var2)))
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
