source("bandits.R")
source("optimization.R")
source("random_forest_opt.R")
library(ggplot2)
library(pryr)
library(dplyr)
####Optimization experiments

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
#' @param obj_name Name of the objective function
run_opt_exp <- function(resources, n_per_resource, objective,
                        noise_model, bounds, true_max, algo,
                        algo_name, obj_name) {
    print(algo_name)
    partial_algo <- pryr::partial(algo,
                                  objective=objective,
                                  noise_model=noise_model,
                                  bounds=bounds)
    
    results <- vapply(resources,
                      function(x) replicate(n_per_resource,
                                            true_max -
                                            objective(partial_algo(x)[[1]])),
                      double(n_per_resource))
    results <- data.frame(t(results))
    if(identical(algo, hyperband_3) || identical(algo, hyperband_4)) {
        results$n_resources <- sapply(resources, function(x) partial_algo(x)[[2]])
    } else {
        results$n_resources <- resources
    }
    results$algorithm <- algo_name
    results$noise_model <- noise_model
    results$objective <- obj_name
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
#' @param obj_name Name of the objective function
run_mult_opt_exp <- function(resources, n_per_resource, objective,
                             noise_model, bounds, true_max, algos,
                             obj_name) {
    return(plyr::ldply(algos,
          function(a) run_opt_exp(resources,
                                  n_per_resource,
                                  objective,
                                  noise_model,
                                  bounds,
                                  true_max,
                                  match.fun(a),
                                  a,
                                  obj_name)))
}

summary_exp <- function(results) {
    melted_results <- reshape2::melt(results,
                                     id.vars=c("n_resources", "algorithm",
                                               "noise_model", "objective"))
    melted_results$log.val <- log10(melted_results$value)
    return(summarize_results(melted_results))
}

summarize_results <- function(melted_results) {
    return(plyr::ddply(melted_results, c("n_resources","algorithm",
                                         "noise_model", "objective"),
                       summarize,
                       N = length(!is.na(value)),
                       mean = mean(value[!is.na(value)]),                       
                       median = median(value[!is.na(value)]),
                       sd = sd(value[!is.na(value)]),
                       se = sd / sqrt(N),
                       log.mean = mean(log.val[!is.na(log.val)]),
                       log.sd = sd(log.val[!is.na(log.val)]),
                       log.se = log.sd / sqrt(N),
                       upper.quantile = quantile(value[!is.na(value)], 3/4),
                       lower.quantile = quantile(value[!is.na(value)], 1/4)))
}

plot_summary_exp <- function(results) {
    summary_res <- summary_exp(results)
    plot_summarized_results(summary_res)
}

plot_summarized_results <- function(summary_res) {
    return(ggplot(summary_res, aes(x=log10(n_resources), y=median, color=algorithm)) +
        geom_line() +
        geom_errorbar(aes(ymin=lower.quantile,
                          ymax=upper.quantile),
                      width=.1) +
        geom_point() +
        scale_y_log10() + 
        xlab("Total Number of Samples (Log Scale)") +
        ylab("Function Error (Log Scale)") +
        facet_grid(objective ~ noise_model))

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
    return(b_max)
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
    return(list(bayes_opt(objective, noise_model, n_samples, n_values, bounds)))
}

## Bayesian optimization with a growing number of experiments and samples
## Use twice the number of rounds as maximum exploration sequential halving
bayes_opt_growing_halving <- function(objective, noise_model, bounds, limit) {
    n_values <- 2 * floor(emdbook::lambertW_base(limit * log(2)) - .1)
    n_samples <- limit / (n_values + 2)
    print(n_values)
    return(list(bayes_opt(objective, noise_model, n_samples,
                          n_values, bounds, n_init=2)))
}

## Bayesian optimization with a growing number of experiments and samples
## Use the number of rounds Hyperband_3 gets
bayes_opt_growing_hyper <- function(objective, noise_model, bounds, limit) {
    s_max <- floor(logb(floor(1 /3 * exp(2 * emdbook::lambertW_base(sqrt(3) *
                                                      sqrt(limit) * log(3) / 2) +
                               .27)), 3))
    n_values <- 0
    for(s in s_max:0) { for(i in 0:s) { n_values <- n_values + 1}}
    n_samples <- limit / n_values
    print(n_values)
    return(list(bayes_opt(objective, noise_model, n_samples, n_values,
                          bounds, n_init=2))    )
}

seq_tree_fixed_prop <- function(objective, noise_model, bounds, limit) {
    max_nodes <- limit / 20
    return(sequential_tree(objective,
                           noise_model,
                           bounds,
                           box_runif,
                           limit,
                           2,
                           max_nodes,
                           1))
    
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
                               1))
}

seq_tree_10_tree <- function(objective, noise_model, bounds, limit) {
    return(seq_tree_eta_n_tree(objective, noise_model, bounds, limit, 2,
                               10))
}

seq_tree_3 <- function(objective, noise_model, bounds, limit) {
    return(seq_tree_eta_n_tree(objective, noise_model, bounds, limit, 3,
                               1))
}

seq_tree_4 <- function(objective, noise_model, bounds, limit) {
    return(seq_tree_eta_n_tree(objective, noise_model, bounds, limit, 4,
                               1))
}

hyperband_eta <- function(objective, noise_model, bounds, limit, eta) {
    # solve for the number of resources to give the required total budget
    r = floor(1 /eta * exp(2 * emdbook::lambertW_base(sqrt(eta) *
                                                      sqrt(limit) * log(eta) / 2) +
                           .27))
    return(hyperband(objective, noise_model, bounds, box_runif,
                     r, eta))
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

