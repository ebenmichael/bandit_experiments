###Code to create the figures for the paper and poster
source("experiments.R")
source("random_forest_opt.R")
source("optimization.R")
library(ggplot2)

### Berkeley color palette
berk_palette <- c("#003262", "#FDB515", "#C4820E", "#3B7EA1",
                  "#D9661F", "#ED4E33", "#00A598", "#46535E",
                  "#859438", "#EE1F60", "#584F29")


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
bran_max <- -0.397887


## Hartmann 3d function, from https://www.sfu.ca/~ssurjano/hart3.html
hart3 <- function(xx) {  
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
hart3_max <- 3.86278

## Hartmann 6d function, from https://www.sfu.ca/~ssurjano/hart6.html
hart6 <- function(xx) {
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(10, 3, 17, 3.5, 1.7, 8,
         0.05, 10, 17, 0.1, 8, 14,
         3, 3.5, 1.7, 10, 17, 8,
         17, 8, 0.05, 10, 0.1, 14)
  A <- matrix(A, 4, 6, byrow=TRUE)
  P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
                   2329, 4135, 8307, 3736, 1004, 9991,
                   2348, 1451, 3522, 2883, 3047, 6650,
                   4047, 8828, 8732, 5743, 1091, 381)
  P <- matrix(P, 4, 6, byrow=TRUE)
  
  xxmat <- matrix(rep(xx,times=4), 4, 6, byrow=TRUE)
  inner <- rowSums(A[,1:6]*(xxmat-P[,1:6])^2)
  outer <- sum(alpha * exp(-inner))
  
  y <- -(2.58 + outer) / 1.94
  return(y)
}
neg_hart6 <- function(x) -hart6(x)
hart6_bd <- rbind(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1))
hart6_arg_max <- c(.20169, .150011, .476874, .275332, 0.311652, 0.6573)
hart6_max <- 3.32237


### Empirical Evaluation

bayes_opt_branin_gaus <- function(n_per_resource) {
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_branin,
                            "gaussian", bran_bd, bran_max,
                            c("bayes_opt_growing_halving",
                              "bayes_opt_growing_hyper"),
                            "negative.branin")
    return(res)
}

bandit_opt_branin_gaus <- function(n_per_resource) {
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_branin,
                            "gaussian", bran_bd, bran_max,
                            c("hyperband_3",
                              "hyperband_4",
                              "seq_halving_max_rand",
                              "seq_halving_100_rand"),
                            "negative.branin")
    return(res)
    
}

tree_opt_branin_gaus <- function(n_per_resource) {
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_branin,
                            "gaussian", bran_bd, bran_max,
                            c("seq_tree_1_tree",
                              "seq_tree_fixed_prop",
                              "partition_tree_growing_4",
                              "partition_tree_fixed_4",
                              "partition_tree_growing_3"),
                            "negative.branin")
    
}

bayes_opt_hart3_gaus <- function(n_per_resource) {
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_hart3,
                            "gaussian", hart3_bd, hart3_max,
                            c("bayes_opt_growing_halving",
                              "bayes_opt_growing_hyper"),
                            "negative.hartmann.3")
    return(res)
}

bandit_opt_hart3_gaus <- function(n_per_resource) {
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_hart3,
                            "gaussian", hart3_bd, hart3_max,
                            c("hyperband_3",
                              "hyperband_4",
                              "seq_halving_max_rand",
                              "seq_halving_100_rand"),
                            "negative.hartmann.3")
    return(res)
    
}

tree_opt_hart3_gaus <- function(n_per_resource) {
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_hart3,
                            "gaussian", hart3_bd, hart3_max,
                            c("seq_tree_1_tree",
                              "seq_tree_fixed_prop",
                              "partition_tree_growing_4",
                              "partition_tree_fixed_4",
                              "partition_tree_growing_3"),
                            "negative.hartmann.3")
    
}

bayes_opt_hart6_gaus <- function(n_per_resource) {
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_hart6,
                            "gaussian", hart6_bd, hart6_max,
                            c("bayes_opt_growing_halving",
                              "bayes_opt_growing_hyper"),
                            "negative.hartmann.6")
    return(res)
}

bandit_opt_hart6_gaus <- function(n_per_resource) {
    
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_hart6,
                            "gaussian", hart6_bd, hart6_max,
                            c("hyperband_3",
                              "hyperband_4",
                              "seq_halving_max_rand",
                              "seq_halving_100_rand"),
                            "negative.hartmann.6")
    return(res)
    
}

tree_opt_hart6_gaus <- function(n_per_resource) {
    res <- run_mult_opt_exp(2.5 * 10^seq(3, 5, .5), n_per_resource, neg_hart6,
                            "gaussian", hart6_bd, hart6_max,
                            c("seq_tree_fixed_prop",
                              "partition_tree_growing_4"),
                            "negative.hartmann.6")
    
}


### plots for the results

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
    # only keep some algorithms for clarity
    keep_algos <- c("seq_tree_fixed_prop", "partition_tree_growing_4",
                    "hyperband_4", "seq_halving_max_rand",
                    "bayes_opt_growing_hyper")
    summary_res <- summary_res[summary_res$algorithm %in% keep_algos,]

    noise_names <- c(gaussian = "Gaussian: Variance = 1/4")
    objective_names <- c(negative.branin = "Branin",
                         negative.hartmann.3 = "Hartmann 3D",
                         negative.hartmann.6 = "Hartmann 6D")
    
    plt <- ggplot(summary_res, aes(x=n_resources, y=mean, color=algorithm)) +
        geom_line(size=2) +
        geom_errorbar(aes(ymin=mean+1.96*se,
                          ymax=mean-1.96*se),
                      width=.1) +
        geom_point() +
        scale_y_log10() +
        scale_x_log10() + 
        scale_color_manual("", values=berk_palette,
                           labels=c("Bayesian Optimization",
                                    "Hyperband",
                                    "SequentialHalving",
                                    "PartitionTree",
                                    "SequentialTree")) +
        xlab("Total Number of Samples (Log Scale)") +
        ylab("Function Error (Log Scale)") +
        facet_grid(noise_model ~ objective,
                   labeller = labeller(noise_model = noise_names,
                                       objective = objective_names)) +
        theme_bw() + 
        theme(text=element_text(size=30), plot.title = element_text(hjust=0.5),
              legend.position = "top")
    return(plt)

}


### More fine grained experiments for sequential tree

#' Count the proportion of times that the final partition in sequential
#' tree contains an arg max
#'
#' @param algorithm The algorithm to run, {sequential, partition}_tree
#' @param resources Vector of resources levels to try
#' @param n_per_resources Number of times to optimize at each resource level
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_arg_max A matrix of true arg maxes to compare to
#' @param eta Amount to cut the number of boxes by
#' @param nodes_or_rounds The maximum number of leaf nodes per tree or rounds
#' @param n_tree Number of trees in forest
arg_max_in_partition_exp <- function(algorithm, resources, n_per_resource,
                                     objective, noise_model, bounds,
                                     true_arg_max, eta, nodes_or_rounds, n_tree) {

    print(c(eta, n_tree, nodes_or_rounds))
    
    results <- vapply(resources,
                      function(x)
                          rowMeans(
                              replicate(n_per_resource,
                                        contains_arg_max(
                                            objective,
                                            true_arg_max,
                                            algorithm(objective,
                                                      noise_model,
                                                      bounds,
                                                      box_runif,
                                                      x,
                                                      eta,
                                                      nodes_or_rounds,
                                                      n_tree)[[3]]))),
                      double(3))
    #results <- data.frame(t(results))
    #names(results) <-  c("percent.true", "max.side.length")
    #results$n_resources <- resources
    #results$n_tree <- n_tree
    #results$eta <- eta
    #results$max_nodes <- max_nodes
    return(t(results))
}

contains_arg_max <- function(objective, true_arg_max, box) {
    # get the maximum side length of the box
    max_len <- max(box[,2] - box[,1])
    # get the function difference of the midpoint
    f_opt <- objective(true_arg_max[1,])
    f_final <- objective(rowMeans(box))
    return(c( 1* any(sapply(1:dim(true_arg_max)[1] ,
                      function(i)
                          any(sapply(1:dim(true_arg_max)[2],
                                     function(j)
                                         true_arg_max[i,j] >= box[j,1] &&
                                         true_arg_max[i,j] <= box[j,2])))),
             max_len, f_opt - f_final))
}

#' Count the proportion of times that the final partition in sequential
#' tree contains an arg max
#'
#' @param algorithm The algorithm to run, {sequential, partition}_tree
#' @param resources Vector of resources levels to try
#' @param n_per_resources Number of times to optimize at each resource level
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_arg_max A matrix of true arg maxes to compare to
#' @param etas Vector of amounts to cut the number of boxes by
#' @param nodes_or_rounds Vector of the maximum number of leaf nodes or rounds
#' @param n_trees Vector of number of trees in forest
arg_max_in_partition_mult_exp <- function(algorithm, resources,
                                          n_per_resource, objective,
                                          noise_model, bounds, true_arg_max, etas,
                                          nodes_or_rounds, n_trees) {

    grid <- expand.grid(etas, nodes_or_rounds, n_trees)
    dat <- as.data.frame(t(mapply(
        function(eta, node_or_round, n_tree)
            arg_max_in_partition_exp(algorithm,
                                     resources,
                                     n_per_resource,
                                     objective,
                                     noise_model,
                                     bounds,
                                     true_arg_max,
                                     eta,
                                     node_or_round,
                                     n_tree),
        grid$Var1,
        grid$Var2,
        grid$Var3)))
    names(dat) <- c("percent.true", "max.side.length", "f.error")
    dat$eta <- grid$Var1
    dat$max_nodes <- grid$Var2
    dat$n_tree <- grid$Var3
    return(dat)
}

#' Count the proportion of times that seq tree doesn't throw away the arg max
#' 
#' @param algorithm The algorithm to run, {sequential, partition}_tree
#' @param limit Number of resources to use
#' @param n_per_resources Number of times to try optimizing
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_max The true maximum
#' @param eta Amount to cut the number of boxes by
#' @param nodes_or_rounds The maximum number of leaf nodes per tree or rounds
#' @param n_tree Number of trees in forest
error_per_round_exp <- function(algorithm, limit, n_per_resource, objective,
                                noise_model, bounds, true_max, eta,
                                nodes_or_rounds, n_tree) {

    print(c(eta, n_tree))
    
    results <- replicate(n_per_resource, apply(
                                             algorithm(objective,
                                                       noise_model,
                                                       bounds,
                                                       box_runif,
                                                       limit,
                                                       eta,
                                                       nodes_or_rounds,
                                                       n_tree)[[4]],
                                             1,
                                             function(x) true_max - objective(x)))
                                                        
    results <- data.frame(results)
    results$round <- 1:dim(results)[1]
    results$eta <- rep(eta, dim(results)[1])
    return(results)
}

#' Count the proportion of times that seq tree doesn't throw away the arg max,
#' for various etas and number of trees
#' 
#' @param algorithm The algorithm to run, {sequential, partition}_tree
#' @param limit Number of resources to use
#' @param n_per_resources Number of times to try optimizing
#' @param objective The objective function to maximize
#' @param noise_model The distribution of the noise
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param true_max The true maximum
#' @param etas Vector of amounts to cut the number of boxes by
#' @param nodes_or_rounds The maximum number of leaf nodes per tree or rounds
#' @param n_tree Number of trees in forest
error_per_round_mult_exp <- function(algorithm, limit, n_per_resource, objective,
                                       noise_model, bounds, true_max, etas,
                                       nodes_or_rounds, n_tree) {
    return(plyr::ldply(
        etas,
        function(eta)
            error_per_round_exp(algorithm,
                                limit,
                                n_per_resource,
                                objective,
                                noise_model,
                                bounds,
                                true_max,
                                eta,
                                nodes_or_rounds,
                                n_tree)))
}


### Empirical results and plots

seq_tree_arg_max <- function(budget, n_exps) {
    max_max_nodes <- floor(exp(emdbook::lambertW_base(2 * budget * log(2)) - .1))
    res <- arg_max_in_partition_mult_exp(sequential_tree,
                                         c(budget),
                                         50,
                                         neg_branin,
                                         "gaussian",
                                         bran_bd,
                                         bran_arg_max,
                                         c(4),
                                         round(seq(10,max_max_nodes,
                                                   length.out=n_exps)),
                                         c(1))
    res$log.error <- log10(res$f.error)
    res$log.side.length <- log10(res$max.side.length)
    res$eta <- as.factor(res$eta)
    return(res)
    # make a plot
    melt_res <- reshape2::melt(res[,c("percent.true", "log.side.length", "eta",
                                      "max_nodes", "log.error")],
                               id.vars=c("eta","max_nodes"))
    
    plt <- ggplot(melt_res, aes(x=max_nodes, y=value, color=eta)) +
        geom_point() + geom_line(size=1.5) +
        facet_grid(variable ~ . , scales="free_y",
                   labeller = as_labeller(
                       c("percent.true" = "Proportion Containing Arg Max",
                         "log.side.length" =
                             "Average Side Length (Log Scale)",
                         "log.error" = "Function Error (Log Scale)"))) +
        theme_minimal() +
        xlab("Number of Elements in Initial Partition") +
        ylab("") +
        scale_color_manual(values=berk_palette) +
        ggtitle("SequentialTree Performance with 10,000 Samples on the Branin Function")
                       
    return(list(res, plt))
}

part_tree_arg_max <- function(budget, n_exps) {
    max_rounds <- floor(sqrt(budget / 20))
    res <- arg_max_in_partition_mult_exp(partition_tree,
                                         c(budget),
                                         50,
                                         neg_branin,
                                         "gaussian",
                                         bran_bd,
                                         bran_arg_max,
                                         c(2,3,4),
                                         round(seq(5, max_rounds,
                                                   length.out=n_exps)),
                                         c(10))
    res$log.error <- log10(res$f.error)
    res$log.side.length <- log10(res$max.side.length)
    res$eta <- as.factor(res$eta)
    return(res)
    # make a plot
    melt_res <- reshape2::melt(res[,c("percent.true", "log.side.length", "eta",
                                      "max_nodes", "log.error")],
                               id.vars=c("eta","max_nodes"))
    
    plt <- ggplot(melt_res, aes(x=max_nodes, y=value, color=eta)) +
        geom_point() + geom_line(size=1.5) +
        facet_grid(variable ~ . , scales="free_y",
                   labeller = as_labeller(
                       c("percent.true" = "Proportion Containing Arg Max",
                         "log.side.length" =
                             "Average Side Length (Log Scale)",
                         "log.error" = "Function Error (Log Scale)"))) +
        theme_minimal() +
        xlab("Number of Rounds") +
        ylab("") +
        scale_color_manual(values=berk_palette) +
        ggtitle("PartitionTree Performance with 10,000 Samples on the Branin Function")
                       
    return(list(res, plt))
}


seq_tree_error_per_round <- function(budget, n_exps) {
    max_nodes <- budget / 20
    results <- error_per_round_mult_exp(sequential_tree,
                                       budget,
                                       n_exps,
                                       neg_branin,
                                       "gaussian",
                                       bran_bd,
                                       bran_max,
                                       c(2,3,4),
                                       max_nodes,
                                       1)
    results$eta <- as.factor(results$eta)
    melted_results <- reshape2::melt(results,
                                     id.vars=c("eta", "round"))
    melted_results$log.val <- log10(melted_results$value)

    summary <- plyr::ddply(melted_results, c("eta", "round"),
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
                       lower.quantile = quantile(value[!is.na(value)], 1/4))

    plt <- ggplot(summary, aes(x=round, y=mean, color=eta)) +
        geom_point() +
        geom_line(size=2) +
        geom_errorbar(aes(ymin=mean+1.96*se,
                          ymax=mean-1.96*se),
                      width=.3) +
        theme_bw() +
        xlab("Round") +
        ylab("Function Error (Log Scale)") +
        scale_color_manual(expression(eta), values=berk_palette) +
        scale_y_log10() +
        scale_x_continuous(breaks=1:max(results$round)) +
        theme(text=element_text(size=30), plot.title = element_text(hjust=0.5)) + 
        ggtitle("SequentialTree Accuracy with 10,000 Samples on the Branin Function")
    
    return(result)
}

part_tree_error_per_round <- function(budget, n_exps) {
    max_rounds <- floor(sqrt(budget / 10))
    results <- error_per_round_mult_exp(partition_tree,
                                       budget,
                                       n_exps,
                                       neg_branin,
                                       "gaussian",
                                       bran_bd,
                                       bran_max,
                                       c(2,3,4,5),
                                       max_rounds,
                                       10)
    results$eta <- as.factor(results$eta)

    melted_results <- reshape2::melt(results,
                                     id.vars=c("eta", "round"))
    melted_results$log.val <- log10(melted_results$value)

    summary <- plyr::ddply(melted_results, c("eta", "round"),
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
                       lower.quantile = quantile(value[!is.na(value)], 1/4))

    plt <- ggplot(summary, aes(x=round, y=mean, color=eta)) +
        geom_point() +
        geom_line(size=2) +
        geom_errorbar(aes(ymin=mean+1.96*se,
                          ymax=mean-1.96*se),
                      width=.3) +
        theme_bw() +
        xlab("Round") +
        ylab("Function Error (Log Scale)") +
        scale_color_manual(expression(eta), values=berk_palette) +
        scale_y_log10() +
        scale_x_continuous(breaks=1:max(results$round)) +
        theme(text=element_text(size=30), plot.title = element_text(hjust=0.5)) +
        ggtitle("PartitionTree Accuracy with 10,000 Samples on the Branin Function")
    
    return(list(results, plt))
}
