source("bandits.R")
library(Rsolnp)
library(gaussianProcess)

#' Uniformly sample points in a box
#'
#' @param n_samples Number of points to sample in that box
#' @param bounds A d x 2 matrix of box constraints for each variable
box_runif <- function(n_samples, bounds) {
    samples <- apply(bounds, 1, function(x) runif(n_samples, x[1], x[2]))
    return(matrix(samples, ncol=dim(bounds)[1]))
}

#' Create a bandit object given a function, values, and a noise model
#'
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param values
create_bandit <- function(objective, noise_model, values) {
    n_values <- length(values)
    obj_values <- apply(values, 1, objective)
    if(noise_model == "gaussian") {
        # create a normal bandit with variance=1/4 
        return(normal_bandit(obj_values, rep(1/4, n_values)))
    } else {
        stop("Only gaussian noise is supported")
    }
}



#' Optimize a function given noisy evaluations using a bandit method
#'
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param get_values A function to propose a set of candidate values
#' @param n_values The number of values to try
#' @param limit The confidence or budget for the best arm selection algorithms
#' @param mab_algo The best arm selection algorithm to use
#'
#' @return The best parameter setting found and the number of times each
#'         setting was tried
bandit_opt <- function(objective, noise_model, get_values, n_values, limit,
                           mab_algo, ...) {
    # get the candidate values
    values <- get_values(n_values)
    # create bandit object for these values
    function_bandit <- create_bandit(objective, noise_model, values)
    # run the best arm selection algorithm
    output <- mab_algo(function_bandit, limit, ...)
    return(list(values[output$best_idx,],
                cbind(apply(values[output$arm_pull_mat[,1],], 1, objective), output$arm_pull_mat[,2])))
                #output$arm_pull_mat))
}

#' Perform the Hyperband routine (Li 2017)
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param get_values A function to propose a set of candidate values
#' @param limit The maximum amount of pulls for any parameter setting
#' @param eta The proportion of configurations discarded at each round
#'
#' @return The best parameter setting found and the number of times each
#'         setting was tried
hyperband <- function(objective, noise_model, bounds, get_values, limit, eta) {
    s_max <-floor(logb(limit, eta))
    B <- (s_max +1) * limit
    ## book keeping
    n_samples <- 0
    best_val <- -Inf
    best_arg_max <- NA
    for(s in s_max:0) {
        # get some parameters
        n_values <- ceiling( B / limit / (s+1) * eta^s)
        resources <- floor(limit * eta^(-s))
        values <- get_values(n_values, bounds)
        # sequential having inner loop
        for(i in 0:s) {
            n_i <- floor(n_values * eta^(-i))
            r_i <- resources * eta^i 
            ## book keeping
            n_samples <- n_samples + n_i * r_i
            # sample and get the mean
            samples <- apply(values,
                             1,
                             function(x) sample_function(objective,
                                                         noise_model,
                                                         x, r_i))
            # if only one sample is drawn, samples will be a vector
            if(is.null(dim(samples))) {
                emp_avgs <- samples
            } else {
                emp_avgs <- colMeans(samples)
            }
            
            n_avgs <- length(emp_avgs)
            if(n_avgs > eta) {
                ## if there is more than one value keep the top n_i / eta
                num_top <- floor(n_i / eta)
                comparitor <- sort(emp_avgs,
                                   partial=(n_avgs -
                                            num_top + 1))[n_avgs -
                                                          num_top + 1]
                keep_bool <- emp_avgs >= comparitor
                values <- values[keep_bool, , drop=FALSE]
                emp_avgs <- emp_avgs[keep_bool]
            } 
        }
        if(max(emp_avgs) > best_val) {

            best_arg_max <- values[which.max(emp_avgs),]
            best_val <- max(emp_avgs)
        }
    }
    return(list(best_arg_max, n_samples))
}

#' Sample from the value of a function under a given noise model
#'
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param value The value of the objective to sample at
#' @param n_samples The number of times to sample each point
sample_function <- function(objective, noise_model, value, n_samples) {
    obj_value <- objective(value)
    if(noise_model == "gaussian") {
        return(t(obj_value + rnorm(n_samples) / 4))
    } else {
        stop("Only gaussian noise is supported")
    }
    
}


#' Compute the expected improvement at x for a given gp
#'
#' @param gp trained gaussianProcess
#' @param value The value to compute at
#' @param max_mean The maximum of the mean function
expected_improvement <- function(gp, value, max_mean) {

    pred <- predict(gp, value)
    mean <- pred$mean
    sd <- sqrt(pred$covariance)
    z <- (mean - max_mean) / sd
    return(sd * z * pnorm(z) + sd * dnorm(z))
}

#' Perform Bayesian optimization using a GP with RBF kernel and
#' hyper parameters chosen by optimizing the marginal log likelihood
#'
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param n_samples The number of times to sample each point
#' @param n_values The number of values to try
#' @param bounds A d x 2 matrix of box constraints for each variable
#' @param n_init The number of initial samples to feed the GP, defaults to 5
bayes_opt <- function(objective, noise_model, n_samples, n_values,
                      bounds, n_init=10) {


    # get starting values randomly
    values <- apply(bounds, 1, function(x) runif(n_init,x[1],x[2]))
    values <- rbind(values, matrix(NA, nrow=(n_values - n_init), ncol=dim(values)[2]))
    targets <- apply(values, 1, function(x) mean(sample_function(objective,
                                                            noise_model,
                                                            x, n_samples)))

    targets <- c(targets, rep(NA, n_values - n_init))
    arg_max_mu <- as.vector(apply(bounds, 1, function(x) runif(1,x[1],x[2])))
    next_x <- as.vector(apply(bounds, 1, function(x) runif(1,x[1],x[2])))
    for(i in (1 + n_init):n_values) {
        print(i)
        # fit the gp
        print("fitting gp")
        tryCatch({
            gp <<- gaussianProcess(values[1:(i-1),],
                                   targets[1:(i-1)],
                                   noise.var=1 / (n_samples * 4))
            ## return(gp_tmp)
        }, error = function(e) {
            gp <<- NA
            })
        if(is.na(gp)) return(NA)
        # compute the maximum of the mean function
        print("Finding max of mean function")
        tryCatch({ arg_max_mu <<- rgenoud::genoud(function(x) predict(gp,
                                                                      t(x))$mean,
                                                  dim(values)[2],
                                                  max=TRUE,
                                                  Domains=bounds,
                                                  boundary.enforcement = 1,
                                                  print.level=0,
                                                  gradient.check=FALSE,
                                                  pop.size=10)$par
        }, error = function(e) {
            arg_max_mu <<- NA
        })
        if(is.na(arg_max_mu)) return(NA)
        max_mu <- predict(gp, t(arg_max_mu))$mean
        
        # optimize expected improvement
        print("Optimizing EI")
        tryCatch({next_x <<- t(rgenoud::genoud(function(x)
            expected_improvement(gp,
                                 t(x),
                                 max_mu),
            dim(values)[2], max=TRUE,
            Domains=bounds,
            boundary.enforcement = 1,
            gradient.check=FALSE,
            print.level=0,
            pop.size=100)$par)
        }, error = function(e) {
            next_x <<- NA
        })
        if(is.na(next_x)) return(NA)
        values[i,] <- next_x
        targets[i] <- mean(sample_function(objective,
                                           noise_model, next_x, n_samples))
    }

    return(next_x)
}
