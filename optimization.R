source("bandits.R")


#' Create a bandit object given a function, values, and a noise model
#'
#' @param objective The objective function to minimize (maximize)
#' @param noise_model The type of noise model
#' @param values
create_bandit <- function(objective, noise_model, values) {
    if(noise_model == "gaussian") {
        # create a normal bandit with variance=1/4 
        return(normal_bandit(objective(values), rep(1/4, length(values))))
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
bandit_opt <- function(objective, noise_model, get_values, n_values, limit,
                           mab_algo, ...) {
    # get the candidate values
    values <- get_values(n_values)
    # create bandit object for these values
    function_bandit <- create_bandit(objective, noise_model, values)
    # run the best arm selection algorithm
    output <- mab_algo(function_bandit, limit, ...)
    return(values[output$best_idx])
}
