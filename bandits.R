### Best Arm algorithms for MAB



## Create an arm with normally distributed rewards
normal_arm <- function(mu, sig2) {
    arm <- list(mean=mu, var=sig2)
    class(arm) <- "normal_arm"
    return(arm)
}

## Override the sample method to be able to have an S3 generic (funky stuff)
sample.default = sample
sample <- function(obj, ...) {
    UseMethod("sample")
}

sample.normal_arm <- function(arm, n) {
    return(arm$mean + sqrt(arm$var) * rnorm(n))
}

## Bandit class
bandit <- function(arms) {
    bandit <- arms
    class(bandit) <- "bandit"
    return(bandit)
}

## Create a bandit with normal arms
normal_bandit <- function(means, vars) {
    arms <- list()
    for(i in 1:length(means)) {
        arms[[i]] <- normal_arm(means[i], vars[i])
    }
    n_bandit <- arms
    class(n_bandit) <- c("normal_bandit", "bandit")
    return(n_bandit)
}


##S3 generic method
pull_arm <- function(obj, ...) {
    UseMethod("pull_arm")
}

## Method to pull arm from a bandit
pull_arm.bandit <- function(bandit, i, n) {
    return(sample(bandit[[i]], n))
}

## S3 generic method
pull_arms <- function(obj, ...) {
    UseMethod("pull_arms")
}

## Pull all the arms from a bandit
pull_arms.bandit <- function(bandit, n_samples) {
    samples <- matrix(0, nrow = length(bandit), ncol = n_samples)
    for(i in 1:length(bandit)) {
        samples[i,] <- pull_arm(bandit, i, n_samples)
    }
    return(samples)
}

## Method to pull all the arms from a normal bandit at once
pull_arms.normal_bandit <- function(bandit, n_samples) {
    means <- integer(length(bandit))
    vars <- integer(length(bandit))
    for(i in 1:length(bandit)) {
        means[i] <- bandit[[i]]$mean
        vars[i] <- bandit[[i]]$var
    }
    # get random iid standard gaussians 
    iids <- matrix(rnorm(length(bandit) * n_samples), nrow = length(bandit),
                   ncol = n_samples)
    samples <- means + vars * iids
    return(samples)
}

## Naive PAC(epsilon, delta) algo (Even-Dar 2006)
## Need to update this for general sub-Gaussian RVs
naive_best <- function(bandit, epsilon, delta)  {
    # pull each arm the same number of times
    n_arms <- length(bandit)
    n_pulls <- ceiling((4 / epsilon ^ 2) * log( 2 * n_arms / delta))
    print(n_pulls)
    sample_avgs <- integer(n_arms)
    for(i in 1:n_arms) {
        sample_avgs[i] <- mean(pull_arm(bandit, i, n_pulls))
    }
    # return arm with best sample average
    return(which.max(sample_avgs))
}

## Median Eliminiation Algorithm (Even-Dar 2006)
median_elimination <- function(bandit, epsilon, delta) {
    # initialize set of accepted arms to be all of them
    accepted <- bandit
    eps_l <- epsilon / 4
    delta_l <- delta / 2
    l <- 1
    # keep track of which arms remain
    idxs <- 1:length(bandit)
    while(length(accepted) > 1) {
        # number of times to sample 
        n_pulls <- ceiling(1 / (eps_l / 2) ^ 2 * log(3 / delta_l))
        print("----------")
        print(n_pulls)
        print(l)
        
        #sample the arms
        emp_avgs <- rowMeans(pull_arms(accepted, n_pulls))

        # get the median     
        med <- median(emp_avgs)
        ## keep everything above the median
        accepted <- bandit(accepted[emp_avgs >= med])
        idxs <- idxs[emp_avgs >= med]
        eps_l <- 3 * eps_l / 4
        delta_l < delta_l / 2
        l <- l + 1
    }
    return(idxs[1])
}



## Exponential gap algorithm for fixed confidence, Karnin 2013
exponential_gap <- function(bandit, delta) {
    # initialize set of accepted arms to be all of them    
    accepted <- bandit
    r <- 1
    idxs <- 1:length(bandit)
    while(length(accepted) > 1) {
        eps_r <- 2^(-r) / 4
        delta_r <- delta / (50 * r^3)
        # sample each arm 
        n_pulls <- ceiling((2 / eps_r^2) * log (2 / delta_r))
        print(sprintf("NUMBER OF PULLS PER ARM: %f", n_pulls))
        #sample the arms
        emp_avgs <- rowMeans(pull_arms(accepted, n_pulls))
        # call median_elimination to find an arm within eps_r/2
        # with prob delta_r
        print("STARTING MEDIAN ELIMINATION")
        med_elim_idx <- median_elimination(accepted, eps_r / 2,
                                           delta_r / 2)
        print("FINISHED MEDIAN ELIMINATION")
        # throw out everything more than eps_r away
        keep_bool <- emp_avgs >= emp_avgs[med_elim_idx] - eps_r
        accepted <- bandit(accepted[keep_bool])
        idxs <- idxs[keep_bool]
        r <= r + 1
    }
    return(idxs[1])
}


#' Perform the successive elimnination algorithm
#'
#' @param bandit The multi-armed bandit object
#' @param delta The confidence that we have the right mean
#' @param n_min The minimum number of samples we can take
#'
#' @return index of the chosen best arm
successive_elimination <- function(bandit, delta, n_min) {

    accepted = bandit
    n_arms <- length(bandit)
    t <- n_min
    idxs <- 1:n_arms
    emp_avgs <- integer(n_arms)
    while(length(accepted) > 1) {
        # pull n_min times and update current averages
        emp_avgs <- (rowSums(pull_arms(accepted, n_min)) +
                     (t- n_min) * emp_avgs) / t
        max_avg <- max(emp_avgs)
        bound <- sqrt(log(pi^2 / 3 * n_arms * t^2 / delta) / t)
        # throw out anything we can confidently reject
        keep_bool <- (max_avg - emp_avgs) < 2 * bound
        accepted <- bandit(accepted[keep_bool])
        idxs <- idxs[keep_bool]
        emp_avgs <- emp_avgs[keep_bool]
        t <- t + n_min
    }
    return(idxs[1])
}
