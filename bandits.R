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
pull_arm.bandit <- function(band, i, n) {
    return(sample(band[[i]], n))
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

#' Run the naive PAC(epsilon, delta) algo (Even-Dar 2006)
#'
#' @param bandit The bandit object
#' @param epsilon How close we are to the best with high probability
#' @param delta The confidence
#'
#' @return The chosen arm and a matrix keeping track of each arm
#'         and how much it was pulled
naive_best <- function(bandit, epsilon, delta)  {
    # pull each arm the same number of times
    n_arms <- length(bandit)
    n_pulls <- ceiling((4 / epsilon ^ 2) * log( 2 * n_arms / delta))
    emp_avgs <- rowMeans(pull_arms(bandit, n_pulls))
    # keep track of how many times each arm was sampled
    arm_pull_mat <- cbind(seq(1,n_arms), rep(n_pulls, n_arms))
    # return arm with best sample average
    output <- list(best_idx = which.max(emp_avgs),
                   arm_pull_mat = arm_pull_mat)
    return(output)
}

## Median Eliminiation Algorithm (Even-Dar 2006)
#'
#' @param bandit The bandit object
#' @param epsilon How close we are to the best with high probability
#' @param delta The confidence
#'
#' @return The chosen arm and a matrix keeping track of each arm
#'         and how much it was pulled
median_elimination <- function(bandit, epsilon, delta) {
    # initialize set of accepted arms to be all of them
    accepted <- bandit
    eps_l <- epsilon / 4
    delta_l <- delta / 2
    l <- 1
    n_arms <- length(bandit)
    # keep track of which arms remain
    idxs <- 1:n_arms
    emp_avgs <- integer(n_arms)
    total_pulls <- 0
    # start book-keeping which arms are pulled how many times
    arm_pull_mat <- matrix(ncol = 2)

    while(length(accepted) > 1) {
        # number of times to sample 
        n_pulls <- ceiling(1 / (eps_l / 2) ^ 2 * log(3 / delta_l))
        total_pulls <- n_pulls + total_pulls
        
        #sample the arms
        emp_avgs <- rowMeans(pull_arms(accepted, n_pulls))

        # book-keeping, the arms pulled and how much they were pulled
        arm_pull_mat <- rbind(arm_pull_mat,
                              cbind(idxs, rep(n_pulls, length(idxs))))
        # get the median     
        med <- median(emp_avgs)
        ## keep everything above the median
        accepted <- bandit(accepted[emp_avgs >= med])
        idxs <- idxs[emp_avgs >= med]
        # update values
        eps_l <- 3 * eps_l / 4
        delta_l < delta_l / 2
        l <- l + 1
    }
    output <- list(best_idx = idxs[1],
                   arm_pull_mat = arm_pull_mat[-1,])
    return(output)
}



#' Exponential gap algorithm for fixed confidence, Karnin 2013
#'
#' @param bandit The bandit object
#' @param delta The confidence
#'
#' @return The chosen arm and a matrix keeping track of each arm
#'         and how much it was pulled
exponential_gap <- function(bandit, delta) {
    # initialize set of accepted arms to be all of them    
    accepted <- bandit
    r <- 1
    n_arms <- length(bandit)
    idxs <- 1:n_arms
    emp_avgs <- integer(n_arms)
    total_pulls <- 0
    # start book-keeping which arms are pulled how many times
    arm_pull_mat <- matrix(ncol = 2)
    
    while(length(accepted) > 1) {
        eps_r <- 2^(-r) / 4
        delta_r <- delta / (50 * r^3)
        # number of times to sample each arm 
        n_pulls <- ceiling((2 / eps_r^2) * log (2 / delta_r))
        total_pulls <- n_pulls + total_pulls
        
        # book-keeping, the arms pulled and how much they were pulled
        arm_pull_mat <- rbind(arm_pull_mat,
                              cbind(idxs, rep(n_pulls, length(idxs))))
        
        #sample the arms
        emp_avgs <- (rowSums(pull_arms(accepted, n_pulls)) +
                     (total_pulls - n_pulls) * emp_avgs) / total_pulls

        # call median_elimination to find an arm within eps_r/2
        # with prob delta_r
        print("STARTING MEDIAN ELIMINATION")
        me_output <- median_elimination(accepted, eps_r / 2,
                                        delta_r / 2)
        med_elim_idx <- me_output$best_idx
        # more book-keeping on which arms are pulled how many times
        arm_pull_mat <- rbind(arm_pull_mat, me_output$arm_pull_mat)
        print("FINISHED MEDIAN ELIMINATION")
        # throw out everything more than eps_r away
        keep_bool <- emp_avgs >= emp_avgs[med_elim_idx] - eps_r
        accepted <- bandit(accepted[keep_bool])
        idxs <- idxs[keep_bool]
        emp_avgs <- emp_avgs[keep_bool]
        r <= r + 1
    }
    output <- list(best_idx = idxs[1],
                   arm_pull_mat = arm_pull_mat[-1,])
    return(output)
}


#' Perform the successive elimnination algorithm (Even-Dar 2006)
#'
#' @param bandit The multi-armed bandit object
#' @param delta The confidence that we have the right mean
#' @param n_min The minimum number of samples we can take
#'
#' @return The chosen arm and a matrix keeping track of each arm
#'         and how much it was pulled
successive_elimination <- function(bandit, delta, n_min) {

    accepted <- bandit
    n_arms <- length(bandit)
    t <- n_min
    idxs <- 1:n_arms
    emp_avgs <- integer(n_arms)
    # start book-keeping which arms are pulled how many times
    all_idxs <- rep(NA, 1000)
    n_samples <- rep(NA, 1000)
    j <- 1
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
        # book keeping
        if(j + length(idxs) > length(all_idxs)) {
            all_idxs <- c(all_idxs, rep(NA, length(all_idxs) * 2))
            n_samples <- c(n_samples, rep(NA, length(n_samples) * 2))
        }
        all_idxs[j:(j + length(idxs) - 1)] <- idxs
        n_samples[j:(j + length(idxs) - 1)] <-
            rep(n_min, length(idxs))
        j <- j + length(idxs)
    }
    
    output <- list(best_idx = idxs[1],
                   arm_pull_mat = cbind(all_idxs[!is.na(all_idxs)],
                                        n_samples[!is.na(n_samples)]))
    return(output)
}


#' Perform the lil'LUCB algorithm (Jamieson review paper)
#'
#' @param bandit The multi-armed bandit object
#' @param delta The confidence that we have the right mean
#' @param epsilon The LIL parameters
#' @param n_min The minimum number of samples we can take
#'
#' @return The chosen arm and a matrix keeping track of each arm
#'         and how much it was pulled
lil_lucb <- function(bandit, delta, epsilon, n_min) {
    accepted <- bandit
    n_arms <- length(bandit)
    t <- n_min
    idxs <- 1:n_arms
    # sample every arm n_min times to initialize the averages
    emp_avgs <- rowMeans(pull_arms(bandit, n_min))
    # keep track of how many times an arm has been sampled
    total_pulls <- rep(n_min, n_arms)
    # start book-keeping which arms are pulled how many times
    all_idxs <- c(idxs, rep(NA, 1000))
    n_samples <- c(rep(n_min, n_arms), rep(NA, 1000))
    j <- n_arms + 1
    # start pulling arms sequentially
    finished <- FALSE

    while(! finished) {
        # get the current best
        h_t <- which.max(emp_avgs)


        bound <- (1 + sqrt(epsilon)) * sqrt(2 * (1+epsilon) *
                                            log( log((1 + epsilon) *
                                                     total_pulls + 2)
                                                / (delta/n_arms)) /
                                            total_pulls)
        avg_plus_bound <- emp_avgs + bound
        # set the value at h_t to -infinity so we don't choose it
        avg_plus_bound[h_t] <- -Inf
        l_t <- which.max(avg_plus_bound)
        # the stopping condition
        finished <- emp_avgs[h_t] - bound[h_t] > emp_avgs[l_t] + bound[l_t]

        # sample from the arms defined by h_t and l_t

        h_t_pulls <- pull_arm(bandit, h_t, n_min)
        l_t_pulls <- pull_arm(bandit, l_t, n_min)

        emp_avgs[h_t] <- (sum(h_t_pulls) +
                          emp_avgs[h_t] *
                          total_pulls[h_t]) / (total_pulls[h_t] + n_min)
        emp_avgs[l_t] <- (sum(l_t_pulls) +
                          emp_avgs[l_t] *
                          total_pulls[l_t]) / (total_pulls[l_t] + n_min)
        # keep track of the pulls
        total_pulls[h_t] <- total_pulls[h_t] + n_min
        total_pulls[l_t] <- total_pulls[l_t] + n_min
        # book keeping
        if(j + 1 > length(all_idxs)) {
            all_idxs <- c(all_idxs, rep(NA, length(all_idxs) * 2))
            n_samples <- c(n_samples, rep(NA, length(n_samples) * 2))
        }
        all_idxs[j] <- h_t
        all_idxs[j + 1] <- l_t
        n_samples[j] <- n_min
        n_samples[j+1] <- n_min
        j <- j + 2
    }
    output <- list(best_idx=h_t,
                   arm_pull_mat = cbind(all_idxs[!is.na(all_idxs)],
                                        n_samples[!is.na(n_samples)]))
    return(output)
}


#' Perform lil'UCB + LS stopping condition (Jamieson 2013)
#'
#' @param bandit The multi-armed bandit object
#' @param conf The confidence that we have the right mean
#' @param epsilon The LIL parameters
#' @param alpha Multiple of total times any arm is sampled in order to stop
#' @param beta Hyperparameter for the algorithm
#' @param n_min The minimum number of samples we can take
#'
#' @return The chosen arm and a matrix keeping track of each arm
#'         and how much it was pulled

lil_ucb <- function(bandit, conf, epsilon, alpha, beta, n_min) {
    # set delta to account for the LIL (and divide by 2 for the LS stopping criterion)
    delta <- (conf * epsilon / (5 * (2 + epsilon))) ^ (1 / (1 + epsilon)) / 2
    accepted <- bandit
    n_arms <- length(bandit)
    t <- n_min
    idxs <- 1:n_arms
    # sample every arm n_min times to initialize the averages
    emp_avgs <- rowMeans(pull_arms(bandit, n_min))
    # keep track of how many times an arm has been sampled
    total_pulls <- rep(n_min, n_arms)
    # start book-keeping which arms are pulled how many times
    all_idxs <- c(idxs, rep(NA, 1000))
    n_samples <- c(rep(n_min, n_arms), rep(NA, 1000))
    j <- n_arms + 1
    # start pulling arms sequentially
    finished <- FALSE
    while(! finished) {
        # count the number of times every arm but a certain arm was pulled
        sum_pulls <- rep(sum(total_pulls), n_arms) - total_pulls
        # stop if any arm has been pulled more than alpha percent of the time        
        finished <- any(total_pulls >= 1 + alpha * sum_pulls)
        # get the current best
        h_t <- which.max(emp_avgs)
        # compute the finite sample any time lil bound for each arm
        # for choosing the best
        choose_bound <- (1 + beta) * (1 + sqrt(epsilon)) * sqrt(2 * (1+epsilon) *
                                            log( log((1 + epsilon) *
                                                     total_pulls + 2)
                                                / delta) /
                                            total_pulls)
        # take the arm which maximizies the average plus the choose bound
        curr_best <- which.max(emp_avgs + choose_bound)

        # compute the bound for the LIL stopping condition
        stopping_bound <- (1 + sqrt(epsilon)) * sqrt(2 * (1+epsilon) *
                                            log( log((1 + epsilon) *
                                                     total_pulls + 2)
                                                / delta/(2 * n_arms)) /
                                            total_pulls)
        avg_plus_bound <- emp_avgs + stopping_bound        
        # set the value at h_t to -infinity so we don't choose it
        avg_plus_bound[h_t] <- -Inf
        # get the next best arm  
        l_t <- which.max(avg_plus_bound)      
        # the LIL stopping condition
        finished <- emp_avgs[h_t] - stopping_bound[h_t] >
            emp_avgs[l_t] + stopping_bound[l_t]
        # sample from the current best arm
        emp_avgs[curr_best] <- (sum(pull_arm(bandit, curr_best, n_min)) +
                          emp_avgs[curr_best] *
                          total_pulls[curr_best]) / (total_pulls[curr_best] + n_min)
        # book keeping
        total_pulls[curr_best] <- total_pulls[curr_best] + n_min

        if(j > length(all_idxs)) {
            all_idxs <- c(all_idxs, rep(NA, length(all_idxs) * 2))
            n_samples <- c(n_samples, rep(NA, length(n_samples) * 2))
        }
        all_idxs[j] <- curr_best
        n_samples[j] <- n_min
        j <- j + 1
        
    }
    output <- list(best_idx=curr_best,
                   arm_pull_mat = cbind(all_idxs[!is.na(all_idxs)],
                                        n_samples[!is.na(n_samples)]))
    return(output)
}



#' Perform the BatchRacing procedure (Jun 2016)
#'
#' @param bandit The bandit object
#' @param delta The confidence
#' @param num_top The number of top arms we want to find
#' @param batch_size The batch size
#' @param pull_limit The repeated pull limit
batch_racing <- function(bandit, delta, num_top, batch_size, pull_limit) {

    ## Define a helper function round_robin
    round_robin <- function(arms, total_pulls, b, r) {

        n_arms <- length(arms)
        num_new_pulls <- integer(n_arms)
        # get the number of pulls to make
        for(i in 1:min(b, n_arms * r)) {
            j <- which.min(total_pulls + num_new_pulls)
            num_new_pulls[j] <- num_new_pulls[j] + 1
        }
        return(num_new_pulls)
    }

    # initialize the sets used in the algorithm
    n_arms <- length(bandit)
    idxs <- 1:n_arms
    remaining <- idxs
    accepted <- c()
    rejected <- c()
    total_pulls <- integer(n_arms)
    emp_avgs <- integer(n_arms)
    # book keeping for number of times an arm was pulled
    all_idxs <- rep(NA, 1000)
    n_samples <- rep(NA, 1000)
    j <- 1
    
    while(length(remaining) >= 1) {
        
        # use round_robin to see how many pulls to do 
        new_pulls <- round_robin(remaining, total_pulls,
                                 batch_size, pull_limit)

        # book keeping
        if(j  + length(new_pulls[new_pulls!=0])> length(all_idxs)) {

            all_idxs <- c(all_idxs, rep(NA, length(all_idxs) * 2))
            n_samples <- c(n_samples, rep(NA, length(n_samples) * 2))
        }


        #print(total_pulls)        
        for(k in 1:length(new_pulls)) {
            if(new_pulls[k] != 0) {
                arm_k <- remaining[k]
                arm_k_pulls <- pull_arm(bandit, arm_k, new_pulls[k])
                                        # update empirical average
                emp_avgs[k] <- ( sum(arm_k_pulls) +
                                 emp_avgs[k] * total_pulls[k]) /
                    (new_pulls[k] + total_pulls[k])
                total_pulls[k] <- total_pulls[k] + new_pulls[k]

                # keep track of how much each arm was pulled
                all_idxs[j] <- arm_k
                n_samples[j] <- new_pulls[k]
                j <- j + 1
            }
        }

        # get upper and lower confidence bounds from the LIL
        non_zero_pulls <- total_pulls[total_pulls!=0]
        # if we didn't pull an arm yet, set bound to infinity
        bound <- rep(Inf, length(emp_avgs))
        bound[total_pulls != 0] <- sqrt(4 * log(log2(2 * non_zero_pulls) /
                              (sqrt(delta / (6 * n_arms))))
                      / non_zero_pulls)
        upper <- emp_avgs + bound
        lower <- emp_avgs - bound

        # accept and reject those for which we can do it confidently
        k_t <- num_top - length(accepted)
        accept_bool <- rep(FALSE, length(remaining))
        reject_bool <- rep(FALSE, length(remaining))
        for(i in 1:length(remaining)) {
            num_greater <- sum(lower[i] > upper)
            accept_bool[i] <- num_greater >= length(remaining) - k_t
            num_less <- sum(upper[i] < lower)
            reject_bool[i] <- num_less >= k_t
        }
        
        accepted <- c(accepted, remaining[accept_bool])
        rejected <- c(rejected, remaining[reject_bool])
        remain_bool <- !(accept_bool | reject_bool)
        remaining <- remaining[remain_bool]

        emp_avgs <- emp_avgs[remain_bool]
        total_pulls <- total_pulls[remain_bool]

    }


    output <- list(best_idxs=accepted,
                   arm_pull_mat = cbind(all_idxs[!is.na(all_idxs)],
                                        n_samples[!is.na(n_samples)]))
    return(output)
    
    
}
