source("bandits.R")


## Ran an experiment with a given bandit
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
