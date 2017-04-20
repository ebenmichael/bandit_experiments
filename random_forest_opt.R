## Use Random Forests to learn a partition over feature space and
## then use bandit algorithms to choose the best partition


rf_partition <- function(data, targets, n_tree) {
    # fit a random forest with ntree trees
    rf <- randomForest::randomForest(data, targets, ntree=n_tree)
    # get the split points for the trees
    splits <- plyr::ldply(1:n_tree,
                     function(x) randomForest::getTree(rf, x)[,c(3,4)])
    
}
