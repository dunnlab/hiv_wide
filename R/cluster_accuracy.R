#' Given a set of thresholds, a list of "true" clusterings for each threshold, and a list of other tree clusterings for each threshold,
#' provides true positive rate, true negative rate, false positive rate, false negative rate. Tree clusterings must match up in label indices.
#'

#' Given a "true" tree clustering and a set of labels for other tree clusterings, match up the indices for that tree.
#'
#' @param true_tree true tree as ID'd in cluster
#' @params labels labels for tree clusterings
match_indices <- function(true_tree, labels) {
    true_labels <- true_tree[[2]]
    true_indices <- sapply(true_labels, function(x) which(sapply(labels, identical, x)))
    return(true_indices)
}

#' Given a "true" tree clustering and a list of other tree clusterings, provides true positive rate, true negative rate,
#' false positive rate, and false negative rate
#'
#' @param true_tree true tree
cluster_accuracy <- function(true_tree, clusterlist) {

    # number of common clusters = tp
    sharedclust <- lapply(clusterlist, function(y) sapply(y, function(x) intersect(x[[3]],y$wgs[[3]]) %>% length))
tp <- data.frame(sharedclust)
names(tp) <- c("70","80","85","90","95","99")
tp$tree <- rownames(tp)
tp <- pivot_longer(tp, cols=c("70","80","85","90","95","99"), names_to="bootstrap", values_to="tp")
# number of total clusters
totalclust <- lapply(allclusters, function(y) sapply(y, function(x) length(x[[3]])))

# number of total - common = clusters specific to given tree = fp
fp <- data.frame(totalclust)
names(fp) <- c("70","80","85","90","95","99")
fp$tree <- rownames(fp)

# number of wgs - common = clusters specific to wgs = fn
# https://stackoverflow.com/questions/50957998/mutate-inside-group-by-based-on-particular-row-inside-group

# number of possible clusters is number of nodes in tree which here we know is 1194
# 1194 - fn = tn
cluster_compare <- pivot_longer(fp, cols=c("70","80","85","90","95","99"), names_to="bootstrap",values_to="totalclust") %>% left_join(tp) %>%
    mutate(fp=totalclust-tp) %>%
    group_by(bootstrap) %>%
    mutate(fn=totalclust[tree=="wgs"]-tp, tn=1194-fn, precision=tp/(tp+fp),recall=tp/(tp+fn))
}
