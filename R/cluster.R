#' Computes clusters with a given bootstrap value on a given tree
#'
#' @param t treeio object, must have both @phylo and @data component
#' @param b bootstrap value
#' @param labels Set of labels from which to get correspondance with
cluster <- function(t, b, labels) {
    bootstraps <- filter(t@data, UFboot >= b)$node
    copy_bootstraps <- bootstraps

    # to cluster, for a given node with UFboot > b, are its children also UFboot > b?
    # we start at the beginning of the vector since those are the largest subtrees / closest to the root and delete the children from the list until we hit a leaf / there are no more children
    i <- 1
    node <- bootstraps[i]
    tip_length <- length(t@phylo$tip.label)
    #tips <- 1:3860 # this is specific to our phylogenetic trees...?
    tips <- 1:tip_length
    clusters <- c()
    while(!is.na(node)) {
        subnodes <- treeio::offspring(t@phylo, node)
        subnodes <- subnodes[-which(subnodes %in% tips)] # remove tips
        # if all subnodes are in the set of nodes with bootstrap > b, delete the subnodes and add node
        # to set of clusters
        if (length(subnodes)==0) {
            clusters <- c(clusters, node)
        } else if (all(subnodes %in% bootstraps)) {
            # this should delete them...
            copy_bootstraps <- copy_bootstraps[-which(copy_bootstraps %in% subnodes)]
            clusters <- c(clusters, node)
        }
        i <- i + 1
        node <- copy_bootstraps[i]
    }

    cluster_tips <- lapply(clusters, function(x) offspring(t@phylo, x, type="tips"))
    cluster_labels <- lapply(cluster_tips, function(x) t@phylo$tip.label[x] %>% sort)
    cluster_index <- sapply(cluster_labels, function(x) which(sapply(labels, identical, x)))

    return(list(cluster_tips, cluster_labels, cluster_index))
}
