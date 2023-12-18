#' Loads trees, creates phylos object as well as clade labels of all of the trees
#'
#' @param treepath Path to the tree files
loadtree <- function(treepath, names) {
    treefiles <- list.files(treepath,".fa.treefile",full.names=TRUE)
    trees <- lapply(treefiles, read.iqtree)
    names(trees) <- names
    phylos <- lapply(trees, function(x) x@phylo)
    bipartitions <- prop.part(phylos)
    labels <- lapply(1:length(bipartitions), function(x) attributes(bipartitions)$labels[bipartitions[x][[1]]] %>% sort)
    return(list(trees, labels))
}
