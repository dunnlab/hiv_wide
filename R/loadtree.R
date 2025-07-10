#' Loads trees, creates phylos object as well as clade labels of all of the trees
#'
#' @param treepath Path to the tree files
loadtree <- function(treepath, names, masks, wgs=NULL) {
    treefiles <- list.files(treepath,".fa.treefile",full.names=TRUE)
    trees <- lapply(treefiles, function(x)
        tryCatch(read.iqtree(x), error=function(e) NULL))
    if(!is.null(wgs)) {
        trees$wgs <- wgs
    }

    # tree naming
    for(i in 1:length(masks)) {
        mask=masks[i]
        #print(paste("mask:",mask))
        tree_index <- grep(mask, treefiles)
        #print(paste("tree_index: ", tree_index))
        if(length(tree_index)==1) {
            #print(paste("i:",i))
            names(trees)[tree_index] <- names[i]
            #print(i)
            #print(names(trees))
        }
        #print("next round")
    }

    trees <- Filter(Negate(is.null), trees) # filter out nulls
    phylos <- lapply(trees, function(x) x@phylo)
    bipartitions <- prop.part(phylos)
    labels <- lapply(1:length(bipartitions), function(x) attributes(bipartitions)$labels[bipartitions[x][[1]]] %>% sort)
    return(list(trees, labels))
}
