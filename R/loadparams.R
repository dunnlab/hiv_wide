library(stringr)

#' Loads trees, creates phylos object as well as clade labels of all of the trees
#'
#' @param treepath Path to the tree files
loadparams <- function(treepath) {
    iqtreefiles <- list.files(treepath,".iqtree",full.names=TRUE)
    dfs <- lapply(iqtreefiles, paramsdf)

    # names
    maskpercent <- str_extract_all(iqtreefiles, "mask[0-9]+") %>% str_extract_all("[0-9]+") %>% as.numeric()
    wgspercent <- 100 - maskpercent
    names <- paste0("wgs",wgspercent)
    if("wgs0" %in% names) {
        names[which(names=="wgs0")] <- "pol"
    }
    names(dfs) <- names
    return(dfs)
}

# Qij = Rij*pij
# https://groups.google.com/g/iqtree/c/dZjx6sf4-OA?pli=1
# sum of rows of Q matrix need to be 0, so Qii is always negative
# Q matrix is normalized so that -sum(k=1)^4 pi_i Q_ii = 1
# https://en.wikipedia.org/wiki/Substitution_model
paramsdf <- function(file) {
    lines <- readLines(file)
    params <- str_extract_all(lines[c(40:45,49:52,63)], "[0-9]+.[0-9]+") %>%
        flatten() %>% as.numeric()
    relative_rates <- lapply(lines[66:70] %>% str_split("\\s+"), function(x) x[3:4]) %>%
        flatten() %>% as.numeric()
    type <- c(rep("Rate parameter R",6),rep("State freq",4),"Gamma shape",rep(c("relative rate","rr proportion"),5))
    paramname <- c("A-C","A-G","A-T","C-G","C-T","G-T",
                   "pi(A)","pi(C)","pi(G)","pi(T)",
                   "alpha",
                   "0","0",
                   "1", "1",
                   "2", "2",
                   "3", "3",
                   "4", "4")
    params_df <- data.frame("Type"=type,"Param"=paramname,"Value"=c(params,relative_rates))
    return(params_df)
}
