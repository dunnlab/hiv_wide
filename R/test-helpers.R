testtree <- function(dir = file_temp(), env=parent.frame()) {
    treestring <- "(((a:0.08,b:0.08)95/95:0.009,c:0.05)87/87:0.01,((d:0.09,e:0.05)100/100:0.006,(f:0.09,g:0.05)100/100:0.005)95/95:0.004)100/100:0.002;"
    f <- tempfile()
    write(treestring, file=f)
    testtree <- treeio::read.iqtree(f)
    return(testtree)
}