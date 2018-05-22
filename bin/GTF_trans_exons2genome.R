#! /usr/bin/env Rscript

coor_convert <- function(x) {
    end <- cumsum(x)
    start <- c(0, end[-length(end)])
    coor <- cbind(start, end)
    return(coor)
}

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

tab <- read.table(input, sep='\t')
coor <- do.call(rbind, tapply(tab$V5, tab$V4, coor_convert))
tab2 <- cbind(tab, coor)

tab2 <- tab2[, c(4, 10, 11, 5, 7, 6, 8, 9, 1, 2, 3)]

write.table(tab2, output, quote=F, col.names=F, row.names=F, sep='\t')
