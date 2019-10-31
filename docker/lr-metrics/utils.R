suppressPackageStartupMessages(library(tidyverse))
library(tidyverse)
library(ggridges)

options(scipen=8)

# Simple logging function
logger <- function(text) {
    cat(timestamp(prefix = "[", suffix = "]", quiet=TRUE), text, "\n")
}

# Compute NX (N50, N90, etc.).  The vector of lengths, l, must be sorted in decreasing order.
nX <- function(l, X) {
    v = l[last(which(cumsum(as.numeric(l)) <= (X/100)*sum(as.numeric(l))))]
    return(v)
}

# Extremely simple argument parser
parse_args <- function(default_args) {
    if (interactive()) {
        args = default_args
    } else {
        args = commandArgs(trailingOnly=TRUE)
    }
    
    l = list()
    l[["files"]] = c()
    
    for (i in 1:length(args)) {
        if (i > 1 && grepl("^--", args[i-1])) {
            next
        } else if (grepl("^--", args[i])) {
            key = gsub("^--", "", args[i])
            val = args[i+1]
            
            if (key %in% names(l)) {
                l[[key]] = c(l[[key]], val)
            } else {
                l[[key]] = c(val)
            }
            
            i = i+1
        } else {
            l[["files"]] = c(l[["files"]], args[i])
        }
    }
    
    return(l)
}

colors = c( "#DF5A49", "#EFC94C", "#45B29D", "#334D5C", "#E27A3F", "#DF5A49", "#EFC94C", "#45B29D", "#334D5C", "#E27A3F" )
ltys = c( rep(1, 5), rep(2, 5) )
pchs = c( 16, 18, 17, 15, 8, 7, 6, 5, 4, 3 )

devs = c("pdf", "png")
if (interactive()) { devs = c("screen", devs) }
