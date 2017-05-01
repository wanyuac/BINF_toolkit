#!/usr/bin/env Rscript
# Merge a list of genomic regions and find out complementary regions afterwards.
# Commandline: Rscript mergeGenomicRegions.R [input file] [genome length] [output prefix]
# The input file is a CSV file comprised of two columns ("from" and "to") without a header.
# This script assumes that there are always >= 2 regions in any coordinate tables.
# A typical input is the coordinate file produced by the script filterCoords.py of RedDog.
# Example:
#   Rscript mergeGenomicRegions.R coord.txt 5248520
# Outputs:
#   (1) [input filename]__merged.csv, (2) [input filename]__comple.csv
#   Saved under the current working directory.
#
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# Edition: 30 Apr 2017

mergeRegions <- function(x) {
    y <- x[1, ]  # take the first row of x to start with
    j <- 1  # row pointer of y
    ub1 <- y$to[1]  # upper bound of the current region
    for (i in 2 : nrow(x)) {  # must guarantee there are >= 2 rows in the table
        z <- x[i, ]
        lb2 <- z$from[1]  # lower bound of the new region
        ub2 <- z$to[1]
        # There are only two behavious: either merge two regions or adding a separate region.
        if (lb2 <= (ub1 + 1)) {  # two regions overlap or adjacent: merge them into a single one. lb1 <= lb2 because the data frame is sorted in an ascending order.
            if (ub2 > ub1) {  # extend the previous region
                y$to[j] <- ub2
                ub1 <- ub2
            }  # else, do nothing as the second range is a subset of the first one
        } else {  # push a new and non-overlapping region into the stack of regions
            y <- rbind(y, z)
            j <- j + 1  # move the point to the new row
            ub1 <- ub2
        }
    }
    
    return(y)
}

findComplementaryRegions <- function(x, L) {  # L: genome size
    n <- nrow(x)  # number of predefined regions
    r <- x[1, ]
    lb <- r$from[1]
    ub <- r$to[1]
    
    # initialise z
    if (lb > 1) {  # if the first region is not at the start of the genome
        z <- data.frame(from = 1, to = lb - 1)  # lb - 1 may equal s
    } else {
        z <- data.frame(from = integer(0), to = integer(0))
    }
    s <- ub + 1

    for (i in 2 : n) {  # must guarantee there are >= 2 rows in the table
        r <- x[i, ]
        lb <- r$from[1]
        ub <- r$to[1]
        # Notice the function mergeRegions guarantees that lb > s and lb - s >= 1.
        # So the width of any gaps >= 1.
        z <- rbind(z, data.frame(from = s, to = lb - 1))  # Notice lb - 1 may equal s
        s <- ub + 1
    }
    
    # Are there any bases left beyond the last predefined region?
    if (s <= L) {
        z <- rbind(z, data.frame(from = s, to = L))
    }
    
    return(z)
}

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
genome.len <- args[2]
prefix <- ifelse(length(args) >= 3, args[3], "coords")

if (file.exists(input)) {
    x <- read.csv(input, header = FALSE)
} else {
    stop(paste("The input file", input, "is not found.", sep = " "))
}

names(x) <- c("from", "to")

# Sort beginnings of regions so that their heads do not go backwards
x <- x[order(x$from, decreasing = FALSE), ]

# Merge regions
y <- mergeRegions(x)

# Find out complementary regions
z <- findComplementaryRegions(y, genome.len)

# write results
write.table(y, file = paste0(prefix, "__merged.csv"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ",")
write.table(z, file = paste0(prefix, "__comple.csv"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ",")

# sumamrise regions
y$base <- y$to - y$from + 1
z$base <- z$to - z$from + 1

print(paste(nrow(x), "regions have been merged into", nrow(y), "regions of", sum(y$base), "bases.", sep = " "))
print(paste("There are", nrow(z), "regions of", sum(z$base), "bases outside of the merged regions.", sep = " "))
