args <- commandArgs(trailingOnly = TRUE);

if (length(args) != 1) {
    cat("usage: cbs_in file");
}

cbs_in <- read.table(args, header=TRUE);
