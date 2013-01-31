args <- commandArgs(trailingOnly = TRUE);

if (length(args) != 1) {
    cat("usage: cbs_in file");
}

cbs_in <- read.table(args, header=TRUE);

library('DNAcopy');

CNA.object = CNA(cbs_in$signal, cbs_in$chr, cbs_in$pos, data.type="logratio", sampleid=args);
smoothed.CNA.object <- smooth.CNA(CNA.object);
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1);

write(segment.smoothed.CNA.object, file=paste("cbs/out/", args, sep=""), sep = "\t")
