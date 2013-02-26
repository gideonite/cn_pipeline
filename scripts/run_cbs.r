#R CMD BATCH --no-save "--args path/to/DNAcopy in_file out_file" /path/to/run_cbs.r cbs.log

args <- commandArgs(trailingOnly = TRUE);

if (length(args) != 3) {
    cat("usage: run_cbr.r path/to/DNAcopy cbs_in_file cbs_out_file");
    q(-1);
}

dna_copy <- args[1];
in_file <- args[2];
out_file <- args[3];

cbs_in <- read.table(in_file, header=TRUE);

library('DNAcopy', lib.loc=dna_copy);

CNA.object = CNA(cbs_in$signal, cbs_in$chr, cbs_in$pos, data.type="logratio", sampleid=basename(in_file));
smoothed.CNA.object <- smooth.CNA(CNA.object);
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1);

write.table(segment.smoothed.CNA.object$out, row.names=FALSE, quote=FALSE,
sep="\t", file=paste(out_file, sep=""));
