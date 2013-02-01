#R CMD BATCH --no-save "--args ../data/AMAZE_p_TCGASNP_b86_87_88_N_GenomeWideSNP_6_A11_735446.tangent.copynumber.data.txt.cbs_in.clean /Users/dresdneg/Dropbox/dev/pipeline/cbs/out/" ../scripts/run_cbs.r cbs.log

args <- commandArgs(trailingOnly = TRUE);

if (length(args) != 2) {
    cat("usage: run_cbr.r cbs_in_file cbs_out_dir");
    q(-1);
}

in_file <- args[1];
out_dir <- args[2];

cbs_in <- read.table(in_file, header=TRUE);

library('DNAcopy');

CNA.object = CNA(cbs_in$signal, cbs_in$chr, cbs_in$pos, data.type="logratio", sampleid=basename(in_file));
smoothed.CNA.object <- smooth.CNA(CNA.object);
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1);

write.table(segment.smoothed.CNA.object$out, row.names=FALSE, quote=FALSE, file=paste(out_dir, basename(in_file), ".cbs_out", sep=""), sep = "\t")
