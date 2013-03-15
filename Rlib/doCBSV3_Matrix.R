##
## $Source$
## $Date$
##
if(!exists("revision")) revision=list()
revision$doCBS = '$Id$'
cat(revision$gcNorm,"\n")

library(DNAcopy)
source("~/lWork/Rlib/getOpts.R")
source("~/lWork/Rlib/fileTools.R")
source("cbsTools.R")

args=getopts()

if(!is.null(args$taskNo)) {
  ##
  ## Running on cluster
  ##
  colNum=as.numeric(args$taskNo)
} else {
  colNum=1
}
file="ML_121208.Rdata"
cat("File:", file, "\n")

##################
##################

alpha=0.01
nperm=10000

####
#### SNP file from YL pipeline
####

load(file)

dd=result
rm(result)

##dd=dd[sample(nrow(dd),10000),]

dd$Pos=as.numeric(as.character(dd$Pos))
dd$Chr=as.numeric(as.character(dd$Chr))

pp=order(dd$Chr,dd$Pos)
dd=dd[pp,]

cna=CNA(dd[,colNum+3],dd$Chr,dd$Pos,"logratio",colnames(dd)[colNum+3])

smCNA=smooth.CNA(cna)

d=segment(smCNA,verbose=T,alpha=alpha,nperm=nperm,
  undo.splits='sdundo',
  undo.SD=1)

rfile=paste(file,colNum,"CBS_out.Rdata",sep="__")
save(d,file=rfile,compress=T)

source("format.cbs.R")
p.out=segments.p(d)
ffile=paste(file,colNum,"CBS_out.txt",sep="__")
format.cbs(d,p.out,ffile)

d$output=p.out
ffile=paste(file,colNum,"CBS_iCNA.Rdata",sep="_")
save(d,file=ffile,compress=T)

##
## Map segments onto genome
## normalize diploid peak
##
seg = double(nrow(d$data))
ssd = double(nrow(d$data))
d$output$seg.sd=double(nrow(d$output))
for(j in seq(nrow(d$output))) {
  rr = (
        d$data$chrom==d$output$chrom[j]
        & d$data$maploc>=d$output$loc.start[j]
        & d$data$maploc<=d$output$loc.end[j]
        )
  seg[rr] = d$output$seg.mean[j]
  ssd[rr] = sd(d$data[rr,3],na.rm=T)
  d$output$seg.sd[j]=sd(d$data[rr,3],na.rm=T)
}

DN=median(abs(diff(dd[,4])),na.rm=T)
P.norm=density(seg,bw=DN/8,from=-DN/2,to=DN/2,na.rm=T)
pk=get.diploid.peak(P.norm)
seg=seg-P.norm$x[pk]
d$output$seg.mean.orig=d$output$seg.mean
d$output$seg.mean=d$output$seg.mean-P.norm$x[pk]

save(d,p.out,file=rfile,compress=T)

##################
##################
cat(paste("\n__","Finished__",sep=""),args$taskNo,"\t",file,"\n\n")
