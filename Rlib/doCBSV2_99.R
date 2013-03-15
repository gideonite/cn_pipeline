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
  files=scan("files","",sep="\n")
  cat("Task No:", args$taskNo, "\n")
  file=files[as.numeric(args$taskNo)]
} else {
  file=args$file
}

cat("File:", file, "\n")
pos=regexpr(".txt",file)
nameTmpl=paste(substr(file,1,pos-1),"A1",sep="_")

##################
##################

alpha=0.01
nperm=10000

####
####

path=strsplit(file,"/")[[1]]
samp=strsplit(path[length(path)],"_")[[1]]
samp=paste(samp[1:3],sep="",collapse="_")
print(samp)

source("agilent_070421.R")

if(is.bzfile(file)) {
  dd=read.Agilent(bzfile(file))
} else {
  dd=read.Agilent(file)
}

jj=dd$Chr %in% c("X","Y","M")
dd$Ch=dd$Chr
dd$Ch[!jj]=as.numeric(dd$Chr[!jj])
dd$Pos=floor(dd$Pos)

##
## Get rid of probes mapped to random chromosomes
##

ir=grep("_random",dd$Ch)
if(length(ir)>0) {
  dd=dd[-ir,]
}

if(max(as.numeric(unique(dd$Ch)),na.rm=T)==19) {
  ## Mouse
  dd$Ch[dd$Ch=="X"]="20.1"
  dd$Ch[dd$Ch=="Y"]="20.2"
  dd$Ch[dd$Ch=="M"]="20.3"
  dd$Ch=as.numeric(dd$Ch)
} else {
  dd$Ch=as.numeric(factor(((dd$Ch)),levels=c(1:22,"X","Y")))
}

dd=dd[!is.na(dd$Ch),]

##
## Merge Probes that are probing the same gene
## 
tag=paste(dd[,2],dd[,3])
jj=which(duplicated(tag))
cnt=length(jj)
while(cnt>0) {
  cat(cnt,"\n")
  ii=which(tag==tag[jj[1]])
  dd[ii[1],4]=median(dd[ii,4],na.rm=T)
  dd[ii[1],1]=paste(dd[ii,1],sep="",collapse=':')
  dd=dd[-ii[2:length(ii)],]
  tag=paste(dd[,2],dd[,3])
  jj=which(duplicated(tag))
  cnt=length(jj)  
}

cna=CNA(dd$R,dd$Ch,dd$Pos,"logratio",samp)

smCNA=smooth.CNA(cna)

smOut=as.matrix(smCNA)
rownames(smOut)=dd[,1]
ofile=paste(nameTmpl,"Smooth.txt",sep="__")
write.xls(smOut,ofile)

d=segment(smCNA,verbose=T,alpha=alpha,nperm=nperm,
  undo.splits='sdundo',
  undo.SD=1)

##fileName=paste(file,"Np",nperm,"alpha",alpha,"cbs.Rdata",sep="_")
rfile=paste(nameTmpl,"CBS_out.Rdata",sep="__")
save(d,file=rfile,compress=T)

source("format.cbs.R")
p.out=segments.p(d)
ffile=paste(nameTmpl,"CBS_out.txt",sep="__")
format.cbs(d,p.out,ffile)

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

ans=cbind(dd[,2:3],seg)
colnames(ans)[3]=samp
rownames(ans)=dd[,1]
sfile=paste(nameTmpl,"NormSegs.txt",sep="__")
write.xls(ans,sfile)

plotFile=paste(nameTmpl,"CBSPlot.png",sep="__")
bitmap(plotFile,width=11,height=8.5,res=100)
plot(d,ylim=c(-2,2))
dev.off()

##################
##################
cat(paste("\n__","Finished__",sep=""),args$taskNo,"\t",file,"\n\n")
