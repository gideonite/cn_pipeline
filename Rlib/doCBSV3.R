##
## $Source: /home/socci/lWork/CGHPipe/Rlib/doCBSV3.R $
## $Date: 2010/01/07 23:31:53 $
## $Tag: doCBSV3.R 2010/01/07 23:31:53 72fc2c193ef1 $
##

if(!exists("revision")) revision=list()
revision$doCBS = '$Id: doCBSV3.R 72fc2c193ef1 2010/01/07 23:31:53 socci $'
cat(revision$gcNorm,"\n")

library(DNAcopy)
source("Rlib/nick.R")
source("Rlib/getOpts.R")
source("Rlib/fileTools.R")
source("Rlib/cbsTools.R")

args=getopts()

#Sys.setenv("R_GSCMD"="/home/liang/Work/CGHPipe/Rlib/RGS")


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

if(is.bzfile(file)) {
  dd=read.delim(bzfile(file),as.is=T,check.names=F)
} else {
  dd=read.delim(file,as.is=T,check.names=F)
}
samp=colnames(dd)[4]
print(samp)

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
dd[,3]=as.numeric(dd[,3])
dd[,4]=as.numeric(dd[,4])

##
## Merge Probes that are probing the same gene
## 
cat("Starting Probe Merge ...")
tag=paste(dd[,2],dd[,3])
jj=which(duplicated(tag))
cnt=length(jj)
if(cnt>0) {
  dup.i=which(tag %in% tag[jj])
  dup.med=tapply(dd[dup.i,4],tag[dup.i],median,na.rm=T)
  dd[dup.i,4]=dup.med[tag[dup.i]]
  dd=dd[!duplicated(tag),]
}
cat("Finished probe merge\n")

##
## Use the correct names for the probes/rows
##
rownames(dd)=dd[,1]

cna=CNA(dd[,4],dd$Ch,dd$Pos,"logratio",samp)
rownames(cna)=dd[,1]

##
## Fix R's mangling of names
##
colnames(cna)[3]=samp

smCNA=smooth.CNA(cna)

smOut=as.matrix(smCNA)
rownames(smOut)=dd[,1]
ofile=paste(nameTmpl,"Smooth.txt",sep="__")
write.xls(smOut,ofile)

d=segment(smCNA,verbose=T,alpha=alpha,nperm=nperm,
  undo.splits='sdundo',
  undo.SD=1)

##fileName=paste(file,"Np",nperm,"alpha",alpha,"cbs.Rdata",sep="_")
d$params=list()
d$params$revision = '$Id: doCBSV3.R 72fc2c193ef1 2010/01/07 23:31:53 socci $'
d$params$segment.opts=list(alpha=alpha,nperm=nperm,undo.splits='sdundo',undo.SD=1)


source("Rlib/format.cbs.R")
p.out=segments.p(d)
ffile=paste(nameTmpl,"CBS_out.txt",sep="__")
format.cbs(d,p.out,ffile)

##
## Map segments onto genome
##
seg = rep(NA,nrow(d$data))
ssd = rep(NA,nrow(d$data))
nn  = rep(NA,nrow(d$data))
d$output$seg.sd=rep(NA,nrow(d$output))
for(j in seq(nrow(d$output))) {
  rr = (
        d$data$chrom==d$output$chrom[j]
        & d$data$maploc>=d$output$loc.start[j]
        & d$data$maploc<=d$output$loc.end[j]
        )
  seg[rr] = d$output$seg.mean[j]
  ssd[rr] = sd(d$data[rr,3],na.rm=T)
  Srr=sum(rr)
  nn[rr] = Srr
  d$output$seg.sd[j]=sd(d$data[rr,3],na.rm=T)
}

d$output$start.idx=match(
  paste(d$output$chrom,d$output$loc.start),
  paste(d$data[,1],d$data[,2]))
d$output$end.idx=c(d$output$start.idx[2:nrow(d$output)]-1,nrow(d$data))

d$map=data.frame(seg=seg,seg.sd=ssd,num=nn)

d$QC=list()
d$QC$dn=median(abs(diff(d$data[,3])),na.rm=TRUE)
d$segments.p=p.out

rfile=paste(nameTmpl,"CBS_out.Rdata",sep="__")
save(d,file=rfile,compress=T)

##
## normalize diploid peak
##
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
bitmap(plotFile,width=11,height=8.5,res=150)
###png(plotFile,width=1100,height=850)
plot(d,ylim=c(-2,2))
dev.off()

##################
##################
cat(paste("\n__","Finished__",sep=""),args$taskNo,"\t",file,"\n\n")
