##
## $Source: /home/socci/lWork/CGHPipe/Rlib/doCBSV2_999.R $
## $Date: 2009/10/13 20:43:31 $
## $Tag: doCBSV2_999.R 2009/10/13 20:43:31 46093441262b $
##

if(!exists("revision")) revision=list()
revision$doCBS = '$Id: doCBSV2_999.R 46093441262b 2009/10/13 20:43:31 socci $'
cat(revision$gcNorm,"\n")

library(DNAcopy)
source("~/lWork/Rlib/getOpts.R")
source("~/lWork/Rlib/fileTools.R")
source("Rlib/cbsTools.R")

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

source("Rlib/agilent_070421.R")

if(is.bzfile(file)) {
  dd=read.Agilent(bzfile(file))
} else {
  dd=read.Agilent(file)
}

samp=path[len(path)]
print(samp)

##
## Get rid of probes mapped to random chromosomes
##

ir=grep("_random",dd$Chr)
if(length(ir)>0) {
  dd=dd[-ir,]
}

if(max(as.numeric(unique(dd$Chr)),na.rm=T)==19) {
  ## Mouse
  dd$Chr[dd$Chr=="X"]="20.1"
  dd$Chr[dd$Chr=="Y"]="20.2"
  dd$Chr[dd$Chr=="M"]="20.3"
  dd$Chr=as.numeric(dd$Chr)
} else {
  dd$Chr=as.numeric(factor(((gsub("[^0-9XY]","",dd$Chr))),levels=c(1:22,"X","Y")))
}

dd=dd[!is.na(dd$Chr),]
dd$R=as.numeric(dd$R)


cna=CNA(dd$R,dd$Chr,dd$Pos,"logratio",samp)

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
d$params$revision = '$Id: doCBSV2_999.R 46093441262b 2009/10/13 20:43:31 socci $'
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
##png(plotFile,width=1100,height=850)
plot(d,ylim=c(-2,2))
dev.off()

##################
##################
cat(paste("\n__","Finished__",sep=""),args$taskNo,"\t",file,"\n\n")
