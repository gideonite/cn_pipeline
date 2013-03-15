##
## $Source$
## $Date$
##
if(!exists("revision")) revision=list()
revision$agilentNorm = '$Id$'
cat(revision$agilentNorm,"\n")

library(matlab)

dprint <- function(name, klass, plotfun, ...) {

  ##fn = paste(paste(name,klass,sep="_"),".png",sep="")
  ##bitmap(fn,width=8.5,height=11,res=300)

  ##fn = paste(paste(name,klass,sep="_"),".png",sep="")
  ##png(fn,width=850,height=1100)
  
  plotfun(...)
  dev.off()
}

cgh.colors <- function(levels) {
  levels = levels%/%2
  sup = (1:(levels))/levels
  sdn = 1-(0:(levels-1))/levels
  return(c(rgb(sup,1,sup),rgb(1,sdn,sdn)))
}

is.bzfile <- function(x) {
  length(grep(".bz2$",x))>0
}

agilentNorm <- function(fname, do.GC=TRUE, verbose=TRUE,
                        do.plot=TRUE, write.gpr=FALSE) {

  ## do.GC=TRUE; verbose=TRUE; do.plot=TRUE; write.gpr=FALSE

  if (verbose) { cat("normalizing\n") }  
  qa=c()

  filename=basename(fname)
  pos=regexpr(".txt",filename)
  samplename=substr(filename,1,pos-1)
  pos=regexpr(".txt",fname)
  nameTmpl=substr(fname,1,pos-1)
  
  if(is.bzfile(fname)) {
    header=readLines(bzfile(fname),9)
  } else {
    header=readLines((fname),9)
  }
  
  fields=strsplit(header[2],"\t")[[1]]
  xmlfield=which(fields=="FeatureExtractor_DesignFileName")
    
  temp=strsplit(header[3],"\t")[[1]]
  xmlfile=temp[xmlfield]
  root=sub(".xml","",xmlfile)
  temp=paste("Maps/",root,".map.Rdata",sep="")
  if (!exists("map")) { ## do this part once in the beginning
    #
    # Hard code HYBRID 244k mapfile to force 1M to be used as 244k
    #
    mapfile="Maps/021529_D_HYBRID_244k.map.Rdata"
    if (verbose) { cat("\n loading ",mapfile,"... ") }
    if (file.exists(mapfile)) {
      load(mapfile)
    } else { 
      cat("unable to find ",mapfile,"\n")
      print(dir("Maps"))
      stop()
    }
  }

  ##
  ## First make sure the map does not have duplicate probes (fatal error now)
  ##

  if(any(duplicated(map$ProbeID))) {
    cat("\n\nDuplicated probes in map not ALLOWED\n\n")
    stop("agilentNorm::L94")
  }

  if (verbose) { cat(" reading... ") }
  
  if(is.bzfile(fname)) {
    dd=read.delim(bzfile(fname),skip=9,as.is=T,comment.char="")
  } else {
    dd=read.delim(fname,skip=9,as.is=T,comment.char="")
  }

  ##
  ## For pseudoimage
  ## Compute these on raw data before we remove unmapped and merge duplicated probes
  ##
  max(dd$Row)->nr; max(dd$Col)->nc

  ##
  ## First only use probes that are actually in file (in case we are processing
  ## less than the entire file)
  ##
  ## Also the probes must be mapped; throw out unmapped probes
  ##
  ##

  pi = intersect(map$ProbeID, dd$ProbeName)
  dd = dd[dd$ProbeName %in% pi,]
  map = map[map$ProbeID %in% pi,]

  results= matrix(NA,nrow=nrow(map),ncol=1)
  colnames(results)=samplename;
  rownames(results)=map$ProbeID

  di=match(map$ProbeID, dd$ProbeName)

  d=cbind(map,R=dd$LogRatio[di]/log10(2))
  d$A=log2(sqrt(dd$gProcessedSignal[di]*dd$rProcessedSignal[di]))
  d$GS=dd$gProcessedSignal[di]
  d$RS=dd$rProcessedSignal[di]
  d$Col=dd$Col[di]
  d$Row=dd$Row[di]
  
  ## combine flags
  t1=grep("OL$",colnames(dd))
  t2=grep("IsManualFlag",colnames(dd))
  t1=c(t1,t2)
  ddFlags=(rowSums(dd[,t1]))
  d$flags=ddFlags[di]
  

  ## deal with replicate probes
  if (verbose) { cat(" resolving replicate probes... ") }
  dupProbes=unique(dd$ProbeName[duplicated(dd$ProbeName)])
  cnt=length(dupProbes)
  cat("\n")
  for (pid in as.character(dupProbes)) {
    cat(sprintf("%6d",cnt),"\r")
    cnt=cnt-1
    
    ii=which(pid==dd$ProbeName)
    jj=which(d$ProbeID==pid)

    d$R[jj]=median(dd$LogRatio[ii])/log10(2)
    d$GS[jj]=exp(median(log(dd$gProcessedSignal[ii])))
    d$RS[jj]=exp(median(log(dd$rProcessedSignal[ii])))
    d$A[jj]=log2(sqrt(d$RS[jj]*d$GS[jj]))
    d$flags[jj]=sum(ddFlags[ii])
  }
  print("DEBUG agilentNorm::147")
  
  ## if green is always ref channel,
  ## can exclude based on green-only to preserve
  ## sensitivity to homozygous deletions in tumor
  
  ##d$dim=d$RS<100 | d$GS<100
  d$dim= d$GS<50
  
  d$R.all=d$R
  d$R[d$dim]=NA  ## clear out weak spots
  d$R[d$flags>0]=NA ## clear out flagged spots

  ##
  ## Pseduo image plot
  ##
  if(do.plot) {
    tmp=tmp.all=array(NA,dim=c(nr,nc))
    for (ii in 1:nrow(d)) {
      tmp[d$Row[ii],d$Col[ii]]=d$R[ii]
      ##        tmp.all[d$Row[ii],d$Col[ii]]=d$R.all[ii]
      if(ii %% 1000 == 0) cat(ii,"\r")
    }
    tmp[is.na(tmp)]=0
    tmp.all[is.na(tmp.all)]=0
    z=1;
    tmp[tmp>z]=z; tmp[tmp<(-z)]=-z
    ov= z+z/(51/2)
    tmp[is.na(tmp)]=ov
    print("DEBUG 134")
    
    cols<- c(cgh.colors(51),"#0000FF")
    dprint(nameTmpl,"_psimage_flags",imagesc,t(tmp),col=cols)
    cols<- c(cgh.colors(51),cgh.colors(51)[26])
    dprint(nameTmpl,"_psimage",imagesc,t(tmp),col=cols)
  }
  
  
  ## find a subset of probes for %GC normalization (secondary)
  mdi.gc=sample(1:nrow(d),min(100000,0.40*nrow(d))) 
  print("DEBUG 143")
  
  if (do.GC) {
    if (verbose) { cat("GC loess... ") } 
    R.model.gc=loess(R ~ PGC.50KB + PGC.2KB + PGC.02KB,
      data=d[mdi.gc,], na.action="na.exclude",
      control=loess.control(surface="interpolate",statistics="approximate",
        trace.hat="approximate", cell=0.2, iterations=4)
      );print("DEBUG 150")
    d$R.gcnorm=d$R-predict(R.model.gc,d);print("DEBUG 151")
    
    
  } else { d$R.gcnorm=d$R }
  print("DEBUG 155")
  
  ## find probes in "unaltered" regions for R vs I loess (primary)
  m=runmed(na.omit(d$R.gcnorm),k=51);print("DEBUG 158")
  md=density(m,bw=0.01);print("DEBUG 159")
  mx=md$x[which.max(md$y)];print("DEBUG 160") ## find peak mode of median-filtered log2 ratio
  tmp=which(abs(m-mx)<0.15);print("DEBUG 161") ## choose probes within +/- threshold=0.15 of peak
  mdi=tmp[which( abs(
    (d$R.gcnorm[tmp]-
     mean(d$R.gcnorm[tmp],na.rm=TRUE))/sd(d$R.gcnorm[tmp],na.rm=TRUE)
    )<4)];print("DEBUG 165")  ## exclude outliers where raw log2 ratio is >(n=4)SD
  
  if (verbose) { cat(length(mdi)," invariant points... ") }      
  if (verbose) { cat("primary loess... ") }      
  
  R.model.pri=loess(R.gcnorm ~ A , data=d[mdi,], na.action="na.exclude",
    control=loess.control(surface="interpolate",statistics="approximate",
      trace.hat="approximate", cell=0.2, iterations=4)
    )
  d$R.norm=d$R.gcnorm-predict(R.model.pri,d)
  results[,1]=d$R.norm
  
  ## QA
  if (verbose) { cat("QA... ") } 
  if (0 & do.plot) {
    
    newfile=paste(nameTmpl,"__Rplots.png",sep="")
    ##newfile=sub(" ","_",newfile)
    if (verbose) { cat("(plotting) ") } 

    bitmap(file=newfile,width=8.5,height=11,res=300)
    ##png(file=newfile,width=850,height=1100)

    par(mfrow=c(4,2))
    
    plot(d$PGC.50KB,d$R,main="R (unnorm) vs GC% 50KB (genomic)",     
         ylab="log2(R)",pch=".")
    abline(h=0,col=2)
    
    plot(d$PGC.50KB,d$R.gcnorm,main="R (norm) vs GC% 50KB",
         ylab="log2(R)",pch=".")
    abline(h=0,col=2)
    
    plot(d$PGC.2KB,d$R,main="R (unnorm) vs GC% 2KB (genomic)",
         ylab="log2(R)",pch=".")
    abline(h=0,col=2)
    
    plot(d$PGC.2KB,d$R.gcnorm,main="R (norm) vs GC% 2KB",
         ylab="log2(R)",pch=".")
    abline(h=0,col=2)
    
    plot(d$PGC.02KB,d$R,main="R (unnorm) vs GC% 200BP (genomic)",
         ylab="log2(R)",pch=".")
    abline(h=0,col=2)
    
    plot(d$PGC.02KB,d$R.gcnorm,main="R (norm) vs GC% 200BP",
         ylab="log2(R)",pch=".")
    abline(h=0,col=2)
    
    plot( d$A,d$R.gcnorm,main="R vs I (unnorm)",
         ylab="log2(R)",xlab="log2(intensity)",pch=".",col=1)
    points(d$A[mdi],d$R.gcnorm[mdi],pch=".",col="darkgray")
    abline(h=0,col=2)
    points(R.model.pri$fitted ~ R.model.pri$x,col=3,pch="-")
    
    plot(d$A,d$R.norm,main="R vs I (norm)",
         ylab="log2(R)",xlab="log2(intensity)",pch=".")
    abline(h=0,col=2)
    
    dev.off()
  }

  if(verbose) {cat("Gather QC_stats...")}
  RMS=sqrt(sum((d$R.norm-d$R)^2,na.rm=T) ) 
  ##  print(RMS)
  
  ii=grep("IsSat",colnames(dd))
  per.sat = 100*sum(dd[di,ii]>0,na.rm=T)/(2*length(di))
  
  ii=grep("FeatNon",colnames(dd)) 
  feat.nonu = 100*sum(dd[di,ii]>0,na.rm=T)/(2*length(di))
  
  ii=grep("BGNon",colnames(dd)) 
  bg.nonu = 100*sum(dd[di,ii]>0,na.rm=T)/(2*length(di))
  
  sddr = sd(diff(d$R.norm),na.rm=T)/sqrt(2)
  
  per.lt100 = 100*sum(d$A<log(100,2),na.rm=T)/length(d$A)
  per.lt300 = 100*sum(d$A<log(300,2),na.rm=T)/length(d$A)
  
  PGC.50KB.cor=cor(d$R.norm,map$PGC.50KB,method="spearman",use="complete.obs")
  PGC.2KB.cor=cor(d$R.norm,map$PGC.2KB,method="spearman",use="complete.obs")
  PGC.02KB.cor=cor(d$R.norm,map$PGC.02KB,method="spearman",use="complete.obs")
  
  stats = c(fileName=fname,
    date=date(),round(c(per.sat=per.sat,feat.nonu=feat.nonu,
      bg.nonu=bg.nonu,RMS=RMS,PGC.50KB.cor=PGC.50KB.cor,
      PGC.2KB.cor=PGC.2KB.cor,PGC.02KB.cor=PGC.02KB.cor,
      sd.dv=sddr,per.lt100=per.lt100,per.lt300=per.lt300),3)
    )
  
  if(0) {
    ##
    ## No longer dumping stats file
    ##
    newfile=paste(nameTmpl,"__STATS.txt",sep="")
    cat(stats,"\n",file=newfile,sep="\t")
  }
  qa=cbind(qa,stats)
  rownames(qa)=  c("fileName","date","per.sat.G","per.sat.R",
            "feat.nonu.G","feat.nonu.R",
            "bg.nonu.G","bg.nonu.R",
            "RMS","PGC.50KB.cor","PGC.2KB.cor","PGC.02KB.cor","sd.dv",
            "per.lt100","per.lt300"
            )    
  
  if (write.gpr) {
    if(verbose) {cat("Writing GPR File ...")}
    dii=which(!is.na(di))
    dd$LogRatio[di[dii]]=d$R.norm[dii]/log2(10)
    
    newfile=paste(nameTmpl,"__GCNormV4b.txt",sep="")
    
    ff=file(newfile,"w")
    
    for(line in header) {
      cat(line,"\n",file=ff,sep="")
    }
    
    write.table(dd,file=ff,sep="\t",quote=F,col.names=T,row.names=F)
    
    close(ff)
  }
  
  if (verbose) { cat("\n") } 
  
  
  colnames(qa)<-samplename
  
  list(call=match.call(),
       func=agilentNorm,map=map,
       results=results,qa=qa)
}
