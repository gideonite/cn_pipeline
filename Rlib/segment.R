colMedians <- function(m) {
  apply(m,2,median)
}

  
segment.normalized<-function(norm.out,resultFolder, append=NA, verbose=TRUE) {

   map=norm.out$map
   rml=norm.out$results
   
   library("DNAcopy")
   chr=as.character(map$"Ch")
   chn=match(chr,c(1:22,"X","Y"))

   pos=map$"Pos"
   
   if (!is.na(append)) {
if (verbose) { cat("appending\n") }  
      load(paste(resultFolder,"segment.out.D",sep=""))
      segment.pre=segment.out
      nm=colnames(segment.pre$rml)
      ai=which(!colnames(norm.out$results)%in%nm)
      rml=norm.out$results[,ai]
   }
   
   names=colnames(rml); nnames=make.names(names)
   
      
   # combine replicate probes
   dupProbes=unique(map$ProbeID[which(duplicated(map$ProbeID))])

   for (pid in as.character(dupProbes)) {
   	dmi=which(map$ProbeID==pid)
   	rml[dmi,]=NA
   	rml[dmi[1],]=colMedians(rml[dmi,])
   }
   results.CNA <- CNA(rml,chr,pos,data.type="logratio",sampleid=names)
   
   print("Segment: smoothing...")
   results.CNA.sm <- smooth.CNA(results.CNA)
   print("Segment: segmenting...")
   results.segment<- segment(results.CNA.sm,verbose=1,nperm=1000)
#   results.segment$chrom=chr
#   qa.seg=segmentQA(results.segment$output)   
   save(results.segment,file=paste(resultFolder,"results.segment.D",sep=""))

   mml=make.mml(rml,chr,k=3)

   d=results.segment$output
   sml=make.sml(rml,d,chr,pos)
   smlgf=fill.sml(sml)
   
   if (!is.na(append)) {
      rml=cbind(segment.pre$rml,rml)
      mml=cbind(segment.pre$mml,mml)
      sml=cbind(segment.pre$sml,sml)
      smlgf=cbind(segment.pre$smlgf,smlgf)
      
      results.segment$data=cbind(segment.pre$segment.results$data,results.segment$data)
      results.segment$output=rbind(segment.pre$segment.results$output,results.segment$output)
   }
   
   save(rml,file=paste(resultFolder,"rml.D",sep=""))
   save(mml,file=paste(resultFolder,"mml.D",sep=""))
   save(sml,file=paste(resultFolder,"sml.D",sep=""))
   save(smlgf,file=paste(resultFolder,"smlgf.D",sep=""))
   
   save(chr,file=paste(resultFolder,"chr.D",sep=""))
   save(pos,file=paste(resultFolder,"pos.D",sep=""))
   segment.list=results.segment$output
   save(segment.list,file=paste(resultFolder,"segment.list.D",sep=""))
   
   list(call=match.call(),func=segment.normalized,segment.results=results.segment,
    rml=rml,mml=mml,sml=sml,smlgf=smlgf) #,qa.seg)
}


segmentQA<- function(r) {
   r$width=r$loc.end-r$loc.start

   sc=as.data.frame(array(NA,dim=c(length(unique(r$chrom)),length(unique(r$ID)))))
   colnames(sc)=as.character(unique(r$ID))
   rownames(sc)=as.character(unique(r$chrom))
   ns=sc
   ml=c()
   for (id in unique(r$ID)) {
#   cat(id,"\n")
      r1=subset(r,ID==id)
      ml[id]=quantile(r1$width,.98)
      for (ch in as.character(unique(r1$chrom))) {
         r2=subset(r1,chrom==ch)
         n=nrow(r2)
         if (n>1) {
            d=sum(abs(diff(r2$seg.mean)))
         } else { d=1 }
         sc[ch,id]=(n-1)/d
         ns[ch,id]=(n-1)
      }
   }
   sc.qu=apply(sc,2,quantile,.1)
   ns.qu=apply(ns,2,quantile,.1)
   ns.auto=colSums(ns[1:22,])
   
#   plot(log(ml),ns.qu)
#   oi=order(ml)
#   cbind(log(ml),ns.qu,sc.qu)[rev(oi),]
   qa=rbind(frag.score=sc.qu,num.seg=ns.qu,max.length=ml,num.seg.auto=ns.auto)
   list(func=segmentQA,data=qa)
}
