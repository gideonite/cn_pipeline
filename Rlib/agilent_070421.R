read.Agilent <- function(file,nrows=-1,rm.dups=TRUE) {

  dd=read.delim(file,skip=9,as.is=T,nrows=nrows)
  
  ##
  ## Non-control spots
  ##
  ii=dd$ControlType==0
  dd=dd[ii,]

  posInfo="SystematicName"

  i.mapped=grep("^chr([0-9]*|[XYxy]):",dd[,posInfo])

  dd=dd[i.mapped,]
    
  tag=strsplit(dd[,posInfo],":")
  
  chr=substr(sapply(tag,'[',1),4,6)
  
  rr=sapply(tag,'[',2)
  pp=regexpr("-",rr[1])
  start=as.numeric(substr(rr,1,pp-1))
  end=as.numeric(substr(rr,pp+1,99))
  pos=(start+end)/2
  
  d=dd[,c("ProbeName","LogRatio","gProcessedSignal","rProcessedSignal",
    "gIsSaturated","rIsSaturated","gIsFeatNonUnifOL","rIsFeatNonUnifOL",
    "gIsBGNonUnifOL","rIsBGNonUnifOL")]
  
  d$Chr=chr
  ii=nchar(d$Chr)==1 & !(d$Chr %in% c("X","Y"))
  d$Chr[ii]=paste(" ",d$Chr[ii],sep="")

  d$Start=start
  d$End=end
  d$Pos=pos
  
  colnames(d)[2]="R"
  d$R=d$R/log10(2)
  d$A=log2(d$gProcessedSignal*d$rProcessedSignal)/2

  chr.n=chr
  chr.n[chr.n=="X"]=23
  chr.n[chr.n=="Y"]=24
  chr.n=as.numeric(chr.n)

  pp=order(d$Start)
  pj=order(chr.n[pp])
  d = (d[pp,])[pj,]


  if(rm.dups) {
    dup.probes=names(which(table(d$ProbeName)>1))
    
    if(!is.null(dup.probes)) {
      for(i in seq(dup.probes)) {
        
        ii=which(d$ProbeName==dup.probes[i])
        
        feat1=c("R","gIsSaturated","rIsSaturated",
          "gIsFeatNonUnifOL","rIsFeatNonUnifOL",
          "gIsBGNonUnifOL","rIsBGNonUnifOL","A")
        
        feat2=c("gProcessedSignal","rProcessedSignal")
        
        d[ii[1],feat1]= mean(d[ii,feat1])
        d[ii[1],feat2]= 2^mean(log2(d[ii,feat2]))
        d=d[-ii[-1],]
      }
    }
  }

  return(d)

}



