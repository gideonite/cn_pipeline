write.xls <- function(dd,filename,row.names=T,col.names=NA) {
  if (!is.data.frame(dd)) {
    dd <- data.frame(dd,check.names=F)
  }
  if(!row.names) {
    col.names=T
  }
  write.table(dd,file=filename,sep="\t",quote=F,
              col.names=col.names,row.names=row.names)
}

readAffy <- function(filename,nrows=-1,del.affy=T) {
  header=readLines(filename,n=1)
  desc.col=length(strsplit(header,"\t")[[1]])
  
  dd=read.delim(filename,row.names=1,check.names=F,as.is=desc.col,
    nrows=nrows,comment.char="")

  if(del.affy) {
    ig = -grep("AFFX",rownames(dd))
  } else {
    ig = T
  }
  
  names=colnames(dd)
  
  ii=grep("Signal$",names)
  ds=as.matrix(dd[ig,ii,drop=F])
  ##
  ## Get rid of zeros before taking logs
  ## set zeros to 1/2 smallest value
  minVal = range(ds[ds>0])[1]/2
  ds[ds==0]=minVal
  ds=log2(ds)
  
  pos=regexpr("_Signal",colnames(ds))
  colnames(ds)=substr(colnames(ds),1,pos-1)
  
  ii=grep("Detection$",names)
  dp=as.matrix(dd[ig,ii,drop=F])
  pos=regexpr("_Detection",colnames(dp))
  colnames(dp)=substr(colnames(dp),1,pos-1)
  
  ii=grep("Signal Log",names)
  if(length(ii)>0) {
    dl=dd[ig,ii,drop=F]
    pos=regexpr("_Signal",colnames(dl))
    colnames(dl)=substr(colnames(dl),1,pos-1)
    dl=as.matrix(dl[,!is.na(dl[1,])])
    
    ii=grep("Change$",names)
    dc=dd[ig,ii,drop=F]
    pos=regexpr("_Change",colnames(dc))
    colnames(dc)=substr(colnames(dc),1,pos-1)
    dc=as.matrix(dc[,!is.na(dc[1,])])
    
    return(list(ds=ds,dp=dp,dl=dl,dc=dc))
  }

  return(list(ds=ds,dp=dp))
  
}

readAffy.SNP <- function(filename,nrows=-1) {

  dd=read.delim(filename,row.names=1,as.is=2,check.names=F,nrows=nrows,skip=1,comment.char="")
  names=colnames(dd)

  probeID = dd[,1]

  chr=dd[,2]
  names(chr)=probeID


  ##
  ## Get the chromosome to sort correctly
  ##
  chr.orig=chr
  chr=as.numeric(factor(chr,as.character(unique(chr)),ordered=T))

  cpos=dd[,3]
  names(cpos)=probeID

  cname="_Call"
  ii=grep(cname,names)
  call=as.matrix(dd[,ii])
  pos=regexpr(cname,colnames(call))
  colnames(call)=substr(colnames(call),1,pos-1)
  rownames(call)=probeID

  cname="_SPA_CN"
  ii=grep(cname,names)
  spa.cn=as.matrix(dd[,ii])
  pos=regexpr(cname,colnames(spa.cn))
  colnames(spa.cn)=substr(colnames(spa.cn),1,pos-1)
  rownames(spa.cn)=probeID

  cname="_SPA_pVal"
  ii=grep(cname,names)
  spa.pv=as.matrix(dd[,ii])
  pos=regexpr(cname,colnames(spa.pv))
  colnames(spa.pv)=substr(colnames(spa.pv),1,pos-1)
  rownames(spa.pv)=probeID

  cname="_GSA_CN"
  ii=grep(cname,names)
  gsa.cn=as.matrix(dd[,ii])
  pos=regexpr(cname,colnames(gsa.cn))
  colnames(gsa.cn)=substr(colnames(gsa.cn),1,pos-1)
  rownames(gsa.cn)=probeID

  cname="_GSA_pVal"
  ii=grep(cname,names)
  gsa.pv=as.matrix(dd[,ii])
  pos=regexpr(cname,colnames(gsa.pv))
  colnames(gsa.pv)=substr(colnames(gsa.pv),1,pos-1)
  rownames(gsa.pv)=probeID

  cname="_LOH"
  ii=grep(cname,names)
  loh=as.matrix(dd[,ii])
  pos=regexpr(cname,colnames(loh))
  colnames(loh)=substr(colnames(loh),1,pos-1)
  rownames(loh)=probeID

  return(list(chr=chr,chr.orig=chr.orig,pos=cpos,call=call,
              spa.cn=spa.cn,spa.pv=spa.pv,
              gsa.cn=gsa.cn,gsa.pv=gsa.pv,
              loh=loh))
}


read.matrix <- function(file,sep="\t",quote='"',comment.char="#",
                        dec=".",fill=TRUE,skip=0,nrows=-1) {
  
  if (is.character(file)) {
    file <- file(file, "r")
    on.exit(close(file))
  }
  if (!inherits(file, "connection")) 
    stop("'file' must be a character string or connection")
  if (!isOpen(file)) {
    open(file, "r")
    on.exit(close(file))
  }
  
  if (skip > 0) 
    readLines(file, skip)
  
  ##
  ##
  ## Get Header
  ##
  ##
  headerLine="#"
  while(substr(headerLine,1,1)=="#") {
    headerLine = readLines(file,1)
  }
  
  header=strsplit(headerLine,"\t")[[1]]
  cols=length(header)
  
  data=scan(file,what=0,dec=dec,quiet=T,comment.char=comment.char)

  data=matrix(data,ncol=cols,byrow=T)
  rownames(data)=data[,1]
  data=data[,-1]
  colnames(data)=header[2:length(header)]

  return(data)

}

