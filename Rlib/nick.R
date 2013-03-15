cc <-
function(...) {
    paste(...,sep='_')
}

cq <-
function(...) {
    paste(...,sep='')
}

DATE <-
function() {
    gsub("-","",Sys.Date())
}

gexp2 <-
function(x){
  return(sinh(log(2)*x+log(2)))
}

glog2 <-
function(x,p0=0,p1=1) {
  return((asinh(p0+p1*x)-log(2*p1))/log(2))
}

ldensity <-
function(x){
 dx=density(log(x))
 dx$x=exp(dx$x)
 dx$y=dx$y/sum(dx$y)
 return(dx)
}

len <-
function(x) {
    length(x)
}

tty.width <-
function() {
   con=pipe("~/bin/getTTYWidth")
   dat=readLines(con,n=1)
   close(con)
   width=as.numeric(dat)
   if(len(width)>0 && width>80) {
     return(width)
   } else {
     return(80)
   }
}

write.xls <- function(dd,filename,row.names=T,col.names=NA) {
  if (!is.data.frame(dd)) {
    dd <- data.frame(dd,check.names=F)
  }
  if(!row.names) {
    col.names=T
  }
  write.table(dd,file=filename,sep="\t",quote=FALSE,
              col.names=col.names,row.names=row.names)
}
