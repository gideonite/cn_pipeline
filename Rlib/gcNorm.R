##
## $Source: /home/socci/lWork/CGHPipe/Rlib/gcNorm.R $
## $Date: 2010/03/29 19:51:02 $
##
options(stringsAsFactors=FALSE)

if(!exists("revision")) revision=list()
revision$agilentNorm = '$Id: gcNorm.R d3fb5b1859c5 2010/03/29 19:51:02 socci $'
cat(revision$gcNorm,"\n")

source("Rlib/nick.R")
source("Rlib/getOpts.R")
source("Rlib/agilentNorm.R")
source("Rlib/fileTools.R")

source("Rlib/parseINI.R")

glb=parse.INI("gcNorm.ini")

cat("gbl$Options$write.gpr =",as.logical(glb$Options$write.gpr),"\n")
cat("gbl$Options$default.grid =",(glb$Options$default.grid),"\n")

args=getopts()

if(0) {
slp=runif(1)*60*5
cat("Sleeping for ",slp/60," min ...")
Sys.sleep(slp)
cat("Ok\n")

x=as.numeric(system("ssh tobiko w | fgrep load | awk '{print $10}'|tr -d ','",intern=T))
cat("X==",x,"\n")
while(x>10) {
  cat("NSF Overload <",x,">\t")
  slp=runif(1)*30*x+x*10
  cat("Sleeping for ",slp," sec ...")
  Sys.sleep(slp)
  x=as.numeric(system("ssh tobiko w | fgrep load | awk '{print $10}'|tr -d ','",intern=T))
  cat("X2==",x,"\n")
}

progName="gcNorm.R"

count.active <- function(progName) {
  xx=system(paste("ps aux |fgrep -v grep|fgrep exec| fgrep",progName),intern=T)
  cnt=0
  for(x in xx) {
    ss=strsplit(x," ")
    si=ss[[1]][ss[[1]]!=""]
    if(si[8]=="R") cnt = cnt + 1
  }
  return(cnt)
}

active=count.active(progName)
while(active>0) {
  slp=runif(1)*60
  cat("     Collison ... sleeping for ",slp," sec ...\n")
  Sys.sleep(slp)
  active=count.active(progName) 
}
}


cat("Starting...\n")


if(!is.null(args$taskNo)) {
  ##
  ## Running on cluster
  ##
  files=scan(args$tasks,"",sep="\n")
  cat("Task No:", args$taskNo, "\n")
  file=files[as.numeric(args$taskNo)]
} else {
  file=args$file
}
cat("File:", file, "\n")

res=agilentNorm(file,do.plot=TRUE,write.gpr=as.logical(glb$Options$write.gpr),default.grid=glb$Options$default.grid)

pos=regexpr(".txt",file)
nameTmpl=substr(file,1,pos-1)
nfile=paste(nameTmpl,"GCN_V3.txt",sep="__")

ans=data.frame(Ch=res$map$Ch,Pos=res$map$Pos,res$results[,1,drop=FALSE])
rownames(ans)=res$map$ProbeID

ii=grep("_random",ans$Ch)
if(length(ii)>0) {
  ans=ans[-ii,]
}

Ch.Num=ans$Ch
if(max(as.numeric(unique(Ch.Num)),na.rm=T)==19) {
  ## Mouse
  Ch.Num[Ch.Num=="X"]="20.1"
  Ch.Num[Ch.Num=="Y"]="20.2"
  Ch.Num[Ch.Num=="M"]="20.3"
  Ch.Num=as.numeric(Ch.Num)
} else {
  Ch.Num=as.numeric(factor(((Ch.Num)),levels=c(1:22,"X","Y","M")))
}

ii=order(ans[,"Pos"])
jj=order(Ch.Num[ii])
ans=ans[ii[jj],]

write.xls(ans,nfile)

cat("\n###Finished###\t",nfile,"\n\n")
