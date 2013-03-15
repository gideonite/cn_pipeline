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
cat("Starting...\n")
