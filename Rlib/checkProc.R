progName="gcNorm.R"
xx=system(paste("ps aux |fgrep -v grep| fgrep",progName),intern=T)
while(length(xx)>0) {
  slp=runif(1)*60
  cat("     Collison ... sleeping for ",slp," sec ...\n")
  Sys.sleep(slp)
  xx=system(paste("ps aux |fgrep -v grep| fgrep",progName),intern=T)
}
cat("Starting...\n")
