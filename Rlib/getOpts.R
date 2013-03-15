##
## $Source$
## $Date$
##

##
## Process R's twisted cmd line
##



getopts <- function() {
  print(commandArgs())
  raw.args=commandArgs()
  args = list(argv=raw.args[1])
  i.arg=grep("--args",raw.args)
  if(length(i.arg)>0) {    
    args$ARGV = substr(raw.args[i.arg:length(raw.args)],3,999)
    ii=grep("=",args$ARGV)
    args.split=strsplit(args$ARGV[ii],"=")
    for(pair in args.split) {
      print(pair)
      args[[pair[1]]]=pair[2]
    }
  }
  return(args)
}

