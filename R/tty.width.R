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

