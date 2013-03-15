glog2 <-
function(x,p0=0,p1=1) {
  return((asinh(p0+p1*x)-log(2*p1))/log(2))
}

