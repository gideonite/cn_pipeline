is.bzfile <- function(x) {
  length(grep(".bz2$",x))>0
}

get.diploid.peak <- function(pd,Nf=8) {
  ####
  ## Density inside the bandwidth of the provided distribution.
  ii=which(pd$x>=-Nf*pd$bw & pd$x<=Nf*pd$bw); 

  #### Detect peaks inside provided density.
  ##
  pk=try(peaksign(pd$y[ii],11),silent=TRUE); 
  ## Fewer points of density inside bandwidth than the span provided
  if (class(pk)=="try-error") {
    return(ii[which.max(pd$y[ii])])
  }
  pc=sum(pk==1); ## Count of peaks detected (1 is good, >1 is bad, 0 is conditional).
  ## Simple/best case, single peak detected, return it.
  if (pc==1) {
        return(ii[which(pk==1)])
  }
  ## Worst case, multiple peaks, pick max and return it, but flag
  else if (pc>1) {
    cat("-*-*-*-*-- Diploid peak ambiguity, review! --*-*-*-*-\n")	
    return(ii[which(pk==1)[which.max(pd$y[ii[which(pk==1)]])]])
  }
  ## No peak detected (either too high a span or flattened by
  ## bandwidth, pick max and return it.
  else {
    return(ii[which.max(pd$y[ii])])
  }	
}

