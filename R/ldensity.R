ldensity <-
function(x){
 dx=density(log(x))
 dx$x=exp(dx$x)
 dx$y=dx$y/sum(dx$y)
 return(dx)
}

