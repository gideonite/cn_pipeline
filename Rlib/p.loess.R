function (object, newdata = NULL, se = FALSE, ...) 
{
    if (!inherits(object, "loess")) 
        stop("first argument must be a \"loess\" object")
    if (is.null(newdata) & (se == FALSE)) 
        return(fitted(object))
    if (is.null(newdata)) 
        newx <- object$x
    else {
        vars <- as.character(attr(delete.response(terms(object)), 
            "variables"))[-1]
        newx <- if (length(vars) > 1 || NCOL(newdata) > 1) {
            if (any(!match(vars, colnames(newdata), FALSE))) 
                stop("'newdata' does not contain the variables needed")
            as.matrix(newdata[, vars, drop = FALSE])
        }
        else as.matrix(newdata)
    }
    res <- stats:::predLoess(object$y, object$x, newx, object$s, object$weights, 
        object$pars$robust, object$pars$span, object$pars$degree, 
        object$pars$normalize, object$pars$parametric, object$pars$drop.square, 
        object$pars$surface, object$pars$cell, object$pars$family, 
        object$kd, object$divisor, se = se)
    if (!is.null(out.attrs <- attr(newdata, "out.attrs"))) {
        if (se) {
            res$fit <- array(res$fit, out.attrs$dim, out.attrs$dimnames)
            res$se.fit <- array(res$se.fit, out.attrs$dim, out.attrs$dimnames)
        }
        else res <- array(res, out.attrs$dim, out.attrs$dimnames)
    }
    if (se) 
        res$df <- object$one.delta^2/object$two.delta
    res
}
