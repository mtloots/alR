#' (Trimmed) vector.
#'
#' Find a vector excluding trimmed values, on which the calculation of a statistic will result in a trimmed statistic.
#'
#' @param x A numeric vector of real values.
#' @param trim The fraction (0 to 0.5) of observations to be trimmed from each end of \code{x}.
#' @param na.rm A logical value indicating whether NA values should be stripped.
#' @param ... Further arguments passed to or from other methods.
#' @return If \code{trim} is zero (the default), the values in \code{x} is returned. If \code{x} is not logical (coerced to numeric), numeric (including integer) or complex, NA_real_ is returned, with a warning.  If \code{trim} is non-zero, a symmetrically trimmed vector is returned with a fraction of \code{trim} observations deleted from each end.
#'
#' @examples
#' x <- c(rnorm(100), 100)
#' mean(x)
#' mean(x, trim=0.01)
#' mean(trim(x, trim=0.01))
#' sd(x)
#' sd(trim(x, trim=0.01))
#'
#' @export
trim <- function(x, trim=0, na.rm=FALSE, ...)
{
if(!is.numeric(x) && !is.complex(x) && !is.logical(x))
{
warning("argument is not numeric or logical: returning NA")
return(NA_real_)
}
if(na.rm) x <- x[!is.na(x)]
if(!is.numeric(trim) || length(trim) != 1)
stop("'trim' must be numeric of length one")

n <- length(x)
if(trim > 0 && n > 0)
{
if(is.complex(x)) stop("trim is not defined for complex data")
if(trim >= 0.5) return(0)
lo <- floor(n * trim) + 1
hi <- n + 1 - lo
x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
}
x
}