#' @export
safe <- function(expr, errval = NA) {
  tryCatch( {
    expr
  }, error = function(e) {message(paste0(e, "\n"));errval})
}

#' @export
makeID <- function(...) {
  x <- paste0(...)
  as.integer(factor(x, levels = unique(x)))
}

#' @export
samp <- function(samples, match, partial = F) {
  if(!partial)
    samples[,which(dimnames(samples)[[2]] == match)]
  else
    samples[,which(substr(dimnames(samples)[[2]],1,nchar(match)) == match)]
}

#' @export
samp3d <- function(samples, match, partial = F) {
  if(!partial)
    samples[,,which(dimnames(samples)[[3]] == match)]
  else
    samples[,,which(substr(dimnames(samples)[[3]],1,nchar(match)) == match)]
}





