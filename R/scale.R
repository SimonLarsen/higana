scale_eigenstrat <- function(x, mult) {
  sum2 <- colSums(apply(x, 2, function(r) !is.na(r)))
  p <- colSums(x, na.rm=TRUE) / sum2 / 2
  xsd <- sqrt(mult * p * (1-p))
  scale(x, center=2*p, scale=xsd)
}

scale2 <- function(x, type="sd") {
  if(type == "sd") x <- scale(x, center=TRUE, scale=TRUE)
  else if(type == "center") x <- scale(x, center=TRUE, scale=FALSE)
  else if(type == "binom") x <- scale_eigenstrat(x, 1)
  else if(type == "binom2") x <- scale_eigenstrat(x, 2)
  else stop("Unrecognized standardization method \"", type, "\".")
  
  x[is.na(x)] <- 0
  x
}

