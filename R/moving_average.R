

cladeMatrix <- function(x, rooted = FALSE) {
  if (!rooted) x <- unroot(x)
  pp <- prop.part(x)
  pplabel <- attr(pp, "labels")
  if (!rooted) pp <- SHORTwise(pp)
  x <- .uncompressTipLabel(x)
  nnodes <- Nnode(x)
  class(x) <- NULL
  # nnodes <- sapply(x, Nnode)
  l <- length(x)
  from <- cumsum(c(1, nnodes[-l]))
  to <- cumsum(nnodes)
  ivec <- integer(to[l])
  pvec <- c(0, to)
  res <- vector("list", l)
  k <- 1
  for (i in 1:l) {
    ppi <- prop.part(x[[i]])
    if (!rooted) ppi <- SHORTwise(ppi)
    indi <- sort(fmatch(ppi, pp))
    ivec[from[i]:to[i]] <- indi
  }
  X <- sparseMatrix(i = ivec, p = pvec, dims = c(length(pp), l))
  list(X = X, prop.part = pp)
}


moving_average <- function(obj, window = 50) {
  fun <- function(x) {
    cx <- c(0, cumsum(x))
    (cx[(window + 1):length(cx)] - cx[1:(length(cx) - window)]) / (window)
  }
  res <- apply(obj$X, 1, fun)
  rownames(res) <- c()
}


