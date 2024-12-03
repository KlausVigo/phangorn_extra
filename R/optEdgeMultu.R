
# needs to go in phangorn extra
optEdgeMulti <- function(object, control = pml.control(epsilon = 1e-8,
                                                       maxit = 10, trace = 1, tau = 1e-8), ...) {
  tree <- object$tree
  theta <- object$tree$edge.length
  weight <- attr(object$data, "weight")
  ll0 <- object$logLik
  eps <- 1
  iter <- 0
  iter2 <- 0
  scale <- 1
  # l <- length(theta)
  while (abs(eps) > control$eps && iter < control$maxit) {
    dl <- score(object)
    thetaNew <- log(theta) + scale * solve(dl[[2]], dl[[1]]) # + diag(l) * 1e-10
    newtheta <- exp(thetaNew)
    tree$edge.length <- as.numeric(newtheta)
    object <- update(object, tree = tree)
    ll1 <- object$logLik
    eps <- (ll0 - ll1) / ll1
    if (eps < 0) {
      newtheta <- theta
      scale <- scale / 2
      tree$edge.length <- as.numeric(theta)
      ll1 <- ll0
      iter2 <- iter2 + 1
    }
    else {
      scale <- 1
      iter2 <- 0
    }
    theta <- newtheta
    if (iter2 == 0 && control$trace > 0) cat("loglik: ", ll1, "\n")
    ll0 <- ll1
    if (iter2 == 10) iter2 <- 0
    if (iter2 == 0) iter <- iter + 1
  }
  object <- update(object, tree = tree)
  object
}
