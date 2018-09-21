smahal <-
function (X) 
{
  units <- rownames(X)
  X <- as.matrix(X)
  s<-apply(X,2,sd)
  stopifnot(!all(s==0))
  X<-X[,s>0,drop=FALSE]
  m <- dim(X)[1]
  rownames(X) <- 1:m
  k <- dim(X)[2]
  for (j in 1:k) X[, j] <- rank(X[, j])
  cv <- cov(X)
  vuntied <- var(1:m)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat) %*% cv %*% diag(rat)
  out <- matrix(NA, m, m)
  rownames(out) <- units#rownames(X)
  colnames(out) <- units
  if (!requireNamespace("MASS", quietly = TRUE)) {
		stop("Please load package 'MASS'")
	}
    icov<-MASS::ginv(cv)
  for (i in 1:m) {
    out[i, ] <- sqrt(mahalanobis(X, X[i, ], icov, inverted = T))
  }
  out
}
