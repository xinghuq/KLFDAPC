
kmatrixGauss=function (x, sigma = 1)
{
  x <- t(as.matrix(x))
  d <- nrow(x)
  n <- ncol(x)
  X2 <- t(as.matrix(colSums(x^2)))
  distance2 <- repmat(X2, n, 1) + repmat(t(X2), 1, n) - 2 *
    t(x) %*% x
  K <- exp(-distance2/(2 * sigma^2))
  return(K)
}



KLFDA=function(kdata, y, r, sigma=0.5, metric = c("weighted", "orthonormalized",
                                    "plain"),tol=1e-5, knn = 6, reg = 0.001)
{
  require(lfda)
  k=kdata
  metric <- match.arg(metric)
  y <- t(as.matrix(y))
  n <- nrow(k)
  if (is.null(r))
    r <- n
  tSb <- mat.or.vec(n, n)
  tSw <- mat.or.vec(n, n)
  for (i in unique(as.vector(t(y)))) {
    Kcc <- k[y == i, y == i]
    Kc <- k[, y == i]
    nc <- nrow(Kcc)
    Kccdiag <- diag(Kcc)
    distance2 <- repmat(Kccdiag, 1, nc) + repmat(t(Kccdiag),
                                                 nc, 1) - 2 * Kcc
    A <- getAffinityMatrix(distance2, knn, nc)
    Kc1 <- as.matrix(rowSums(Kc))
    Z <- Kc %*% (repmat(as.matrix(colSums(A)), 1, n) * t(Kc)) -
      Kc %*% A %*% t(Kc)
    tSb <- tSb + (Z/n) + Kc %*% t(Kc) * (1 - nc/n) + Kc1 %*%
      (t(Kc1)/n)
    tSw <- tSw + Z/nc
  }
  K1 <- as.matrix(rowSums(k))
  tSb <- tSb - K1 %*% t(K1)/n - tSw
  tSb <- (tSb + t(tSb))/2
  tSw <- (tSw + t(tSw))/2
  F=tSb/tSw
  require(lfda)
  eigTmp <- suppressWarnings(rARPACK::eigs(A = solve(tSw + reg * diag(1, nrow(tSw), ncol(tSw)),tol=tol) %*% tSb, k = r,
                                        which = "LM"))
  eigVec <- Re(eigTmp$vectors)
  eigVal <- as.matrix(Re(eigTmp$values))

  Tr <- getMetricOfType(metric, eigVec, eigVal, n)
  Z <- t(t(Tr) %*% k)
  out <- list(T = Tr, Z = Z)
  class(out) <- "KLFDA"
  return(out)
}


