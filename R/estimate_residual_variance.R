# @title Estimate residual variance
# @param X an n by p matrix of covariantes
# @param y an n vector of data
# @param s a susie fit
estimate_residual_variance = function (X, y, s) {
  n = nrow(X)
  return((1/n)*get_ER2(X,y,s))
}

# @title Estimate residual variance for summary statistics
# @param XtX a p by p matrix
# @param Xty a p vector
# @param s a susie fit
# @param yty a scaler, y'y, where y is centered to have mean 0
# @param n sample size
estimate_residual_variance_ss = function (XtX, Xty, s, yty, n)
  (1/n)*get_ER2_ss(XtX,Xty,s,yty)


MoM = function(alpha, mu, mu2, sigmasq, X, y, verbose) {
  n = nrow(X)
  p = ncol(X)
  L = nrow(mu)

  # Compute XtX and its trace
  XtX = t(X) %*% X
  trXtX = sum(diag(XtX))

  # Compute Xty and yty
  Xty = t(X) %*% y
  yty = sum(y^2)

  # Compute E[beta]
  b = colSums(alpha * mu)

  # Compute E[beta beta^T]
  E_beta_betaT = matrix(0, p, p)
  for (l in 1:L) {
    # Use mu2 to get the second moment directly
    E_beta_betaT = E_beta_betaT + diag(alpha[l,] * mu2[l,]) +
      outer(alpha[l,] * mu[l,], alpha[l,] * mu[l,])
  }

  # Compute m1
  m1 = yty - 2 * sum(b * Xty) + sum(diag(XtX %*% E_beta_betaT))

  # Update sigmasq
  sigmasq = m1 / n

  if (verbose) cat(sprintf("Update sigma^2 to %f\n", sigmasq))

  return(sigmasq)
}
