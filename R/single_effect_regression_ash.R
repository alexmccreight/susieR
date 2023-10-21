#' @rdname single_effect_regression_ash
#'
#' @title Bayesian single-effect linear regression using ASH prior
#'
#' @description These methods fit the regression model \eqn{y = Xb +
#'   e}, where elements of e are \emph{i.i.d.}  \eqn{N(0,s^2)}, and b is
#'   a p-vector of effects to be estimated. The assumption is that b has
#'   exactly one non-zero element, with all elements equally likely to
#'   be non-zero. The prior on the coefficient of the non-zero element
#'   is \eqn{N(0,V)}.
#'
#' @details \code{single_effect_regression_ss} performs single-effect
#' linear regression with summary data, in which only the statistcs
#' \eqn{X^Ty} and diagonal elements of \eqn{X^TX} are provided to the
#' method.
#'
#' \code{single_effect_regression_rss} performs single-effect linear
#' regression with z scores. That is, this function fits the
#' regression model \eqn{z = R*b + e}, where e is \eqn{N(0,Sigma)},
#' \eqn{Sigma = residual_var*R + lambda*I}, and the b is a p-vector of
#' effects to be estimated. The assumption is that b has exactly one
#' non-zero element, with all elements equally likely to be non-zero.
#' The prior on the non-zero element is \eqn{N(0,V)}. The required
#' summary data are the p-vector \code{z} and the p by p matrix
#' \code{Sigma}. The summary statistics should come from the same
#' individuals.
#'
#' @param y An n-vector.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param V A scalar giving the (initial) prior variance
#'
#' @param residual_variance The residual variance.
#'
#' @param prior_weights A p-vector of prior weights.
#'
#' @param optimize_V The optimization method to use for fitting the
#'   prior variance.
#'
#' @param check_null_threshold Scalar specifying threshold on the
#'   log-scale to compare likelihood between current estimate and zero
#'   the null.
#'
#' @return A list with the following elements:
#'
#' \item{alpha}{Vector of posterior inclusion probabilities;
#'   \code{alpha[i]} is posterior probability that the ith coefficient
#'   is non-zero.}
#'
#' \item{mu}{Vector of posterior means (conditional on inclusion).}
#'
#' \item{mu2}{Vector of posterior second moments (conditional on
#'   inclusion).}
#'
#' \item{lbf}{Vector of log-Bayes factors for each variable.}
#'
#' \item{lbf_model}{Log-Bayes factor for the single effect regression.}
#'
#' \code{single_effect_regression} and \code{single_effect_regression_ss}
#' additionally output:
#'
#' \item{V}{Prior variance (after optimization if \code{optimize_V !=
#'   "none"}).}
#'
#' \item{loglik}{The log-likelihood, \eqn{\log p(y | X, V)}.}
#'
#' @importFrom mixsqp mixsqp
#' @importFrom stats dnorm 
#' @importFrom Matrix rowSums
#' @importFrom ashr ash
#'
#' @keywords internal
#'
single_effect_regression_ash =
  function (y, X, V, residual_variance = 1, prior_weights = NULL,
            optimize_V = c("ash"),
            check_null_threshold = 0) {
  Xty = compute_Xty(X,y)
  betahat = (1/attr(X,"d")) * Xty
  shat2 = residual_variance/attr(X,"d")
  sdhat = sqrt(shat2)
  if (is.null(prior_weights))
    prior_weights = rep(1/ncol(X),ncol(X))

  ash_obj = ash(bhat, sdhat, mixcompdist = "normal")
  lbf = log_BF_ash(ash_obj, bhat, sdhat)

  #### Update it to "favor" most likely causal effect
  loglik <- get_loglik(ash_obj, bhat, sdhat)
  new_pi <- m_step(L = loglik, zeta = cal_zeta(lbf))
  new_ash_obj <- update_ash_obj(ash_obj, new_pi)

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0
  maxlbf = max(lbf)

  # w is proportional to BF, but subtract max for numerical stability.
  w = exp(lbf - maxlbf)

  # Posterior prob for each SNP.
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  post_var = (1/V + attr(X,"d")/residual_variance)^(-1) # Posterior variance.
  post_mean = (1/residual_variance) * post_var * Xty
  post_mean2 = post_var + post_mean^2 # Second moment.

  # BF for single effect model.
  lbf_model = maxlbf + log(weighted_sum_w)
  loglik = lbf_model + sum(dnorm(y,0,sqrt(residual_variance),log = TRUE))

  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
              lbf_model = lbf_model,V = V,loglik = loglik))
}

log_BF_ash <- function(m, bhat, sdhat) {
  
  lBF <- numeric(length(sdhat))
  tt <- numeric(length(sdhat))
  pi_k <- m$fitted_g$pi
  sd_k <- m$fitted_g$sd
  
  sd_combined <- outer(sdhat, sd_k, function(s1, s2) sqrt(s1^2 + s2^2))
  dnorm_vals <- matrix(0, nrow = length(bhat), ncol = length(pi_k))
  
  for (k in seq_along(pi_k)) {
    dnorm_vals[, k] <- pi_k[k] * dnorm(bhat, sd = sd_combined[, k])
  }
  
  tt <- rowSums(dnorm_vals)
  lBF <- log(tt) - dnorm(bhat, sd = sdhat, log = TRUE)
  
  return(lBF)
}

get_loglik_mat <- function(ash_obj, bhat, sdhat) {
  
  sdmat <- sqrt(outer(c(sdhat^2), ash_obj$fitted_g$sd^2, "+"))
  
  L <- dnorm(outer(bhat, rep(0, length(ash_obj$fitted_g$sd)), `-`) / sdmat, log = TRUE) - log(sdmat)
  
  # dealing in case of due to small sd due to small sample size
  L <- apply(L, 2, function(x) replace(x, which(is.na(x)), median(x, na.rm = TRUE)))
  
  return(L)
}

m_step <- function(L, zeta, init_pi0_w = 1, control_mixsqp = list(verbose = FALSE)) {
  
  tlength <- ncol(L) - 1
  mixsqp_out <- mixsqp(
    L,
    zeta,
    log = TRUE,
    x0 = c(init_pi0_w, rep(1e-6, tlength))
  )
  
  return(mixsqp_out$x)
}

cal_zeta <- function(lBF) {
  max_lBF <- max(lBF)
  out <- exp(lBF - max_lBF) / sum(exp(lBF - max_lBF))
  return(out)
}

update_ash_obj <- function(ash_obj, new_pi) {
  ash_obj$fitted_g$pi <- new_pi
  return(ash_obj)
}

##### create ash object ---
bhat <- rnorm(100)
sdhat <- runif(100, min = 0.1, max = 2)


print(ash_obj$fitted_g$pi)
print(new_ash_obj$fitted_g$pi)