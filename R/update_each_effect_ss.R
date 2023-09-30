# @title update each effect once
# @param XtX a p by p matrix, X'X
# @param Xty a p vector
# @param s_init a list with elements sigma2, V, alpha, mu, Xr
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param estimate_prior_method The method used for estimating prior
#   variance, 'optim' or 'EM'.
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
# 
#' @importFrom Matrix diag
update_each_effect_ss = function (XtX, Xty, s_init,
                                  estimate_prior_variance = FALSE,
                                  estimate_prior_method = "optim",
                                  check_null_threshold = 0) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"
  
  # Repeat for each effect to update.
  s = s_init
  L = nrow(s$alpha)
  if (L > 0) {
    for (l in 1:L) {
        
      # Remove lth effect from fitted values.
      s$XtXr = s$XtXr - XtX %*% (s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      XtR = Xty - s$XtXr

      # Correct zR discrepancy. For motivation see DENTIST paper Chen et al (2021)
      if (s$correct_zR_discrepancy$to_correct) {
        if (s$correct_zR_discrepancy$is_init) {
          # skip the first iteration but turn on the zR correction mode
          s$force_iterate = TRUE
          s$correct_zR_discrepancy$is_init = FALSE
        } else {
          # Get the current existing non-zero effect variables
          c_index = get_non_zero_effects_proxy(s$alpha, s$V, s$pi) 
          if (length(c_index)>0) {
            # Detect outlier against existing non-zero effect variables
            outlier_index = detect_zR_discrepancy(c_index, XtR, XtX, r2=0.6, p=1E-4)
            # Apply correction
            if (outlier_index>0) {
              s$correct_zR_discrepancy$outlier_index = union(s$correct_zR_discrepancy$outlier_index, outlier_index)
              s$pi[s$correct_zR_discrepancy$outlier_index] = 0
              s$pi = s$pi / sum(s$pi)
            }
          }
          # Check if corrections are done, and if so let IBSS proceed as usual with convergence check
          if (l==L && all(s$correct_zR_discrepancy$outlier_index == s_init$correct_zR_discrepancy$outlier_index)) {
            s$correct_zR_discrepancy$to_correct = FALSE
            s$force_iterate = FALSE
          }
        }
      }

      res = single_effect_regression_ss(as.matrix(XtR),attr(XtX,"d"),s$V[l],
              s$sigma2,s$pi,estimate_prior_method,check_null_threshold)
      
      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l]      = res$V
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$lbf_model +
        SER_posterior_e_loglik_ss(attr(XtX,"d"),XtR,s$sigma2,
                                  res$alpha * res$mu,res$alpha * res$mu2)
      s$XtXr = s$XtXr + XtX %*% (s$alpha[l,] * s$mu[l,])
    }
  }
  s$XtXr = unname(as.matrix(s$XtXr))
  return(s)
}

detect_zR_discrepancy <- function(c_index, z, XtX, r2=0.6,p=1E-4) {
  chisq_cutoff = qchisq(1-p, df = 1)
  max_index = which.max(abs(z))
  if (max_index %in% c_index) {
    return(-1)
  } else {
    R = cor(XtX[c(max_index, c_index), c(max_index, c_index)])
    z_test = z[c_index]
    z_max = z[max_index]
    # DENTIST-S test, $S(\hat{z}_1, \hat{z}_2, r_{12}) = \frac{(\hat{z}_1 - r_{12}\hat{z}_2)^2}{1-r_{12}^2} \sim \chi^2_{(1)}$
    stats_filter = sapply(1:length(z_test), function(i) ((z_max - R[1, i+1] * z_test[i])^2 / 1 - R[1, i+1]^2))
    stats_filter = any(stats_filter > chisq_cutoff)
    r2_filter = sapply(1:length(z_test), function(i) R[1, i+1]^2)
    stats_filter = any(stats_filter > r2)
    if(stats_filter && r2_filter) {
      return (max_index)
    } else {
      return (-1)
    }
  }
}

# KL Divergence
kl_divergence_p_q <- function(p, q, epsilon=1E-10) {
  q <- q + epsilon
  return(sum(ifelse(p == 0, 0, p * (log(p+epsilon)-log(q)))))
}

# Jensen-Shannon Divergence
js_divergence_p_q <- function(p, q) {
  m <- (p + q) / 2
  return((kl_divergence_p_q(p, m) + kl_divergence_p_q(q, m)) / 2)
}

get_non_zero_effects_proxy = function(alpha, V, p, tol=5E-5) {
  # effects to drop if V is too small
  to_drop = which(V<=tol)
  # single effects prior
  p = p/sum(p)
  # effects to drop if alpha too close to single effects prior --- in case uses don't update V
  for (i in 1:nrow(alpha)) {
    js = js_divergence_p_q(alpha[i,], p)
    if (js<tol) {
      # two vectors very similiar
      to_drop = union(to_drop, i)
    }
  }
  max_indices <- apply(alpha, MARGIN = 1, FUN = which.max)
  # only keep non-zero effects based on V
  max_indices[-to_drop]
}
