library(MASS)
library(stats)
library(glue)

#One-sample
#' @import MASS
T_k2 <- function(alpha, n, k, P_k_T, sample_mean1, mu0, S_n) {
  T_k2_value <- ((n - k) / k) * t(P_k_T %*% (sample_mean1 - mu0)) %*% MASS::ginv(P_k_T %*% S_n %*% t(P_k_T)) %*% P_k_T %*% (sample_mean1 - mu0)
  return(T_k2_value)
}

#' @import stats
t_alpha <- function(alpha, n, k) {
  t_alpha_value <- stats::qf(1 - alpha, df1 = k, df2 = n - k)
  return(t_alpha_value)
}

difference_function <- function(alpha, n, k, P_k_T, sample_mean1, mu0, S_n) {
  return(T_k2(alpha, n, k, P_k_T, sample_mean1, mu0, S_n) - t_alpha(alpha, n, k))
}

#' @import stats
find_alpha <- function(n, k, P_k_T, sample_mean1, mu0, S_n) {
  result <- stats::uniroot(difference_function, interval = c(0, 1), n = n, k = k, P_k_T = P_k_T, sample_mean1 = sample_mean1, mu0 = mu0, S_n = S_n)

  if (is.na(result$root)) {
    return(NA)
  } else {
    alpha <- result$root
    alpha_adjusted <- ifelse(alpha < 0, 0, ifelse(alpha > 1, 1, alpha))
    return(alpha_adjusted)
  }
}

# Two-sample
#' @import MASS
T_k2_2 <- function(n1, n2, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0) {
  T_k2_value <- (n1*n2/(n1+n2)) * t(P_k_T %*% (sample_mean1 - sample_mean2 - mu0)) %*% MASS::ginv(P_k_T %*% S_n %*% t(P_k_T)) %*% P_k_T %*% (sample_mean1 - sample_mean2 - mu0)
  return(T_k2_value)
}

#' @import stats
t_alpha_2 <- function(alpha, n, k) {
  t_alpha_value <- (k*n/(n-k+1)) * stats::qf(1 - alpha, df1 = k, df2 = n - k + 1)
  return(t_alpha_value)
}

difference_function_2 <- function(alpha, n1, n2, n, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0) {
  return(T_k2_2(n1, n2, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0) - t_alpha_2(alpha, n, k))
}

#' @import stats
find_alpha_2 <- function(n1, n2, n, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0) {
  result <- stats::uniroot(difference_function_2, interval = c(0, 1), n1 = n1, n2 = n2, n = n, k = k, P_k_T = P_k_T, sample_mean1 = sample_mean1, sample_mean2 = sample_mean2, S_n = S_n, mu0 = mu0)

  if (is.na(result$root)) {
    return(NA)
  } else {
    alpha <- result$root
    alpha_adjusted <- ifelse(alpha < 0, 0, ifelse(alpha > 1, 1, alpha))
    return(alpha_adjusted)
  }
}
