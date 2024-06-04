library(MASS)
library(stats)
library(glue)

###########################################################################################################################
# FUNCION OBJETIVO
###########################################################################################################################

#' @import MASS
#' @import stats
#' @import glue
random_projection_test = function(X, Y = NULL, mu0 = NULL, proj_dimension = NULL){
  p <- dim(X)[2]
  sample_mean1 <- colMeans(X)
  if (is.null(Y) && is.null(mu0) && is.null(proj_dimension)){
    n <- dim(X)[1]
    k <- floor(n/2)
    P_k_T <- matrix(stats::rnorm(k * p, mean = 0, sd = 1), nrow = k, ncol = p)
    mu0 <- rep(0, p)
    S_n <- stats::cov(X)
    T_k2 <- T_k2(0, n, k, P_k_T, sample_mean1, mu0, S_n)

    return_list <- list(statistic = T_k2, p_value = find_alpha(n, k, P_k_T, sample_mean1, mu0, S_n), degrees_freedom = list(k, n-k), null_value = mu0, method = "One-sample projection test")
    return(return_list)
  }
  else if (is.null(Y) && is.null(mu0) && !(is.null(proj_dimension))){
    n <- dim(X)[1]
    k <- proj_dimension
    if (k>=min(n,p)){stop("Error: Projection dimension must be less the minimum value between sample size and sample dimension")}
    P_k_T <- matrix(stats::rnorm(k * p, mean = 0, sd = 1), nrow = k, ncol = p)
    mu0 <- rep(0, p)
    S_n <- stats::cov(X)
    T_k2 <- T_k2(0, n, k, P_k_T, sample_mean1, mu0, S_n)

    return_list <- list(statistic = T_k2, p_value = find_alpha(n, k, P_k_T, sample_mean1, mu0, S_n), degrees_freedom = list(k, n-k), null_value = mu0, method = glue::glue("One-sample projection test with projection dimension {k}"))
    return(return_list)
  }
  else if (is.null(Y) && !(is.null(mu0)) && !(is.null(proj_dimension))){
    n <- dim(X)[1]
    k <- proj_dimension
    if (k>=min(n,p)){stop("Error: Projection dimension must be less the minimum value between sample size and sample dimension")}
    P_k_T <- matrix(stats::rnorm(k * p, mean = 0, sd = 1), nrow = k, ncol = p)
    mu0 <- mu0
    S_n <- stats::cov(X)
    T_k2 <- T_k2(0, n, k, P_k_T, sample_mean1, mu0, S_n)

    return_list <- list(statistic = T_k2, p_value = find_alpha(n, k, P_k_T, sample_mean1, mu0, S_n), degrees_freedom = list(k, n-k), null_value = mu0, method = glue::glue("One-sample projection test with chosen mean under null hypothesis and projection dimension {k}"))
    return(return_list)
  }
  else if (is.null(Y) && !(is.null(mu0)) && is.null(proj_dimension)){
    n <- dim(X)[1]
    k <- floor(n/2)
    P_k_T <- matrix(stats::rnorm(k * p, mean = 0, sd = 1), nrow = k, ncol = p)
    mu0 <- mu0
    S_n <- stats::cov(X)
    T_k2 <- T_k2(0, n, k, P_k_T, sample_mean1, mu0, S_n)

    return_list <- list(statistic = T_k2, p_value = find_alpha(n, k, P_k_T, sample_mean1, mu0, S_n), degress_freedom = list(k, n-k), null_value = mu0, method = glue::glue("One-sample projection test with chosen mean under null hypothesis"))
    return(return_list)
  }
  else if (!(is.null(Y)) && is.null(mu0) && is.null(proj_dimension)){
    n1 <- dim(X)[1]
    n2 <- dim(Y)[1]
    n <- n1+n2-2
    k <- floor(n/2)
    P_k_T <- matrix(stats::rnorm(k * p, mean = 0, sd = 1), nrow = k, ncol = p)
    mu0 <- rep(0,p)
    sample_mean2 <- colMeans(Y)
    S_n <- ((n1-1)/n)*stats::cov(X)+((n2-1)/n)*stats::cov(Y)
    T_k <- T_k2_2(n1, n2, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0)

    return_list <- list(statistic = T_k, p_value = find_alpha_2(n1, n2, n, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0), degrees_freedom = list(k, n-k+1), null_value = mu0, method = glue::glue("Two-sample projection test"))
    return(return_list)
  }
  else if (!(is.null(Y)) && !(is.null(mu0)) && is.null(proj_dimension)){
    n1 <- dim(X)[1]
    n2 <- dim(Y)[1]
    n <- n1+n2-2
    k <- floor(n/2)
    P_k_T <- matrix(stats::rnorm(k * p, mean = 0, sd = 1), nrow = k, ncol = p)
    mu0 <- mu0
    sample_mean2 <- colMeans(Y)
    S_n <- ((n1-1)/n)*stats::cov(X)+((n2-1)/n)*stats::cov(Y)
    T_k2 <- T_k2_2(n1, n2, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0 = mu0)

    return_list <- list(statistic = T_k2, p_value = find_alpha_2(n1, n2, n, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0), degrees_freedom = list(k, n-k+1), null_value = mu0, method = glue::glue("Two-sample projection test with chosen mean under null hypothesis"))
    return(return_list)
  }
  else if (!(is.null(Y)) && is.null(mu0) && !(is.null(proj_dimension))){
    n1 <- dim(X)[1]
    n2 <- dim(Y)[1]
    n <- n1+n2-2
    k <- proj_dimension
    if (k>=min(n,p)){stop("Error: Projection dimension must be less the minimum value between sample size and sample dimension")}
    P_k_T <- matrix(stats::rnorm(k * p, mean = 0, sd = 1), nrow = k, ncol = p)
    mu0 <- rep(0,p)
    sample_mean2 <- colMeans(Y)
    S_n <- ((n1-1)/n)*stats::cov(X)+((n2-1)/n)*stats::cov(Y)
    T_k2 <- T_k2_2(n1, n2, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0)

    return_list <- list(statistic = T_k2, p_value = find_alpha_2(n1, n2, n, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0), degrees_freedom = list(k, n-k+1), null_value = mu0, method = glue::glue("Two-sample projection test with k = {k}"))
    return(return_list)
  }
  else if (!(is.null(Y)) && !(is.null(mu0)) && !(is.null(proj_dimension))){
    n1 <- dim(X)[1]
    n2 <- dim(Y)[1]
    n <- n1+n2-2
    k <- proj_dimension
    if (k>=min(n,p)){stop("Error: Projection dimension must be less the minimum value between sample size and sample dimension")}
    P_k_T <- matrix(stats::rnorm(k * p, mean = 0, sd = 1), nrow = k, ncol = p)
    mu0 <- mu0
    sample_mean2 <- colMeans(Y)
    S_n <- ((n1-1)/n)*stats::cov(X)+((n2-1)/n)*stats::cov(Y)
    T_k2 <- T_k2_2(n1, n2, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0 = mu0)

    return_list <- list(statistic = T_k2, p_value = find_alpha_2(n1, n2, n, k, P_k_T, sample_mean1, sample_mean2, S_n, mu0), degrees_freedom = list(k, n-k+1), null_value = mu0, method = glue::glue("Two-sample projection test with chosen mean under null hypothesis and projection dimension {k}"))
    return(return_list)
  }
}
