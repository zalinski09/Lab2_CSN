# This script just contains the minus log-likelihood functions and how to call
# them.
# It does not work as a standalone.
# A value for N, M, Mp, C must be provided.



#
# 1 - Poisson
#
minus_log_likelihood_poisson <- function(lambda) {
  -M * log(lambda) + N * (lambda + log(1 - exp(-lambda))) + C
}

#mle_poisson <- mle(minus_log_likelihood_poisson, start = list(lambda = M / N),
#                   method = "L-BFGS-B", lower = c(1.0000001))




#
# 2 - Geometric
#
minus_log_likelihood_geometric <- function(q) {
  -(M - N) * log(1-q) - N * log(q)
}

#mle_geometric <- mle(minus_log_likelihood_geometric, start = list(q = N / M),
#                     method = "L-BFGS-B", lower = c(0.0000001), upper = c(0.9999999))




#
# 3 - Zeta with gamma = 2 (not a function, just a computation)
#
#mle_zeta_2 <- N * log((pi^2) / 6) + 2 * Mp




#
# 4 - Zeta
#
minus_log_likelihood_zeta <- function(gamma) {
  N * log(zeta(gamma)) + gamma * Mp
}

#mle_zeta <- mle(minus_log_likelihood_zeta, start = list(gamma = 2),
#                method = "L-BFGS-B", lower = c(1.0000001))




#
# 5 - Right-truncated zeta
#
rtrunc_zeta <- function(gamma, kmax) {
  x <- c(1:kmax)
  sum(x ^ - gamma)
}

minus_log_likelihood_rtrunc_zeta <- function(gamma, kmax) {
  N * log( rtrunc_zeta(gamma, kmax) ) + gamma * Mp
}

#mle_right_truncated_zeta <- mle(minus_log_likelihood_rtrunc_zeta,
#                                start = list(gamma = 2, kmax = N),
#                                method = "L-BFGS-B", lower = c(1.0000001, 1))





# (not used)
# minus_log_likelihood <- function(N, M, C, lambda, q, gamma1, gamma2, kmax)



