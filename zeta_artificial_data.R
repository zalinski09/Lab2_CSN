# test of the correct functioning of mle for zeta distribution
# (just parameter estimation, not model selection)

require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function


# initial value for the parameters
# q0 = N/M
# lambda0 = M/N
# gamma0 = 2      --> this is needed for the zeta
# kmax0 = N

degree_sequence = read.table("./samples_from_discrete_distributions/data/sample_of_zeta_with_parameter_1.5.txt",
                             header = FALSE)
degree_spectrum = table(degree_sequence)

N <- dim(degree_sequence)[1]
x <- degree_sequence$V1

minus_log_likelihood_zeta <- function(gamma) {
  N * log(zeta(gamma)) + gamma * sum(log(x))
}

mle_zeta <- mle(minus_log_likelihood_zeta, start = list(gamma = 2),
                method = "L-BFGS-B", lower = c(1.0000001))

res <- attributes(summary(mle_zeta))$coef[1]
cat("gamma: ", res)

