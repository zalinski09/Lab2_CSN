# test of the correct functioning of mle for geometric distribution
# (just parameter estimation, not model selection)

require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function


# initial value for the parameters
# q0 = N/M          --> this is needed for the geometric
# lambda0 = M/N
# gamma0 = 2
# kmax0 = N

degree_sequence = read.table("./samples_from_discrete_distributions/data/sample_of_geometric_with_parameter_0.05.txt",
                             header = FALSE)
degree_spectrum = table(degree_sequence)

N <- dim(degree_sequence)[1]
x <- degree_sequence$V1
M <- sum(x)

minus_log_likelihood_geometric <- function(q) {
  -(M - N) * log(1-q) - N * log(q)
}

mle_geometric <- mle(minus_log_likelihood_geometric, start = list(q = N / M),
                method = "L-BFGS-B", lower = c(0.0000001), upper = c(0.9999999))

res <- attributes(summary(mle_geometric))$coef[1]
cat("q: ", res)

