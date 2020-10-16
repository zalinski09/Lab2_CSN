# procedure to obtain the best parameters by maximum likelihood
require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function

degree_sequence = read.table("./data/English_out-degree_sequence.txt", header = FALSE)
degree_spectrum = table(degree_sequence)
print(degree_spectrum)

N <- dim(degree_sequence)[1]

x <- degree_sequence$V1

minus_log_likelihood_zeta <- function(gamma) {
  N * log(zeta(gamma)) + gamma * sum(log(x))
}

mle_zeta <- mle(minus_log_likelihood_zeta, start = list(gamma = 2),
                method = "L-BFGS-B", lower = c(1.0000001))

summary(mle_zeta)
attributes(summary(mle_zeta))
res <- attributes(summary(mle_zeta))$coef[1]

# initial value for the parameters
# q0 = N/M
# lambda0 = M/N
# gamma0 = 2
# kmax0 = N