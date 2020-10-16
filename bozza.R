require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function

degree_sequence = read.table("./data/English_out-degree_sequence.txt", header = FALSE)
degree_spectrum = table(degree_sequence)

N <- dim(degree_sequence)[1]

barplot(degree_spectrum, main = "English", xlab = "degree", ylab = "number of vertices")
barplot(degree_spectrum, main = "English", xlab = "degree", ylab = "number of vertices", log = "xy")

x <- degree_sequence$V1
minus_log_likelihood_zeta <- function(gamma) {
  length(x) * log(zeta(gamma)) + gamma * sum(log(x))
}

mle_zeta <- mle(minus_log_likelihood_zeta, start = list(gamma = 2),
                method = "L-BFGS-B", lower = c(1.0000001))

summary(mle_zeta)
attributes(summary(mle_zeta))
gamma <- attributes(summary(mle_zeta))$coef[1]

attributes(summary(mle_zeta))$m2logL


get_AIC <- function(m2logL,K,N) {
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

AIC <- get_AIC(attributes(summary(mle_zeta))$m2logL, 1, N)



