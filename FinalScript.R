# 1) for each out-degree dataset (English,...):
#       read the data: e.g. degree_sequence = read.table(...)
#       remove unlinked nodes (those with ki = 0)
#       (plot the degree spectrum to visualize the data, e.g. in scale log-log)
# 2) produce a summary table as Table [1]
#

degree_sequence <- data.frame("language" = character(), "N" = integer(), "Max degree" = double(),
                              "Mean degree (M/N)" = double(), "Inverse mean degree (N/M)" = double(), check.names = FALSE)

source = read.table("list_out.txt", header = TRUE, as.is = c("language","file"))

for (x in 1:nrow(source)) {
  degree_sequence[x,] <- write_summary(source$language[x], source$file[x])
}



# 3) for each distribution (from 1 to 5):
#       estimate the free parameter(s) (except for model 3, the zeta(2))
#       produce a summary table as Table [3]
#       compute AIC = -2L + 2K * (N / N-K-1), where L is the log-likelihood
N <- dim(degree_sequence)[1]
x <- degree_sequence$V1
M <- sum(x)


mle_poisson <- mle(minus_log_likelihood_poisson, start = list(lambda = M / N),
                   method = "L-BFGS-B", lower = c(1.0000001))
poisson_sequence <- data.frame("language" = character(), "N" = integer(), "Max degree" = double(),
                              "Mean degree (M/N)" = double(), "Inverse mean degree (N/M)" = double(), check.names = FALSE)

source_poisson = read.table("list_poisson.pdf", header = TRUE, as.is = c("language","file"))

for (x in 1:nrow(source_poisson)) {
  poisson_sequence[x,] <- write_summary(source_poisson$language[x], source_poisson$file[x])
}

#mle_geometric <- mle(minus_log_likelihood_geometric, start = list(q = N / M),
#                     method = "L-BFGS-B", lower = c(0.0000001), upper = c(0.9999999))

#mle_zeta <- mle(minus_log_likelihood_zeta, start = list(gamma = 2),
#                method = "L-BFGS-B", lower = c(1.0000001))

#mle_right_truncated_zeta <- mle(minus_log_likelihood_right_truncated_zeta, start = list(gamma = 2, kmax = N),
#                                method = "L-BFGS-B", lower = c(1.0000001, 1))







get_AIC <- function(m2logL,K,N) {
  m2logL + 2*K*N/(N-K-1) 
}

AIC <- get_AIC(attributes(summary(mle_zeta))$m2logL, 1, N)
# 4) compute AICbest = min(AICi)
#
# 5) for each distribution (from 1 to 5):
#       compute DELTA = AIC - AICbest
# 6) produce a summary table as Table [4]
#
# 7) evaluate the data with a new probability distribution (e.g. Altmann function)