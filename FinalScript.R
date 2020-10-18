# Final script for lab 2

#
# ausiliary functions
#
write_summary <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  return (c(language, length(degree_sequence$V1), max(degree_sequence$V1), sum(degree_sequence$V1)/length(degree_sequence$V1),
            length(degree_sequence$V1)/sum(degree_sequence$V1)))
}

# 1 - Poisson
minus_log_likelihood_poisson <- function(lambda) {
  -M * log(lambda) + N * (lambda + log(1 - exp(-lambda))) + C
}

# 2 - Geometric
minus_log_likelihood_geometric <- function(q) {
  -(M - N) * log(1-q) - N * log(q)
}

# 4 - Zeta
minus_log_likelihood_zeta <- function(gamma) {
  N * log(zeta(gamma)) + gamma * Mp
}

# 5 - Right-truncated zeta
rtrunc_zeta <- function(gamma, kmax) {
  x <- c(1:kmax)
  sum(x ^ - gamma)
}

minus_log_likelihood_rtrunc_zeta <- function(gamma, kmax) {
  N * log( rtrunc_zeta(gamma, kmax) ) + gamma * Mp
}

get_AIC <- function(m2logL,K,N) {
  m2logL + 2*K*N/(N-K-1) 
}
#
# end auxiliary functions
#


# 1) for each out-degree dataset (English,...):
#       read the data: e.g. degree_sequence = read.table(...)
#       remove unlinked nodes (those with ki = 0)
#       (plot the degree spectrum to visualize the data, e.g. in scale log-log)
# 2) produce a summary table as Table [1]
#

properties <- data.frame("language" = character(), "N" = integer(),
                              "max_degree" = double(),
                              "Mean degree (M/N)" = double(),
                              "Inverse mean degree (N/M)" = double(),
                              check.names = FALSE)

source = read.table("list_out.txt", header = TRUE, as.is = c("language","file"))

for (x in 1:nrow(source)) {
  properties[x,] <- write_summary(source$language[x], source$file[x])
}



# 3) for each distribution (from 1 to 5):
#       estimate the free parameter(s) (except for model 3, the zeta(2))
#       produce a summary table as Table [3]
#       compute AIC = -2L + 2K * (N / N-K-1), where L is the log-likelihood
# 4) compute AICbest = min(AICi)
#

mle_params <- data.frame("language" = character(), "lambda" = double(),
                         "q" = double(), "gamma1" = double(),
                         "gamma2" = double(), "kmax" = double(),
                         check.names = FALSE)

AIC <- data.frame("language" = character(), "AIC_poisson" = double(),
                  "AIC_geometric" = double(), "AIC_zeta_2" = double(),
                  "AIC_zeta" = double(), "AIC_trunc_zeta" = double(),
                  "AIC_best" = double())



for (i in 1:nrow(source)) {
  # each iteration is a different language
  degree_sequence <- read.table(source$file[i], header = FALSE)
  lang <- source$language[i]
  x <- degree_sequence$V1
  
  N <- dim(degree_sequence)[1]
  M <- sum(x)
  Mp <- sum(log(x))
  C <- sum(sapply(x, function(k) sum(2:k)))
  
  
  # 1: Poisson
  mle_poisson <- mle(minus_log_likelihood_poisson, start = list(lambda = M / N),
                     method = "L-BFGS-B", lower = c(1.0000001))
  
  # 2: geometric
  mle_geometric <- mle(minus_log_likelihood_geometric, start = list(q = N / M),
                       method = "L-BFGS-B", lower = c(0.0000001),
                       upper = c(0.9999999))
  
  # 4: zeta
  mle_zeta <- mle(minus_log_likelihood_zeta, start = list(gamma = 2),
                  method = "L-BFGS-B", lower = c(1.0000001))
  
  # 5: truncated-zeta
  mle_trunc_zeta <- mle(minus_log_likelihood_rtrunc_zeta,
                        start = list(gamma = 2, kmax = properties$max_degree[i]),
                        method = "L-BFGS-B", lower = c(1.0000001, 1),
                        upper = c(7 * attributes(summary(mle_zeta))$coef[1],
                                  properties$max_degree[i]))
  
  
  mle_params[i,] <- (c(lang, attributes(summary(mle_poisson))$coef[1],
                       attributes(summary(mle_geometric))$coef[1],
                       attributes(summary(mle_zeta))$coef[1],
                       attributes(summary(mle_trunc_zeta))$coef[1],
                       attributes(summary(mle_trunc_zeta))$coef[2]))
  
  AIC[i,] <- c(lang, get_AIC(attributes(summary(mle_poisson))$m2logL, 1, N),
               get_AIC(attributes(summary(mle_geometric))$m2logL, 1, N),
               2 * (N * log((pi^2)/6) + 2 * Mp),
               get_AIC(attributes(summary(mle_zeta))$m2logL, 1, N),
               get_AIC(attributes(summary(mle_trunc_zeta))$m2logL, 2, N), 0)
  
  AIC$AIC_best[i] <- min(AIC$AIC_poisson[i],
                         AIC$AIC_geometric[i],
                         AIC$AIC_zeta_2[i],
                         AIC$AIC_zeta[i],
                         AIC$AIC_trunc_zeta[i])
  
}





# 5) for each distribution (from 1 to 5):
#       compute DELTA = AIC - AICbest
# 6) produce a summary table as Table [4]
#
# 7) evaluate the data with a new probability distribution (e.g. Altmann function)