# Final script for lab 2

require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function

#
# auxiliary functions
#
write_summary <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  return (list(language, length(degree_sequence$V1), max(degree_sequence$V1),
            sum(degree_sequence$V1) / length(degree_sequence$V1),
            length(degree_sequence$V1) / sum(degree_sequence$V1)))
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

# Altmann
c_const <- function(gamma, delta) {
  j <- c(1:N)
  1 / sum(j ^ -gamma * exp(-delta * j))
}

minus_log_likelihood_altmann <- function(gamma, delta) {
  gamma * Mp + delta * M - N * log(c_const(gamma, delta))
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
                         "max_degree" = integer(),
                         "Mean degree (M/N)" = double(),
                         "Inverse mean degree (N/M)" = double(),
                         check.names = FALSE)

source = read.table("list_out.txt", header = TRUE, as.is = c("language","file"))

for (i in 1:nrow(source)) {
  properties[i,] <- write_summary(source$language[i], source$file[i])
}



# 3) for each distribution (from 1 to 5):
#       estimate the free parameter(s) (except for model 3, the zeta(2))
#       produce a summary table as Table [3]
#       compute AIC = -2L + 2K * (N / N-K-1), where L is the log-likelihood
# 4) compute AICbest = min(AICi)
#

mle_params <- data.frame("language" = character(), "lambda" = double(),
                         "q" = double(), "gamma1" = double(),
                         "gamma2" = double(), "kmax" = integer(),
                         check.names = FALSE)

AIC <- data.frame("language" = character(), "AIC_poisson" = double(),
                  "AIC_geometric" = double(), "AIC_zeta_2" = double(),
                  "AIC_zeta" = double(), "AIC_trunc_zeta" = double(),
                  "AIC_best" = double(), check.names = FALSE)



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
                        method = "L-BFGS-B",
                        lower = c(1.0000001, properties$max_degree[i])
                        )
  
  
  
  mle_params[i,] <- list(lang, attributes(summary(mle_poisson))$coef[1],
                      attributes(summary(mle_geometric))$coef[1],
                      attributes(summary(mle_zeta))$coef[1],
                      attributes(summary(mle_trunc_zeta))$coef[1],
                      attributes(summary(mle_trunc_zeta))$coef[2])
  
  AIC[i,] <- list(lang, get_AIC(attributes(summary(mle_poisson))$m2logL, 1, N),
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

DELTA <- data.frame("language" = character(), "Model 1" = double(),
                    "Model 2" = double(), "Model 3" = double(),
                    "Model 4" = double(), "Model 5" = double(),
                    check.names = FALSE)

for (i in 1:nrow(source)) {
  DELTA[i,] <- list(source$language[i],
                 AIC$AIC_poisson[i] - AIC$AIC_best[i],
                 AIC$AIC_geometric[i] - AIC$AIC_best[i],
                 AIC$AIC_zeta_2[i] - AIC$AIC_best[i],
                 AIC$AIC_zeta[i] - AIC$AIC_best[i],
                 AIC$AIC_trunc_zeta[i] - AIC$AIC_best[i])
}



# 7) evaluate the data with a new probability distribution (e.g. Altmann function)
#

altmann_params <- data.frame("language" = character(),
                             "gamma" = double(),
                             "delta" = double(),
                             "AIC" = double(),
                             check.names = FALSE)

AIC_new <- data.frame("language" = character(), "AIC_poisson" = double(),
                      "AIC_geometric" = double(), "AIC_zeta_2" = double(),
                      "AIC_zeta" = double(), "AIC_trunc_zeta" = double(),
                      "AIC_altmann" = double(), "AIC_best" = double(),
                      check.names = FALSE)


for (i in 1:nrow(source)) {
  degree_sequence <- read.table(source$file[i], header = FALSE)
  lang <- source$language[i]
  x <- degree_sequence$V1
  
  N <- dim(degree_sequence)[1]
  M <- sum(x)
  Mp <- sum(log(x))
  
  
  mle_altmann <- mle(minus_log_likelihood_altmann, start = list(gamma = 2, delta = 1),
                     method = "L-BFGS-B", lower = c(1.0000001, 0.0000001))
  
  altmann_params[i,] <- list(lang, attributes(summary(mle_altmann))$coef[1],
                             attributes(summary(mle_altmann))$coef[2],
                             get_AIC(attributes(summary(mle_altmann))$m2logL, 2, N))
  
  AIC_new[i,] <- list(lang, AIC$AIC_poisson[i], AIC$AIC_geometric[i],
                      AIC$AIC_zeta_2[i], AIC$AIC_zeta[i], AIC$AIC_trunc_zeta[i],
                      get_AIC(attributes(summary(mle_altmann))$m2logL, 2, N), 0)
  
  AIC_new$AIC_best[i] <- min(AIC_new$AIC_poisson[i],
                         AIC_new$AIC_geometric[i],
                         AIC_new$AIC_zeta_2[i],
                         AIC_new$AIC_zeta[i],
                         AIC_new$AIC_trunc_zeta[i],
                         AIC_new$AIC_altmann[i])
}

DELTA_new <- data.frame("language" = character(), "Model 1" = double(),
                    "Model 2" = double(), "Model 3" = double(),
                    "Model 4" = double(), "Model 5" = double(),
                    "Model 6" = double(), check.names = FALSE)

for (i in 1:nrow(source)) {
  DELTA_new[i,] <- list(source$language[i],
                        AIC_new$AIC_poisson[i] - AIC_new$AIC_best[i],
                        AIC_new$AIC_geometric[i] - AIC_new$AIC_best[i],
                        AIC_new$AIC_zeta_2[i] - AIC_new$AIC_best[i],
                        AIC_new$AIC_zeta[i] - AIC_new$AIC_best[i],
                        AIC_new$AIC_trunc_zeta[i] - AIC_new$AIC_best[i],
                        AIC_new$AIC_altmann[i] - AIC_new$AIC_best[i])
}


