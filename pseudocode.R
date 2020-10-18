# pseudocode of the final script (Lab2_CSN)
# _____________________________________________________________________________
#
#
#
# 1) for each out-degree dataset (English,...):
#       read the data: e.g. degree_sequence = read.table(...)
#       remove unlinked nodes (those with ki = 0)
#       (plot the degree spectrum to visualize the data, e.g. in scale log-log)
# 2) produce a summary table as Table [1]
#
# 3) for each distribution (from 1 to 5):
#       estimate the free parameter(s) (except for model 3, the zeta(2))
#       produce a summary table as Table [3]
#       compute AIC = -2L + 2K * (N / N-K-1), where L is the log-likelihood
# 4) compute AICbest = min(AICi)
#
# 5) for each distribution (from 1 to 5):
#       compute DELTA = AIC - AICbest
# 6) produce a summary table as Table [4]
#
# 7) evaluate the data with a new probability distribution (e.g. Altmann function)
#