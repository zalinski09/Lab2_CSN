# Script for the model selection


write_summary <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  return (c(language, length(degree_sequence$V1), max(degree_sequence$V1), sum(degree_sequence$V1)/length(degree_sequence$V1),
            length(degree_sequence$V1)/sum(degree_sequence$V1)))
}


# 1) for each out-degree dataset (English,...):
#       read the data: e.g. degree_sequence = read.table(...)
#       remove unlinked nodes (those with ki = 0)
#       (plot the degree spectrum to visualize the data, e.g. in scale log-log)
#
# 2) produce a summary table as Table [1]
#
degree_sequence <- data.frame("language" = character(), "N" = integer(), "Max degree" = double(),
                              "Mean degree (M/N)" = double(), "Inverse mean degree (N/M)" = double())

source = read.table("list_out.txt", header = TRUE, as.is = c("language","file"))

for (x in 1:nrow(source)) {
  degree_sequence[x,] <- write_summary(source$language[x], source$file[x])
}

degree_sequence