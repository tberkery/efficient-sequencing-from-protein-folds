library(tidyverse)
library(ggplot2)
library(gridExtra)

df_tc = read_csv("benchmarking_results.csv") %>%
  filter(iteration <= 9)
# assuming protein-fold-pipeline is your working directory

plots = vector(mode = "list", length = 4)
data_structures = as.vector(unique(df_tc$data_structure))
sample_types = as.vector(unique(df_tc$sample))
colors = c("blue", "blue", "blue", "blue")
binwidths = c(15, 0.01, 0.1, 0.005)
counter = 1

for (ds in data_structures) {
  for (s in sample_types) {
    df_tc_sub = df_tc %>%
      filter(data_structure == ds, sample == s)
    p = ggplot(df_tc_sub, aes(x = runtime)) +
      geom_histogram(binwidth = binwidths[counter])
    plots[[counter]] = p
    counter = counter + 1
  }
}

grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2)
ggsave("time_complexity_histograms.png")

df_tc_100k = df_tc %>%
  filter(sample == "hypothetical_protein_test")
df_tc_comp = df_tc %>%
  filter(sample == "compressed_hypothetical_protein_test")
p1 = ggplot(df_tc_100k, aes(x = runtime, fill = data_structure)) +
  geom_density() +
  labs(title = "Runtime for Example Protein",
       x = "Runtime",
       y = "Density")
p2 = ggplot(df_tc_comp, aes(x = runtime, fill = data_structure)) +
  geom_density() +
  labs(title = "Runtime for Compressed Protein",
       x = "Runtime",
       y = "Density")

ggsave("runtime_density_test_protein.png", p1, width = 12, height = 8, units = "in")
ggsave("runtime_density_compressed_protein.png", p2, width = 12, height = 8, units = "in")

