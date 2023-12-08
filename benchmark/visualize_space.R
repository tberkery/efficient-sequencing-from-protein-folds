library(tidyverse)
library(ggplot2)
library(gridExtra)

df_tc = read_csv("memory_usage_updated.csv")
# assuming protein-fold-pipeline is your working directory

plots = vector(mode = "list", length = 4)
data_structures = as.vector(unique(df_tc$data_structure))
sample_types = as.vector(unique(df_tc$sample))
colors = c("blue", "blue", "blue", "blue")
binwidths = c(0.5, 0.1, 0.5, 0.1)
counter = 1
sample_titles = c("Full", "Comp", "Full", "Comp")

for (ds in data_structures) {
  for (s in sample_types) {
    df_tc_sub = df_tc %>%
      filter(data_structure == ds, sample == s)
    p = ggplot(df_tc_sub, aes(x = memory_usage)) +
      geom_histogram(binwidth = binwidths[counter]) +
      labs(title = paste0(ds, ", ", sample_titles[counter]),
           x = "Memory Used (MiB)",
           y = "Density")
    plots[[counter]] = p
    counter = counter + 1
  }
}

grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2, widths = rep(7, 2), heights = c(5, 5))
#ggsave("space_complexity_histograms.png", width = 12, height = 8, units = "in") # cannot save a grid arrange object... this line just saves last plot in grid

df_tc_100k = df_tc %>%
  filter(sample == "hypothetical_protein_test")
df_tc_comp = df_tc %>%
  filter(sample == "compressed_hypothetical_protein_test")
p1 = ggplot(df_tc_100k, aes(x = memory_usage, fill = data_structure)) +
  geom_density() +
  labs(title = "Space Complexity for Example Protein",
       x = "Memory Used (MiB)",
       y = "Density")
p2 = ggplot(df_tc_comp, aes(x = memory_usage, fill = data_structure)) +
  geom_density() +
  labs(title = "Space Complexity for Compressed Protein",
       x = "Memory Used (MiB)",
       y = "Density")

ggsave("space_density_test_protein_updated.png", p1, width = 12, height = 8, units = "in")
ggsave("space_density_compressed_protein_updated.png", p2, width = 12, height = 8, units = "in")

