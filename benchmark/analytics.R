library(tidyverse)

# File paths assume you run from home directory of project repo 
# (not the benchmarking directory that contains this R Script)

memory_usage = read_csv("memory_usage_updated.csv")
runtime = read_csv("benchmarking_results_updated.csv")

# Get statistics for regular simulated protein

memory_usage_stats = memory_usage %>%
  filter(sample == "hypothetical_protein_test") %>%
  group_by(data_structure) %>%
  summarize(mean = mean(memory_usage, na.rm = TRUE),
            median = median(memory_usage, na.rm = TRUE),
            stdev = sd(memory_usage, na.rm = TRUE),
            count = n(),
            .groups = 'keep')

memory_usage_stats %>% write_csv("memory_usage_statistics_updated.csv")

runtime_stats = runtime %>%
  filter(sample == "hypothetical_protein_test") %>%
  group_by(data_structure) %>%
  summarize(mean = mean(runtime, na.rm = TRUE),
            median = median(runtime, na.rm = TRUE),
            stdev = sd(runtime, na.rm = TRUE),
            count = n(),
            .groups = 'keep')

runtime_stats %>% write_csv("runtime_statistics_updated.csv")

# t-Test for difference between means across algorithms for simulated protein
group1 = "MLSE-Viterbi"
non_viterbi_data_structures = setdiff(unique(runtime$data_structure), group1)
group2a = non_viterbi_data_structures[[1]]
group2b = non_viterbi_data_structures[[2]]

memory_usages_group1 = (memory_usage %>% 
                   filter(sample == "hypothetical_protein_test",
                          data_structure == group1))$memory_usage
runtimes_group1 = (runtime %>% 
              filter(sample == "hypothetical_protein_test",
                     data_structure == group1))$runtime

memory_usages_group2a = (memory_usage %>% 
                          filter(sample == "hypothetical_protein_test",
                                 data_structure == group2a))$memory_usage
runtimes_group2a = (runtime %>% 
                     filter(sample == "hypothetical_protein_test",
                            data_structure == group2a))$runtime

memory_usages_group2b = (memory_usage %>% 
                           filter(sample == "hypothetical_protein_test",
                                  data_structure == group2b))$memory_usage
runtimes_group2b = (runtime %>% 
                      filter(sample == "hypothetical_protein_test",
                             data_structure == group2b))$runtime

t_test_memory_viterbi_and_bwt_hmm = t.test(memory_usages_group1, memory_usages_group2a, alternative = 'less')
t_test_runtime_viterbi_and_bwt_hmm = t.test(runtimes_group1, runtimes_group2a, alternative = 'less')

t_test_memory_viterbi_and_bwt_hmm_with_merge = t.test(memory_usages_group1, memory_usages_group2b, alternative = 'less')
t_test_runtime_viterbi_and_bwt_hmm_with_merge = t.test(runtimes_group1, runtimes_group2b, alternative = 'less')

# Remark: to view t-test results (they are listed in the report), physically open the t_test_... objects in e.g. RStudio.
# They are nested with many fields, and it will be easiest to view them this way.