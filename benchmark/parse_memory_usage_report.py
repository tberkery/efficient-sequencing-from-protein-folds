import sys
import pandas as pd
mem_usage_report_file = sys.argv[1]
use_next = False
df = pd.DataFrame(columns = ["data_structure", "sample", "memory_usage"])
data_structure_names = ["BWT-HMM", "MLSE-Viterbi"]
sample_names = ["hypothetical_protein_test", "compressed_hypothetical_protein_test"]
i = 0
with open(mem_usage_report_file, "r") as fp:
    for line in fp.readlines():
        if use_next:
            local_count += 1
            use_next = False
            line_split = line.split("    ") # this is not a tab, but rather a copy of a specific number of spaces
            print(line_split)
            mem_usage = line_split[2].split()
            memory_used = mem_usage[0]
            print(memory_used)
        if line[0:1] == "=":
            print("detected =")
            use_next = True
            local_count = 0
        else: