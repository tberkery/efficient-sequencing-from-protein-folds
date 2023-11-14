import time
from memory_profiler import profile
import os
import sys
import pandas as pd

# Make sure to install all of the above packages

# Once all fields are complete, run the following line in the command line: python .\benchmark\profile.py > space_complexity_info.txt

# Import any functions that you wish to profile in terms of space and time complexity
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from main_method.compress_runs import compress_runs
from main_method.propose_paths import query_seq
from BWT_HMM.bwthmm import bwt_hmm

#user = "Tad Berkery" # specify who you are so we now whose computer the profiling corresponds to
user = 'Richard Hu'

# Define functions that return desired sample sequence to profile with data structures here
def get_sample_sequence():
    seq = ""
    with open("BWT_HMM/test.txt", "r") as fp:
        for line in fp.readlines():
            seq += line.rstrip("\n")
    # when running on test.txt sequence input, gets killed on Amazon EC2 t3.2xlarge instance
    seq = seq[0:100000] # try on just first 100,000 characters to see if we can benchmark
    return seq

def get_smaller_sequence():
    seq, runs, len_seq = compress_runs.compress_runs("BWT_HMM/test.txt")
    return seq, runs, len_seq

# For every data structure function you import, define custom function with @profile annotation to enable space complexity analysis
# Make sure to pass sequence as argument to each defined function and to call data structure function defined in other file with this sequence passed as well.
@profile
def bwt_hmm_with_space_annotation(seq):
    bwt_hmm(seq)

@profile
def make_proposals_mlse(seq, num_proposals=10, max_path_len=10, sep='_'):
    len_seq = len(seq)
    query_seq(seq, len_seq, num_proposals, max_path_len, sep=sep, threshold=0.01, verbose=True)

# It is highly recommended that you make sure any data structure functions you call don't print anything to output. Otherwise, it will be hard to read space complexity info printed by memory_profiler to stdout.

ds_functions = [bwt_hmm_with_space_annotation, make_proposals_mlse] # add any functions for data structures you wish to benchmark (make sure you have made a copy and used the @profile annotation for memory_profiler space complexity anlaysis.)
sample_functions = [get_sample_sequence, get_smaller_sequence] # now define any functions that will provide samples you wish to profile

sample_function_descriptors = ["hypothetical_protein_test", "query_test"] # write appropriate names describing example sequences to test here
data_structure_descriptors = ["BWT-HMM", "MLSE-Viterbi"] # write appropriate names describing data structures to profile here
num_iterations = 100 # set number of times to run each unique sample sequence/data structure combination

results = pd.DataFrame(columns = ["data_structure", "sample", "iteration", "user", "runtime"])
sample_counter = 0
for func in sample_functions:
    ds_counter = 0
    sample_name = sample_function_descriptors[sample_counter]
    for ds in ds_functions:
        ds_name = data_structure_descriptors[ds_counter]
        for i in range(num_iterations):
            start = time.time()
            ds(func())
            end = time.time()
            runtime = end - start
            results.loc[len(results.index)] = [ds_name, sample_name, i, user, runtime]
        ds_counter += 1
    sample_counter += 1

results.to_csv("benchmarking_results.csv", index = False)
