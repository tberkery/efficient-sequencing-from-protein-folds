import time
from memory_profiler import profile
import os
import sys

# import any functions that you wish to profile in terms of space and time complexity
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from BWT_HMM.bwthmm import bwt_hmm

print("foo")