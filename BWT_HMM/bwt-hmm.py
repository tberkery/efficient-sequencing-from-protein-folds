from hmmlearn import hmm
import numpy as np
import sys

seq = ""
with open(sys.argv[1], "r") as fp:
    for line in fp.readlines():
        seq += line.rstrip("\n")

def burrows_wheeler_transform(s):
    """ Apply Burrows-Wheeler Transform to a given string. """
    # Append a unique symbol to the end of the string to mark its end
    s = s + '$'
    # Create a table of rotations of the string
    table = sorted(s[i:] + s[:i] for i in range(len(s)))
    # Extract the last column
    last_column = ''.join(row[-1] for row in table)
    return last_column

# Example protein conformation data as a string of numbers from 1 to 6
example_data = seq #"123456123456123456123456123456123456123456123456"
# Apply BWT to the example data
bwt_result = burrows_wheeler_transform(example_data)
print(bwt_result)
print(type(bwt_result))

# Since we're dealing with discrete data (numbers 1 to 6), we use a Multinomial HMM
# Number of states in the HMM needs to be defined. This requires domain knowledge or experimentation.
n_states = 4  # Example value

# Initialize the HMM
model = hmm.MultinomialHMM(n_components=n_states)
start_probs = np.full(n_states, 1.0 / n_states)  # initially assume equal probability of each state
transmat_probs = np.full((n_states, n_states), 1.0 / n_states**2)
model.set_params(startprob_prior = start_probs, transmat_prior = transmat_probs)

# The data needs to be in a specific format for hmmlearn, usually as a 2D numpy array
# Here, each number is considered a separate observation
# For example, '66$1122334455' becomes [[6], [6], [1], [1], [2], [2], [3], [3], [4], [4], [5], [5]]
observed_data = np.array([[int(x)] for x in bwt_result if x.isdigit()]) # Note that $ get filtered out
print("Data provided to HMM:", observed_data)

# Training the model (this requires actual data and might need to adjust the parameters)
# For demonstration, I am just using the observed_data for training, but in practice,
# you'd train this on a larger, representative dataset
model.fit(observed_data) # fits behind-the-scenes using EM algorithm.

# After training, this model can be used to predict states or evaluate new observations
print(model)