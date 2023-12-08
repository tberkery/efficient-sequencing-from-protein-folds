from hmmlearn import hmm
import numpy as np
import sys

def burrows_wheeler_transform(s):
    """ Apply Burrows-Wheeler Transform to a given string. """
    # Append a unique symbol to the end of the string to mark its end
    s = s + '$'
    # Create a table of rotations of the string
    table = sorted(s[i:] + s[:i] for i in range(len(s)))
    # Extract the last column
    last_column = ''.join(row[-1] for row in table)
    return last_column

def create_bwts(sequence, chunk_size=30):
    """ Creates BWTs for every chunk of size `chunk_size` in a protein conformation sequence. """
    bwts = []
    for i in range(0, len(sequence), chunk_size):
        chunk = sequence[i:i+chunk_size]

        # Apply the Burrows-Wheeler Transform to the chunk
        bwt = burrows_wheeler_transform(chunk)

        # Append the BWT of the chunk to the list
        bwts.append(bwt)

    return bwts


def bwt_hmm(seq, chunk_size=30):
    """ Applies BWT-HMM to long protein conformation sequence by merging BWTs of smaller chunks """

    # Create BWTs for each chunk
    bwts = create_bwts(seq, chunk_size)

    # Define HMM parameters
    n_states = 4  # Example value
    start_probs = np.full(n_states, 1.0 / n_states)
    transmat_probs = np.full((n_states, n_states), 1.0 / n_states**2)

    # Initialize the HMM
    model = hmm.MultinomialHMM(n_components=n_states)

    # Prepare observed data
    observed_data = []
    for bwt in bwts:
        # Filter out non-numeric characters
        filtered_bwt = ''.join(c for c in bwt if c.isdigit())
        # Convert to integers and store in observed data
        observed_data.extend([int(c) for c in filtered_bwt])

    # Convert observed data to NumPy array
    observed_data = np.array([[x] for x in observed_data])

    # Train the HMM
    model.fit(observed_data)

    return model, bwts

# Example usage
seq = 66666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666665555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555552222222222222222222222222222222222222222222222222222222222222222222222222222222222444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
seq = str(seq)
model, bwts = bwt_hmm(seq)

# Model can be used for prediction or other analysis tasks
print("Trained HMM model:", model)
#print("Transition Matrix:\n", model.transmat_)
#print("Emission Probabilities:\n", model.emissionprob_)
#print("Initial State Probabilities:\n", model.startprob_)


# BWTs can be used for further analysis or indexing
print("BWTs for individual chunks:", bwts)
