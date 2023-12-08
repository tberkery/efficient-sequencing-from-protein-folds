from hmmlearn import hmm
import numpy as np
import sys

def bwt_hmm(seq, chunk_size=30):
    """
    Applies BWT-HMM to a very long protein conformation sequence by merging BWTs of smaller chunks.

    Args:
        seq: A string representing the protein conformation sequence.
        chunk_size: The size of each chunk to create a BWT for.

    Returns:
        A trained HMM model and a list of BWTs for each chunk.
    """

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
seq = "...your very long protein conformation sequence..."
model, bwts = bwt_hmm(seq)

# Model can be used for prediction or other analysis tasks
print("Trained HMM model:", model)
#print("Transition Matrix:\n", model.transmat_)
#print("Emission Probabilities:\n", model.emissionprob_)
#print("Initial State Probabilities:\n", model.startprob_)

# BWTs can be used for further analysis or indexing
print("BWTs for individual chunks:", bwts)
