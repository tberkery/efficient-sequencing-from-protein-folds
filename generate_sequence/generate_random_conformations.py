"""
Generate a randomized test sequence. If indicated, allow for random, probabilistic "patterns" of conformational change. 
"""

from random import random, sample
from copy import deepcopy

"""
ConformationPattern

A designated, sequential folding path, e.g.: A->B->C or A->B->C->A

The path can be deterministic (P=1) if ANY protein in the specified conformation will ALWAYS
continue following the folding path. Otherwise, a Markovian probability should be specified for
each position in the pattern. If a probability is specified, then the pattern will preserved
at index 0, 1, 2, ... so long as the previous index was also preserved based on the randomly
sampled probability.

For instance, if the pattern is ABCDE, then:
    p=[0,0,...]: '' (no sequence insertion)
    p=[1,1,...] or p=None: 'ABCDE'
    p=[1,0.5,0.5,...]: the final sequence will vary. 50% 'A', 25% 'AB', 25% 'ABC...'
"""
class ConformationPattern:
    """
    pattern: List: the ordered list that defines the conformation pattern, 
        e.g.: A->B->C = ['A', 'B', 'C']

    p: probability of each character in the pattern
        None: pattern cannot change
        List: the ordered list corresponding to Markovian probabilities at each sequence
            position

    preserved:
        The previous sequence was preserved, so the pattern can continue
    """
    def __init__(self, pattern, p):
        self.pattern = pattern
        self.length = len(pattern)

        self.p = p
        if len(self.p) != len(self.pattern):
            raise ValueError("A probability must be designated at each position of the ConformationPattern!")

        self.preserved = True
    
    def generate_pattern(self, deterministic=False, idx_seed=0):
        if deterministic:
            return ''.join(self.pattern[idx_seed:])
        
        idx_curr = idx_seed
        while self.preserved and idx_curr < self.length:
            sample = random()
            if sample > self.p[idx_curr]:
                self.preserved = False
            else:
                idx_curr += 1
        return ''.join(self.pattern[idx_seed:idx_curr])

"""
ConformationGenerator

Make a random sequence of conformations within a designated alphabet. 

"""
class ConformationGenerator_TransitionMatrix:
    """
    alphabet: Set: a list of possible characters

    transitions: 
        None: uniform probabilities for all characters
        Dict[str, Dict[float]]: 1D probability distribution transition matrix
            Will be processed into a cumulative distribution transition matrix
    """
    def __init__(self, alphabet, transitions, patterns=[], seed=None, alpha=0.01):
        self.alphabet = sorted(list(set(alphabet)))
        self.len_alphabet = len(self.alphabet)

        self.seed = seed
        if self.seed is None:
            self.seed = self.sample()
        if self.seed not in self.alphabet:
            raise ValueError("Invalid seed")
        
        self.patterns = patterns

        self.transitions_cdf = dict()
        if transitions is None:
            self.transitions_cdf = None
        else:
            # transition probabilities
            for key in self.alphabet:
                p_uniform = 1/self.len_alphabet
                self.transitions_cdf[key] = list() # sorted list of cumulative probabilities
                cum_next = 0
                if key not in transitions:
                    for next in self.alphabet:
                        cum_next += p_uniform
                        self.transitions_cdf[key].append(cum_next)
                else:
                    # Use the transition matrix to make cdf.
                    #   Also allow for sparse transition matrix, in which we assign uniform distribution
                    #   to remaining undefined transitions

                    num_seen = len(transitions[key].keys())

                    p_defined = 0
                    for next in self.alphabet:
                        if next in transitions[key]:
                            p_defined += transitions[key][next]
                    if p_defined > 1 + alpha:
                        raise ValueError("Probability sum exceeds 1 for some transition in the transition matrix")
                    if p_defined < 1 - alpha and num_seen == self.len_alphabet:
                        raise ValueError("Probability does not sum up to 1 despite some transition being fully defined in the transition matrix")
                    n_remaining = self.len_alphabet - num_seen
                    p_remaining = (1-p_defined)/n_remaining

                    for next in self.alphabet:
                        if next in transitions[key]:
                            cum_next += transitions[key][next]
                        else:
                            cum_next += p_remaining
                        self.transitions_cdf[key].append(cum_next)

    def sample(self):
        return sample(self.alphabet, k=1)[0]
    """
    Use cumulative distribution and sample the next state. 
    """
    def transition(self, prev):
        if self.transitions_cdf is None:
            return self.sample()
        rand = random()
        idx_next = 0
        while idx_next < self.len_alphabet and self.transitions_cdf[prev][idx_next] < rand:
            idx_next += 1
        idx_next = min(idx_next, self.len_alphabet)
        return self.alphabet[idx_next]


    def generate(self, n):
        if n <= 0:
            return ''
        
        seq = str(self.seed)
        prev = self.seed
        offset = 1
        while offset < n:
            next = self.transition(prev)
            seq += str(next)
            offset += 1
            prev = next

        return seq
    
    """
    Near-duplicate of generate, but writes sequence instead of saving the whole thing to memory
    """
    def generate_and_write(self, n, outfile):
        with open(outfile, 'w') as fh:
            if n > 0:
                prev = self.seed
                fh.write(str(prev))
                offset = 1
                while offset < n:
                    next = self.transition(prev)
                    offset += 1
                    prev = next
                    fh.write(str(prev))
        return

class ConformationGenerator_Patterns:
    def __init__(self, alphabet, transitions, patterns=[], seed=None, alpha=0.01):
        self.alphabet = sorted(list(set(alphabet)))
        self.len_alphabet = len(self.alphabet)

        self.seed = seed
        if self.seed is None:
            self.seed = self.sample()
        if self.seed not in self.alphabet:
            raise ValueError("Invalid seed")
        
        self.patterns = patterns

def main():
    alphabet = []
    for i in range(65,65+7):
        alphabet.append(chr(i))
    transitions = { # sparse transition matrix (any "unused" probability density will be distributed among remaining classes)
        'A':{'B':0.5},
        'B':{'C':0.5},
        'C':{'D':0.5},
        'D':{'E':0.5}
    }
    generator = ConformationGenerator_TransitionMatrix(alphabet, transitions)
    print(generator.generate(100))
    generator.generate_and_write(10**5, 'sim_ABCDE.txt')
    return

if __name__ == '__main__':
    main()

