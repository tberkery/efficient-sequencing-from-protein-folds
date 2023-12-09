import sys

inputName = "../generate_sequence/simulation.txt" #"compressed_MDseq.txt"
#k = 7
outputName = "output.txt"

with open(inputName, 'r') as ifile:
    lines = ifile.readlines()

# Process sequence
sequence = ''
for i in range(len(lines)):
    sequence = sequence + lines[i].strip()


'''
Takes a sequence and outputs the 10 best matching kmers with all ks between kmin and kmax
'''
def kmerindex(sequence, kmin, kmax):
    
    # Records all kmers with their occurances in the form (occurance * length, kmer)
    commonkmers = []

    # Make kmer of different lengths and find best matches
    for k in range(kmin, kmax):

        # Build kmer index
        kmerIndex = {}
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer not in kmerIndex:
                kmerIndex[kmer] = []
            kmerIndex[kmer].append(i)

        # Append all kmers in the form (occurance * length, kmer)
        for kmer in kmerIndex:
            number = len(kmerIndex[kmer])
            commonkmers.append((len(kmer) * number, kmer))
    
    # Print out the 5 best matches
    commonkmers.sort(reverse=True)
    return commonkmers[:10]


bestkmers = kmerindex(sequence, kmin=5, kmax=10)

for pair in bestkmers:
    print("Sequence: {}, Occurance: {}".format(pair[1], int(pair[0] / len(pair[1]))))
