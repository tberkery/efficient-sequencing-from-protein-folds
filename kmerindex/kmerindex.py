import sys

inputName = "../generate_sequence/simulation.txt" #"compressed_MDseq.txt"
#k = 7
outputName = "output.txt"


def kmerindex(inputName, kmin, kmax):
    # Read the sequence into a list then convert to string
    with open(inputName, 'r') as ifile:
        lines = ifile.readlines()

    # Process sequence
    sequence = ''
    for i in range(len(lines)):
        sequence = sequence + lines[i].strip()

    
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


bestkmers = kmerindex("../generate_sequence/simulation.txt", kmin=5, kmax=15)
print(bestkmers)