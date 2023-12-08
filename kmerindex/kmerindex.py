import sys

inputName = "../generate_sequence/simulation.txt"
k = 10
outputName = "output.txt"

# Read the sequence into a list then convert to string
with open(inputName, 'r') as ifile:
    lines = ifile.readlines()

sequence = ''
for i in range(len(lines)):
    sequence = sequence + lines[i].strip()

kmerIndex = {}
for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
    if kmer not in kmerIndex:
        kmerIndex[kmer] = []
    kmerIndex[kmer].append(i)

# Write the length of the kmerIndex into a file
with open(outputName, 'w') as ofile:
    ofile.write(str(len(kmerIndex)))


# Find most common kmer
commonkmers = []
for kmer in kmerIndex:
    number = len(kmerIndex[kmer])
    commonkmers.append((number, kmer))
    
commonkmers.sort(reverse=True)
for i in range(5):
    print(commonkmers[i]) 

#print(commonKmer, number)
