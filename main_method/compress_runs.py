"""
Compress single-character runs. 

Return:
    shortSeq - the compressed sequence; "AAABBC" --> "ABC"
    runLengths - the runs indexed to shortSeq; "ABC" --> [3,2,1]
"""
def compress_runs(fileName):
    seq = ""
    with open(fileName, 'r') as ifile:
        for line in ifile:
            seq += line.strip()
    
    shortSeq = seq[0]
    runLengths = [1]
    j = 0 # index in shortSeq
    
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            runLengths[j] += 1
        else:
            shortSeq = shortSeq + seq[i]
            runLengths.append(1)
            j += 1
    len_seq = j
    return shortSeq, runLengths, len_seq