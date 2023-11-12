# genomics-sequences-final-project
Final project for Computational Genomics: Sequences

Live goals
1. What are the common folding paths; querying them, counts, etc.?
    a. Find kth order probabilities, build HMM/Viterbi? Wheeler-HMM variation of BWT-HMM?
    b. Use data structure in (1) to "propose" favorable paths of protein folding
    c. Query folding paths against live sequence
    d. Kinetics of a folding path (number of runs, RSA, etc.)
    d. Sequence of folding paths? Longest common substring?
2. Check "equilibrium" (kth order probabilities remain consistent at different points of simulation).

Good article on viterbi: https://www.sciencedirect.com/topics/mathematics/viterbi-algorithm#:~:text=A%20Viterbi%20algorithm%20is%20a,Technology%20(Third%20Edition)%2C%202003


TODO:
1. Try existing: BWT-HMM --> produce the index, run some queries? Also, analyze time and space
2. Viterbi (existing; also look at implementation into Wheeler-HMM or Viterbi design)? --> Make the algorithm and explain it in clear slides; analysis of time and space and why it is better for our sequence
3. Check "equilibrium" (kth order probabilities remain consistent at different points of simulation). Rank-select-access (but how do we improve on existing?)
4. nothing worked and we try some of the small points in "Live goals"
