# Introduction
This repository corresponds to our final project for Computational Genomics: Sequences. The team consists of Tad Berkery, Xinlei (Lily) Chen, Richard Hu, and Yuseong (Nick) Oh. We outline in our report our review of prior work, methodologies, and results. What follows are some auxiliary instructions to complement this on how to run our code with the goal of aiding with reproducibility.

# Protein Fold Clustering
A key source of inspirating for this project is code that generates a sequence of numbers indicating clusters representing protein folds. This code is written in MatLab and is primarily contained in the "msm_code_Richard", "clustering_code", "Lagnevin_clustering" and "Langevin dynamics code" subfolders. This code is notably outside the scope of our project and comes from a biophysics class taken by Richard and Lily. This code is what produced the file text.txt within the BWT_HMM subfolder, which is a sequence of 1,000,000 numbers representing cluster assignments at each position. We included this code for comprehensiveness, but, in relation to reproducing this project, the code associated with our deliverables begins at the level of treating this sequence file as already existent and working from there. We then created the "generate_sequence" directory and code inside it to provide a mechanism for generating additional such sequences.

# Several Algorithms
With the above context, we then proceeded to implement several architectures, headlined by the Burrows-Wheeler-Transform Hidden-Markov-Model (BWT-HMM folder), Burrows-Wheeler-Transform Hidden-Markov-Model with merging (BWT-HMM_Merging folder), and MLSE Viterbi, which is spread across the "Compress_sequence" and "main_method" subdirectories. The best way to see our underlying implementations is to look at the code in these directories. The code in these directories is a mix of code that can be executed by running a given script and code that works together across many scripts.

# Benchmarking
A key deliverable for our project beyond any architecture implementation itself is a mini suite of software for benchmarking performance in terms of time and space associated with different algorithms when used on different sequences. This software is built out in the "benchmark" subdirectory, which contains Python scripts that methodically run all of the algorithm implementations referenced above in a very methodical, trackable manner. The most important file, which is executed to obtain the benchmarking data, is the `profile.py` file. Note that the general structure of this file at a high-level is that it imports functions for each architecture of interest that take only the sequence to run the architecture on as its only argument, wraps these functions in wrapper functions critically annotated with `@profile` toe neable use of the `memory_profiler` Python package (which you will need to install)

# Miscellaneous
* For the code in the "BWT_HMM" and "BWT_HMM_Merging" directories, much of the Python code (such as `bwthmm.py`) relies on the [`hmmlearn`](https://pypi.org/project/hmmlearn/) package. This package proved to be quite difficult to install on several occasions, perhaps likely because it is self-described by its creators as being "under limited-maintenance mode." If you have issues with it, please ensure that you have Python >= 3.6, NumPy >= 1.10, and scikit-learn >= 0.16. Of note, `hmmlearn` is a direct offshoot of the origin scikit-learn hidden markov model package, which is now deprecated "due to it no longer matching the scope and the API of the project" and "scheduled for removal in the 0.17 release of the project". scikit-learn recommends using hmmlearn as we did, and we did not want our code to depend on a package set to soon go out of existence, so we used hmmlearn and endured its often painful installation difficulties. If this poses issues on a Mac, try [here](https://github.com/hmmlearn/hmmlearn/issues/475). If this poses issues on a PC, try [here](https://stackoverflow.com/questions/51002441/unable-to-install-hmmlearn-in-python-3).
* The "kmerindex" directory contains code that we used at times to contextualize how different architectures fare in relation to a generic, typical kmer index.

# List of Python Package Dependencies
(beyond standard packages)
* `pandas`
* `numpy`
* `memory_profiler`
* `hmmlearn`
* `random`
* `copy`
* `matplotlib`
* `bisect`
* `networkx`
* `collections`

# List of R Package Dependencies
(beyond standard packages)
* `tidyverse`
* `ggplot2`
* `gridExtra`
