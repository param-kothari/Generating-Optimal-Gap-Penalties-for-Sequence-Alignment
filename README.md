# Generating Optimal Gap Penalties for Sequence Alignment
We built a mathematical model for measuring the quality of Protein/DNA sequence alignments.
We implemented an algorithm to determine the best open and closed gap penalty values for given sequences which result in the best possible alignment.
Thereby, we are eliminating the need for user-defined values which are usually random with no basis to be used on the sequence alignment DP algorithm by Needleman and Wunsch.

The 'needleman_wunsch_affine.py' file can be executed with 2 strings as inputs which will output the best possible alignment with parameter values for the particular strings based on BLOSUM62 matrix.
