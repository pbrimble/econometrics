#------------------------------------------------------------------------------#
These two R files apply the Generic Machine Learning for Heterogeneous Treatment Effects
methodology from Chernozhukhov et al. 2019. The code has been modified from previous code written by Demirer.
Several modifications have been made to make the code more flexible, including the introduction of
k-fold cross-validation and separating the machine learning code and analysis code into two files.
#------------------------------------------------------------------------------#
Authors:  Paul Brimble - Mind and Behaviour Research Group,
          Centre for the Study of African Economies (CSAE), University of Oxford
Email:    paul.brimble@bsg.ox.ac.uk
#------------------------------------------------------------------------------#
Preamble for original code:
This program obtains empirical results for the paper "Generic Machine Learning Discovery
and Classification Analysis of Heterogenous Treatment Effects in Randomized Experiments"
by V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, I. FERNANDEZ-VAL
Authors:  V. CHERNOZHUKOV, M. DEMIRER, E. DUFLO, I. FERNANDEZ-VAL
#------------------------------------------------------------------------------#

The code is split into two files:
1) "gml-hte_ml" - machine learning estimation.
2) "gml-hte_an" - analysis using inputs from "gml-hte_ml".

The "gml-hte_ml" code returns an RData file to be used in the "gml_hte_an" code

The "gml_hte_an" code returns five outputs:
1) PLOTS: Two plots for each outcome variable reporting ATE and GATES (first plot is top 2 methods, second plot is all methods)
2) BLP: Latex table reporting ATE and heterogeneity loading factor (HET).
3) GATES: Latex table reporting the GATE for the most and least affected groups, and the difference.
4) CLAN: Latex table reporting the most affected and least affected averages for the CLAN variables and the difference.
5) BEST: Latex table reporting best machine learning methods.

ML_Functions.R code is also required.
