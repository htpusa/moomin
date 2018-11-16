# MOOMIN
MOOMIN is a tool for analysing differential expression data. It takes as its input a metabolic network and a the results of
a DE analysis: a posterior probability of differential expression and a (logarithm of a) fold change for a list of genes.
It then forms a hypothesis of a metabolic shift, determining for each reaction its status as "increased flux",
"decreased flux", or "no change". These are expressed as colours: green for an increase, red for a decrease, and grey for no
change. A paper is in preparation describing the theoretical framework.

# Dependencies
MOOMIN runs in Matlab (developped in R2016a) and relies on the COBRA Toolbox. Additionally, a MILP-solver compatible with
COBRA TB is needed. To obtain COBRA TB and how to use it:
https://opencobra.github.io/cobratoolbox/stable/

# Usage
