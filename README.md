# MOOMIN
MOOMIN (Mathematical explOration of Omics data on a MetabolIc Network) is a tool for analysing differential expression data. It takes as its input a metabolic network and a the results of a DE analysis: a posterior probability of differential expression and a (logarithm of a) fold change for a list of genes.
It then forms a hypothesis of a metabolic shift, determining for each reaction its status as "increased flux",
"decreased flux", or "no change". These are expressed as colours: green for an increase, red for a decrease, and grey for no
change. A paper is in preparation describing the theoretical framework.

# Dependencies
MOOMIN runs in Matlab (developped in R2016a) and relies on the COBRA Toolbox. Additionally, a MILP-solver compatible with
COBRA is needed. To obtain COBRA and how to use it:
https://opencobra.github.io/cobratoolbox/stable/

# Usage
In order to use MOOMIN, you need DE results obtained using Bayesian methods. In other words, a posterior probability of differential expression (PPDE) is needed instead of the more common p-value.
You also need a metabolic network of the organism under study. You can read an SBML-file (.xml) using the COBRA-function
"readSBML" or you can download a Matlab-structure containing a COBRA model directly if one is available.

COBRA Toolbox needs to be in the Matlab path.

A MILP-solver needs to be set for COBRA.

The file "example.m" is a Matlab script that contains an example pipeline. It can be run with simply "example". Note that the script will clear your Matlab workspace.

The E.coli core network used in the example was obtained from
http://systemsbiology.ucsd.edu/Downloads/EcoliCore
