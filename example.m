clear;

% A MILP-solver needs to be set for COBRA TB
% e.g. changeCobraSolver('ibm_cplex','milp');

% load the model and the expression data
model = readSBML('example/Ec_core_flux1.xml',10);
load('example/expression');

% read a list of inactive genes
inactiveGenes = tdfread('example/inactiveGenes.txt');
inactiveGenes = cellstr(inactiveGenes.Gene);

% remove inactive reactions
inactiveReactions = findRxnsInActiveWithGenes(model, inactiveGenes);
model = removeRxns(model, inactiveReactions);

% solve for one solution using stoichiometric contraints
[model_stoich, MILPsolution_stoich, MILPproblem_stoich] = moomin(model, expression);

% solve for up to 10 solutions with only topological contraints
[model_topo, MILPsolution_topo, MILPproblem_topo] = moomin(model, expression, 'stoichiometry', 0, 'enumerate', 10);

% write the standard output file
writeMoominOutput(model_stoich, 'example/output_stoich.txt');

% write the standard output for input colours
writeMoominOutput(model_topo, 'example/input_topo.txt','type','input');