function [updatedModel inactiveReactions] = removeInactive(model, inactiveGenes)

% Removes reactions from a COBRA model based on a list of inactive genes.
%
% USAGE:
%
%	[updatedModel inactiveReactions] = removeInactive(model, inactiveGenes)
%
% INPUTS:
%	model				the metabolic model (COBRA model structure)
%	inactiveGenes		a cell array containing a list of inactive genes
%
% OUTPUTS:
%	updatedModel		the input model with inactive reactions removed
%							the IDs of the removed rxns are stored in an
%							additional field 'inactiveRxns'
%	inactiveReactions	a cell array that contains the IDs of the removed rxns

	GPR = GPRparser(model);
	
	inactiveReactions = {};
	
	for reacInd=1:size(model.rxns,1)
		geneAss = GPR{reacInd};
		inactive = [];
		if ~isempty(geneAss)
			for col=1:size(geneAss,2)
				assGenes = intersect(model.geneNames, geneAss{1,col});
				inactive = [inactive any(ismember(assGenes,inactiveGenes))];
			end
		end
		if ~isempty(inactive) && all(inactive)
			inactiveReactions = [inactiveReactions; model.rxns(reacInd)];
		end
	end
	
	updatedModel = removeRxns(model, inactiveReactions);
	updatedModel.inactiveRxns = inactiveReactions;