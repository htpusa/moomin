function geneSets = getGeneSets(model)

% auxiliary function for NSDE.m to get gene association as lists of genes from a COBRA model

	nReactions = size(model.rxns,1);
	geneSets = cell(nReactions,1);
	for reacInd=1:nReactions
		indsAsCell = regexp(model.rules{reacInd,1},'\d+','match');
		inds = [];
		if ~isempty(indsAsCell)
			for i=1:size(indsAsCell,2)
				inds = [inds; str2double(indsAsCell{1,i})];
			end
		end
		geneSets{reacInd,1} = inds;
	end