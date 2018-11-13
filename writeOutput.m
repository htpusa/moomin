function writeOutput(model, fileName, varargin)

% Print the output of moomin.m in various formats
%
% USAGE:
%
%	writeOutput(model, fileName, varargin)
%
% INPUTS:
%	model:			input model (COBRA model structure) as given by the function NSDE.m
%	filename:		name of the output file
%
% OPTIONAL INPUTS
%	nSolution:		number of the solution that is output (default 1)
%	type:			customise output format
%					The options are
%						'full': a table with additional information about the genes 
%						'json': a .json file is produced to be used with apps such as Escher
%						'EC': output is written using EC numbers
%						'input': only input colours are written
%
% .. Author: - T.P.

% default parameters
	nSolution = 1;
	full = 0;
	json = 0;
	EC = 0;
	input = 0;
	
% read optional inputs
	if ~isempty(varargin)
		if rem(size(varargin,2),2) ~= 0
			error('Check optional inputs.');
		else
			for i=1:2:size(varargin,2)
				switch varargin{1,i}
					case 'nSolution'
						nSolution = varargin{1,i+1};
						if nSolution > size(model.outputColours,2)
							error('There are only %d solutions.', size(model.outputColours,2));
						end
					case 'type'
						switch varargin{1,i+1}
							case 'full'
								full = 1;
							case 'json'
								json = 1;
							case 'input'
								input = 1;
							case 'EC'
								EC = 1;
							otherwise
								error('Could not recognise type.\nNo type named "%s"', varargin{1,i+1})
						end
					otherwise
						error('Could not recognise optional input names.\nNo input named "%s"', varargin{1,i});
				end
			end
		end
	end
	
	if EC
		IDlist = getECnumbers(model);
	else
		IDlist = cellfun(@(x) strcat('R_',x), model.rxns,'UniformOutput',false);
	end

	inputAsString = coloursAsString(model.inputColours);
	outputAsString = coloursAsString(model.outputColours(:,nSolution));
	
	if full
		[weights,sortByWeight] = sort(model.weights,'descend');
		ID = IDlist(sortByWeight);
		output=model.outputColours(sortByWeight,:);
		input = model.inputColours(sortByWeight);
		leadingGenes = model.leadingGenes(sortByWeight);
		leadingGenesNames = {};
		FC = {};
		for i=1:size(leadingGenes,1)
			if leadingGenes(i) > 0
				leadingGenesNames = [leadingGenesNames; model.expression.GeneID(leadingGenes(i),:)];
				FC = [FC; model.expression.FC(leadingGenes(i),1)];
			else
				leadingGenesNames = [leadingGenesNames; 'NA'];
				FC = [FC; 'NA'];
			end
		end
		reactionName = model.rxnNames(sortByWeight);
		model = creategrRulesField(model);
		GPR = model.grRules;
		subsystem = model.subSystems;
		%EC = getECnumbers(model);
		%EC = EC(sortByWeight);
		outputTable = table(ID,reactionName,output,input,weights,GPR,subsystem,leadingGenesNames,FC);
		writetable(outputTable,fileName,'Delimiter','\t');
	elseif json
		jsonStr = '{';
		for reacInd=1:size(model.rxns,1)-1
			if input
    			jsonStr = [jsonStr '"' model.rxns{reacInd,1} '": ' num2str(model.inputColours(reacInd)) ', '];
    		else
    			jsonStr = [jsonStr '"' model.rxns{reacInd,1} '": ' num2str(model.outputColours(reacInd,nSolution)) ', '];
    		end
		end
		if input
			jsonStr = [jsonStr '"' model.rxns{end,1} '": ' num2str(model.inputColours(end)) '}'];
		else
			jsonStr = [jsonStr '"' model.rxns{end,1} '": ' num2str(model.outputColours(end,nSolution)) '}'];
		end
		fileID = fopen(fileName,'w');
		fprintf(fileID,jsonStr);
		fclose(fileID);
	else
		ID = {};
		outputColour = {};
		inputColour = {};
		for reacInd=1:size(model.rxns,1)
			if input && model.inputColours(reacInd) ~= 0
				ID = [ID; IDlist{reacInd,1}];
				inputColour = [inputColour; inputAsString{reacInd,1}];
			elseif ~input && model.outputColours(reacInd,nSolution) ~= 0
				ID = [ID; IDlist{reacInd,1}];
				outputColour = [outputColour; outputAsString{reacInd,1}];
			end
		end
		if input
			outputTable = table(ID,inputColour);
		else
			outputTable = table(ID,outputColour);
		end
		writetable(outputTable,fileName,'Delimiter','\t');
	end
	
	
	
	