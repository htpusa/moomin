function [model, MILPsolutions, MILPproblem] = moomin(model, expression, varargin)

% Main function of the MOOMIN method. Formulates a MILP-problem and solves it to provide
% a hypothesis of a metabolic shift
%
% USAGE:
%
%	[model, MILPsolutions, MILPproblem] = moomin(model, expression, varargin)
%
% INPUTS:
%	model:				the metabolic model (COBRA model structure)
%	expression:			gene expression data (a structure with the fields 
%						'GeneID', 'PPDE', and 'FC')
%
% OPTIONAL INPUTS:
%   pThresh:			threshold for differential expression (default 0.9)
%	alpha:				alpha parameter of the weight function (default 3)
%	beta:				beta parameter of the weight function (default 2)
%	stoichiometry		boolean variable to choose if stoichiometry is considered (default 1)
%	removeReactions		a cell array containing a list of reactions to be removed from the model
%	enumerate			integer to determine the maximum number of alternative optimal
%						solutions that are enumerated  (default 1 ie only one solution)
%	precision			integer to determine up to how many significant numbers the weights are
%						calculated. In practice will influence the uniqueness of weight values.
%						(default 7)
%	subSystemThresh		if this option is given, subsystem information is used to further guide
%						the choice of reactions. Value given is the threshold: if a subsystem
%						is coloured above this percentage, grey reactions in that subsystem
%						are considered to have PPDE='pThresh'-0.1
%	solverParameters	a structure containing solver parameters to be passed on to
%						'solveCobraMILP'
%	solverTimelimit		time limit for the MILP solver (default 1000)
%
% OUTPUT:
%    model:				the metabolic model with additional fields containing inputs and
%						outputs of the algorithm
%						Solutions are in a matrix in a field called 'outputColours' where
%						every column is an optimal solution and rows correspond to reactions
%						The colours are coded as
%							2	reverse red
%							1	red
%							0	grey
%							-1	blue
%							-2	reverse blue
%							6	yellow (a priori grey reversible reaction)
%						Additional fields are:
%						'inputColours': the colours obtained without any inference
%						'weights': the reaction weights
%						'leadingGenes': the gene "responsible" for the colour and score of a reaction
%						'frequency': how often a reaction appears in an optimal solution
%						'combined': an attempted consensus between all optimal solutions
%									if the status differs between solutions, a 7 is entered
%						'expression': the expression data
%						
%	MILPsolutions		outputs of 'solveCobraMILP'
%	MILPproblem			the final MILP-problem solved (or was attempted to be solved)
%
%
% .. Author: - T.P.

	tic;

% default parameters
	pThresh = 0.9;
	alpha = 3;
	beta = 2;
	epsilon = 1;
	useStoichiometry = 1;
	enumerate = 1;
	precision = 7;
	useSubsystems = 0;
    solverParameters.emphasis.numerical = 0;
	tolerance = 1e-6;
	solverParameters.mip.tolerances.absmipgap = tolerance;
	solverParameters.mip.tolerances.mipgap = tolerance;
	solverParameters.mip.tolerances.integrality = tolerance;
	solverParameters.simplex.tolerances.feasibility = tolerance;
	solverParameters.simplex.tolerances.optimality = tolerance;
	solverTimelimit = 1000;
	
	model.ub(:) = 100;
    model.lb(model.lb~=0) = -100;
	
% read optional inputs

	if ~isempty(varargin)
		if rem(size(varargin,2),2) ~= 0
			error('Check optional inputs.');
		else
			for i=1:2:size(varargin,2)
				switch varargin{1,i}
					case 'pThresh'
						pThresh = varargin{1,i+1};
					case 'alpha'
						alpha = varargin{1,i+1};
					case 'beta'
						beta = varargin{1,i+1};
					case 'stoichiometry'
						useStoichiometry = varargin{1,i+1};
					case 'removeReactions'
						model = removeRxns(model, varargin{1,i+1});
					case 'enumerate'
						enumerate = varargin{1,i+1};
					case 'precision'
						precision = varargin{1,i+1};
					case 'subsystemThresh'
						useSubsystems = 1;
						subsystemThresh = varargin{1,i+1};
					case 'solverParameters'
						solverParameters = varargin{1,i+1};
					case 'solverTimelimit'
						solverTimelimit = varargin{1,i+1};
					otherwise
						error('Could not recognise optional input names.\nNo input named "%s"', varargin{1,i});
				end
			end
		end
	end
	
	disp('Determining gene colours and weights...');
	geneColours = (expression.PPDE > pThresh) .* sign(expression.FC);
	geneWeights = arrayfun(@(x)geneWeight(x,alpha,beta,pThresh),expression.PPDE);
	
	disp('Determining reaction colours and weights...');
	% get the indices of ass. genes from model.rules
	geneSets = getGeneSets(model);
	
	nReactions = size(model.rxns,1);
	reactionColours = zeros(nReactions,1);
	reactionWeights = zeros(nReactions,1);
	leadingGenes = zeros(nReactions,1);
	for reacInd=1:nReactions
		[~,assGeneInds,~] = intersect(expression.GeneID,model.genes(geneSets{reacInd,1}));
		if isempty(assGeneInds)
			reactionWeights(reacInd) = geneWeight(0,alpha,beta,pThresh);
			leadingGenes(reacInd,1) = 0;
		else
			assColours = geneColours(assGeneInds);
			assWeights = geneWeights(assGeneInds);
			% a contradiction
			if sum(assColours==1)>0 && sum(assColours==-1)>0
				reactionWeights(reacInd) = geneWeight(0.5,alpha,beta,pThresh);
				leadingGenes(reacInd,1) = 0;
			else
				reactionColours(reacInd) = sign(sum(assColours));
				[reactionWeights(reacInd), ind] = max(assWeights);
				if reactionColours(reacInd) ~= 0
					leadingGenes(reacInd,1) = assGeneInds(ind);
				else
					leadingGenes(reacInd,1)= 0;
				end
			end
		end
	end
	
	model.leadingGenes = leadingGenes;
	reactionWeights = round(reactionWeights,precision,'significant');
	
% Subsystem heuristic
	if useSubsystems
		disp('Using subsystem information...');
		subsystemsUnique = {};
		nMembers = [];
		for i=1:nReactions
		   index = find(strcmp(subsystemsUnique,model.subSystems{i,1}));
		   if isempty(index)
			  subsystemsUnique = [subsystemsUnique; model.subSystems{i,1}];
			  nMembers = [nMembers; 1];
		   else
			   nMembers(index) = nMembers(index) + 1;
		   end
		end
		if size(subsystemsUnique,1)<2
			warning('Using subsystem heuristic but only found 0 or 1 subsystems in the model.');
		end
		SSred = zeros(size(subsystemsUnique,1),1);
		SSblue = zeros(size(subsystemsUnique,1),1);
		for reac=1:nReactions
		   index = find(strcmp(subsystemsUnique,model.subSystems{reac,1}));
		   if reactionColours(reac)==1
			   SSred(index) = SSgreen(index) + 1;
		   elseif reactionColours(reac)==-1
			   SSblue(index) = SSred(index) + 1;
		   end          
		end
		SScol = SSred+SSblue;
		for reac=1:nReactions
			if reactionColours(reac)==0
				index = find(strcmp(subsystemsUnique,model.subSystems{reac,1}));
		      	if SScol(index)/nMembers(index) > subsystemThresh
			    	reactionWeights(reac) = geneWeight(pThresh-0.1,alpha,beta,pThresh);
				end
			end
		end
	end
	
	disp('Creating the original MILP...');
	
	nMetabs = size(model.S,1);
	optimum = -sum(abs(reactionWeights));

% With stoichiometric constraints
	if useStoichiometry
		ub = repmat(max(model.ub),nReactions,1);
		lb = repmat(-max(model.ub),nReactions,1);
		% impose a priori colours
		for i=1:nReactions
			if reactionColours(i)==1 && model.lb(i)==0
				lb(i) = 0;
			elseif reactionColours(i)==-1 && model.lb(i)==0
				ub(i) = 0;
			end
		end

		[i,j,v] = find([model.S;repmat(eye(nReactions),4,1)]);
	% Stoichiometry
		A = sparse(i,j,v,nMetabs+4*nReactions,3*nReactions);
	% x+=1 -> v>=epsilon
		A(nMetabs+1:nMetabs+nReactions,nReactions+1:2*nReactions) = diag(lb-epsilon);
	% x+=0 -> v<=0
		A(nMetabs+nReactions+1:nMetabs+2*nReactions,nReactions+1:2*nReactions) = diag(-ub);
	% x-=1 -> v<=-epsilon
		A(nMetabs+2*nReactions+1:nMetabs+3*nReactions,2*nReactions+1:3*nReactions) = diag(ub+epsilon);
	% x-=0 -> v>=0
		A(nMetabs+3*nReactions+1:nMetabs+4*nReactions,2*nReactions+1:3*nReactions) = diag(-lb);
	% x+ and x- cannot be 1 at the same time
		%A = [A; zeros(nReactions) eye(nReactions) eye(nReactions)];
	% place holder for optimality constraint
		A = [A; zeros(1,nReactions) reactionWeights' reactionWeights'];

		csense(1:nMetabs) = 'E';
		csense(nMetabs+1:nMetabs+nReactions) = 'G';
		csense(nMetabs+nReactions+1:nMetabs+2*nReactions) = 'L';
		csense(nMetabs+2*nReactions+1:nMetabs+3*nReactions) = 'L';
		csense(nMetabs+3*nReactions+1:nMetabs+4*nReactions) = 'G';
		%csense(nMetabs+4*nReactions+1:nMetabs+5*nReactions) = 'L';
		csense = [csense 'G'];

		c = [zeros(nReactions,1); reactionWeights; reactionWeights];

		b = [zeros(nMetabs,1); lb; zeros(nReactions,1); ub; zeros(nReactions,1)];
		%b = [b; ones(nReactions,1)];
		b = [b; optimum];

		ub = [ub; ones(2*nReactions,1)];
		lb = [lb; zeros(2*nReactions,1)];

		vartype(1:nReactions) = 'C';
		vartype(nReactions+1:3*nReactions) = 'B';
		
% Without stoichiometric constraints
	else
		ub = ones(nMetabs+2*nReactions,1);
		lb = ub-1;
		% impose a priori colours
		for i=1:nReactions
			if reactionColours(i)==1 && model.lb(i)==0
				ub(nMetabs+nReactions+i,1) = 0;
			elseif reactionColours(i)==-1 && model.lb(i)==0
				ub(nMetabs+i) = 0;
			end
		end   

		A = sparse(nReactions+3*nMetabs, nMetabs+nReactions*2);
	% x+ and x- cannot be 1 at the same time
		A(1:nReactions,nMetabs+1:end) = [eye(nReactions) eye(nReactions)];
	% if a connected arc is included, a node is included
		A(nReactions+1:nReactions+nMetabs,1:nMetabs) = -diag(sum(model.S~=0,2));
		A(nReactions+1:nReactions+nMetabs, nMetabs+1:end) = [model.S~=0 model.S~=0];
	% if a node is included, it has to have an outgoing arc
		A(nReactions+nMetabs+1:nReactions+2*nMetabs, 1:nMetabs) = -eye(nMetabs);
		A(nReactions+nMetabs+1:nReactions+2*nMetabs, nMetabs+1:nMetabs+nReactions) = model.S<0;
		A(nReactions+nMetabs+1:nReactions+2*nMetabs, nMetabs+nReactions+1:end) = model.S>0;
	% if a node is included, it has to have an incoming arc
		A(nReactions+2*nMetabs+1:end, 1:nMetabs) = -eye(nMetabs);
		A(nReactions+2*nMetabs+1:end, nMetabs+1:nMetabs+nReactions) = model.S>0;
		A(nReactions+2*nMetabs+1:end, nMetabs+nReactions+1:end) = model.S<0;
	% place holder for optimum
		A = [A; zeros(1,nMetabs) reactionWeights' reactionWeights'];

		csense(1:nReactions) = 'L';
		csense(nReactions+1:nReactions+nMetabs) = 'L';
		csense(nReactions+nMetabs+1:nReactions+2*nMetabs) = 'G';
		csense(nReactions+2*nMetabs+1:nReactions+3*nMetabs) = 'G';
		csense = [csense 'G'];

		c = [zeros(nMetabs,1); reactionWeights; reactionWeights];

		b = [ones(nReactions,1); zeros(3*nMetabs,1)];
		b = [b; optimum];

		vartype(1:nMetabs+2*nReactions) = 'B';
	end
	
	MILPproblem.A = A;
	MILPproblem.b = b;
	MILPproblem.c = c;
	MILPproblem.lb = lb;
	MILPproblem.ub = ub;
	MILPproblem.csense = csense;
	MILPproblem.vartype = vartype;
	MILPproblem.osense = -1;
	MILPproblem.x0 = [];
	
% Solve the MILP

	cont = 1;
	counter = 1;
	model.outputColours = [];
	MILPsolutions = {};

	while cont
	
	% With stoichiometric constraints
		if useStoichiometry
			fprintf('\nSolving MILP #%d with stoichiometric constraints...\n',counter);
			solution = solveCobraMILP(MILPproblem, 'timeLimit', solverTimelimit, 'printLevel', 3, solverParameters);
	% Without stoichiometric constraints
		else
			fprintf('\nSolving MILP #%d without stoichiometric constraints...\n',counter);
			solution = solveCobraMILP(MILPproblem, 'timeLimit', solverTimelimit, 'printLevel', 3, solverParameters);
		end
		
		toc
	
	% "Interpret" solution	
		if solution.stat==1
			outputColours = zeros(nReactions,1);
			if useStoichiometry
				outputColours(solution.int(1:nReactions)>1e-4) = 1;
				outputColours(solution.int(nReactions+1:end)>1e-4) = -1;
			else
				outputColours(solution.int(nMetabs+1:nMetabs+nReactions)>1e-4) = 1;
				outputColours(solution.int(nMetabs+nReactions+1:end)>1e-4) = -1;
			end
			outputColours = sign(outputColours);
			for i=1:nReactions
				if outputColours(i)==1 && reactionColours(i)==-1
					outputColours(i) = -2;
				elseif outputColours(i)==-1 && reactionColours(i)==1
					outputColours(i) = 2;
				elseif outputColours(i)~=0 && reactionColours(i)==0 && model.lb(i)<0
					outputColours(i) = 6;
				end
			end
		else
			disp('No optimal solution was found for the MILP.');
			outputColours = [];
		end
	
		model.outputColours = [model.outputColours outputColours];
		
		MILPsolutions = [MILPsolutions; solution];
		
		cont = solution.stat==1 && counter < enumerate;
		if counter==1
			if useStoichiometry
				MILPproblem.b(end,1) = solution.obj;
			else
				MILPproblem.b(end,1) = solution.obj;
			end
		end
	% add constraints for enumeration
		previousSol = outputColours~=0;
		if cont
			if useStoichiometry
				MILPproblem.A = [MILPproblem.A; zeros(1,nReactions) (2*previousSol-1)' (2*previousSol-1)'];
			else
				MILPproblem.A = [MILPproblem.A; zeros(1,nMetabs) (2*previousSol-1)' (2*previousSol-1)'];
			end
			MILPproblem.b = [MILPproblem.b; sum(previousSol)-1];
			MILPproblem.csense = [MILPproblem.csense 'L'];
		end
		counter = counter+1;
	end
	
	model.inputColours = reactionColours;
	model.weights = reactionWeights;
	
% count how often an arc appears in a solution
    model.frequency = sum(model.outputColours~=0,2)/size(model.outputColours,2);
    
% combine solutions
    combined = zeros(nReactions,1);
    for i=1:nReactions
        row = model.outputColours(i,:);
        if any(row)
            colours = row(find(row));
            if all(colours(1)==colours)
                combined(i,1) = colours(1);
            else
                combined(i,1) = 7;
            end
        end
    end
    model.combined = combined;
    
    model.expression = expression;
    