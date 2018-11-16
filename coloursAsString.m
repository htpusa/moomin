function colours = coloursAsString(colourVector)

% auxiliary function for writeOutput.m to turn reaction colour numbers into string format

	colours = {};
	
	for ind=1:size(colourVector,1)
		switch colourVector(ind)
		case 2
			colours = [colours;{'r.green'}];
		case 1
			colours = [colours;{'green'}];
		case 0
			colours = [colours;{'grey'}];
		case -1
			colours = [colours;{'red'}];
		case -2
			colours = [colours;{'r.red'}];
		case 6
			colours = [colours;{'yellow'}];
		end
	end