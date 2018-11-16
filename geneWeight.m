function weight = geneWeight(p,alpha,beta,t)

% auxiliary function for moomin.m to generate gene weights

	if p==1
		weight = -alpha*beta*log(1-t);
	else
		weight = min(beta*(-log(1-p)+log(1-t)), -alpha*beta*log(1-t));
	end