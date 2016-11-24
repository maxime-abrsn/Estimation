%  GMM  function library.
%
%	Written by:  Mike Cliff,  UNC Finance,  mcliff@unc.edu
%
%	Uses the MINZ libarary for minimization
%
%	gmm		Main program controlling estimation
%	gmmS		Creates GMM weighting matrix
%	andmon		Andrews-Monahan routine for gmmS
%
%	lingmmm		Moment conditions for standard linear model
%	lingmmj		Jacobian for linear model (deriv of moments)
%	lingmmh		Hessian for linear model (second deriv of obj)
%
%       msdm		Moments for mean/cov estimates of a sample
%	msdj		Jacobian accompanying msdm
%
%	prt_gmm		Prints GMM results
%	drawobj		Plots objective function
%	diagw		Plots W matrix for diagnostics
