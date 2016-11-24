function prob = stdn_prb(x)
% PURPOSE: computes the standard Normal probability function
%---------------------------------------------------
% USAGE: prob = stdn_prb(x,v)
% where: x = the value to test
%            (may be a matrix size(v), or a scalar)
%---------------------------------------------------
% RETURNS:
%        prob = the probability of observing a standard 
%               normal value > |x|, i.e., prob(value > |x|).  
%                
% --------------------------------------------------
% SEE ALSO:stdn_pdf, stdn_cdf, stdn_inv, stdn_rnd
%---------------------------------------------------

% Written by: Mike Cliff,   UNC Finance   mcliff@unc.edu
% Created  12/23/98

if nargin ~= 1
error('Wrong # of arguments to stdn_prb');
end;

rtail = ones(size(x)) - normcdf(abs(x));   % 1 - Prob(stdn < +|x|)
ltail = normcdf(-abs(x));                  % Prob(stdn < -|x|)

prob = rtail + ltail;

