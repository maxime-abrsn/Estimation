function [m,e] = msdm(b,infoz,stat,y,x,z,w)
% PURPOSE:  Provide moment conditions and error term for Mean/Covar Estimates
%-------------------------------------------------------------------------
% USAGE:  [m,e] = msdm(b,infoz,stat,y,x,z,w)
%  b      parameters: k-means, then vech(cov(y)), 
%  infoz   MINZ infoz structure
%  stat   MINZ status structure
%  y      data of interest
%  x      not used
%  z      vector of ones, same # rows as y
%  w      GMM weighting matrix
%-------------------------------------------------------------------------
% RETURNS:
%  m      vector of moment conditions
%  e      Model errors  (Nobs x Neq)
%-------------------------------------------------------------------------
% VERSION: 1.1 (6/6/00)
  
% written by:
% Mike Cliff,  Purdue Finance  mcliff@mgmt.purdue.edu
% Created: 6/6/00

[T,k] = size(y);
e = y - repmat(b(1:k)',T,1);
ct = k;
for i = 1:k
  for j = i:k
    ct = ct + 1;
    e1 = e(:,i).*e(:,j) - b(ct);
    e = [e e1];
  end
end
m = e'*z/T;