function [m,e] = lingmmm(b,infoz,stat,y,x,z)
% PURPOSE:  Provide moment conditions and error term for 
%             linear GMM estimation
%-------------------------------------------------------------------------
% USAGE:  [m,e] = lingmmm(b,infoz,stat,y,x,z,w)
%  b      model parameters
%  infoz   MINZ infoz structure
%  stat   MINZ status structure
%  y,x,z  Data: dependent, independent, and instruments
%  w      GMM weighting matrix
%-------------------------------------------------------------------------
% RETURNS:
%  m      vector of moment conditions
%  e      Model errors  (Nobs x Neq)
%-------------------------------------------------------------------------
% VERSION: 1.1.2

% written by:
% Mike Cliff,  Purdue Finance  mcliff@mgmt.purdue.edu
% Created: 12/10/98
% Modified 9/26/00 (1.1.1 Does system of Eqs)
%          11//13/00 (1.1.2 No W as input argument)
  
k = rows(b);
nx = cols(x);
neq = k/nx;
e = [];
if mod(k,nx) ~= 0
  error('Problem determining number of equations')
end

for i = 1:neq
  ei = y(:,i) - x*b((i-1)*nx+1:i*nx);
  e = [e ei];
end
       
m = vec(z'*e/rows(e));

