function M = lingmmj(b,infoz,stat,y,x,z)
% PURPOSE:  Provide derivative of moment conditions (Jacobian) for
%             linear GMM estimation
%-------------------------------------------------------------------------
% USAGE:  M = lingmmj(b,infoz,stat,y,x,z)
%  b      k-vector of model parameters
%  infoz   MINZ infoz structure
%  stat   MINZ status structure
%  y,x,z  Data: dependent, independent, and instruments
%-------------------------------------------------------------------------
% RETURNS:
%  M      Jacobian  (rows(m) x k)
%-------------------------------------------------------------------------
%VERSION: 1.1.2

% written by:
% Mike Cliff,  Purdue Finance  mcliff@mgmt.purdue.edu
% Created: 12/10/98
% Updated: 9/26/00 (1.1.1 Does System of Eqs)
%          11//13/00 (1.1.2 No W as input argument)
  
k = rows(b);
nx = cols(x);
neq = k/nx;
e = [];
if mod(k,nx) ~= 0
  error('Problem determining number of equations')
end

M = kron(eye(neq),-z'*x/rows(x));
