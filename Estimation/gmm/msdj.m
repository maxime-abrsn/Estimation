function [M,e] = msdj(b,infoz,stat,y,x,z,w)
% PURPOSE:  Provide Jacobian for Mean/StdDev Est
%-------------------------------------------------------------------------
% USAGE:  M = msdj(b,infoz,stat,y,x,z,w)
%  b      model parameters
%  infoz  MINZ infoz structure
%  stat   MINZ status structure
%  y,x,z  Data: dependent, independent, and instruments
%  w      GMM weighting matrix
%-------------------------------------------------------------------------
% RETURNS:
%  M      Jacobian  (k+k(k+1)/2 by k+k(k+1)/2) (k is cols(y)
%-------------------------------------------------------------------------
% VERSION: 1.2 (5/20/05)
  
% written by:
% Mike Cliff,  Virginia Tech Finance  mcliff@vt.edu
% Created: 6/6/00
% Updated: 5/20/05

[T,k] = size(y);
kk = k*(k+1)/2;
M = zeros(k+kk,k+kk);
M(1:k,1:k) = -eye(k);
M(k+1:k+kk,k+1:k+kk) = -eye(kk);
e = y - repmat(b(1:k)',T,1);
ct = k;
for i = 1:k
  for j = i:k
    ct = ct + 1;
    if i == j
      M(ct,i) = -2*e(:,i)'*z/T;
    else
      M(ct,i) = -e(:,j)'*z/T;
      M(ct,j) = -e(:,i)'*z/T;
    end
  end
end
