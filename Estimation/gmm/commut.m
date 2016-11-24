function K = commut(m,n);
%function K = commut(m,n);
%
% K is the mn x mn commutation matrix that transforms 
% vec(A) into vec(A'):
%   K*vec(A) = vec(A')


% Written by Mike Cliff,  mcliff@unc.edu, 10/12/99

K = zeros(m*n,m*n);
H0 = zeros(m,n);

for i = 1:m
  for j = 1:n
    H = H0;
    H(i,j) = 1;
    K = K + kron(H,H');
  end
end

