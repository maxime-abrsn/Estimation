function Y=makelags(X,a,b,hush)
%function makelags(X,a,b)
%
% MAKELAGSS     Create lags a through b of a matrix X 
%
%       X       (T x k) Matrix
%       a       first lag, a >= 0 
%       b       last lag,  a <=b <= T
%	hush	Dummy argument for error checking
%
%       Example, for X = [1:10 10:10:100]', makelags(X,0,2) gives
%
%
%		     3    30     2    20     1    10
%		     4    40     3    30     2    20
%		     5    50     4    40     3    30
%		     6    60     5    50     4    40
%		     7    70     6    60     5    50
%		     8    80     7    70     6    60
%		     9    90     8    80     7    70
%		    10   100     9    90     8    80


if nargin == 3
  disp('Please check use of MAKELAGS')
end

if b < a disp('ERROR: INCREASE LAGS. '); end

Y=[];

for i = a:b
  Y = [X(i+1:end+i-b,:) Y];
end
