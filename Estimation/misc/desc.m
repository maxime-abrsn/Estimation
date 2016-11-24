function results = desc(X)
% PURPOSE: calculate sample statistics for X
%---------------------------------------------------
% USAGE: results = desc(X,prtdum)
% where: X        data, organized in columns
%        prtdum   set to 0 to suppress printing 
%---------------------------------------------------
% RETURNS:
%  results   structure with fields
%    .mean
%    .std
%    .skew
%    .kurt
%    .rho
%    .n
%    .min
%    .max
%    .pct   percentiles: 1, 5, 10, 25, 50, 75, 90, 95, 99
% --------------------------------------------------
% SEE ALSO: prt_desc
%---------------------------------------------------

% WRITTEN BY: Mike Cliff,  UNC Finance,  mcliff@unc.edu
 
if nargin ~= 1
error('Wrong # of arguments to desc');
end;

pcts = [1;5;10;25;50;75;90;95;99];
k = cols(X);

for i = 1:k
   temp = nonan(X(:,i));
   results(i).mean = mean(temp);
   results(i).std  = std(temp);
   results(i).skew = skewness(temp);
   results(i).kurt = kurtosis(temp)-3;
   results(i).rho    = spacf(temp,1,1);   
   results(i).n    = rows(temp);
   results(i).min  = min(temp);
   results(i).max  = max(temp);
   results(i).pct  = prctile(temp,pcts);
end





































































































