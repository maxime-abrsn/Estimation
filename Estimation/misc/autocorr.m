function Rho = autocorr(data,order)
%  rho = autocorr(data,order)
% data is a (T x k) matrix of data
% order is the highest number of autocorrelations [1]

if nargin == 1
  order = 1;
end

for i = 1:cols(data)
  temp = makelags(data(:,i),0,order,0);		% last 0 stops err check
  rho = corrcoef(temp);
  Rho(1:order,i) =   rho(2:order+1,1);
end
