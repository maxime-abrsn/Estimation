function result = arma(y,arp,maq,cflag,maxit)
% PURPOSE: computes Box-Jenkins ARMA models
%---------------------------------------------------
% USAGE: result = arma(y,arp,maq,cflag,maxit)
%     or result = arma(y,arp,maq,cflag)
% where: y = time series vector (nobs x 1)
%      arp = autoregressive order
%      maq = moving average order
%    cflag = 0,1 flag for constant term, 0 = no constant, 1 = constant
%            (default = 0)
%    maxit = maximum # of iterations (default = 50)
%---------------------------------------------------
% RETURNS: a structure
%        result.meth   = 'arma'
%        result.beta   = (ar, cterm, ma) estimates
%        result.tstat  = t-stats (order = ar, cterm, ma)
%        result.yhat   = yhat (1-step-ahead predictions)
%        result.resid  = residuals (1-step-ahead errors)
%        result.sige   = e'*e/(n-k)
%        result.rsqr   = r-squared
%        result.rbar   = r-bar squared
%        result.like   = log-likelihood function value
%        result.iter   = # of iterations
%        result.nobs   = (nobs unadjusted for AR term lags)
%        result.nvar   = nvars (#ar + #ma + 1)
%        result.y      = y data vector
%        result.ar     = # ar parameters
%        result.ma     = # ma parameters
%        result.cflag  = 0,1 flag value for constant term
%
%        result.eflag  = 1 if an error, else 0
% --------------------------------------------------
% NOTES: robust estimation is used with errors > 1.6*std(resid)
%           converted to linear loss rather than quadratic loss
%	AR done via OLS
% --------------------------------------------------
% SEE ALSO: prt_reg(result), plt_reg(result), armaf
%---------------------------------------------------
% REQUIRES: System Identification Toolbox

%
% Written by J. LeSage
% Modified by Mike Cliff,  Purdue Finance,  mcliff@mgmt.purdue.edu
% UPDATED: 11/6/99 (ols for AR models)
%          8/17/01 (new armax)
%          9/14/01 (don't use se in ols to avoid using my version of OLS.M)
  
if nargin > 5
  error('Wrong # of arguments to arma');
elseif nargin == 4
  maxit = 50; % default maximum iterations
  freq=1;
elseif nargin == 3
  maxit = 50;
  cflag = 0;
elseif nargin <= 2
  error('Wrong # of arguments to arma');
end;

nobs = length(y);

result.eflag = 0;
result.meth = 'arma';
result.nobs = nobs;
result.y = y;
result.nvar = arp+maq+cflag;
result.ma = zeros(maq,1);
result.ar = zeros(arp,1);
result.cterm = 0;

na = arp;
nc = maq;

if cflag == 1
  nb = 1;
  iota = ones(nobs,1);
  z = [y iota];
  nn = [na nb nc 0];    
else
  z = y;
  nn = [na nc];
end;

% a call to the system identification routine armax or ols
if maq == 0
  Y = y(1+arp:end);
  X = [];
  for i = 1:arp
    X = [X y(1+arp-i:end-i)];
  end
  X = [X ones(rows(Y),cflag)];
  theta = ols(Y,X);
  result.beta = -theta.beta;
  result.beta(end) = theta.beta(end);
  result.se = sqrt(theta.sige*diag(inv(X'*X)));
  result.tstat = -theta.tstat;
  result.tstat(end) = theta.tstat(end);
  sige = theta.sige;
  itinfo = 1;
else
%[theta, itinfo] = armax(z,nn,maxit);  
  M = armax(z,nn,'MaxIter',maxit);
  itinfo = get(M,'EstimationInfo');
  itinfo = itinfo.Iterations;
%  [par, p, lam] = th2par(theta);
  result.beta = get(M,'ParameterVector');
  result.se = sqrt(diag(get(M,'CovarianceMatrix')));
  result.tstat = result.beta./result.se;
  sige = get(M,'NoiseVariance');
end

result.iter = itinfo(1);

if itinfo(1) == maxit
  warning('No convergence in arma --- increase maxit');
  result.eflag = 1;
  return
end;

if arp > 0
arpar = zeros(arp,1);
result.ar = arpar;
end;
if maq > 0
mapar = zeros(maq,1);
result.ma = mapar;
end;

for i=1:arp
  result.beta(i,1) = -1*result.beta(i,1);
  result.tstat(i,1) = -1*result.tstat(i,1);
  arpar(i,1) = result.beta(i,1);
  result.ar(i) = arpar(i,1);
end;

if cflag > 0
cterm = result.beta(arp+1,1);
result.cterm = cterm;
        for i=1:maq
        mapar(i,1) = result.beta(arp+1+i,1);
        result.ma(i) = mapar(i,1);
        end;
else
cterm = 0;
result.cterm = cterm;
        for i=1:maq
        mapar(i,1) = result.beta(arp+i,1);
        result.ma(i) = mapar(i);
        end;
end;


% find yhat, residuals, rsqred, sige
if maq == 0
  e = theta.resid;
  e = [y(1:arp) - mean(y); e];
else
%  e = pe(z,theta);        % residuals
  e = pe(M,z);        % residuals
end
yhat = y - e;

skip = max(arp,maq);

% compute log likelihood function value
epe = e'*e;
nskip = nobs-skip;

like = -(nskip/2)*log(2*pi) - (nskip/2)*log(sige) - epe/(2*sige);

result.sige = sige;
result.resid = e;
result.yhat = yhat;
result.like = like;

ym = y(skip+1:nobs,1) - ones(nobs-skip,1)*mean(y(skip+1:nobs,1));
rsqr1 = epe;
rsqr2 = ym'*ym;
result.rsqr = 1.0 - rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(nobs-result.nvar);
rsqr2 = rsqr2/(nobs-1.0);
result.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
