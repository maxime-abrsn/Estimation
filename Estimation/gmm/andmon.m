function [band,uout,D,eflag] = andmon(gmmopt,e,M,Z)
%function [band,uout,D,eflag] = andmon(gmmopt,e,M,Z)
%
% ANDREWS-MONAHAN HAC ESTIMATOR
%
% Calculate objects needed by gmmS for Andrews-Monahan HAC estimate.
% Basic idea is to take residuals from normal estimation then fit a
% simple time-series model such as ARMA(1,1), AR(1), or MA(q), then use
% these 'whitened' residuals to calculate the Spectral Density matrix.
% The parameters of the auxillary time-series model are used to 'recolor'
% the Spectral Density matrix of the whitened residuals.
%
% Another key feature is the use of the Quadratic-Spectral kernel with
% 'automatic' bandwidth selection.  This procedure still requires choosing
% a weighting matrix to the relative importance of each of the disturbance
% vectors.  I use the inverse of the standard deviation of the residuals
% as weights in attempt to make the bandwidth selection scale-invariant.
% This effort seems to be partially successful.
%
% The procedure can also be used to generate the automatic bandwidth,
% but without the pre-whitening. This is the Andrews procedure, and is
% implemented by setting gmmopt.aminfo.nowhite = 1.
%
% INPUTS:
%  gmmopt	structure from GMM package.  Relevant fields are:
%    gmmopt.prt      Printing (0=none, 1=screen,else file)	[1]
%    gmmopt.aminfo   another structure with fields:
%      aminfo.p        autoregressive lag (p<=1)		[1]
%      aminfo.q        moving average lags			[0]
%      aminfo.vardum   1 for VAR(1) (overrides p and q)		[0]
%      aminfo.kernel   Weighting scheme for lags		['QS']
%      aminfo.nowhite  1 suppresses pre-whitening		[0]
%      aminfo.diagdum  1 for diagnostic graphs			[0]
%
%  e		matrix of residuals -- we want its Spectral Density
%  M		Jacobian (derivative of moments wrt parmeters)
%  Z		Instruments
%
% OUTPUTS:
%  band		bandwidth determined from the data
%  uout		matrix of whitened residuals
%  D		Matrix to 're-color' Spectral Density of whitened resids
%  eflag	=1 if an error in estimating ARMA model, else 0

% WRITTEN BY:  Mike Cliff, Purdue Finance,  mcliff@mgmt.purdue.edu
% VERSION 1.1  
% CREATED: 11/6/99

%========================================================================
%   INITIALIZATIONS
%========================================================================

% --- Set Defaults ------------------------------------------------------
if ~isfield(gmmopt,'prt'),  gmmopt.prt = 1; end
if ~isfield(gmmopt,'aminfo'),  gmmopt.aminfo.p = 1; end
if ~isfield(gmmopt.aminfo,'kernel'), gmmopt.aminfo.kernel = 'QS'; end
if ~isfield(gmmopt.aminfo,'p'), gmmopt.aminfo.p = 1; end
if ~isfield(gmmopt.aminfo,'q'), gmmopt.aminfo.q = 0; end
if ~isfield(gmmopt.aminfo,'vardum'), gmmopt.aminfo.vardum = 0; end
if ~isfield(gmmopt.aminfo,'nowhite'), gmmopt.aminfo.nowhite = 0; end
if ~isfield(gmmopt.aminfo,'diagdum'), gmmopt.aminfo.diagdum = 0; end

eflag = 0;
aminfo = gmmopt.aminfo;
p = aminfo.p;
q = aminfo.q;
vardum = aminfo.vardum;
diagdum = aminfo.diagdum;
nobs = rows(e);

if p > 0
  if p > 1, error('ANDMON written for p<=1 in ARMA(p,q)/VAR(p)');  end
  if q > 1, error('ANDMON written for ARMA(p,q), q <=1 when p = 1'); end
end

switch aminfo.kernel
  case 'H',  const = 0.6611; ktype = 2; kname = 'Truncated (Hansen)';
  case 'NW', const = 1.1447; ktype = 1; kname = 'Bartlett (Newey-West)';
  case 'G',  const = 2.6614; ktype = 2; kname = 'Parzen (Gallant)';
  case 'QS', const = 1.3221; ktype = 2; kname = 'Quadratic-Spectral';
  case 'TH', const = 1.7462; ktype = 2; kname = 'Tukey-Hanning';
  otherwise, error('Need kernel in ANDMON');
end


% --- Get Residuals times Instruments ------------------------------------
u = [];
for i = 1:cols(e)
  u = [u repmat(e(:,i),1,cols(Z)).*Z];
end
north = cols(u);  
In = eye(north);


%========================================================================
%   VAR(1) to WHITEN RESIDUALS
%========================================================================
if vardum == 1
  error('VAR(1) Not Correctly Implemented')
  [varout,Phi,Omega] = varx(u,1);
  Phi = Phi(:,2:end);				% Drop intercept
  PHI = Phi;
  THETA = zeros(north,north);
  term = (In - Phi)\In;
  f = term*Omega*term;
  if ktype == 1
    H = [];
    for j = 1:128				% Truncate infinite sum
      H = H + Phi^j*Omega*(Phi')^j;
    end
    H = term^2*Phi*H;
    fq = H + H';
  else
    fq = Phi*Omega + Phi^2*Omega + Phi^2*Omega*Phi';
    fq = fq + fq' - 6*Phi*Omega*Phi';
    fq = term^3*fq*(term')^3;
  end
  f = f/(2*pi);
  fq = fq/(2*pi);
  K = commut(north,north);
  nvar = cols(M);
  w = 1./sqrt(diag(Omega));			% Calc weighting matrix
  W = diag(vec(w*w'));

%  W = ones(nvar^2,1);
%  W([1:nvar+1:nvar^2]) = 2;
%  W = diag(W);  
%  W = kron(M',M')'*W*kron(M',M');
%   w = diag(W);
  alpha = 2*vec(fq)'*W*vec(fq);
  alpha = alpha/trace(W*(eye(north^2)+K)*kron(f,f));
  for i = 1:north
    uout(:,i) = varout(i).resid;
  end
  uout = [zeros(1,north); uout];		% Initial values
  msg = '\n VAR(1) ';
   
else
%========================================================================
%   ARMA(1,1), AR(1), or MA(q) for WHITENING
%========================================================================

  THETA = zeros(north,north);
  PHI = zeros(north,north);
  w = 1./std(u)';				% Wt to standardize
  for i = 1:north
    tempout = arma(u(:,i),p,q,1,100);		% Estimate the ARMA model
    if tempout.eflag == 1
      eflag = 1; band = 0; uout = 0; D = 0;
      return
    end
    if p == 0					% then package estimates
      phi(i,1) = 0;
    elseif p == 1
      phi(i,1) = tempout.beta(1);
    end
    if q == 0
      theta(i,1) = 0;
    else
      theta(i,:) = tempout.beta(p+2:end)';
    end
    sigma(i,1) = sqrt(tempout.sige);
    uout(:,i) = tempout.resid;
    uout(1:max(p,q),i) = 0;			% Set initial cond to 0
    PHI(i,i) = sum(phi(i,:)');
    THETA(i,i) = sum(theta(i,:)');
        
% --- ARMA(1,1) or AR(1) Case --------------------------------------------
    if p == 1
      num(i,1) = 4*(1+phi(i)*theta(i))^2*(phi(i)+theta(i))^2;
      if ktype == 1
        term = (1-phi(i)^6)*(1+phi(i)^2);
      else
        term = (1-phi(i))^8;
      end
      num(i,1) = num(i)/term;
      den(i,1) = (1+theta(i))^4/(1-phi(i))^4;
      if q == 0
        msg = ['\n AR(' int2str(p) ') ']; 
      else
        msg = ['\n ARMA(' int2str(p) ',' int2str(q) ') '];
      end
% --- MA(q) Case ---------------------------------------------------------
    else
      tempnum = 0;
      tempden = 0;
      for j = 1:q
        temp = theta(i,j) + theta(i,1:q-j)*theta(i,1+j:q)';
        tempnum = tempnum + temp*j^ktype;
        tempden = tempden + temp;
      end
      num(i,1) = (2*tempnum)^2;
      den(i,1) = (2*tempden + 1 + theta(i,:)*theta(i,:)')^2;
      msg = ['\n MA(' int2str(q) ') '];
    end
    num(i,1) = w(i)*sigma(i)^4*num(i,1);
    den(i,1) = w(i)*sigma(i)^4*den(i,1);
  end
  alpha = sum(num)/sum(den);
end

band = const*(alpha*nobs)^(1/(2*ktype+1));


%========================================================================
%   FINAL CALCULATIONS
%========================================================================
%%in.fmt='%5.4f';in.hspc=2;
%%abs(eig(Phi))
%%mprint1(PHI,in);
% --- Adjust Eigenvalues of PHI to avoid near-singularity ----------------
[B,lambda] = eig(PHI*PHI');
[C,lambda] = eig(PHI'*PHI);
Delta = diag(B'*PHI*C);
if max(abs(Delta)) >.97
  fprintf(gmmopt.prt,'Modified Large Eigenvalues in Andrews-Monahan'); 
  Delta = min(Delta,.97);
  Delta = max(Delta,-.97);
  Delta = diag(Delta);
%%%  PHI = B*Delta*C';
end
%%fprintf(1,'\n');
%%mprint1(PHI,in);

% --- The Matrix D for 'Recoloring' --------------------------------------
D = (In - PHI)\(In + THETA);


% --- Andrews case, No Pre-whitening -------------------------------------
if aminfo.nowhite == 1
  D = eye(size(D));
  uout = u;
  msg = [msg 'bandwidth selection in Andrews HAC\n'];
else
  msg = [msg 'Prewhitening in Andrews-Monahan HAC\n'];
end

% --- Print some Info ---------------------------------------------------
fprintf(gmmopt.prt,msg);
fprintf(gmmopt.prt,' %s Kernel, Automatic bandwidth = %5.2f\n',kname,band);
info.rnames = [repmat(' M',north,1) int2str([1:north]')];
info.fmt = '%5.2f';  info.hspc = 2;

if vardum == 1					% VAR case
  mprint1(PHI,info,gmmopt.prt);
else						% AR/ARMA models
  info.cnames = info.rnames;
  if q > 0 
    out = theta; info.rnames = [repmat('  Theta',q,1) int2str([1:q]')];
  else
    out = []; info.rnames = [];
  end
  if p == 1
    out = [phi out]; info.rnames = strvcat('  Phi1',info.rnames); 
  end
  mprint1(out',info,gmmopt.prt);
end


% --- Diagnostics -------------------------------------------------------
if diagdum == 1
  lags = 12;				% # autocorr to plot
  uterm = input(['Enter index for residual of interest (1:' ...
    int2str(north) '), other exits ']);
  while uterm > 0
    tittext = ['Sample Partial Autocorr: u' int2str(uterm)];
    clf
    spacf(u(:,uterm),lags);		% Sample partial autocorr for u
    title([tittext ' Before Whitening']);
    disp('Strike a key for the next graph'), pause
    spacf(uout(:,uterm),lags);		% same for uout
    title([tittext ' After Whitening']);
    disp('Strike a key to continue'), pause
    uterm = input(['Enter index for residual of interest (1:' ...
      int2str(north) '), other exits ']);
  end
end
