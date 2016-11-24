function [S,eflag,gmmopt] = gmmS(b,gmmopt,Y,X,Z)
% PURPOSE: Calculates the spectral density matrix S in GMM
%--------------------------------------------------------------------
% USAGE: S=gmmS(b,gmmopt,Y,X,Z)
% where
%  b       parameter values
%  gmmopt  options structure from GMM
%    .W      Set to 'S' to use gmmS for inverse Spectral Density
%		as GMM weighting matrix
%    .W0     Initial weighting matrix:
%		'Z'	inv(I kron Z'Z) DEFAULT
%		'I'	Identity
%		'C'	Calculated from parms in b
%               'U'     User-defined (pass to gmm as Win)
%
%    .S      Type of spectral density matrix:
%		'I'	Identity
%		'P'	Plain (like OLS,  e'e.*.Z'Z)
%		'W'	White (hetero-consistent)
%		'H'	Hansen (truncated, wts = 1 or 0)
%		'NW'	Newey-West (Bartlett weights)
%		'G'	Gallant (Parzen weights)
%		'AM'	Andews or Andrews-Monahan (see ANDMON)
%
%    .lags   Lags in H, NW, or G kernels
%    .wtvec  User-defined weights for Hansen matrix
%		Allows for seasonals, etc (eg. wtvec = [1 0 1])
%    .Strim  Controls demeaning of moments in calculating S
%		0	None
%		1	Demean model errors e DEFAULT
%		2	Demean moments, Z'e
%
%  Y,X,Z   data: dependent, independent, and instruments
%--------------------------------------------------------------------
% RETURNS: S  the spectral density matrix of the GMM moments
%          eflag = 1 if matrix is neg. definite
%
%--------------------------------------------------------------------
% NOTES:
%    gmmopt.Seval  Lets momt know it is eval S   (used for special 
%              application where added a penalty in momt to impose bound,
%              don't want penalty in calc of S)  
%--------------------------------------------------------------------
% VERSION: 1.1.7 (2/13/03)

% written by:
% Mike Cliff,  Purdue Finance  mcliff@mgmt.purdue.edu
% CREATED:  12/10/98 (1.1.1 converted from gmmW)
% UPDATED:  11/19/99 (1.1.2 inv(Z'Z); Seval)
%           5/11/00  (1.1.3 Seval) 
%           8/8/00   (1.1.4 error checking for User-defined W)
%           9/23/00  (1.1.5 gmmopt.W and fcnchk)
%           7/11/02  (1.1.6) Added Plain option
%           2/13/03  (1.1.7) Fixed if-stmt for plotting
  
%====================================================================
%   INITIALIZATIONS
%====================================================================
gmmopt.infoz.Seval = 1;              % Let momt know it is calc S
eflag = 0;

momt = fcnchk(gmmopt.infoz.momt);
jake = fcnchk(gmmopt.infoz.jake);
if nargout(momt) == 1
  error(['Function ' gmmopt.infoz.momt ' requires two outputs.'])
end  

% --- Get info on size of problem ---------------------------------------
[m,e] = feval(momt,b,gmmopt.infoz,[],Y,X,Z);
M = feval(jake,b,gmmopt.infoz,[],Y,X,Z); 
nobs = rows(e);
nz = cols(Z);
north = rows(m);
neq = north/nz;
Stype = gmmopt.S;
gmmopt.infoz.Seval = 0;

% --- Give only Identity if specified ------------------------------
if strcmp(Stype,'I')
  S = eye(north);
  return
end

% --- Determine number of lags to use ----------------------------------
switch Stype
  case {'W' 'P'}
    maxlags = 0; amdum = 0; band = 1;
  case {'NW','H','G'}
    maxlags = gmmopt.lags; amdum = 0; band = maxlags+1;
  case 'AM'
    maxlags = nobs-1; amdum = 1; Stype = gmmopt.aminfo.kernel;
  case {'QS','TH'}
    error('Set gmmopt.S=''AM'' and gmmopt.aminfo.kernel=''QS'' or ''TH''');
  otherwise
    error('Incorrect Spectral Matrix choice.  Check gmmopt.S')
end

% --- For Truncated Kernel with User-defined Wts -----------------------
if (strcmp(Stype,'H') & isfield(gmmopt,'wtvec'))
  wtvec = [1; gmmopt.wtvec];			% Wt of 1 on lag 0
  maxlags = rows(wtvec)-1;
end

%==========================================================================
%   CALCULATE INITIAL WEIGHTING MATRIX
%==========================================================================
%%if (gmmopt.i == 1 & strcmp(gmmopt.W,'S'))  % If 1st GMM iter & 'optimal' W
if gmmopt.i == 1                        % If 1st GMM iter & 'optimal' W
  if strcmp(gmmopt.W0,'I')
    S = eye(north);			% Identity 			
    return
  elseif strcmp(gmmopt.W0,'Z')		% Instruments
    S = kron(eye(neq),(Z'*Z));
    return
  elseif strcmp(gmmopt.W0,'C')		% Calculate from Starting Values
    S=zeros(north,north);
  elseif strcmp(gmmopt.W0,'U')
    error('GMMS: Shouldn''t be here to evaluate W0.')
  end
else
  S=zeros(north,north);			% 2nd or more GMM iter
end

if gmmopt.north ~= north
  error('Mismatch of matrix dimensions')
end

if gmmopt.Strim == 1
  e = e - repmat(mean(e),nobs,1);	% Remove sample mean
end


%===========================================================================
%   ANDREWS-MONAHAN PROCEDURE
%===========================================================================
if amdum == 1
  [band,u,D,amerr] = andmon(gmmopt,e,M,Z);
  if amerr == 1
    fprintf(gmmopt.prt,'Trouble with ARMA, Switching to AR(1)\n');
    gmmopt.aminfo.p=1; gmmopt.aminfo.q=0; gmmopt.aminfo.vardum=0;
    [band,u,D,amerr] = andmon(gmmopt,e,M,Z);
    disp('ANDMON Error')
  end
  maxlags = nobs-1;
else
  u = [];
  for i = 1:cols(e)
    u = [u repmat(e(:,i),1,nz).*Z];
  end
end

if gmmopt.Strim == 2
  u = u - repmat(mean(u),nobs,1);	% Remove sample mean
end

%corrcoef(u)				% Uncomment to help debug
%save uout u				% Finds linear dependencies


%====================================================================
%   LOOP THROUGH TO BUILD UP S   Follows Greene (1997, p.528)
%====================================================================

for lag = 0:maxlags
  
  Rho = u(1:nobs-lag,:)'*u(1+lag:nobs,:);	% Vectorized (fast!!)

%  Rho = zeros(north,north);			% Loop (~25 x slower)
%  for i = lag+1:nobs				% Old way
%    e1 = e(i,:);  e2 = e(i-lag,:);
%    Z1 = Z(i,:);  Z2 = Z(i-lag,:);
%    Rho  = Rho + kron(e1'*e2,Z1'*Z2);
%  end

  Rho = Rho/nobs;
  if lag >= 1, Rho = Rho + Rho'; end

%in.fmt = '%5.2f'; mprint1(Rho,in,1); pause	% Printing for debugging

  x = lag/band;
  switch Stype
   case 'P'		% --- Plain Standard Errors --------------------
    wt = 1;
    Rho = kron(e'*e/nobs,Z'*Z/nobs);
   case {'I','W','H'}	% --- Identity, White, or Hansen ---------------
    plotname = 'Truncated';
    if isfield(gmmopt,'wtvec')
      wt = wtvec(lag+1);
    else
      wt = 1;
    end
  case 'NW'		% --- Newey West (Bartlett) --------------------
    plotname = 'Bartlett (Newey-West)';
    wt = 1 - x;
    if x > 1, wt = 0; end
  case 'G'		% --- Gallant (Parzen) -------------------------
    plotname = 'Parzen (Gallant)';
    if x < 0.5  
      wt = 1 - 6*x^2 + 6*x^3;
    elseif x < 1
      wt = 2*(1-x)^3; 
    else 
      wt = 0;
    end
  case 'QS'		% --- Quadratic Spectral -----------------------
    plotname = 'Quadratic-Spectral';
    term = 6*pi*x/5;
    if lag == 0
      wt = 1;
    else
      wt = 25*(sin(term)/term - cos(term))/(12*pi^2*x^2);
    end
  case 'TH'		% --- Tukey-Hanning ----------------------------
    plotname = 'Tukey-Hanning';
    if x < 1
      wt = (1+cos(pi*x))/2;
    else
      wt = 0;
    end
  end
  wt = wt*nobs/(nobs-rows(b));			% Degrees of freedom adj
  S = S + wt*Rho;
  wtout(lag+1,1) = wt;
end



%====================================================================
%   FINISHING STUFF
%====================================================================


% --- Plot Kernel Weights if Desired --------------------------------
if (gmmopt.plot == 1 & ~(strcmp(Stype,'W') | strcmp(Stype,'P') ))
  maxplot = min(maxlags,50);
  xdata = [0:maxplot]';
  figure(200+gmmopt.i)
  plot(xdata,wtout(1:maxplot+1),'-',xdata,0*xdata,'r:')
  title([plotname ' Kernel Weights in HAC Estimate'])
  xlabel('Lag')
  ylabel('Weight')
  axis([0 maxplot -.1 1.1])
end

% --- Recolor S if Using Andrews-Monahan Procedure ------------------
if amdum == 1
  S = D*S*D';
  gmmopt.lags = band;			% Send bandwidth for printing
end

% --- Check that S is Positive Definite -----------------------------
try 					% Safeguard against errors
  mineig = min(eig(S));
catch					% If problem with eig(S)
  mineig = -Inf;
end
if mineig < max(size(S))*norm(S)*eps	% If S no P.D. give info
  save sout S
  eflag=1;
  fprintf(gmmopt.prt,...
    'Spectral Density Matrix not Pos. Def.  Min(eig(S)) = %10.8f\n',...
    mineig);
end

