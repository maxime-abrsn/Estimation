function [gmmout, gmmopt]=gmm(b,gmmopt,Y,X,Z,Win)
% PURPOSE: Estimate model parameters using GMM
%--------------------------------------------------------------------------
% USAGE: [gmmout, gmmopt]=gmm(b,gmmopt,Y,X,Z,Win)
% 
% b  Vector of starting values for parameters
% gmmopt  structure of gmm options                              [default]
%   gmmopt.infoz   Nested structure of infoz needed in MINZ
%   gmmopt.infoz.momt  Filename of moment conditions            REQUIRED
%   gmmopt.infoz.jake  Filename of Jacobian of moment cond      ['numz']
%   gmmopt.infoz.hess  Hessian updating (see gmmS function)     ['gn']
%
%   gmmopt.gmmit  Number of GMM iterations (NaN is iterated)	[2]
%   gmmopt.maxit  Cap on number of GMM iterations		[25]
%   gmmopt.tol    Convergence criteria for iter. GMM		[1e-7]
%   gmmopt.W0     Initial GMM weighting matrix			['Z']
%	'I' = Identity, 'Z' = Instruments (Z'Z), 'C' = Calculate from b, 
%       'Win' = Fixed passed as Win, myfile = user's own m-file
%   gmmopt.W      Subsequent GMM Weighting Matrix               ['S']
%       'S' = inverse Spectral Density from gmmS
%       myfile = user's m-file name
%   gmmopt.S      Type of Spectral Density matrix		['NW']
%	'W'=White, 'NW'=Newey-West (Bartlett), 'G'=Gallant (Parzen)
%	'H'=Hansen (Truncated), 'AM'=Andrews-Monahan, 'P'=Plain (OLS)
%	myfile = user's m-file
%   gmmopt.aminfo structure if gmmopt.S='AM'.   See ANDMON.M
%   gmmopt.lags   Lags used in truncated kernel for S	        [nobs^(1/3)]
%   gmmopt.wtvec  User-provided vector of wts for truncated kernel
%   gmmopt.Strim  Contols demeaning of moments in calc of S	[1]
%	0 = none, 1 = demean e, 2 = demean Z'e
%   gmmopt.Slast  1 to recalc S at final param est. 2 updates W	[1]
%   gmmopt.null   Vector of null hypotheses for t-stats		[0]
%   gmmopt.prt    Fid for printing (0=none,1=screen,else file)	[1]
%   gmmopt.plot   1 does some plots, else suppress		[1]
%   gmmopt.vname  Optional k-vector of parameter names
%
% Y  "Dependent" variables
% X  "Independent" variables
% Z  Instruments (can be same as X)
% Win User-defined Initial weighting matrix                     OPTIONAL
%    To use a function to calculate W0, don't use Win, but set
%    gmmopt.W0 to 'U' and give the m-file name in gmmopt.W
%
% See Also:
%	gmmS.m     more info on the spectral density matrix
%	hessz.m    more info on Hessian methods
%
%--------------------------------------------------------------------------
% RETURNS: 
% gmmout  results structure
%  gmmout.f     function value
%  gmmout.J     chi-square stat for model fit
%  gmmout.p     p-value for model fit
%  gmmout.b     coefficient estimates
%  gmmout.se    standard errors for each parameter
%  gmmout.bcov  cov matrix of parameter estimates
%  gmmout.t     t-stats for parms = null
%  gmmout.pb    p-values for coefficients
%  gmmout.m     moments
%  gmmout.mse   standard errors of moments
%  gmmout.varm  covariance matrix of moments
%  gmmout.mt    t-stats for moments = 0
%  gmmout.mp    p-vals for moments
%  gmmout.nobs  number of observations
%  gmmout.north number of orthogonality conditions
%  gmmout.neq   number of equations
%  gmmout.nz    number of instruments
%  gmmout.nvar  number of parameters
%  gmmout.df    degrees of freedom for model
%  gmmout.stat  stat structure from MINZ
%  gmmout.null  vector of null hypotheses for parameter values
%  gmmout.W     weighting matrix
%  gmmout.S     spectral density matrix
%  gmmout.eflag error flag for spectral density matrix
%  gmmout.ithist History of MINZ iterations
%
% gmmopt  updated options structure (also includes nobs, neq, etc.)
%--------------------------------------------------------------------------
% VERSION: 1.3.5 (12/27/04)

% written by:
% Mike Cliff,  Purdue Finance  mcliff@mgmt.purdue.edu
% CREATED:   12/10/98
% MODIFIED:  11/14/99 (1.2.x modified calc of f with user's W; pass M)
%            7/21/00  (1.2.1 |dof| in iterated GMM, NW lags if error, GN dflt) 
%            8/8/00   (1.2.2 iter after using user's W0; print momt wts)
%            9/23/00  (1.3.x Win can be fixed matrix or m-file)
%            11/30/00 (1.3.1) improved usage of sub-optimal W
%            5/11/01  (1.3.2) Fixed Wuse to properly do suboptimal W
%            7/11/02  (1.3.3) Added Plain S, mod user's own S
%            4/3/03  (1.3.4) Don't call print functions unless needed
%            12/27/04  (1.3.5) Set MINZ default printing to match GMM
  
%===================================================================
%   INITIALIZATIONS
%===================================================================

if ~isstruct(gmmopt)
  error('GMM options should be in a structure variable');
end;

nobs = rows(Y);
nz = cols(Z);
k = rows(b);

% --- Basic Iterations, etc -----------------------------------------
if ~isfield(gmmopt,'gmmit'), 	gmmopt.gmmit = 2;		end
if ~isfield(gmmopt,'maxit'),	gmmopt.maxit = 25;		end
if ~isfield(gmmopt,'tol'),	gmmopt.tol   = 1e-7;		end
if ~isfield(gmmopt,'lags'),	gmmopt.lags = floor(nobs^(1/3));end
if ~isfield(gmmopt,'null'),	gmmopt.null  = zeros(k,1);	end

% --- File references and Optimization Options ----------------------
if ~isfield(gmmopt.infoz,'momt'),
  error('I Need Some Moment Conditions')
else
  momt = fcnchk(gmmopt.infoz.momt);
  m = feval(momt,b,gmmopt.infoz,[],Y,X,Z);
  north = rows(m);
  if north < k
    error(sprintf('Model is not Identified. %1d moments, %1d parms',north,k)); 
  end
end
if ~isfield(gmmopt.infoz,'jake'),  gmmopt.infoz.jake='numz';	end
if ~isfield(gmmopt.infoz,'hess'),  gmmopt.infoz.hess = 'gn';	end
if ~isfield(gmmopt.infoz,'maxit'), gmmopt.infoz.maxit = 100;	end
if (~isfield(gmmopt.infoz,'lambda') & strcmp(gmmopt.infoz.hess,'marq'))
  gmmopt.infoz.lambda = .01;  
end

if ~isfield(gmmopt,'W'),	gmmopt.W = 'S'; 	end
if ~isfield(gmmopt,'Slast'),	gmmopt.Slast=1;		end
if ~isfield(gmmopt,'W0'),	gmmopt.W0 = 'Z';	end  
if ~isfield(gmmopt,'S'),	gmmopt.S = 'NW';	end
if ~isfield(gmmopt,'Strim'),	gmmopt.Strim = 1;	end
if strcmp(gmmopt.S,'AM')
  if ~isfield(gmmopt,'aminfo'), gmmopt.aminfo.p = 1; end
  if ~isfield(gmmopt.aminfo,'kernel'), gmmopt.aminfo.kernel = 'QS'; end
  if ~isfield(gmmopt.aminfo,'p'), gmmopt.aminfo.p = 1; end
  if ~isfield(gmmopt.aminfo,'q'), gmmopt.aminfo.q = 0; end
  if ~isfield(gmmopt.aminfo,'vardum'), gmmopt.aminfo.vardum = 0; end
  if ~isfield(gmmopt.aminfo,'nowhite'), gmmopt.aminfo.nowhite = 0; end
end

if strcmp(gmmopt.W0,'Win')  % Mod 8/8/00 to allow iter after initial user's W
  if nargin ~= 6, error('You Specified a Fixed W but didn''t Pass One.'); end
%  gmmopt.gmmit = 1;
end

% --- Some Printing Controls ---------------------------------------
if ~isfield(gmmopt,'prt')
  gmmopt.prt = 1;
else
  if ~isfield(gmmopt.infoz,'prt')
    gmmopt.infoz.prt = gmmopt.prt;
  end
end
if ~isfield(gmmopt,'plot'),	gmmopt.plot=1;		end
if ~isfield(gmmopt,'vname'),	gmmopt.vname=[];	end

gmmopt.infoz.call = 'gmm';
gmmopt.infoz.func = 'lsfunc';	gmmopt.infoz.grad = 'lsgrad';
gmmopt.nobs = nobs;		gmmopt.north = north;
gmmopt.neq = north/nz;		gmmopt.nvar = k;
gmmopt.nz = nz;
gmmout.eflag = 0;		% Set to zero for error checking

%===================================================================
%   DISPLAY SOME HEADING INFOZ
%===================================================================
gmmopt.hprt = 0; gmmopt.eprt = 0;
if gmmopt.prt ~= 0
  gmmopt.hprt = 1; 
  prt_gmm(gmmopt,gmmopt.vname,gmmopt.prt)
end


%===================================================================
%   LOOP FOR ITGMM
%===================================================================

% --- Get Ready for Loop -------------------------------------------
if isnan(gmmopt.gmmit)
  maxiter = gmmopt.maxit;
  loopdum = 1;
else
  maxiter = gmmopt.gmmit;
  loopdum = 0;
end

bold = b;
of = 1/eps;			% Initialize objective function
Wuse = gmmopt.W0;		% Initialize Weighting Matrix choice
ithist = [];			% Store history of MINZ iterations
i = 1;

% --- Do the Loop --------------------------------------------------
while i <= maxiter
  if i == 1                     % Mod 5/11/01 to fix sub-optimal W
    Wuse = gmmopt.W0;
  else
    Wuse = gmmopt.W;
  end  
  gmmopt.i = i;
  fprintf(gmmopt.prt,[blanks(20) 'STARTING GMM ITERATION %2d\n'],i);

% --- Determine the Weighting Matrix -------------------------------
  if strmatch(Wuse,strvcat('I','Z','S','C'),'exact')
    if strmatch(gmmopt.S,strvcat('P','W','NW','G','H','AM'),'exact')
      [S,eflag,gmmopt] = gmmS(b,gmmopt,Y,X,Z);
      gmmout.eflag = max(gmmout.eflag,eflag);
      if eflag == 1
	gmmopt.S = 'NW';
	gmmopt.lags = min(gmmopt.lags,floor(nobs^(1/3)));
	S = gmmS(b,gmmopt,Y,X,Z);
	fprintf(gmmopt.prt,'Switching to Newey-West (%d lags)\n',...
		gmmopt.lags);
      end
    else
      S = feval(fcnchk(gmmopt.S),b,gmmopt.infoz,Y,X,Z);
    end
    W = S\eye(north);
  else
% --- Next little bit is new for suboptimal W --------------------------
    if strcmp(Wuse,'Win')
      W = Win;
    else
      W = feval(fcnchk(Wuse),b,gmmopt,Y,X,Z);
    end
  end
  if gmmopt.plot == 1
    figure(100+i)
    diagw(W,i);
  end
  
% --- Print Weights: M'Wm = w'm = 0 (Added 8/8/00) -------------------
% Now only print if fairly small.  Increase limit of 10 if desired
  if max(north,k) <= 10 & gmmopt.prt > 0
    fprintf(gmmopt.prt,'  Weights Attached to Moments\n');
    jake = fcnchk(gmmopt.infoz.jake);
    M = feval(jake,b,gmmopt.infoz,[],Y,X,Z);
    wt = M'*W./repmat(sum(M'*W,2),1,north);
    wtinfo.cnames = repmat('Moment ',north,1);  
    wtinfo.cnames = [wtinfo.cnames int2str([1:north]')];
    if isempty(gmmopt.vname)
      wtinfo.rnames = strvcat(' ',[repmat('Var',k,1) int2str([1:k]')]);
    else
      wtinfo.rnames = strvcat(' ',gmmopt.vname); 
    end  
    mprint1(wt,wtinfo,gmmopt.prt);
  end  
    
% --- Do the Minimization ------------------------------------------
  bold = b;
  [b,gmmopt.infoz,stat] = minz(b,gmmopt.infoz.func,gmmopt.infoz,Y,X,Z,W);
  dof = (of-stat.f)/of;
  of = stat.f;
  if (abs(dof) <= gmmopt.tol & loopdum == 1)
    i = maxiter+1;
  end

% --- Print Est from this iteration --------------------------------
  if i < maxiter & gmmopt.prt > 0
    bupdate.rnames = strvcat(' ',[' b' int2str(i)]);
    if isempty(gmmopt.vname)
      bupdate.cnames = [repmat('Var',k,1) int2str([1:k]')];
    else
      bupdate.cnames = gmmopt.vname; 
    end
    mprint1(b',bupdate,gmmopt.prt);  
  end

% --- Some Housekeeping --------------------------------------------
  ithist = [ithist; stat.iter];
  i = i + 1;
end



%===================================================================
%   NOW SOME CALCULATIONS.  MAKE STD ERRs, CREATE OUTPUT STUFF
%===================================================================

% --- Re-evaluate moments and Jacobian -----------------------------
gmmopt.infoz.Seval = 1;                  % Can use in momt cond
momt = fcnchk(gmmopt.infoz.momt);
jake = fcnchk(gmmopt.infoz.jake);
m = feval(momt,b,gmmopt.infoz,stat,Y,X,Z);
M = feval(jake,b,gmmopt.infoz,stat,Y,X,Z);

% --- Re-evaluate Spectral Density if Needed -----------------------
if gmmopt.Slast >= 1
  fprintf(gmmopt.prt,...
    '\n            EVALUATING S at FINAL PARAMETER ESTIMATES\n');
  gmmopt.i = i;                        % Want SEs, not initial W
  if strmatch(gmmopt.S,strvcat('P','W','NW','G','H','AM'),'exact')
    S = gmmS(b,gmmopt,Y,X,Z);
  else
    S = feval(fcnchk(gmmopt.S),b,gmmopt,Y,X,Z);
  end
  gmmopt.i = i - 1;                    % Restore actual iteration #
  if gmmopt.Slast == 2
    W = S\eye(north);
  end
end

% --- Calculate Covariance matrices --------------------------------
term = (M'*W*M)\eye(k);
bcov = term*(M'*W*S*W*M)*term/nobs;	% Cov(b)
term = (eye(north) - M*term*M'*W);
varm = term*S*term'/nobs;		% Cov(m)

% --- The J-stat ---------------------------------------------------
if strcmp(Wuse,'S')
  gmmout.J = nobs*stat.f;		% T x Obj Function from min
else
  gmmout.J = nobs*m'*pinv(nobs*varm)*m;	% Sub-optimal W, Cov(m) singular
end

% --- Assign Results to output structure --------------------------
gmmout.m = m;                   gmmout.M = M;
if north > k
  gmmout.mse = diag(sqrt(varm));
  gmmout.mt = gmmout.m./gmmout.mse;
  gmmout.mp = stdn_prb(gmmout.mt);
end
gmmout.nobs = nobs;		gmmout.north = north;
gmmout.neq = north/nz;		gmmout.nvar = k;
gmmout.nz = nz;
gmmout.null = gmmopt.null;	gmmout.df = north - k;
gmmout.b = b;
gmmout.bcov = bcov;		gmmout.varm = varm;
gmmout.se = sqrt(diag(bcov));	gmmout.t = (b-gmmopt.null)./gmmout.se;
gmmout.f = stat.f;              gmmout.p = 1 - chi2cdf(gmmout.J,north-k);
gmmout.pb = stdn_prb(gmmout.t); 
gmmout.W = W;			gmmout.S = S;
gmmout.stat = stat;		gmmout.ithist = ithist;


%===================================================================
%   PRINT THE OUTPUT
%===================================================================

if isfield(gmmopt.infoz,'ftol')			% This checks if
  gmmout.ftol = gmmopt.infoz.ftol;		% moments = 0 for
else						% Just-identified model
  gmmout.ftol = 1e-7; 
end
gmmout.hprt = 0; gmmout.eprt = 0;
if gmmopt.prt ~= 0
  gmmout.eprt = 1; 
  prt_gmm(gmmout,gmmopt.vname,gmmopt.prt)
end
gmmout = rmfield(gmmout,'ftol');
