function prt_gmm(results,vnames,fid)
% PURPOSE: Prints output using GMM results structure
%---------------------------------------------------
% USAGE: prt_gmm(results,vnames,fid)
% Where: results = a structure returned by GMM
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_reg(results,[],fid) to print to a file with no vnames   
% --------------------------------------------------
%  RETURNS: nothing, just prints the regression results
% --------------------------------------------------
% SEE ALSO: prt, plt
%---------------------------------------------------   
% VERSION: 1.1.9 (7/11//02)

% written by:
% Mike Cliff, Purdue Finance  mcliff@mgmt.purdue.edu
% CREATED: 12/9/98
% UPDATED: 11/1/99  (1.1.2 Max # minz iterations)
%	   5/4/00   (1.1.3 fixed fid as arg to mprint1)
%          7/12/00  (1.1.4 Fixed Max # minz iterations)
%          8/8/00   (1.1.5 Label assoc w/ W0)
%          9/23/00  (1.1.6 gmmopt.W and fcnchk)  
%          11/30/00 (1.1.7 tidy up sub-optimal W) 
%          7/5/01   (1.1.8 label for optimal W)
%	   7/11/02  (1.1.9 minor label adj: hess, jake, S)
%
% Adapted from PRT_REG in Jim LeSage's Econometrics Toolbox

%===================================================================
%   INITIALIZATIONS
%===================================================================

if ~isstruct(results)
  error('prt_reg requires structure argument');
elseif nargin == 1
  nflag = 0; fid = 1; vsize = 0; vnames = [];
elseif nargin == 2
  fid = 1; nflag = 1; vsize = rows(vnames);
elseif nargin == 3
  nflag = 0;
  [vsize, junk] = size(vnames); % user may supply a blank argument
   if vsize > 0
   nflag = 1;          
   end;
else
 error('Wrong # of arguments to prt_reg');
end;

nobs = results.nobs;
north = results.north;
neq = results.neq;
k = results.nvar;
nz = results.nz;

if isfield(results,'hprt'),  hprt = results.hprt;
else, hprt =0;
end
if isfield(results,'eprt'),  eprt = results.eprt;
else, eprt = 0;
end
if isfield(results,'infoz')
  infoz = results.infoz;
end

%===================================================================
%   MAKE PARAMETER NAMES IF NEEDED
%===================================================================

Vname = '  Parameter';
if vsize == 0, nflag = 0; end
if (rows(vnames) ~= k & nflag == 1)
  warning('Wrong # of parameter names in prt_gmm');
  fprintf(fid,'Will use generic parameter names \n');
  nflag = 0;
end

for i=1:k
  if nflag == 1
    Vname = strvcat(Vname,[blanks(2) vnames(i,:)]);
  else
    tmp = ['  parameter ',num2str(i)];
    Vname = strvcat(Vname,tmp);
  end
end

eqline = ...
' ===============================================================\n';


%===================================================================
%   OUTPUT HEADING    Do if hprt == 1
%===================================================================

if hprt == 1

  if strcmp(infoz.hess,'gn'), dirname='Gauss-Newton';
  elseif strcmp(infoz.hess,'marq')
    dirname=['Marquardt (lambda >= ' num2str(infoz.lambda),')'];
  elseif strcmp(infoz.hess,'dfp') dirname='DFP';
  elseif strcmp(infoz.hess,'bfgs') dirname='BFGS';
  elseif strcmp(infoz.hess,'sd') dirname='Steepest Descent';
  else dirname=['User''s Hessian (' infoz.hess ')'];
  end
  if strcmp(infoz.jake,'numz')   
    dertype='Numerical';
  else   
    dertype=['Analytical (' infoz.jake ')'];
  end
  if strcmp(results.S,'I'), stype='Identity';
  elseif strcmp(results.S,'P'), stype='Plain';  
  elseif strcmp(results.S,'W'), stype='White';
  elseif strcmp(results.S,'H')
    if isfield(results,'wtvec')
      stype=['User-defined weights to ' ...
	     int2str(rows(results.wtvec)) ' lags'];
    else
      stype=['Hansen (',num2str(results.lags),' lags)'];
    end
  elseif strcmp(results.S,'NW') 
    stype=['Newey-West (',num2str(results.lags),' lags)'];
  elseif strcmp(results.S,'G') 
    stype=['Gallant (',num2str(results.lags),' lags)'];
  elseif strcmp(results.S,'AM') 
    if results.aminfo.nowhite == 1
      stype='Andrews (automatic bandwidth)';
    else
      stype='Andrews-Monahan (automatic bandwidth)';
    end
  else
    stype = ['User''s (' results.S ')'];
  end
  if strcmp(results.W0,'Z'), W0type ='inv(Z''Z)';
  elseif strcmp(results.W0,'I'), W0type='I';
  elseif strcmp(results.W0,'Win'), W0type='Fixed';
  elseif strcmp(results.W0,'C'), W0type='Calc from b0';
  else, W0type=results.W0;
  end
  if results.gmmit == 1
    Wtype = ['N/A'];
  else
    if strcmp(results.W,'S')
      Wtype = 'Optimal';    
    else
      Wtype = results.W;     
    end
  end
  fprintf(fid,eqline);
  fprintf(fid,...
  '                      GMM ESTIMATION PROGRAM                     \n');
  fprintf(fid,eqline);
  fprintf(fid,' \n');
  fprintf(fid,' %d Parameters, %d Moment Conditions\n',k,north);
  fprintf(fid,' %d Equation Model, %d Instruments\n',neq,nz);
  fprintf(fid,' %d Observations\n',nobs);
  fprintf(fid,' %d Passes, Max., %d Iterations/Pass\n',...
    results.gmmit,results.infoz.maxit);
  fprintf(fid,' Search Direction:         %s\n',dirname);
  fprintf(fid,' Derivatives:              %s\n',dertype);
  fprintf(fid,' Initial Weighting Matrix: %s\n',W0type);
  fprintf(fid,' Weighting Matrix:         %s\n',Wtype);
  fprintf(fid,' Spectral Density Matrix:  %s\n',stype);
  fprintf(fid,' \n\n');

end


%===================================================================
%   ESTIMATES ETC       Do if eprt == 1
%===================================================================


if eprt == 1

  fprintf(fid,...
  '\n  ------------------  GMM PARAMETER ESTIMATES  -----------------\n');

  tmp = [results.b results.se results.null results.t results.pb];  
  cnames = strvcat('Coeff','Std Err','Null','t-stat','p-val');
  in.cnames = cnames;
  in.rnames = Vname;
  if mean(results.b) > 1000
    bfmt = '%12.2f';
  else
    bfmt = '%12.6f';
  end
  if mean(results.se) > 1000
    sefmt = '%12.2f';
  else
    sefmt = '%12.6f';
  end
  in.fmt = strvcat(bfmt,sefmt,'%6.2f','%6.2f','%6.4f');
  mprint1(tmp,in,fid);


% --- Do This Stuff for Over-Identified Models ------------------------
  if  results.df >= 1 

    fprintf(fid,...
    '  -------------------  GMM MOMENT CONDITIONS  ------------------\n');

    tmp2 = [results.m results.mse results.mt results.mp];  
    in2.cnames = strvcat('Moment','Std Err','t-stat','p-val');
    Vname2 = repmat('     Moment ',results.north,1);
    Vname2 = [Vname2 int2str([1:results.north]')];
    Vname2 = strvcat(' ',Vname2);
    in2.rnames = Vname2;
    if mean(results.m) > 1000
      bfmt = '%12.2f';
    else
      bfmt = '%12.6f';
    end
    if mean(results.mse) > 1000
      sefmt = '%12.2f';
    else
      sefmt = '%12.6f';
    end
    in2.fmt = strvcat(bfmt,sefmt,'%6.2f','%6.4f');
    mprint1(tmp2,in2,fid);
    fprintf(fid,...
      '      J-stat = %5.4f    Prob[Chi-sq.(%d) > J] = %5.4f\n',...
      results.J,results.df,results.p);

% --- Otherwise, Make Sure Moments to Zero ----------------------------
  else
    if results.f > results.ftol         % Mod 9/23/00 to .f (was .J)
      fprintf(fid,'  *** Moments <> 0 in Just-Identified Model.');
      fprintf(fid,'  J-stat = %10.8f ***\n',results.J);
      fprintf(fid,'    Moments = ');
      fprintf(fid,'%6.4f ',results.m);
      fprintf(fid,'\n  Check Convergence tolerances \n');
    end
  end

  fprintf(fid,eqline);

end
