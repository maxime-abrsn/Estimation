function [out1, out2] = prt_desc(results,type,vnames,fid)
% PURPOSE: Prints summary stats from DESC
%---------------------------------------------------
% USAGE: prt_desc(results,type,vnames,fid)
% Where: results = a structure returned by desc
%        type    = 0 suppresses percentiles in output
%        vnames  = an optional vector of variable names
%        fid     = file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%---------------------------------------------------               
%  NOTES:   e.g. vnames = ['y    ',
%                          'x1   ',  NOTE: fixed width
%                          'x2   ',        like all MATLAB
%                          'cterm'];       strings
%           e.g. fid = fopen('ols.out','wr');
%  use prt_desc(results,[],fid) to print to a file with no vnames  
% --------------------------------------------------
%  RETURNS: nothing, just prints the regression results
% --------------------------------------------------
 
% written by: Mike Cliff,  UNC Finance,  mcliff@unc.edu
% CREATED  12/23/98

if ~isstruct(results);
  error('prt_desc requires structure argument');
elseif nargin == 1
  nflag = 0; fid = 1; type = 1;
elseif nargin == 2
  nflag = 0; fid = 1;
elseif nargin == 3
  fid = 1; nflag = 1;
elseif nargin == 4
  nflag = 0;
  [vsize junk] = size(vnames); % user may supply a blank argument
  if vsize > 0
    nflag = 1;          
  end;
else
  error('Wrong # of arguments to prt_desc');
end;
if isempty(type), type = 1; end

k = cols(results);

%  make up some variable names if needed
if nflag == 0
  vnames = [];
  for i=1:k
     vnames = strvcat(vnames,['Var ' int2str(i)]);
  end;
end;

xmean = [];
xstd = [];
xskew = [];
xkurt = [];
xrho = [];
xn = [];
xmin = [];
xmax = [];
xpct = [];

for i = 1:k
   xmean = [xmean; results(i).mean];
   xstd  = [xstd; results(i).std];
   xskew = [xskew; results(i).skew];
   xkurt = [xkurt; results(i).kurt];
   xrho = [xrho; results(i).rho];
   xn    = [xn; results(i).n];
   xmin  = [xmin results(i).min];
   xmax  = [xmax results(i).max];
   xpct  = [xpct results(i).pct];
end

out1 = [xmean xstd xskew xkurt xrho xn];
out2 = [xmin; xpct; xmax];
   
fmt1 = [];    fmt2 = [];
for i = 1:5
   if min(out1(:,i)) > 1000
      newfmt = '%8.1f';
   else
      newfmt = '%8.4f';
   end
   fmt1 = strvcat(fmt1,newfmt); 
end
if max(out1(:,6)) > 1000
  in1.fmt = strvcat(fmt1,'%8d');
else
  in1.fmt = strvcat(fmt1,'%8d');
end
for i = 1:k
   if min(out2(:,i)) > 1000
      newfmt = '%8.1f';
   else
      newfmt = '%8.4f';      
   end
   fmt2 = strvcat(fmt2,newfmt); 
end
in2.fmt = fmt2;

in1.rnames = strvcat('Var',vnames);
in1.cnames = strvcat('Mean','Std Dev','Skew','Ex Kurt','Rho','N');
in2.rnames = strvcat(' ','Min',' 1%',' 5%','10%','25%','50%',...
   '75%','90%','95%','99%','Max');
in2.cnames = vnames;

mprint1(out1,in1,fid);
if type ~= 0
  mprint1(out2,in2,fid);
end
