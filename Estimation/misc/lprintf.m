function output = lprintf(prtdat,crtdat,info,file)
%-----------------------------------------------------------------
% PURPOSE:   
% Prints a matrix of data with a criteria-based symbol next
% to each element.  Data are formatted globally or column-wise,
% as specified. Designed for indicating significance levels of
% parameter estimates.
%-----------------------------------------------------------------
% USAGE:  output = lprintf(prtdat,crtdat,info,fid)
%
%  prtdat  (R,C) matrix of data to print
%  crtdat  (R,C) matrix of data to use as criteria for symbols
%  info    structure variable (all optional): [default]
%    info.crit     vector of criteria, largest first  [10%,5%,1%]
%    info.symb     vector of symbols (strings)   [c,b,a]
%    info.fmt      scalar or C-vector of format strings ['%7.4f']
%    info.lnnum    vector indexing rows with horiz. lines [none] 
%    info.rnames    R-vector of strings with row names  [none]
%    info.cnames    C-vector of strings with column names [none]
%    info.head     vector of strings to put as header [none]
%    info.tail     vector of strings to put as footer [none]
%    info.nandum   1 to have NaNs print as blanks [1]
%
% file  'filename' for output. Screen if missing, none if file = 0
%-----------------------------------------------------------------
% RETURNS:  output   a string matrix with the results
%----------------------------------------------------------------
% NOTES: 
%   Symbol(i) assigned when crtdat(r,c) <= info.crit(i) 
%   The number of criteria is flexible.  Default is 10%,5%,1% with
%   corresponding symbols of \superscript c,b,a
%
%   For consistency w/ mprint1, can use rname or rnames and
%   cname or cnames
%-----------------------------------------------------------------

% Written by:  Mike Cliff, Virginia Tech Finance,  mcliff@vt.edu
% CREATED:  12/9/98
% UPDATED: 7/11/03
%          12/28/04	Can use cell structure for row/col names


%==========================================================================
%  INITIALIZATIONS
%==========================================================================

[R,C] = size(prtdat);

pflag = 1;
if nargin == 2, file = 1; info.null = 1; end
if nargin == 3
  if ~isstruct(info)
    file = info; info.null = 1;
  else
    file = 1;
  end
end
if file == 0, pflag = 0; end

if ~isfield(info,'ldum'), info.ldum = 1; end
if info.ldum == 1
  ampstr = ' & ';
else
  ampstr = ' ';
  info.noeol = 1;
end
if ~isfield(info,'crit'), info.crit = [.1;.05;.01]; end
if ~isfield(info,'symb'), info.symb = strvcat('\ ^c','\ ^b','\ ^a'); end
ncrit = rows(info.crit);    % CHECK SORTING OF crit
if info.crit(1) ~= max(info.crit)
  error('Criteria out of order')
end

if ~isfield(info,'fmt'), info.fmt = '%7.4f'; end
if rows(info.fmt) == 1, info.fmt = repmat(info.fmt,C,1); end

if ~isfield(info,'noeol'), info.noeol =  0; end
if ~isfield(info,'lnnum'), info.lnnum =  []; end
if ~isfield(info,'head')
  info.head  = ' '; hdum = 0;
else hdum = 1; end
if ~isfield(info,'tail')
  info.tail  = ' '; tdum = 0;
else, tdum = 1; end
if ~isfield(info,'rname')
  if ~isfield(info,'rnames')
    info.rname = ' '; rdum = 0; 
  else
    info.rname = info.rnames; rdum = 1;
  end
else
  rdum = 1; 
end
if ~isfield(info,'cname')
  if ~isfield(info,'cnames')
    info.cname = ' '; cdum = 0; 
  else
    info.cname = info.cnames; cdum = 1;
  end
else
  cdum = 1; 
end
if iscell(info.cname)
  info.cname = strvcat(info.cname);
end
if iscell(info.rname)
  info.rname = strvcat(info.rname);
end

if ~isfield(info,'nandum'), info.nandum = 1; end

output = [];
amp = repmat(ampstr,R,1);
if info.noeol == 1
  eol = repmat(' ',R,1);
else
  eol = repmat(' \\',R,1);
end

%==========================================================================
%  LOOP TO BUILD UP OUTPUT
%    For each column of data, check the crit matrix row by row,
%    returning the appropriate string
%==========================================================================

for c = 1:C
  symv = repmat(' ',R,cols(info.symb));
  for r = 1:R
    for i = 1:ncrit
      if crtdat(r,c) <= info.crit(i)
        symv(r,:) = info.symb(i,:);
      end
    end
  end
  temp = num2str(prtdat(:,c),info.fmt(c,:));
  output = [output amp temp symv];
end 

if info.nandum == 1
  for r = 1:R
    output(r,:) = strrep(output(r,:),'NaN','   ');
  end
end


% ---- Assemble any header and footer text --------------------------------
head = []; tail = []; rname = []; cname = [];
for i = 1:rows(info.head), head = strvcat(head,info.head(i,:)); end
for i = 1:rows(info.tail), tail = strvcat(tail,info.tail(i,:)); end
if rdum == 1
  for r = 1:R, rname = strvcat(rname,info.rname(r,:)); end
end
if cdum == 1
  if rows(info.cname) == C+1
    cname = info.cname(1,:);
    cstart = 2;
  else
    cstart = 1;
  end
  for c = cstart:C+cstart-1, cname = [cname ampstr info.cname(c,:)]; end
  if info.noeol ~= 1
    cname = [cname ' \\'];
  end
end
if min(info.lnnum) == 0
  cname = [cname ' \hline'];
  info.lnnum(find(info.lnnum==0)) = [];
  if cols(info.lnnum) == 0
    info.lnnum = []
  end
end

% ---- Create hlines to generate horizontal lines -------------------------
hlines = repmat(' ',R,7);
hlines(info.lnnum,:) = repmat(' \hline',rows(info.lnnum),1);


output = [output eol hlines];
if rdum == 1, output = [rname output]; end
if cdum == 1, output = strvcat(cname,output); end
if hdum == 1, output = strvcat(head,output); end
if tdum == 1, output = strvcat(output,tail); end


if pflag == 1
%  if file ~= 1, eval(['diary ',file]); end
%  disp(output)
%  if file ~= 1, eval('diary'); end

  for i = 1:rows(output)
    fprintf(file,'%s\n',output(i,:));
  end
%  if file ~= 1, fclose(file); end
end
