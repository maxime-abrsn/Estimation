function Out=mprint1(x,info,fid)
%---------------------------------------------------------------
% PURPOSE: Pretty-prints a matrix
%---------------------------------------------------------------
% USAGE: Out = mprint1(x,info,fid)
%
%  x         data
%  info      structure variable (all optional)       [default]
%   .fmt     format string for matrix                ['%8.4f']
%   .vspc    number of blank lines between blocks    [1]
%   .hspc    spacing between colums                  [3]
%   .ldum    set to 1 for LaTeX version              [0]
%   .rnames  k-vector of row names                   (optional)
%             can incl. col heading if using cnames
%   .cnames  k-vector of column names                (optional)
%   .swidth  width of the display                    [80]
%   .nandum  set to 1 toreplace NaNs with blanks     [1]
%  fid       set to 0 to suppress display results    [1]
%              1 is screen, higher is fid
%---------------------------------------------------------------

% Written by  Mike Cliff, mcliff@vt.edu
% CREATED  12/15/98
% UPDATED  3/17/08 (minor: fixed when fid = 0)
%          12/23/04 (minor: can pass cell strings for rnames/cnames)
%          10/27/99 (minor: 1st elmt of rnames optional)
%          8/3/01  (minor: handle a blank column name)
%	   6/5/02  (minor: added nandum to print blanks)
%	   6/12/02 (minor: nandum tweaking, comma format (see num2strc))
  
%======================================================================
%  INITIALIZATIONS
%======================================================================

if nargin == 3 && fid == 0
  fid = fopen('/temp/junk.out','w');
end

if nargin == 1,  info.null = 1; end
if nargin < 3, fid = 1; end 

[R,C] = size(x);
if ~isfield(info,'vspc'), info.vspc = 0; end
if ~isfield(info,'hspc'), info.hspc = 3; end
if ~isfield(info,'ldum'), info.ldum = 0; end
if ~isfield(info,'swidth'), info.swidth = 80; end
if ~isfield(info,'nandum'), info.nandum = 1; end

if ~isfield(info,'fmt'), info.fmt = '%7.4f'; end
if rows(info.fmt) == 1, info.fmt = repmat(info.fmt,C,1); end


if isfield(info,'cnames')
  cdum = 1; 
  if iscell(info.cnames)
    info.cnames = strvcat(info.cnames);
  end
else 
  cdum = 0; 
end
if isfield(info,'rnames')
  rdum = 1; 
  if iscell(info.rnames)
    info.rnames = strvcat(info.rnames);
  end
else 
  rdum =0; 
  info.rnames = []; 
end

if (cdum == 1 & rdum == 1)
  if rows(info.rnames) == R
    info.rnames = strvcat(' ',info.rnames);
  end
end

ldum = info.ldum;
if ldum == 1
  amp = repmat(' & ',R*(1 + info.vspc)-info.vspc+cdum,1);
  eol = repmat(' \\',R*(1 + info.vspc)-info.vspc+cdum,1);
else
  amp = repmat(' ',R*(1 + info.vspc)-info.vspc+cdum,info.hspc);
  eol = repmat(' ',R*(1 + info.vspc)-info.vspc+cdum,1);  
end

out = repmat(' ',R*(1 + info.vspc)-info.vspc+cdum,1);
Out = out;

%======================================================================
%  BUILD UP OUTPUT
%======================================================================
for c = 1:C  
  if any(isnan(x(:,c))) & info.nandum == 1
    temp = [];
    for r = 1:rows(x)
      if isnan(x(r,c))
	temp = strvcat(temp,' ');
      else
	if ~isempty(findstr(info.fmt(c,:),'C'))
	  temp = strvcat(temp,cfmt(x(r,c),info.fmt(c,:)));
	else
	  temp = strvcat(temp,num2str(x(r,c),info.fmt(c,:)));	
	end
      end
    end
  else
    if ~isempty(findstr(info.fmt(c,:),'C'))
      temp = cfmt(x(:,c),info.fmt(c,:));
    else
      temp = num2str(x(:,c),info.fmt(c,:));
    end    
  end
  if cdum == 1
    ctemp = deblank(info.cnames(c,:));
    if isempty(ctemp), ctemp = ' '; end
    temp = strvcat(ctemp,temp);
  end
  temp = strjust(temp);
  Out = [Out amp temp];  
  if cols([info.rnames out amp temp]) > info.swidth
    if fid > 0
      prtout = [info.rnames out];
      if ldum == 0
        for r = 1:rows(out)
          fprintf(fid,'%s\n',prtout(r,:));
        end
        fprintf(fid,' \n');
      end
      if c == C
        prtout = [info.rnames amp temp];
        if ldum == 0
          for r = 1:rows(out)
            fprintf(fid,'%s\n',prtout(r,:));
          end
          fprintf(fid,' \n');
        end
      end
    end
    out = [amp temp];
  else
    out = [out amp temp];
    if c == C
      if fid > 0
        prtout = [info.rnames out];
        if ldum == 0
          for r = 1:rows(out)
            fprintf(fid,'%s\n',prtout(r,:));
          end
          fprintf(fid,' \n');
        end
      end
    end
  end
end

Out = [info.rnames Out eol];
if ldum == 1 & fid > 0
  for r = 1:rows(Out)
    fprintf(fid,'%s\n',Out(r,:));
  end
end


function cfmtout = cfmt(x,fmt)
  bfmt = fmt(1:findstr(fmt,'%')-1);
  efmt = fmt(findstr(fmt,'C')+1:end);
  n = fmt(findstr(fmt,'.')+1:findstr(fmt,'C')-1);
  cfmtout = num2strc(x,str2num(n),bfmt,efmt);
