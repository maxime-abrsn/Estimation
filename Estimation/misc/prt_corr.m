function out = prt_corr(in,info)
% function out = prt_corr(in,info)
%
%  PRT_CORR	Prints a correlation matrix.  Drops 1st row and last col
%		and leaves blanks on and above diagonal
%  info contains
%   .names      Names of variables 
%   .fmt	Format of output ['%6.2f']
%   .ldum	Set to 1 for Latex [1]

% WRITTEN BY: Mike Cliff,  Virginia Tech Finance,  mcliff@vt.edu
% Updated 8/31/99
%         12/29/04

if nargin == 1
  info.null = [];
end

if ~isfield(info,'fmt')
  info.fmt = '%6.2f';
end

if ~isfield(info,'names')
  info.names = [repmat('Var',rows(in),1) int2str([1:rows(in)]')];
end

if iscell(info.names)
  info.names = strvcat(info.names);
end

if ~isfield(info,'ldum')
  info.ldum = 1;
end

data = in(2:end,1:end-1);
if info.ldum == 1
  spc = ' & ';
  eol = ' \\ ';
  amp = repmat(spc,rows(data),1);
  EOL = repmat(eol,rows(data),1);
  colhead1 = '\multicolumn{1}{c}{';
  colhead2 = '}';
else
  spc = ' ';
  eol = ' ';
  amp = repmat(spc,rows(data),1);
  EOL = repmat(eol,rows(data),1);
  colhead1 = ' ';
  colhead2 = ' ';
end

tempcols = [];
out = [];
for c = 1:cols(data)
  tempcols = [tempcols spc colhead1 deblank(info.names(c,:)) colhead2];
  temp = [];
  for r = 1:rows(data)
    if r < c
      new = ' ';
    else
      new = num2str(data(r,c),info.fmt);
    end
    temp = strvcat(temp,new);
  end

  out = [out amp strjust(temp)];
end

out = [info.names(2:end,:) out EOL];
out = strvcat([tempcols eol],out);
