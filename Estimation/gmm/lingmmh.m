function stat = lingmmh(b,infoz,stat,y,x,z,w)
% PURPOSE:  Provide analytic Hessian for linear GMM estimation
%-------------------------------------------------------------------------
% USAGE:  [stat,H] = lingmmh(b,infoz,stat,y,x,z,w)
%  b      model parameters
%  infoz   MINZ infoz structure
%  stat   MINZ status structure
%  y,x,z  Data: dependent, independent, and instruments
%  w      GMM weighting matrix
%-------------------------------------------------------------------------
% RETURNS:   INVERSE Hessian!!!
%  stat  updated status structure variable, including new inverse Hessian
%-------------------------------------------------------------------------
%VERSION: 1.2.1 (11/13/00)

% written by:
% Mike Cliff,  Purdue Finance  mcliff@mgmt.purdue.edu
% Created: 12/10/98
% Updated: 7/24/00  (Scaling of H by 2x)
%          11/13/00 (1.2.1 W not used as arg to lingmm anymore)
  
M = lingmmj(b,infoz,stat,y,x,z);
H = 2*M'*w*M;
stat.Hi = (H\eye(size(H)));
stat.Hcond=sqrt(cond(stat.Hi));
