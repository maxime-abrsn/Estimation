function diagw(W,iter)
%function diagw(W,iter)
% Plots GMM Weighting matrix
%  W is the matrix
%  iter indicates iteration number

% WRITTEN BY:	Mike Cliff,  Purdue Finance,  mcliff@mgmt.purdue.edu
% DATE:		11/6/99

msg = 'GMM Weighting Matrix';
if nargin == 2
  msg = [msg ' at Iteration ' int2str(iter)];
end 

N = cols(W);

colormap([1 1 1]);			% Turn off colors
bar3(abs(W))				% Make plot
title(msg)				% Add title
for i = 1:N				% Add X's for negative wts
  for j = i:N
    if W(i,j) < 0
      text(i,j,1.01*abs(W(i,j)),'X')
      text(j,i,1.01*abs(W(j,i)),'X')
    end
  end
end

set(gca,'Xtick',[1:N])			% Add labels

