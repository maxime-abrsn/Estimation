function [a,b,height]=objplot(infoz,Y,X,Z,W,bin,index,r1,r2,labs)
%function [a,b,height]=objplot(infoz,Y,X,Z,W,bin,index,r1,r2,labs)
%
% OBJPLOT       Graphs the quadratic objective function
%               g'Wg as a function of two parameters
%
% INPUT:
% infoz		structure from MINZ/GMM library
% Y		Y data (dep vars)
% X		X data (indep vars)
% Z		Z data (instruments)
% W		GMM Weighting matrix (set to 1 if not desired)
% bin		Parameter estimates 
% index		Pair of integers indexing parms to plot
% r1, r2	Upper and lower bounds on range of deviations from bin
%		  if left blank +/- 10% of bin is used
% labs		Lables for plot if desired
%
% OUTPUT
% a             Values of first paramter 
% b             Values of second parmeter
% height  :     Value of Obj Function at a-b pairs
%===========================================================================
% Example:
%  out = gmm(beta,gmmopt,Y,X,Z);
%  objplot(gmmopt.infoz,Y,X,Z,out.W,out.b,[3 4],...
%    [0 2],[.9 1.1]*out.b(4),myname);
%
%  The structure gmmopt is used in the GMM estimation, along with 
%  data Y, X, and instruments Z.  The results are saved in out.
%  The function is plotted using the same data, and model informatio
%  (gmmopt.infoz) as the estimation.  The weighting matrix from the
%  estimation (out.W) is used to evaluate the objective function.
%  The parameter estimates (out.b) are used to evaluate the function and 
%  is marked as a * on the plot.  The plot varies the 3rd parameter from 0
%  to 2 and the fourth parameter from 90%-110% of its estimated value.
%  The names of these parameters are given in myname.
%
%
%  To do a 1-variable plot, do
%    plot(a,H(:,11)) for the first variable
%    plot(b,H(11,:)) for the second variable  
%
%VERSION: 1.1.3 (5/15/01)

% WRITTEN BY:	Mike Cliff,  Purdue Finance,  mcliff@mgmt.purdue.edu
% CREATED:	 10/19/99
% MODIFIED:	 11/9/99 (1.1.1 raised minimized value on surface)
%                9/23/00 (1.1.2 fcnck)
%                5/15/01 (1.1.3 minor adj for 1-d plots)
  
pts = 21;				% Number of grid points
zoom = .1;				% look at parm +/- zoom %

if nargin <= 9
  labs = strvcat('Var1','Var2');
  if nargin == 7
    r1 = bin(index(1))*(1+[-zoom; zoom]);
    r2 = bin(index(2))*(1+[-zoom; zoom]);
  end
end

% --- Set up Grid Points ------------------------------------------------
nin = rows(bin);	
vec = bin;
b = linspace(r2(1),r2(2),pts)';		% Var on y-axis (2nd NaN)
a = linspace(r1(1),r1(2),pts)';		% Var on x-axis (1st NaN)

na=rows(a);		nb=rows(b);
zv=zeros(na,1);		height=zeros(na,nb);
func = fcnchk('lsfunc');
infoz.call='gmm';

% --- Evaluate Function on Grid ---------------------------------------------
for j=1:nb
  vec(index(2)) = b(j);
  for i=1:na
    vec(index(1)) = a(i);
     zv(i) =  feval(func,vec,infoz,[],Y,X,Z,W);
  end
  height(:,j) = zv(:);
end

htmin = feval(func,bin,infoz,[],Y,X,Z,W);

% --- Plot the Results -----------------------------------------------------
figure(1)
clf
mesh(b,a,height), view(320,30)
xlabel(labs(2,:)), ylabel(labs(1,:)), zlabel('Obj')
axout = axis;
htmin = htmin + .01*(axout(6)-axout(5));
hold on, plot3(bin(index(2)),bin(index(1)),htmin,'r*')

vv = linspace(min(min(height)),max(max(height)),pts)';

figure(2)
clf
[c,h]=contour(b,a,height,vv);
clabel(c,h,'fontsize',6,'color','k','rotation',40)
xlabel(labs(2,:)), ylabel(labs(1,:)),zlabel('Obj')
hold on, plot(bin(index(2)),bin(index(1)),'r*')

