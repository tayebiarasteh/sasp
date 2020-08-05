function [p,xp,yp]=gmix_view2(parm,data,index1,index2,do_ellip, M,iplot);
%function [p,xp,yp]=gmix_view2(parm,data,index1,index2,<do_ellip>,<M>,<iplot>);
%
%  parm:           gsum parameters
%  data:           feature data
%  index1, index2: the feature indexes of the two-dimensional view
%  do_ellip:       (optional) if 1,  causes ellipses to be plotted
%  M:              (optional) specify number of grid points.
%  p,xp,yp:        optional output arguments to generate
%                  the PDF on a grid  (No plotting will be
%                  done if here are any output arguments).
%                  Use: imagesc(xp,yp,p). Note that automatic
%                  plotting uses a log intensity, but output p is
%                  probablity.
% iplot (optional) 0: don't plot
%                  1: plot using subplots
%                  2: plot using figures

%  Dr. Paul M. Baggenstoss
%  Naval Undersea Warfare Center
%  Newport, RI
%  p.m.baggenstoss@ieee.org
%  Date              Reason
%  ------------      --------------
%  Mar 28, 1997      Initial release
%  Jan 25, 1999      To un-normalize the axes
%  Feb 7,  1999      Added outputs
%  May 16, 1999      Added other options for iplot
%  Dec,    1999      Converted to MATLAB 5

if nargin < 5,    do_ellip = 0; end;
if nargin < 6,	  M=60;         end;
if nargin < 7,
    if nargout > 0,   iplot=0;	    else, iplot=1; end;
end;

iplot = 0; % modification SM

nmode = length(parm.modes);
xmin = min(data,[],2);
xmax = max(data,[],2);
del=xmax-xmin;
xmax=xmax+del/3;
xmin=xmin-del/3;

[DIM,N]=size(data);


%--- plot points ---
w1 = (xmax(index1)-xmin(index1));
w2 = (xmax(index2)-xmin(index2));
x = [data(index1,:)' data(index2,:)'];
s1=1;
m1=0;
s2=1;
m2=0;

if(iplot),
	if(iplot ==1), subplot(211); else, figure(1); end;
	plot(x(:,1)*s1+m1, x(:,2)*s2+m2,'b.');
	axis([ xmin(index1)*s1+m1 xmax(index1)*s1+m1 xmin(index2)*s2+m2 xmax(index2)*s2+m2]);
	ax=axis;
end;

% create the grid
idx = 1+ [  floor([0:(M*M-1)]' * (1.0/M) )   rem([0:(M*M-1)]',M) ];
idx = [(idx(:,1)*w1/M + xmin(index1)) (idx(:,2)*w2/M + xmin(index2))];
p = zeros(M,M);

if(iplot),
	xlabel(sprintf('%s',parm.features(index1).name));
	ylabel(sprintf('%s',parm.features(index2).name));
	if(iplot ==1), subplot(212); else, figure(2); end;
end;

for i=1:nmode,
    tmpvar = zeros(DIM,DIM);
    tmpvar = parm.modes(i).cholesky_covar;
    [q,tmpvar] = qr(tmpvar(:,[index1 index2]),0);
    tmpvar = triu(tmpvar);
    tmpvar = tmpvar(1:2,:);
    gtmp = zeros(M,M);
    gtmp(:) = exp(lqr_eval(idx(:,[1:2])',  ...
              parm.modes(i).mean([index1 index2]),tmpvar));
    if(do_ellip==1),
		% contour at 1 times sigma
		mx=max(max(gtmp));
		xp=(([1:M])*w1/M+xmin(index1))*s1+m1;
		yp=(([1:M])*w2/M+xmin(index2))*s2+m2;
	    contour(xp,yp, gtmp/mx,[.3679 .3679],'b-');
		hold on;
     	xlabel(sprintf('%s',parm.features(index1).name ));
     	ylabel(sprintf('%s',parm.features(index2).name ));
        axis(ax);
    end;
    p = p + parm.modes(i).weight*gtmp;
end;

if(do_ellip==1), hold off; end;

xp=[xmin(index1)*s1+m1 xmax(index1)*s1+m1];
yp=[xmax(index2)*s2+m2 xmin(index2)*s2+m2];

%if(iplot>0 & do_ellip==0), & modification SM
	p = 10*log10(p);
	mx=max(max(p));
	p = 2*(p - mx + 32);
    idx = find(p < 0);
    p(idx) = zeros(length(idx),1);
	image( xp,fliplr(yp),p);
	%surfl(p);
	%shading interp;
    axis xy
    colormap('hot');
    xlabel(sprintf('%s',parm.features(index1).name));
    ylabel(sprintf('%s',parm.features(index2).name));
%end;

