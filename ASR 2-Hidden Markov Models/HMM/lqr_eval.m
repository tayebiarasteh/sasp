function l = lqr_eval(x,means,cholesky_covar);
% function l = lqr_eval(x,means,cholesky_covar);
% evaluates the log Gaussian PDF 
%
% This version of lqr_eval is for transposed input data (DIM-by-N)

%  Dr. Paul M. Baggenstoss
%  Naval Undersea Warfare Center
%  Newport, RI
%  p.m.baggenstoss@ieee.org
%  Date              Reason
%  ------------      --------------
%  Mar 28, 1997      Initial release
%  Dec 7, 1999       Did math a better way

[DIM,N]=size(x);
DIM0=length(means);
if(DIM0 ~= DIM),  
     error('lqr_eval: data and means have different dimension'); 
end;

logdet = 2 * sum(log(abs(diag(cholesky_covar)))) ;
tmpidx = x - repmat(means(:),1,N);
tmpidx = 0.5 * sum(( cholesky_covar' \ tmpidx) .^ 2,1);
l = -tmpidx' - (DIM/2)*log(2*pi) - 0.5*logdet;

