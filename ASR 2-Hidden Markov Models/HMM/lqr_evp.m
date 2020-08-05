function lg=lqr_evp(parm,x,flag)
%function lg=lqr_evp(parm,x,flag)
% Version of lqr_evp for transposed input data (DIM-by-nsamp)
% Flag=1,
%   if flag==1, returns log mode PDF in columns, that is the
%   output of lqr_eval only, not including mixing weights
%   for each column.
% Flag=0,
%   Computes the total log-PDF for normalized data in a single 
%   output column (i.e. the log of the weighted sum of the rows of
%   the flag=1 output) plus the Jacobian (it outputs the 
%   PDF value with respect to the raw un-normalized data by 
%   taking into account the jacobian of the normalization operation).

%  Dr. Paul M. Baggenstoss
%  Naval Undersea Warfare Center
%  Newport, RI
%  p.m.baggenstoss@ieee.org
%  Date              Reason
%  ------------      --------------
%  Mar 28, 1997      Initial release
%  Jul 19, 1997      Accounted for data scaling
%  Oct 20, 1998      Fixed bug occuring when nmode=1
%  Dec 6, 1999       MATLAB 5 Upgrade

nmode = length(parm.modes);
wts = [parm.modes.weight]';
DIM0 = length(parm.features);
[DIM,N] = size(x);

if( DIM ~= DIM0) ,
    error('LQREVAL_P: data and parms have different DIM');
end;

lg = zeros(N,nmode);
for k=1:nmode,
       lg(:,k) = lqr_eval(x,parm.modes(k).mean,parm.modes(k).cholesky_covar);
end;

if(flag==1),	return;   end;

if(nmode > 1),
    % find the max across columns
    mx = max(lg,[],2);

    % subtract that from each column
    for i=1:nmode,   lg(:,i)=lg(:,i)-mx;   end;

    % compute the exponent
    g = exp(lg);

    % Compute combined Gaussian sum function at each sample.
    y = g * wts;

    % add back in the scaling
    lg = log(y) + mx;
end;
