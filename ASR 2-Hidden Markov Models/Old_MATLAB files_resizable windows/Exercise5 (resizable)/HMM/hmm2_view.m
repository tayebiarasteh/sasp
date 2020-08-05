function hmm2_view(parm,x,i1,i2)
%function hmm2_view(parm,x,i1,i2);

% File:   see_hmm.m
% Author: Dr. Paul M. Baggenstoss,
%         NUWC, Newport, RI, code 21
%         p.m.baggenstoss@ieee.org
%         401-832-8240
% Date:   Mar 5, 1999
% Subsequent Revisions:
%         June 15, 1999 Added 1-D capability
%         Dec 16, 1999 MATLAB5 Upgrade
% modified by Stefan Meier, 2011

if(nargin  < 3),
	error('Syntax: see_hmm(parm,x,i1,[i2]);');
end;

N=parm.N;

for i=1:N,
	subplot(2,N,i+3);
	p = parm.pdf(i);
	gmix_view2(p,x,i1,i2,0,200);

	title(sprintf('State %d',i));

end;
