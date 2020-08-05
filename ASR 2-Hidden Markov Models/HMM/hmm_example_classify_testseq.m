% Example of HMM estimation 3 states

% File:   hmm_exmp.m
% Author: Dr. Paul M. Baggenstoss,
%         NUWC, Newport, RI, code 21
%         p.m.baggenstoss@ieee.org
%         401-832-8240
% Date:   Mar 5, 1999
% Subsequent Revisions:



% Classification using the trained HMM parameters
% by calculation of the log likelihood
nsamp2=100;
istart2=1;
[q,parm]=hmm2_reest(parm,x2,istart2,nsamp2,0);
q



