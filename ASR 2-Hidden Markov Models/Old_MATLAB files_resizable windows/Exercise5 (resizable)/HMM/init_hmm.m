function hparm=init_hmm(data,NSTATES,NMODES,names,min_std);
% function hparm=init_hmm(data,NSTATES,NMODES,names,min_std);

% Initialize HMM parameters
%names={'MEAN','STDV'};
min_std=[.1 .1];
NSTATES=3;
NMODES=3;
hparm.N=NSTATES;
hparm.Pi=ones(NSTATES,1)/NSTATES;
hparm.A=ones(NSTATES,NSTATES)/NSTATES;
for i=1:NSTATES,
    hparm.pdf(i)= init_gmix(data,NMODES,names,min_std);
end;


