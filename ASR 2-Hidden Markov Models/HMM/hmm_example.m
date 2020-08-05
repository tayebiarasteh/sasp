% Example of HMM estimation (3 states)

echo off
format compact

% State transition and state priors
A=[.8 .1 .1;
   .1 .8 .1;
   .1 .1 .8];

Pi=[ 1 0 0];

% create 'nrecord' records, each with 'nsteps' steps.
nrecord=10;
nsteps=400; % time steps per record
N=16;       % length of each segment
NFEAT=2;

if(~exist('x')),
   [x,istart,nsamp]=hmm_maketestdata(Pi,A,nrecord,nsteps,N,NFEAT);
end;

clf;
plot(x(1,:),x(2,:),'b.');
xlabel('MEAN');
ylabel('STDV');
clc;
echo on
% We've created 10 records of 400 samples each.
% The figure shows the distribution of the data.
% Here are the State transition matrices and
% prior state probabilities:
A
Pi
% Let's see if we can estimate them from the data
% using the Baum-Welch algorithm.  Press a key to start
pause;
clc;

% Initialize HMM parameters
names={'MEAN','STDV'};
min_std=[.1 .1];
NSTATES=3;
NMODE=10;
parm=init_hmm(x,NSTATES,NMODE,names,min_std);

% Press a key to start training
pause;

% Train
NIT=100; % number of iterations
[log_pdf_val, parm] = hmm2_reest(parm, x, istart, nsamp, NIT);


% Observe state PDF's.
%A=parm.A;
%Pi=parm.Pi;
hmm2_view(parm,x,1,2);
% Press a key..
pause
clc

%% Anneal, then do 30 more iterations
%parm=ann_hmm(parm,2,1.2);
%[log_pdf_val, parm] = hmm2_reest(parm, x, istart, nsamp, NIT);
%hmm2_view(parm,x,1,2);
%% That should have fixed the problem..
%pause;
%clc

% Generate 100 synthetic samples
[x2,states]=hmm2_synth_mex(parm,100);
x2=x2';
clf;
plot(x(1,:),x(2,:),'r+',x2(1,:),x2(2,:),'b*');
xlabel('MEAN');
ylabel('STDV');
title('real (red), synthetic (blue)');
pause

% estimate the states of the synthetic data
log_probs = hmm2_get_probs(parm,x2);
est_states=hmm2_viterb(parm,log_probs);
%plot([1:100],states,'ro',[1:100],est_states,'b+');
plot([1:100],states+1,'ro',[0:99],est_states,'b+');


% Classification using the trained HMM parameters
% by calculation of the log likelihood
nsamp2=100;
istart2=1;
[q,parm]=hmm2_reest(parm,x2,istart2,nsamp2,0);
q



