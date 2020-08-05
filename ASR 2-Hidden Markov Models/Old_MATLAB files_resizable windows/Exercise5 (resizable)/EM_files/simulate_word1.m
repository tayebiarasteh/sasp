function [data] = simulate_word1
%SIMULATE_WORD Summary of this function goes here
%   Detailed explanation goes here

% State transition and state priors
A=[.7 .2 .1;
   .2 .7 .1;
   .1 .1 .8];

Pi=[ 1 0 0];

% create 'nrecord' records, each with 'nsteps' steps.
nrecord=10;
nsteps=400; % time steps per record
N=16;       % length of each segment
NFEAT=2;

[data,istart,nsamp]=hmm_maketestdata1(Pi,A,nrecord,nsteps,N,NFEAT);

end

