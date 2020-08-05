% Example of HMM estimation 3 states

% File:   hmm_exmp.m
% Author: Dr. Paul M. Baggenstoss,
%         NUWC, Newport, RI, code 21
%         p.m.baggenstoss@ieee.org
%         401-832-8240
% Date:   Mar 5, 1999
% Subsequent Revisions:


% Generate 100 synthetic samples

nsteps = 100;
N = 16;

states = zeros(nsteps,1);
x2 = zeros(2,nsteps);

if testword==1
	A  = [.8 .1 .1; .1 .8 .1; .1 .1 .8];
	Pi = [ 1 0 0];
else
	A  = [.7 .2 .1; .2 .7 .1; .1 .1 .8];
	Pi = [ 1 0 0];
end


%% create state sequence
states(1) = select_discrete_rv(Pi);
for i=2:nsteps,
	states(i,1) = select_discrete_rv(A(states(i-1),:));
end;

%% create feature sequence
for k=1:nsteps
	if states(k,1)==1
		x=randn(N,1);
	elseif states(k,1)==2
		if testword == 1
			x = randn(N,1)*2;
		else
			x = randn(N,1)*2.5;
		end
	else
		x = randn(N,1)+2;
	end;
	z=[mean(x) std(x)];
	x2(:,k)=z';
end;



