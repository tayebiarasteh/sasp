function [ft,istart,nsamp]=hmm_maketestdata(Pi,A,nrecord,nsteps,N,NFEAT);
%function [ft,istart,nsamp]=hmm_maketestdata(Pi,A,nrecord,nsteps,N,NFEAT);

ft=zeros(NFEAT,nsteps*nrecord);

fprintf('Creating data...\n');
% open the feature file

istart=[0:nrecord-1]'*nsteps+1;
nsamp = ones(nrecord,1)*nsteps;
	
for irec=1:nrecord,
		fprintf('Creating record %d of %d\n',irec,nrecord);
	
		% first create the state sequence
		states=zeros(nsteps,1);
		states(1) = select_discrete_rv(Pi);
		for i=2:nsteps,
			states(i) = select_discrete_rv(A(states(i-1),:));
		end;
	
		% now create the observations
		for i=1:nsteps,
			x=randn(N,1);
			if(states(i)==1),
			     x=randn(N,1);
			elseif(states(i)==2),
				x=randn(N,1)*2;
			else,
				x=randn(N,1)+2;
			end;
			z=[mean(x) std(x)];
                        ft(:,i+(irec-1)*nsteps)=z(:);
		end;
end;
