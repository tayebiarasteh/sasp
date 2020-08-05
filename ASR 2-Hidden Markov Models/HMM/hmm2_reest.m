function [log_pdf_val, hparm] = hmm2_reest(hparm, event_data, istart, nsamp, NIT, data_wts);
% function [log_pdf_val, hparm] = hmm2_reest(hparm, event_data, istart, nsamp, NIT, data_wts);

global data_directory;
if(nargin < 6), data_wts=[]; end;

if( NIT > 0 ),	% do training
	
	BAD= -9e37;
	BIAS=0;
	plast= BAD;
	N=hparm.N;
	[dim,ntot]=size(event_data);
	nrec=length(istart);
	
	pout_increase_per_thousand_threshold = .1;
	n_pout_no_increase = 0;
	max_n_pout_no_increase = 10;
	
	% determine current n_modes
	for istate=1:N,
		current_n_modes_a(istate) = length(hparm.pdf(istate).modes);
	end;
	
	% if older version of parameter without "powr", set field
	if(~isfield(hparm,'powr')), hparm.powr = 2; end;
	
	% The input "powr" is the power of the exponent in the PDF, normally
	% equal to 2 (for Gaussian).  Higher numbers make the kernels roll off faster,
	% with powr=infinity producing a flat elliptically-shaped table kernel.
	% Smaller numbers produce heavier tails. Powr=1 is exponential kernel.
	%
	% The input var_fac is a kind of annealing term which is normally 1,
	% but if higher will expand the variances of all the kernels by that factor.
	%
	
	var_fac = 1.0;
	
	if(0),
		BIAS=1;
		min_std_org = [hparm.pdf(1).features.min_std]';
		min_std = min_std_org*0+.001;
		dim=length(min_std_org);
		for i=1:N,
			for j=1:dim,
				hparm.pdf(i).features(j).min_std = min_std(j);
			end;
		end;
	end;
	
	for it=1:NIT,
		
		% zero accumulators
		A_hat=zeros(N,N);
		Pi_hat=zeros(N,1);
		for istate=1:N,
			nmode = length(hparm.pdf(istate).modes);
			newmean{istate}=zeros(dim,nmode);
			newvar{istate}=zeros(dim*dim,nmode);
			atot{istate}=zeros(nmode,1);
		end;
		
		if(N==1),
			A_hat = 1.0;
			Pi_hat = 1.0;
		end;
		
		% Accumulate the record
		[pout,newmean,newvar,atot,A_hat,Pi_hat] = ...
			hmm2_accum(hparm, event_data, istart, nsamp, ...
			newmean, newvar,atot,A_hat,Pi_hat,data_wts,hparm.powr,var_fac);
		
		
		% sometimes, mex function returns values that are less than realmin
		% but greater than 0, this causes problems for Matlab
		for istate=1:N,
			% in particular we have seen this problem with newvar
			i_weird = find((newvar{istate}(:) > 0.0 & newvar{istate}(:) < realmin));
			newvar{istate}(i_weird) = 0.0;
		end;
		
		% Write a dump file if there is a NaN
		if(any(isnan(pout(:)))),
			fprintf('bad pout value\n')
			logfile_dir = sprintf('%s/logfiles/hmm2_reest', data_directory);
			mkdir_if_needed(logfile_dir);
			logfile_name = sprintf('%s/%f.mat', logfile_dir, now);
			fprintf('see log file: %s\n', logfile_name);
			eval(sprintf('save %s hparm event_data istart nsamp newmean newvar atot A_hat Pi_hat', ...
				logfile_name));
			break;
		end;
		
		
		
		% determine relative weight of each state
		for istate=1:N,
			state_wts(istate) = sum(atot{istate});
		end;
		state_wts=state_wts/sum(state_wts);
		
		% Prune very small states
		if(any(state_wts < 1e-10)),
			fprintf('Warning: pruning some states\n');
			igood = find(state_wts >= 1e-10);
			A_hat = A_hat(igood,igood);
			Pi_hat=Pi_hat(igood);
			atot= atot(igood);
			current_n_modes_a= current_n_modes_a(igood);
			newmean= newmean(igood);
			newvar= newvar(igood);
			state_wts= state_wts(igood);
			N=length(igood);
			hparm.pdf=hparm.pdf(igood);
			hparm.A=hparm.A(igood,igood);
			hparm.Pi=hparm.Pi(igood);
			hparm.N=N;
		end;
		
		% Update Gaussian mixtures
		for istate=1:N,
			gparm = hparm.pdf(istate);
			dim=length(gparm.features);
			gparm = gmix_norm(gparm,newmean{istate}, ...
				newvar{istate},atot{istate}, BIAS);
			wts=[gparm.modes.weight]';
			
			% determine effective number of (segment) samples for each mode
			nsamp_mode = wts * (state_wts(istate) *  ntot);
			
			% Prune the very very small weights
			thresh = min(0.00001,max(wts));
			igood = find(wts >= thresh);
			if(length(igood)<length(wts)),
				gparm.modes=gparm.modes(igood);
				fprintf('hmm2: pruning very small modes from %d to %d\n', ...
					length(wts),length(igood));
				wts=[gparm.modes.weight]';
			end;
			
			% Prune one small weight
			[mw,iw]=min(nsamp_mode);
			if(mw<dim & length(wts)>1 ),
				igood = find([1:length(gparm.modes)] ~= iw);
				gparm.modes=gparm.modes(igood);
				fprintf('hmm2: pruning from %d to %d\n',length(wts),length(igood));
				wts=[gparm.modes.weight]';
			end;
			hparm.pdf(istate)=gparm;
		end;
		
		% Normalize A_hat and Pi_hat
		if(isfield(hparm,'mask')),
			A_hat = A_hat .* hparm.mask;
		end;
		asum=sum(A_hat,2);
		if(any(asum==0)),
			error('zeros in asum');
		end;
		A_hat=A_hat ./repmat(asum,1,N);
		Pi_hat=Pi_hat/sum(Pi_hat);
		hparm.A=A_hat;
		hparm.Pi=Pi_hat;
		
		fprintf('Iteration %g, Pout=%g: %g',it,sum(pout),sum(pout)-plast);
		if( sum(pout) - plast < pout_increase_per_thousand_threshold*(ntot/1000)),
			n_pout_no_increase = n_pout_no_increase + 1;
			if( n_pout_no_increase >= max_n_pout_no_increase ),
				break;
			end;
		end;
		plast=sum(pout);
		
		% check whether any of the n_modes changed
		n_modes_changed = 0;
		for istate=1:N,
			if( current_n_modes_a(istate) ~= length(hparm.pdf(istate).modes) ),
				n_modes_changed = 1;
				current_n_modes_a(istate) = length(hparm.pdf(istate).modes);
			end;
		end;
		if( n_modes_changed ),
			n_pout_no_increase = ceil(n_pout_no_increase / 2);
		end;
		
		fprintf('  n_no_incr=%d\n',  n_pout_no_increase);
		
		plotHMM(event_data,hparm,it,n_pout_no_increase >= max_n_pout_no_increase-1);    % added SM
		pause(0.05); % added SM
		
	end; % it
	fprintf('  n_no_incr=%d\n',  n_pout_no_increase);
	
else	% NIT = 0, eval only
		
		% new SM-----------
		% zero accumulators
		N=hparm.N;
		min_std_org = [hparm.pdf(1).features.min_std]';
		dim=length(min_std_org);
		A_hat=zeros(N,N);
		Pi_hat=zeros(N,1);
		for istate=1:N,
			nmode = length(hparm.pdf(istate).modes);
			newmean{istate}=zeros(dim,nmode);
			newvar{istate}=zeros(dim*dim,nmode);
			atot{istate}=zeros(nmode,1);
		end;
		var_fac = 1.0;		
		if(N==1),
			A_hat = 1.0;
			Pi_hat = 1.0;
		end;
		
		% Accumulate the record
		[pout,~,~,~,~,~] = ...
			hmm2_accum(hparm, event_data, istart, nsamp, ...
			newmean, newvar,atot,A_hat,Pi_hat,data_wts,hparm.powr,var_fac);	
		%---------- new SM
		
% old --------
% 	% Accumulate the record
% 	pout = hmm2_accum(hparm, event_data, istart, nsamp);
% ---------old
	
	% Write a dump file if there is a NaN
	if(any(isnan(pout(:)))),
		fprintf('bad pout value\n')
		logfile_dir = sprintf('%s/logfiles/hmm2_reest', data_directory);
		mkdir_if_needed(logfile_dir);
		logfile_name = sprintf('%s/%f.mat', logfile_dir, now);
		fprintf('see log file: %s\n', logfile_name);
		eval(sprintf('save %s hparm event_data istart nsamp', logfile_name));
	end;
	
end;

log_pdf_val = pout;

return

function parm = gmix_norm(parm,newmean,newvar,atot, bias);
%function parm = gmix_norm(parm,newmean,newvar,atot, [bias]);
%
% Subrouting to complete Gaussian mixture update
% (see gmix_accum.m)
%
% Inputs:
%        parm     Input/Output parameters,
%        newmean  accumulated mean estimates
%        newvar   accumulated covariance estimates
%        atot     total weights used to normalize the above
%        bias     (optional) bias=1 for BIAS method, 0 for CONSTRAINT method
%
%  Dr. Paul M. Baggenstoss
%  Naval Undersea Warfare Center
%  Newport, RI
%  p.m.baggenstoss@ieee.org
%  Date              Reason
%  ------------      --------------
%  Nov 21, 1999      Initial release
%  Dec 7,  1999      MATLAB 5 Upgrade
%  Mar 3,  2000      Fixed bug - not normalizing min_std
%                     and use eig() in place of chol()
%  Mar 13,  2000     Fixed bug - could divide by 0 if atot(k)=0;

global data_directory;

nmode = length(parm.modes);
DIM = length(parm.features);
min_std = [parm.features.min_std]';

if(nargin < 5), bias=0; end;

assert(all(size(newmean) == [DIM,nmode]));
assert(all(size(newvar) == [DIM*DIM,nmode]));
assert(all(size(atot) == [nmode,1]));

tmpvar=zeros(DIM,DIM);
for k=1:nmode,
	
	if(atot(k) > 0) ,
		parm.modes(k).mean=newmean(:,k)/atot(k);
		
		tmpvar(:)=newvar(:,k)/atot(k);
		
		if(bias==0),
			if(any(isinf(tmpvar(:))) | any(isnan(tmpvar(:)))),
				fprintf('bad tmpvar value\n')
				logfile_dir = sprintf('%s/logfiles/gmix_norm', data_directory);
				mkdir_if_needed(logfile_dir);
				logfile_name = sprintf('%s/tmpvar.%f.mat', logfile_dir, now);
				fprintf('see log file: %s\n', logfile_name);
				eval(sprintf('save %s', logfile_name));
			end;
			% better not to use chol - in case of singular matrix
			[V,D]=eig(tmpvar);
			S=sqrt(abs(D));
			
			% constrains on the variances
			% first determine eigen decomp of C=R' * R = V * S^2 * V'
			S = diag(S);
			S = max(S,sqrt(diag( V' * diag(min_std.^2) * V )));
			tmpvar = chol(V * diag(S).^2 * V');
		else,
			
			for i=1:DIM;
				tmpvar(i,i)= tmpvar(i,i)+min_std(i).^2;
			end;
			tmpvar=chol(tmpvar);
		end;
		
		parm.modes(k).cholesky_covar = tmpvar;
	end;
	
	wts(k) = atot(k)/sum(atot);
	parm.modes(k).weight = wts(k) * (wts(k) > 0);
end;

%%%	function [newmean,newvar,atot,A_hat,Pi_hat,pout,gamma,l] = ...
%%%	        hmm_accum_(hparm,x,newmean, ...
%%%	                  newvar,atot,A_hat,Pi_hat,pout);
%%%
%%%	  BAD = -9e37;
%%%	  N=hparm.N;
%%%	  [dim,nsamp]=size(x);
%%%
%%%	  % Compute probabilities
%%%	  pk={};
%%%	  l= BAD*ones(nsamp,N);
%%%	  for istate=1:N,
%%%	            gparm = hparm.pdf(istate);
%%%	            nmode = length(gparm.modes);
%%%	            wts = [gparm.modes.weight]';
%%%
%%%	            % log-PDF for each mode
%%%	            lg=lqr_evp(gparm,x,1);
%%%
%%%	            % exponentiate
%%%	            if(nmode>1),
%%%	               lg = lg + log(repmat(wts(:)',nsamp,1));
%%%	               mx = max(lg')';
%%%	               lg=lg-repmat(mx,1,nmode);
%%%	               g = exp(lg);
%%%	               psum = sum(g')';
%%%	               lg = log(psum) + mx;
%%%	               p = g ./ repmat(psum,1,nmode);
%%%	            else,
%%%	               p = ones(ngood,1);
%%%	            end;
%%%	            pk{istate}=p;
%%%	            l(:,istate)=lg;
%%%	   end;
%%%	   if(N>1),lmax=max(l')'; else, lmax=l; end;
%%%	   igood=find(lmax > BAD/2);
%%%	   gamma=l;
%%%
%%%	   [alphas,alognorm,betas,blognorm,tgamma,A_hat] = ...
%%%	            baumwelsh(l(igood,:),lmax(igood),hparm.A,hparm.Pi,A_hat);
%%%	   gamma(igood,:)=tgamma;
%%%
%%%	   Pi_hat = Pi_hat + gamma(1,:)';
%%%	   pout = pout + log(sum(alphas(end,:)))+alognorm(end);
%%%
%%%	   for istate=1:N,
%%%	            gparm = hparm.pdf(istate);
%%%	            nmode = length(gparm.modes);
%%%	            wts = [gparm.modes.weight]';
%%%	            datawts = pk{istate} .* repmat(gamma(:,istate),1,nmode);
%%%	            for imode=1:nmode,
%%%	               atot{istate}(imode)= atot{istate}(imode)+sum(datawts(:,imode));
%%%	               newmean{istate}(:,imode) =  newmean{istate}(:,imode) ...
%%%	                         + (datawts(:,imode)'*x)';
%%%	               y = x - repmat(gparm.modes(imode).mean',nsamp,1);
%%%	               y = y .* repmat(sqrt( datawts(:,imode)),1,dim);
%%%	               tmp = y'*y;
%%%	               newvar{istate}(:,imode) =  newvar{istate}(:,imode) + tmp(:);
%%%	            end;
%%%	   end;
%%%	return
%%%

function plotHMM(x,parm,itNum,finished)
	figure(111);
	set(111,'position',[100 100 1200 800]);

	hmm2_view(parm,x,1,2);
	subplot(2,3,3)
	cla
	text(0.2,0.2,['Iteration ' num2str(itNum)],'FontSize',18)
	axis off
	subplot(2,3,2)
	plot(x(1,:),x(2,:),'b.');
	xlabel('OFFSET');
	ylabel('POWER');
	title('Training set');

	if itNum==1
		pause(1)
	end
	if finished
		subplot(2,3,3)
		cla
		text(0.2,0.2,'Finished!','FontSize',18)
		axis off
	end
		