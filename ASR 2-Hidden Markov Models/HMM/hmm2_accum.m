function      [pout,newmean,newvar,atot,A_hat,Pi_hat] = ...
	  hmm2_accum(hparm, event_data, starts, nsamp, ...
	      newmean, newvar,atot,A_hat, ...
          Pi_hat,rec_wts,powr,var_fac);

% function      [pout,newmean,newvar,atot,A_hat,Pi_hat] = ...
%	  hmm2_accum(hparm, event_data, istart, nsamp, ...
%	      newmean, newvar,atot,A_hat,Pi_hat,data_wts,hparm.powr,var_fac);
      
if(nargin <= 10 | isempty(powr) ), powr=2; end;
if(nargin <= 11 | isempty(var_fac) ) , var_fac=1; end;
[dim,ntot]=size(event_data);
nrec = length(starts);
m=length(nsamp);
if(m ~= nrec), error('mismatch nsamp vs. starts'); end;
if(nargin <= 9 | isempty(rec_wts) ), rec_wts = ones(nrec,1); end;
nmax = max(nsamp);
if(ntot ~= sum(nsamp)), error('ntot ~= sum(nsamp)'); end;
if(starts(1) ~= 1), error('starts(1) ~= 1'); end;

if(nargin > 4),
end;

pout=zeros(1,nrec);
for irec=1:nrec,
       idx = starts(irec) + [0:nsamp(irec)-1];
       x = event_data(:,idx);
       [newmean,newvar,atot,A_hat,Pi_hat,pout2,gamma,l] = ...
                hmm_accum_rec(hparm,x,newmean, newvar,atot,A_hat,Pi_hat,pout(irec));
       pout(irec)=pout2;
end;

return


function   [newmean,newvar,atot,A_hat,Pi_hat,pout,gamma,l] = ...
                hmm_accum_rec(hparm,x,newmean, newvar,atot,A_hat,Pi_hat,pout);

      BAD = -9e37;
	  N=hparm.N;
	  [dim,nsamp]=size(x);
	
	  % Compute probabilities
	  pk={};
	  l= BAD*ones(nsamp,N);
	  for istate=1:N,
	            gparm = hparm.pdf(istate);
	            nmode = length(gparm.modes);
	            wts = [gparm.modes.weight]';
	
	            % log-PDF for each mode
	            lg=lqr_evp(gparm,x,1);
	
	            % exponentiate
	            if(nmode>1),
	               lg = lg + log(repmat(wts(:)',nsamp,1));
	               mx = max(lg')';
	               lg=lg-repmat(mx,1,nmode);
	               g = exp(lg);
	               psum = sum(g')';
	               lg = log(psum) + mx;
	               p = g ./ repmat(psum,1,nmode);
               else,
%  	              p = ones(ngood,1); % edit SM
              end;
	            pk{istate}=p;
	            l(:,istate)=lg;
            end; 
	   if(N>1),lmax=max(l')'; else, lmax=l; end;
	   igood=find(lmax > BAD/2);
	   gamma=l;

       [alphas,alognorm,betas,blognorm,tgamma,A_hat] = ...
	            baumwelsh(l(igood,:),lmax(igood),hparm.A,hparm.Pi,A_hat);
	   gamma(igood,:)=tgamma;
	
	   Pi_hat = Pi_hat + gamma(1,:)';
	   pout = pout + log(sum(alphas(end,:)))+alognorm(end);
	
	   for istate=1:N,
	            gparm = hparm.pdf(istate);
	            nmode = length(gparm.modes);
	            wts = [gparm.modes.weight]';
	            datawts = pk{istate} .* repmat(gamma(:,istate),1,nmode);
	            for imode=1:nmode,
	               atot{istate}(imode)= atot{istate}(imode)+sum(datawts(:,imode));
	               newmean{istate}(:,imode) =  newmean{istate}(:,imode) ...
	                         + (datawts(:,imode)'*x')';
	               y = x' - repmat(gparm.modes(imode).mean',nsamp,1);
	               y = y .* repmat(sqrt( datawts(:,imode)),1,dim);
	               tmp = y'*y;
	               newvar{istate}(:,imode) =  newvar{istate}(:,imode) + tmp(:);
               end;
           end;
	return
