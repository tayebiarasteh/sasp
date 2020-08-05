function l = get_probs(hparm,x);

    BAD= -9e37;

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
 	              p = ones(ngood,1);
              end;
	            pk{istate}=p;
	            l(:,istate)=lg;
   end; 
return

