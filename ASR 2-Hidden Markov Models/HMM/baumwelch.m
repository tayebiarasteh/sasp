function [alphas,alognorm,betas,blognorm,gamma,Ahat]= baumwelch(lbtj,lmax,A,Pi,Ahat);
% function [alphas,alognorm,betas,blognorm,gamma,Ahat]= baumwelch(lbtj,lmax,A,Pi,Ahat);
%
%  Dr. Paul M Baggenstoss
%  NUWC, Newport, RI 
%  p.m.baggenstoss@ieee.org
%
% This routine inputs the raw observation log-probabilities
% (lbtj - which stands for log b[t]_j - the log observation probability
% at time t for state).  It implements the forward and backward 
% procedures as well as accumulating re-estimated transition matrix 
% (Ahat) and gets state probabilities given the data (gamma).
% It maintains log-scaling quantities alognorm, blognorm
% It is ised by hmm_accum.
%

% MATLAB version of baumwelch.  Not verctorized (very slow)

   [nsamp,N]=size(lbtj);
   alognorm=zeros(nsamp,1);

   aout = lbtj(1,:)';
   max1=max(aout);
   alognorm(1)=max1;
   aout = Pi .* exp( aout -max1 );
   alphas(1,:) = aout';

   for t=2:nsamp,
         tmp = (alphas(t-1,:) * A)';
         aout = lbtj(t,:)';
         max1=max(aout);
         aout = tmp .* exp(aout-max1);
         max2=max(max(aout),eps);
         aout = aout / max2;
         alognorm(t) = max1 + log(max2) + alognorm(t-1);
         alphas(t,:) = aout';
   end;

   if(nargout<3), return; end;

   blognorm=zeros(nsamp,1);
   bout=ones(N,1);
   blognorm(nsamp)=0;
   betas(nsamp,:)=bout';


   for t=nsamp-1 : -1 : 1,
      tmp = lbtj(t+1,:);
      max1=max(tmp);
      tmp = exp(tmp-max1).* betas(t+1,:);
      bout = A * tmp';
      max2=max(max(bout),eps);
      bout = bout / max2;
      betas(t,:)=bout';
      blognorm(t) = log(max2) + blognorm(t+1) + max1;
   end;

   psi=zeros(N,N);
   for t=1:nsamp-1,
      psisum=0;
      for i=1:N,
          for j=1:N,
                 ptmp= alphas(t,i)*betas(t+1,j) *A(i,j)*exp(lbtj(t+1,j)-lmax(t+1));
                 psi(i,j)=ptmp;
                 psisum =psisum + ptmp;
          end;
       end;
       Ahat = Ahat + psi/psisum;
   end;

   gamma = zeros(nsamp,N);
   for t=1:nsamp,
      gammasum=0;
      for i=1:N,
            gamma(t,i)=alphas(t,i)*betas(t,i);
            gammasum = gammasum + gamma(t,i);
      end;
      gamma(t,:) = gamma(t,:) / gammasum;
   end;

return

